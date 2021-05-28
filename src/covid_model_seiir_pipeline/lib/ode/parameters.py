"""Subroutines to normalize alternative parameterizations of the ODE system."""
from typing import Tuple

import numba
import numpy as np

from covid_model_seiir_pipeline.lib.ode.constants import (
    AGGREGATES,
    DEBUG,
    FIT_PARAMETERS,
    FORECAST_PARAMETERS,
    DISTRIBUTION_PARAMETERS,
    NEW_E,
    WANED,
    PARAMETERS,
)
from covid_model_seiir_pipeline.lib.ode import (
    timer,
)


@numba.njit
def normalize_parameters(input_parameters: np.ndarray,
                         dist_parameters: np.ndarray,
                         aggregates: np.ndarray,
                         past_aggregates: np.ndarray,
                         forecast: bool) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Coerces all parameterizations of the ODE model to the same format.

    Parameters
    ----------
    input_parameters
        An array whose first several fields align with the :obj:`PARAMETERS`
        index mapping, whose middle elements align with either the
        :obj:`FIT_PARAMETERS` or :obj:`FORECAST_PARAMETERS`, and whose
        remaining elements are the :obj:`VACCINE_TYPES` for each of the
        :obj:`N_GROUPS`.
    dist_parameters
    aggregates
        An array that aligns with the :obj:`AGGREGATES` index mapping
        representing aggregates of the full ODE system like the total
        infectious population.
    past_aggregates
    forecast
        The `input_parameters` are for the forecast if `True`, otherwise
        `input_parameters` are for the beta fit.

    Returns
    -------
    parameters
        An array aligning with the :obj:`PARAMETERS` index mapping representing
        the standard parameterization of the ODE system.
    new_e
        An array aligning with the :obj:`NEW_E` index mapping representing
        new infections split by variant and source.
    vaccines
        An array with all vaccine types for all groups.

    """
    alpha = input_parameters[PARAMETERS.alpha]
    if forecast:
        param_size = len(PARAMETERS) + len(FORECAST_PARAMETERS)
        params, vaccines = input_parameters[:param_size], input_parameters[param_size:]

        beta_wild = params[FORECAST_PARAMETERS.beta_wild]
        beta_variant = params[FORECAST_PARAMETERS.beta_variant]
        params = params[np.array(PARAMETERS)]

        b_wild = (
            beta_wild * aggregates[AGGREGATES.infectious_wild]**alpha / aggregates[AGGREGATES.n_total]
        )
        b_variant = (
            beta_variant * aggregates[AGGREGATES.infectious_variant]**alpha / aggregates[AGGREGATES.n_total]
        )

        new_e = np.zeros(len(NEW_E))
        new_e[NEW_E.wild] = b_wild * aggregates[AGGREGATES.susceptible_wild]
        new_e[NEW_E.variant_naive] = b_variant * aggregates[AGGREGATES.susceptible_wild]
        new_e[NEW_E.variant_reinf] = b_variant * aggregates[AGGREGATES.susceptible_variant_only]
        new_e[NEW_E.total] = new_e.sum()

    else:
#        start = timer.timenow()
        param_size = len(PARAMETERS) + len(FIT_PARAMETERS)
        params, vaccines = input_parameters[:param_size], input_parameters[param_size:]

        new_e_total = params[FIT_PARAMETERS.new_e]
        kappa = params[FIT_PARAMETERS.kappa]
        rho = params[FIT_PARAMETERS.rho]
        rho_b1617 = params[FIT_PARAMETERS.rho_b1617]
        phi = params[FIT_PARAMETERS.phi]
        psi = params[FIT_PARAMETERS.psi]
        params = params[np.array(PARAMETERS)]

        scale_wild = (1 + kappa * rho)
        scale_variant = (1 + kappa * (phi * (1 - rho_b1617) + rho_b1617 * psi))
        scale = scale_variant / scale_wild

        susceptible_wild, susceptible_variant_only, infectious_wild, infectious_variant, n_total = aggregates[np.array([
            AGGREGATES.susceptible_wild, AGGREGATES.susceptible_variant_only,
            AGGREGATES.infectious_wild, AGGREGATES.infectious_variant,
            AGGREGATES.n_total,
        ])]

        si_wild = susceptible_wild * infectious_wild ** alpha
        si_variant_naive = scale * susceptible_wild * infectious_variant ** alpha
        si_variant_reinf = scale * susceptible_variant_only * infectious_variant ** alpha

        z = si_wild + si_variant_naive + si_variant_reinf

        new_e = np.zeros(len(NEW_E))
        new_e[NEW_E.wild] = si_wild / z * new_e_total
        new_e[NEW_E.variant_naive] = si_variant_naive / z * new_e_total
        new_e[NEW_E.variant_reinf] = si_variant_reinf / z * new_e_total
        new_e[NEW_E.total] = new_e_total
#        end = timer.timenow()
#        print('Normal normalization: ', end - start)

#    start = timer.timenow()
    waning_dist = dist_parameters[DISTRIBUTION_PARAMETERS.waning_immunity_time]
    waned = np.zeros(len(WANED))
    wild_removed = np.append(past_aggregates[AGGREGATES.NewR_wild], aggregates[AGGREGATES.NewR_wild])
    variant_removed = np.append(past_aggregates[AGGREGATES.NewR_variant], aggregates[AGGREGATES.NewR_variant])
    newR_wild = wild_removed[1:] - wild_removed[:-1]
    newR_variant = variant_removed[1:] - variant_removed[:-1]
    waned[WANED.wild] = (newR_wild[::-1] * waning_dist[1:]).sum()
    waned[WANED.variant] = (newR_variant[::-1] * waning_dist[1:]).sum()
#    end = timer.timenow()
#    print('Waning: ', end - start)

    if DEBUG:
        assert np.all(np.isfinite(params))
        assert np.all(np.isfinite(vaccines))
        assert np.all(np.isfinite(new_e))

    return params, vaccines, new_e, waned
