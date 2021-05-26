"""Subroutines to normalize alternative parameterizations of the ODE system."""
from typing import Tuple

import numba
import numpy as np

from covid_model_seiir_pipeline.lib.ode.constants import (
    AGGREGATES,
    AGG_MAP,
    COMPARTMENTS,
    DEBUG,
    FIT_PARAMETERS,
    FORECAST_PARAMETERS,
    DISTRIBUTION_PARAMETERS,
    INFECTIOUS_WILD,
    INFECTIOUS_VARIANT,
    N_GROUPS,
    NEW_E,
    WANED,
    PARAMETERS,
    SUSCEPTIBLE_WILD,
    SUSCEPTIBLE_VARIANT_ONLY,
    REMOVED_WILD,
    REMOVED_VARIANT,
)


@numba.njit
def make_aggregates(y: np.ndarray, y_past: np.ndarray) -> np.ndarray:
    """Make total system aggregates for use in group systems.

    Parameters
    ----------
    y
        The state of the full system on the last ODE iteration.
    y_past

    Returns
    -------
    aggregates
        An array that with values that can be indexed by the AGGREGATES
        index mapping.

    """
    aggregates = np.zeros((len(AGGREGATES), y_past.shape[1] + 1))

    for target, compartments in AGG_MAP:
#        for group_y_past in np.split(y_past, N_GROUPS):
#            aggregates[target, :-1] = group_y_past[compartments, :].sum(axis=0)
        for group_y in np.split(y, N_GROUPS):
            aggregates[target, -1] = group_y[compartments].sum()

    if DEBUG:
        assert np.all(np.isfinite(aggregates))

    return aggregates


@numba.njit
def normalize_parameters(input_parameters: np.ndarray,
                         dist_parameters: np.ndarray,
                         aggregates: np.ndarray,
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
            beta_wild * aggregates[AGGREGATES.infectious_wild, -1]**alpha / aggregates[AGGREGATES.n_total, -1]
        )
        b_variant = (
            beta_variant * aggregates[AGGREGATES.infectious_variant, -1]**alpha / aggregates[AGGREGATES.n_total, -1]
        )

        new_e = np.zeros(len(NEW_E))
        new_e[NEW_E.wild] = b_wild * aggregates[AGGREGATES.susceptible_wild, -1]
        new_e[NEW_E.variant_naive] = b_variant * aggregates[AGGREGATES.susceptible_wild, -1]
        new_e[NEW_E.variant_reinf] = b_variant * aggregates[AGGREGATES.susceptible_variant_only, -1]
        new_e[NEW_E.total] = new_e.sum()

    else:
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
        ]), -1]

        si_wild = susceptible_wild * infectious_wild ** alpha
        si_variant_naive = scale * susceptible_wild * infectious_variant ** alpha
        si_variant_reinf = scale * susceptible_variant_only * infectious_variant ** alpha

        z = si_wild + si_variant_naive + si_variant_reinf

        new_e = np.zeros(len(NEW_E))
        new_e[NEW_E.wild] = si_wild / z * new_e_total
        new_e[NEW_E.variant_naive] = si_variant_naive / z * new_e_total
        new_e[NEW_E.variant_reinf] = si_variant_reinf / z * new_e_total
        new_e[NEW_E.total] = new_e_total

    waning_dist = dist_parameters[DISTRIBUTION_PARAMETERS.waning_immunity_time]
    waned = np.zeros(len(WANED))
    #newR_wild = aggregates[AGGREGATES.removed_wild, 1:] - aggregates[AGGREGATES.removed_wild, :-1]
    #newR_variant = aggregates[AGGREGATES.removed_variant, 1:] - aggregates[AGGREGATES.removed_variant, :-1]
    waned[WANED.wild] = 0. #(newR_wild[::-1] * waning_dist[1:]).sum()
    waned[WANED.variant] = 0. # (newR_variant[::-1] * waning_dist[1:]).sum()

    if DEBUG:
        assert np.all(np.isfinite(params))
        assert np.all(np.isfinite(vaccines))
        assert np.all(np.isfinite(new_e))

    return params, vaccines, new_e, waned
