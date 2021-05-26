"""Subroutines to initialize escape variant invasion."""
import numba
import numpy as np

from covid_model_seiir_pipeline.lib.ode.constants import (
    AGGREGATES,
    COMPARTMENTS,
    DEBUG,
    PARAMETERS,
)


@numba.njit
def maybe_invade(group_y: np.ndarray,
                 aggregates: np.ndarray,
                 params: np.ndarray,
                 transition_map: np.ndarray) -> np.ndarray:
    """Cause escape variants to invade if the criteria is correct.

    Parameters
    ----------
    group_y
        The state of the group system on the last ODE step. The first
        elements in this array align with the :obj:`COMPARTMENTS` index
        map.
    aggregates
        An array that aligns with the :obj:`AGGREGATES` index map representing
        total system aggregates like the total infectious subpopulation.
    params
        An array that aligns with the :obj:`PARAMETERS` index map representing
        the standard set of ODE parameters.
    transition_map

    Returns
    -------
    group_dy
        The change in state of the group system on this ODE step adjusted
        for whether there is escape variant invasion this step.

    """
    # Short circuit if we don't have variant invasion this step.
    no_variant_present = params[PARAMETERS.rho_variant] < 0.01
    already_invaded = aggregates[AGGREGATES.infectious_variant] > 0.0
    if no_variant_present or already_invaded:
        return transition_map

    alpha, pi = params[PARAMETERS.alpha], params[PARAMETERS.pi]

    transition_map = _invade_compartment_subset(
        group_y,
        alpha, pi,
        COMPARTMENTS.S, COMPARTMENTS.E,
        COMPARTMENTS.E_variant, COMPARTMENTS.I1_variant,
        transition_map,
    )
    transition_map = _invade_compartment_subset(
        group_y,
        alpha, pi,
        COMPARTMENTS.S_u, COMPARTMENTS.E_u,
        COMPARTMENTS.E_variant_u, COMPARTMENTS.I1_variant_u,
        transition_map,
    )
    transition_map = _invade_compartment_subset(
        group_y,
        alpha, pi,
        COMPARTMENTS.S_p, COMPARTMENTS.E_p,
        COMPARTMENTS.E_variant_u, COMPARTMENTS.I1_variant_u,
        transition_map,
    )
    transition_map = _invade_compartment_subset(
        group_y,
        alpha, pi,
        COMPARTMENTS.S_pa, COMPARTMENTS.E_pa,
        COMPARTMENTS.E_variant_pa, COMPARTMENTS.I1_variant_pa,
        transition_map,
    )
    return transition_map


@numba.njit
def _invade_compartment_subset(group_y: np.ndarray,
                               alpha: float, pi: float,
                               susceptible: int, exposed: int,
                               exposed_variant: int, infectious1_variant: int,
                               transition_map: np.ndarray) -> np.ndarray:
    """Shift a small set of folks to the escape variant compartments."""
    # Shift at least 1 person to the escape variants if we can
    min_invasion = 1
    # Cap the the invasion so we don't take everyone. The handles corner cases
    # where variants invade just as vaccination is starting.
    max_invasion = 0.5 * group_y[susceptible]
    delta = min(max(pi * group_y[exposed], min_invasion), max_invasion)

    # Set the boundary condition so that the initial beta for the escape
    # variant starts at 5 (for consistency with ancestral type invasion).
    transition_map[susceptible, exposed_variant] += delta
    transition_map[susceptible, infectious1_variant] += (delta / 5)**(1 / alpha)

    if DEBUG:
        assert np.all(np.isfinite(transition_map))

    return transition_map
