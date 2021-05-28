import numba
import numpy as np

from covid_model_seiir_pipeline.lib.ode.constants import (
    AGGREGATES,
    AGG_MAP,
    COMPARTMENTS,
    DEBUG,
    REMOVED_WILD,
    REMOVED_VARIANT,
    TRACKING_COMPARTMENTS,
    VACCINE_TYPES,
)


@numba.njit
def compute_tracking_columns(group_dy: np.ndarray,
                             transition_map: np.ndarray,
                             vaccines_out: np.ndarray) -> np.ndarray:
    """Computes some aggregates that don't directly correspond to compartments.

    We need to do some accounting to capture, e.g., total new wild type
    infections. The only things we're tracking here are quantities that cannot
    be inferred from the time-series of the ODE compartments after the model
    is run. These quantities come about when there are compartments that have
    multiple in-paths.

    Parameters
    ----------
    group_dy
        An array who's first set of elements correspond to the
        :obj:`COMPARTMENTS` index map representing the resulting change in
        state from an ODE step. The remaining set of elements correspond to
        the :obj:`TRACKING_COMPARTMENTS` index map and should all be zero when
        passed to this function.
    transition_map
        An 2-d matrix representation whose rows and columns both correspond
        to the :obj:`COMPARTMENTS` index map. The value at position (i, j)
        in the matrix is the quantity of people moving from compartment i
        to compartment j.
    vaccines_out
        An array with a row for every compartment in `group_dy` representing
        who the vaccine is delivered to and a column for every vaccine
        efficacy in `group_vaccines` representing how effective the vaccines
        were (and therefore what compartment folks will move to on the next
        step).

    """
    # New wild type infections
    group_dy[TRACKING_COMPARTMENTS.NewE_wild] = (
        transition_map[COMPARTMENTS.S, COMPARTMENTS.E]
        + transition_map[COMPARTMENTS.S_u, COMPARTMENTS.E_u]
        + transition_map[COMPARTMENTS.S_p, COMPARTMENTS.E_p]
        + transition_map[COMPARTMENTS.S_pa, COMPARTMENTS.E_pa]
    )
    # New variant type infections
    group_dy[TRACKING_COMPARTMENTS.NewE_variant] = (
        transition_map[COMPARTMENTS.S, COMPARTMENTS.E_variant]
        + transition_map[COMPARTMENTS.S_variant, COMPARTMENTS.E_variant]
        + transition_map[COMPARTMENTS.S_u, COMPARTMENTS.E_variant_u]
        + transition_map[COMPARTMENTS.S_variant_u, COMPARTMENTS.E_variant_u]
        + transition_map[COMPARTMENTS.S_p, COMPARTMENTS.E_variant_u]
        + transition_map[COMPARTMENTS.S_pa, COMPARTMENTS.E_variant_pa]
        + transition_map[COMPARTMENTS.S_variant_pa, COMPARTMENTS.E_variant_pa]
        + transition_map[COMPARTMENTS.S_m, COMPARTMENTS.E_variant_pa]
    )
    # New wild type protected infections
    group_dy[TRACKING_COMPARTMENTS.NewE_p_wild] = (
        transition_map[COMPARTMENTS.S_p, COMPARTMENTS.E_p]
        + transition_map[COMPARTMENTS.S_pa, COMPARTMENTS.E_pa]
    )
    # New variant type protected infections
    group_dy[TRACKING_COMPARTMENTS.NewE_p_variant] = (
        transition_map[COMPARTMENTS.S_pa, COMPARTMENTS.E_variant_pa]
        + transition_map[COMPARTMENTS.S_variant_pa, COMPARTMENTS.E_variant_pa]
        + transition_map[COMPARTMENTS.S_m, COMPARTMENTS.E_variant_pa]
    )
    return group_dy


@numba.njit
def compute_aggregates(transition_map: np.ndarray, vaccines_out: np.ndarray):
    aggregates = np.zeros(len(AGGREGATES))
    # New variant type infections breaking through natural immunity
    aggregates[AGGREGATES.NewE_nbt] = (
        transition_map[COMPARTMENTS.S_variant, COMPARTMENTS.E_variant]
        + transition_map[COMPARTMENTS.S_variant_u, COMPARTMENTS.E_variant_u]
        + transition_map[COMPARTMENTS.S_variant_pa, COMPARTMENTS.E_variant_pa]
    )
    # New variant type infections breaking through vaccine immunity
    aggregates[AGGREGATES.NewE_vbt] = (
        transition_map[COMPARTMENTS.S_m, COMPARTMENTS.E_variant_pa]
    )

    aggregates[AGGREGATES.NewR_wild] = transition_map[:, REMOVED_WILD].sum()
    aggregates[AGGREGATES.NewR_variant] = transition_map[:, REMOVED_VARIANT].sum()

    aggregates[AGGREGATES.Waned_wild] = (
        transition_map[COMPARTMENTS.R, COMPARTMENTS.S]
        + transition_map[COMPARTMENTS.R_u, COMPARTMENTS.S_u]
        + transition_map[COMPARTMENTS.R_p, COMPARTMENTS.S_p]
        + transition_map[COMPARTMENTS.R_pa, COMPARTMENTS.S_pa]
    )

    aggregates[AGGREGATES.Waned_variant] = (
        transition_map[COMPARTMENTS.R_variant, COMPARTMENTS.S_variant]
        + transition_map[COMPARTMENTS.R_variant_u, COMPARTMENTS.S_variant_u]
        + transition_map[COMPARTMENTS.R_variant_pa, COMPARTMENTS.S_variant_pa]
    )

    aggregates[AGGREGATES.effective_vaccines] = (
        vaccines_out[:, VACCINE_TYPES.p].sum()
        + vaccines_out[:, VACCINE_TYPES.m].sum()
        + vaccines_out[:, VACCINE_TYPES.pa].sum()
        + vaccines_out[:, VACCINE_TYPES.ma].sum()
    )

    for target, compartments in AGG_MAP:
        compartments_out = transition_map[compartments, :].sum()
        compartments_in = transition_map[:, compartments].sum()
        aggregates[target] = compartments_in - compartments_out

    if DEBUG:
        assert np.all(np.isfinite(aggregates))

    return aggregates
