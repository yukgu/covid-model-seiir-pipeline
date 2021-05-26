import itertools
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scipy.stats
import tqdm

from covid_model_seiir_pipeline.lib import (
    math,
    ode,
)
from covid_model_seiir_pipeline.pipeline.regression.model.containers import (
    ODEParameters,
)


def prepare_ode_fit_parameters(past_infections: pd.Series,
                               population: pd.DataFrame,
                               rhos: pd.DataFrame,
                               vaccinations: pd.DataFrame,
                               sampled_params: Dict[str, pd.Series]) -> ODEParameters:
    past_index = past_infections.index
    population_low_risk, population_high_risk = split_population(past_index, population)

    vaccinations = math.adjust_vaccinations(vaccinations)
    vaccinations = {k: v.rename(k).reindex(past_index, fill_value=0.) for k, v in vaccinations.items()}

    return ODEParameters(
        new_e=past_infections,
        population_low_risk=population_low_risk,
        population_high_risk=population_high_risk,
        rho=rhos['rho'].reindex(past_index, fill_value=0.0),
        rho_variant=rhos['rho_variant'].reindex(past_index, fill_value=0.0),
        rho_b1617=rhos['rho_b1617'].reindex(past_index, fill_value=0.0),
        **sampled_params,
        **vaccinations,
    )


def split_population(past_index: pd.Index, population: pd.DataFrame):
    population_low_risk = (population[population['age_group_years_start'] < 65]
                           .groupby('location_id')['population']
                           .sum()
                           .reindex(past_index, level='location_id'))
    population_high_risk = (population[population['age_group_years_start'] >= 65]
                            .groupby('location_id')['population']
                            .sum()
                            .reindex(past_index, level='location_id'))
    return population_low_risk, population_high_risk


def sample_params(past_index: pd.Index,
                  param_dict: Dict,
                  params_to_sample: List[str]) -> Dict[str, pd.Series]:
    sampled_params = {}
    for parameter in params_to_sample:
        param_spec = param_dict[parameter]
        if isinstance(param_spec, (int, float)):
            value = param_spec
        else:
            value = np.random.uniform(*param_spec)

        sampled_params[parameter] = pd.Series(
            value,
            index=past_index,
            name=parameter,
        )

    return sampled_params


def clean_infection_data_measure(infection_data: pd.DataFrame, measure: str) -> pd.Series:
    """Extracts measure, drops nulls, adds a leading zero.

    Infections and deaths have a non-overlapping past index due to the way
    the infections ES is built. This function, pulls out a measure, drops
    nulls from the non-overlaping region, and then pads the front of the
    series with a 0 so that the resulting series has the property:

        s == s.groupby('location_id').cumsum().groupby('location_id').diff().fillna(0)

    which is to say we can preserve the counts under conversions between daily
    and cumulative space.

    """
    data = infection_data[measure].dropna()
    min_date = data.reset_index().groupby('location_id').date.min()
    prepend_date = min_date - pd.Timedelta(days=1)
    prepend_idx = prepend_date.reset_index().set_index(['location_id', 'date']).index
    prepend = pd.Series(0., index=prepend_idx, name=measure)
    return data.append(prepend).sort_index()


def run_ode_fit(ode_parameters: ODEParameters, progress_bar: bool) -> Tuple[pd.DataFrame, pd.DataFrame]:
    fit_results = []
    ode_parameter_list = tqdm.tqdm(list(ode_parameters), disable=not progress_bar)
    for location_id, location_params in ode_parameter_list:
        loc_fit_results = run_loc_ode_fit(
            location_params
        )
        loc_fit_results['location_id'] = location_id
        loc_fit_results = loc_fit_results.set_index(['location_id', 'date'])

        fit_results.append(loc_fit_results)
    fit_results = pd.concat(fit_results).sort_index()
    return fit_results[['beta', 'beta_wild', 'beta_variant']], fit_results[[c for c in fit_results if 'beta' not in c]]


def run_loc_ode_fit(ode_parameters: ODEParameters) -> pd.DataFrame:
    # Filter out early dates with few infections
    # to reduce noise in the past fit from leaking into the beta regression.
    infections = ode_parameters.new_e
    full_index = infections.index
    infections = filter_to_epi_threshold(infections)
    ode_parameters = ode_parameters.reindex(infections.index)

    date = pd.Series(infections.index.values)
    t = (date - date.min()).dt.days.values
    obs = infections.values

    pop_groups = {
        'lr': ode_parameters.population_low_risk.iloc[0],
        'hr': ode_parameters.population_high_risk.iloc[0]
    }
    pop = sum(pop_groups.values())

    all_group_compartments = list(ode.COMPARTMENTS._fields) + list(ode.TRACKING_COMPARTMENTS._fields)
    aggregate_compartments = list(ode.AGGREGATES._fields)
    system_size = len(pop_groups)*len(all_group_compartments)
    initial_condition = np.zeros(system_size + len(aggregate_compartments))
    params = [ode_parameters.to_df().loc[:, list(ode.PARAMETERS._fields) + list(ode.FIT_PARAMETERS._fields)].values]

    for i, (risk_group, group_pop) in enumerate(pop_groups.items()):
        offset = i * len(all_group_compartments)
        initial_condition = _distribute_initial_condition(
            infections=obs[0],
            group_pop=group_pop,
            total_pop=pop,
            alpha=ode_parameters.alpha[0],
            initial_condition=initial_condition,
            offset=offset,
        )
        initial_condition[system_size + ode.AGGREGATES.susceptible_wild] += (
            initial_condition[offset + ode.COMPARTMENTS.S]
        )
        initial_condition[system_size + ode.AGGREGATES.infectious_wild] += (
            initial_condition[offset + ode.COMPARTMENTS.I1]
        )
        initial_condition[system_size + ode.AGGREGATES.n_total] += group_pop
        params.append(
            ode_parameters.get_vaccinations(ode.VACCINE_TYPES._fields, risk_group)
        )

    params = np.hstack(params).T
    dist_params = [get_waning_dist(ode_parameters)]
    result, *_ = math.solve_dde(
        system=ode.fit_system,
        t=t,
        init_cond=initial_condition,
        params=params,
        dist_params=dist_params,
    )
    compartments_by_group = [f'{compartment}_{risk_group}' for risk_group, compartment
                             in itertools.product(pop_groups, all_group_compartments)]
    components = pd.DataFrame(
        data=result.T,
        columns=compartments_by_group + aggregate_compartments,
    )
    components['date'] = date
    components = components.set_index('date')

    new_e_wild = components.filter(like='NewE_wild').sum(axis=1)
    new_e_variant = components.filter(like='NewE_variant').sum(axis=1)

    s_wild = components['susceptible_wild']
    s_variant_only = components['susceptible_variant_only']
    s_variant = s_wild + s_variant_only

    i_wild = components['infectious_wild']
    i_variant = components['infectious_variant']

    disease_density_wild = s_wild * i_wild**ode_parameters.alpha.values / pop
    beta_wild = (new_e_wild.diff() / disease_density_wild).reindex(full_index)
    disease_density_variant = s_variant * i_variant**ode_parameters.alpha.values / pop
    beta_variant = (new_e_variant.diff() / disease_density_variant).reindex(full_index)

    components = components.reindex(full_index, fill_value=0.)
    for risk_group, group_pop in pop_groups.items():
        components.loc[components[f'S_{risk_group}'] == 0, f'S_{risk_group}'] = group_pop

    components['beta_wild'] = beta_wild
    components['beta_variant'] = beta_variant

    rho = ode_parameters.rho
    rho_b1617 = ode_parameters.rho_b1617
    kappa = ode_parameters.kappa
    phi = ode_parameters.phi
    psi = ode_parameters.psi
    beta1 = beta_wild / (1 + kappa * rho)
    beta2 = beta_variant / (1 + kappa * (phi * (1 - rho_b1617) + rho_b1617 * psi))
    components['beta'] = (i_wild * beta1 + (i_variant * beta2).fillna(0)) / (i_wild + i_variant)

    return components.reset_index()


def _distribute_initial_condition(infections: float, group_pop: float, total_pop: float, alpha: float,
                                  initial_condition: np.ndarray, offset: int = 0) -> np.ndarray:
    new_e = infections * group_pop / total_pop
    initial_condition[offset + ode.COMPARTMENTS.S] = (
        group_pop - new_e - (new_e / 5)**(1.0 / alpha)
    )
    initial_condition[offset + ode.COMPARTMENTS.E] = new_e
    initial_condition[offset + ode.COMPARTMENTS.I1] = (new_e / 5) ** (1.0 / alpha)
    return initial_condition


def filter_to_epi_threshold(infections: pd.Series,
                            threshold: float = 50.) -> pd.Series:
    # noinspection PyTypeChecker
    start_date = infections.loc[threshold <= infections].index.min()
    while infections.loc[start_date:].count() <= 2:
        threshold *= 0.5
        # noinspection PyTypeChecker
        start_date = infections.loc[threshold <= infections].index.min()
        if threshold < 1e-6:
            start_date = infections.index.min()
            break
    return infections.loc[start_date:]


def get_waning_dist(ode_parameters: ODEParameters):
    start = ode_parameters.waning_start.mean()
    mean = ode_parameters.waning_mean.mean()
    var = ode_parameters.waning_sd.mean()**2
    return scipy.stats.gamma(loc=start, scale=var/mean, a=mean**2/var)
