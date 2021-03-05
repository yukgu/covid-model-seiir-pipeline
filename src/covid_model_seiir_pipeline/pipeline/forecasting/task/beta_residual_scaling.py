import functools
import multiprocessing
from typing import Dict, List
from pathlib import Path

import click
import pandas as pd
import numpy as np

from covid_model_seiir_pipeline.lib import (
    cli_tools,
    static_vars,
)
from covid_model_seiir_pipeline.pipeline.forecasting.specification import (
    ForecastSpecification,
    ScenarioSpecification,
    FORECAST_JOBS,
)
from covid_model_seiir_pipeline.pipeline.forecasting.data import ForecastDataInterface


logger = cli_tools.task_performance_logger


def run_compute_beta_scaling_parameters(forecast_version: str, scenario: str):
    """Pre-compute the parameters for rescaling predicted beta and write out.

    The predicted beta has two issues we're attempting to adjust.

    The first issue is that the predicted beta does not share the same
    y-intercept as the beta we fit from the ODE step. To fix this,
    we compute an initial scale factor `scale_init` that shifts the whole
    time series so that it lines up with the fit beta at the day of transition.

    The second issue is the long-term scaling of the predicted beta. If we use
    our initial scaling computed from the transition day, we bias the long
    range forecast heavily towards whatever unexplained variance appears
    in the regression on the transition day. Instead we use an average of
    the residual from the regression over some period of time in the past as
    the long range scaling.

    For locations with a large number of deaths, the average of the residual
    mean across draws represents concrete unexplained variance that is biased
    in a particular direction. For locations with a small number of deaths,
    the average of the residual mean is susceptible to a lot of noise week to
    week and so we re-center the distribution of the means about zero. For
    locations in between, we linearly scale this shift between the distribution
    of residual means centered  around zero and the distribution of
    residual means centered around the actual average residual mean.
    We call this final shift `scale_final`

    Parameters
    ----------
    forecast_version
        The path to the forecast version to run this process for.
    scenario
        Which scenario in the forecast version to run this process for.

    Notes
    -----
    The last step of average residual re-centering requires information
    from all draws.  This is the only reason this process exists as separate
    from the main forecasting process.

    """
    logger.info(f"Computing beta scaling parameters for forecast "
                f"version {forecast_version} and scenario {scenario}.", context='setup')

    forecast_spec: ForecastSpecification = ForecastSpecification.from_path(
        Path(forecast_version) / static_vars.FORECAST_SPECIFICATION_FILE
    )
    num_cores = forecast_spec.workflow.task_specifications[FORECAST_JOBS.scaling].num_cores
    data_interface = ForecastDataInterface.from_specification(forecast_spec)

    logger.info('Loading input data.', context='read')
    scenario_spec = forecast_spec.scenarios[scenario]

    logger.info('Computing scaling parameters.', context='compute')
    scaling_data = compute_initial_beta_scaling_parameters(scenario_spec, data_interface, num_cores)

    logger.info('Writing scaling parameters to disk.', context='write')
    write_out_beta_scale(scaling_data, scenario, data_interface, num_cores)

    logger.report()


def compute_initial_beta_scaling_parameters(scenario_spec: ScenarioSpecification,
                                            data_interface: ForecastDataInterface,
                                            num_cores: int) -> List[pd.DataFrame]:
    # Serialization is our bottleneck, so we parallelize draw level data
    # ingestion and computation across multiple processes.
    _runner = functools.partial(
        compute_initial_beta_scaling_parameters_by_draw,
        scenario_spec=scenario_spec,
        data_interface=data_interface
    )
    draws = list(range(data_interface.get_n_draws()))
    with multiprocessing.Pool(num_cores) as pool:
        scaling_data = list(pool.imap(_runner, draws))
    return scaling_data


def compute_initial_beta_scaling_parameters_by_draw(draw_id: int,
                                                    scenario_spec: ScenarioSpecification,
                                                    data_interface: ForecastDataInterface) -> pd.DataFrame:

    # Construct a list of pandas Series indexed by location and named
    # as their column will be in the output dataframe. We'll append
    # to this list as we construct the parameters.
    draw_data = []

    beta_scaling = scenario_spec.beta_scaling
    betas = data_interface.load_betas(draw_id)
    covariates = data_interface.load_covariates(scenario_spec)
    variant_prevalence = covariates[['variant_prevalence_B1351', 'variant_prevalence_P1']].sum(axis=1)
    average_over_min_min = variant_prevalence[variant_prevalence > 0].reset_index().groupby('location_id').date.min()

    # Select out the transition day to compute the initial scaling parameter.
    beta_transition = betas.groupby('location_id').last()

    draw_data.append(beta_transition['beta'].rename('fit_final'))
    draw_data.append(beta_transition['beta_hat'].rename('pred_start'))
    draw_data.append((beta_transition['beta'] / beta_transition['beta_hat']).rename('scale_init'))

    # Compute the beta residual mean for our parameterization and hang on
    # to some ancillary information that may be useful for plotting/debugging.
    rs = np.random.RandomState(draw_id)

    a = rs.randint(1, beta_scaling['average_over_min'])
    b = rs.randint(a + 21, beta_scaling['average_over_max'])

    draw_data.append(pd.Series(a, index=beta_transition.index, name='history_days_start'))
    draw_data.append(pd.Series(b, index=beta_transition.index, name='history_days_end'))
    draw_data.append(pd.Series(beta_scaling['window_size'], index=beta_transition.index, name='window_size'))

    log_beta_residual_mean = (np.log(betas['beta'] / betas['beta_hat'])
                              .groupby(level='location_id')
                              .apply(lambda x: x.iloc[-b: -a].mean())
                              .rename('log_beta_residual_mean'))
    draw_data.append(log_beta_residual_mean)
    draw_data.append(pd.Series(draw_id, index=beta_transition.index, name='draw'))

    return pd.concat(draw_data, axis=1)


def write_out_beta_scale(beta_scales: List[pd.DataFrame],
                         scenario: str,
                         data_interface: ForecastDataInterface,
                         num_cores: int) -> None:
    _runner = functools.partial(
        write_out_beta_scales_by_draw,
        data_interface=data_interface,
        scenario=scenario
    )
    with multiprocessing.Pool(num_cores) as pool:
        pool.map(_runner, beta_scales)


def write_out_beta_scales_by_draw(beta_scales: pd.DataFrame,
                                  data_interface: ForecastDataInterface,
                                  scenario: str) -> None:
    # Compute these draw specific parameters now that we have the offset.
    beta_scales['scale_final'] = np.exp(beta_scales['log_beta_residual_mean'])
    draw_id = beta_scales['draw'].iat[0]
    data_interface.save_beta_scales(beta_scales, scenario, draw_id)


@click.command()
@cli_tools.with_task_forecast_version
@cli_tools.with_scenario
@cli_tools.add_verbose_and_with_debugger
def beta_residual_scaling(forecast_version: str, scenario: str,
                          verbose: int, with_debugger: bool):
    cli_tools.configure_logging_to_terminal(verbose)
    run = cli_tools.handle_exceptions(run_compute_beta_scaling_parameters, logger, with_debugger)
    run(forecast_version=forecast_version,
        scenario=scenario)


if __name__ == '__main__':
    beta_residual_scaling()
