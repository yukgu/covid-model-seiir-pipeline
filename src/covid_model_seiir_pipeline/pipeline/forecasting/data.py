from functools import reduce
from itertools import product
from pathlib import Path
from typing import Dict, List, Union

from loguru import logger
import pandas as pd

from covid_model_seiir_pipeline.lib import (
    io,
)
from covid_model_seiir_pipeline.pipeline.regression import (
    RegressionDataInterface,
    RegressionSpecification,
    HospitalParameters,
)
from covid_model_seiir_pipeline.pipeline.forecasting.specification import (
    ForecastSpecification,
    ScenarioSpecification,
)
from covid_model_seiir_pipeline.pipeline.forecasting.model import (
    RatioData,
    HospitalMetrics,
    HospitalCorrectionFactors,
    HospitalCensusData,
    ScenarioData,
    VariantScalars,
)


class ForecastDataInterface:

    def __init__(self,
                 regression_root: io.RegressionRoot,
                 covariate_root: io.CovariateRoot,
                 forecast_root: io.ForecastRoot,
                 fh_subnationals: bool):
        self.regression_root = regression_root
        self.covariate_root = covariate_root
        self.forecast_root = forecast_root
        self.fh_subnationals = fh_subnationals

    @classmethod
    def from_specification(cls, specification: ForecastSpecification) -> 'ForecastDataInterface':
        # TODO: specify input format from config
        regression_root = io.RegressionRoot(specification.data.regression_version)
        covariate_root = io.CovariateRoot(specification.data.covariate_version)
        # TODO: specify output format from config.
        forecast_root = io.ForecastRoot(specification.data.output_root)

        return cls(
            regression_root=regression_root,
            covariate_root=covariate_root,
            forecast_root=forecast_root,
            fh_subnationals=specification.data.fh_subnationals
        )

    def make_dirs(self, **prefix_args):
        io.touch(self.forecast_root, **prefix_args)

    ############################
    # Regression paths loaders #
    ############################

    def get_n_draws(self) -> int:
        return self._get_regression_data_interface().get_n_draws()

    def load_location_ids(self) -> List[int]:
        return self._get_regression_data_interface().load_location_ids()

    def load_regression_coefficients(self, draw_id: int) -> pd.DataFrame:
        return self._get_regression_data_interface().load_regression_coefficients(draw_id=draw_id)

    def load_transition_date(self, draw_id: int) -> pd.Series:
        dates_df = self._get_regression_data_interface().load_date_file(draw_id=draw_id)
        dates_df['end_date'] = pd.to_datetime(dates_df['end_date'])
        transition_date = dates_df.set_index('location_id').sort_index()['end_date'].rename('date')
        return transition_date

    def load_beta_regression(self, draw_id: int) -> pd.DataFrame:
        beta_regression = self._get_regression_data_interface().load_regression_betas(draw_id=draw_id)
        beta_regression['date'] = pd.to_datetime(beta_regression['date'])
        return beta_regression

    def load_infection_data(self, draw_id: int) -> pd.DataFrame:
        infection_data = self._get_regression_data_interface().load_infection_data(draw_id=draw_id)
        infection_data['date'] = pd.to_datetime(infection_data['date'])
        return infection_data

    def load_beta_params(self, draw_id: int) -> Dict[str, float]:
        df = self._get_regression_data_interface().load_beta_param_file(draw_id=draw_id)
        return df.set_index('params')['values'].to_dict()

    def get_hospital_parameters(self) -> HospitalParameters:
        return self._get_regression_data_interface().load_specification().hospital_parameters

    def load_hospital_usage(self) -> HospitalMetrics:
        df = self._get_regression_data_interface().load_hospital_data(measure='usage')
        df['date'] = pd.to_datetime(df['date'])
        df = df.set_index(['location_id', 'date']).sort_index()
        return HospitalMetrics(**{metric: df[metric] for metric in df.columns})

    def load_hospital_correction_factors(self) -> HospitalCorrectionFactors:
        df = self._get_regression_data_interface().load_hospital_data(measure='correction_factors')
        df['date'] = pd.to_datetime(df['date'])
        df = df.set_index(['location_id', 'date']).sort_index()
        return HospitalCorrectionFactors(**{metric: df[metric] for metric in df.columns})

    def load_hospital_census_data(self) -> HospitalCensusData:
        return self._get_regression_data_interface().load_hospital_census_data()

    def load_ratio_data(self, draw_id: int) -> RatioData:
        return self._get_regression_data_interface().load_ratio_data(draw_id=draw_id)

    ##########################
    # Covariate data loaders #
    ##########################

    def check_covariates(self, scenarios: Dict[str, ScenarioSpecification]) -> List[str]:
        regression_spec = self._get_regression_data_interface().load_specification().to_dict()
        # Bit of a hack.
        forecast_version = str(self.covariate_root._root)
        regression_version = regression_spec['data']['covariate_version']
        if not forecast_version == regression_version:
            logger.warning(f'Forecast covariate version {forecast_version} does not match '
                           f'regression covariate version {regression_version}. If the two covariate '
                           f'versions have different data in the past, the regression coefficients '
                           f'used for prediction may not be valid.')

        regression_covariates = set(regression_spec['covariates'])

        for name, scenario in scenarios.items():
            if set(scenario.covariates).symmetric_difference(regression_covariates) > {'intercept'}:
                raise ValueError('Forecast covariates must match the covariates used in regression.\n'
                                 f'Forecast covariates:   {sorted(list(scenario.covariates))}.\n'
                                 f'Regression covariates: {sorted(list(regression_covariates))}.')

            if 'intercept' in scenario.covariates:
                # Shouldn't really be specified, but might be copied over from
                # regression.  No harm really in just deleting it.
                del scenario.covariates['intercept']

            for covariate, covariate_version in scenario.covariates.items():
                if not io.exists(self.covariate_root[covariate](covariate_scenario=covariate_version)):
                    raise FileNotFoundError(f'No {covariate_version} file found for covariate {covariate}.')

        return list(regression_covariates)

    def load_covariate(self, covariate: str, covariate_version: str, with_observed: bool = False) -> pd.DataFrame:
        location_ids = self.load_location_ids()
        covariate_df = io.load(self.covariate_root[covariate](covariate_scenario=covariate_version))
        covariate_df = self._format_covariate_data(covariate_df, location_ids, with_observed)
        covariate_df = (covariate_df
                        .rename(columns={f'{covariate}_{covariate_version}': covariate})
                        .loc[:, [covariate]])
        return covariate_df

    def load_covariates(self, scenario: ScenarioSpecification) -> pd.DataFrame:
        covariate_data = []
        for covariate, covariate_version in scenario.covariates.items():
            if covariate != 'intercept':
                covariate_data.append(self.load_covariate(covariate, covariate_version))
        covariate_data = reduce(lambda x, y: x.merge(y, left_index=True, right_index=True), covariate_data)
        return covariate_data.reset_index()

    #########################
    # Scenario data loaders #
    #########################

    def load_scenario_specific_data(self, scenario_spec: ScenarioSpecification,
                                    variant_prevalence: pd.Series) -> ScenarioData:
        if scenario_spec.system == 'vaccine':
            forecast_scenario = scenario_spec.system_params.get('forecast_version', 'reference')
            vaccinations = self.load_vaccine_info(f'vaccinations_{forecast_scenario}')
            location_ids = vaccinations.reset_index().location_id.tolist()

            variant_start_threshold = pd.Timestamp('2021-05-01')
            bad_variant_entrance_date = (variant_prevalence[variant_prevalence > 1]
                                         .reset_index()
                                         .groupby('location_id')
                                         .date
                                         .min())
            locs_with_bad_variant = (bad_variant_entrance_date[bad_variant_entrance_date < variant_start_threshold]
                                     .reset_index()
                                     .location_id
                                     .tolist())

            locs_without_bad_variant = list(set(location_ids).difference(locs_with_bad_variant))
            bad_variant_vaccinations = vaccinations.loc[locs_with_bad_variant]
            not_bad_variant_vaccinations = vaccinations.loc[locs_without_bad_variant]
            # FIXME: should get from population partition
            risk_groups = ['lr', 'hr']
            vaccination_groups = ['unprotected', 'protected', 'immune']
            out_cols = [f'{vaccination_group}_{risk_group}'
                        for vaccination_group, risk_group in product(vaccination_groups, risk_groups)]
            vaccinations = pd.DataFrame(columns=out_cols, index=vaccinations.index)
            for risk_group in risk_groups:
                vaccinations.loc[bad_variant_vaccinations.index, f'unprotected_{risk_group}'] = (
                    bad_variant_vaccinations[f'unprotected_{risk_group}']
                    + bad_variant_vaccinations[f'effective_protected_wildtype_{risk_group}']
                    + bad_variant_vaccinations[f'effective_wildtype_{risk_group}']
                )
                vaccinations.loc[bad_variant_vaccinations.index, f'protected_{risk_group}'] = (
                    bad_variant_vaccinations[f'effective_protected_variant_{risk_group}']
                )
                vaccinations.loc[bad_variant_vaccinations.index, f'immune_{risk_group}'] = (
                    bad_variant_vaccinations[f'effective_variant_{risk_group}']
                )

                vaccinations.loc[not_bad_variant_vaccinations.index, f'unprotected_{risk_group}'] = (
                    not_bad_variant_vaccinations[f'unprotected_{risk_group}']
                )
                vaccinations.loc[not_bad_variant_vaccinations.index, f'protected_{risk_group}'] = (
                    not_bad_variant_vaccinations[f'effective_protected_wildtype_{risk_group}']
                    + not_bad_variant_vaccinations[f'effective_protected_variant_{risk_group}']
                )
                vaccinations.loc[not_bad_variant_vaccinations.index, f'immune_{risk_group}'] = (
                    not_bad_variant_vaccinations[f'effective_wildtype_{risk_group}']
                    + not_bad_variant_vaccinations[f'effective_variant_{risk_group}']
                )
        else:
            vaccinations = None

        if scenario_spec.algorithm == 'draw_level_mandate_reimposition':
            mobility_scenario = scenario_spec.covariates['mobility']
            percent_mandates = self.load_mobility_info(f'{mobility_scenario}_mandate_lift')
            mandate_effects = self.load_mobility_info(f'effect')
        else:
            percent_mandates = None
            mandate_effects = None

        scenario_data = ScenarioData(
            vaccinations=vaccinations,
            percent_mandates=percent_mandates,
            mandate_effects=mandate_effects
        )
        return scenario_data

    def load_mobility_info(self, info_type: str):
        location_ids = self.load_location_ids()
        info_df = io.load(self.covariate_root.mobility_info(info_type=info_type))
        return self._format_covariate_data(info_df, location_ids)

    def load_vaccine_info(self, info_type: str):
        location_ids = self.load_location_ids()
        info_df = io.load(self.covariate_root.vaccine_info(info_type=info_type))
        return self._format_covariate_data(info_df, location_ids)

    ##############################
    # Miscellaneous data loaders #
    ##############################

    def load_thetas(self, theta_specification: Union[str, int], sigma: float) -> pd.Series:
        location_ids = self.load_location_ids()
        if isinstance(theta_specification, str):
            thetas = pd.read_csv(theta_specification).set_index('location_id')['theta']
            thetas = thetas.reindex(location_ids, fill_value=0)
        else:
            thetas = pd.Series(theta_specification,
                               index=pd.Index(location_ids, name='location_id'),
                               name='theta')

        if ((1 < thetas) | thetas < -1).any():
            raise ValueError('Theta must be between -1 and 1.')
        if (sigma - thetas >= 1).any():
            raise ValueError('Sigma - theta must be smaller than 1')

        return thetas

    def get_infections_metadata(self):
        return self._get_regression_data_interface().get_infections_metadata()

    def get_model_inputs_metadata(self):
        infection_metadata = self.get_infections_metadata()
        return infection_metadata['model_inputs_metadata']

    def load_full_data(self) -> pd.DataFrame:
        metadata = self.get_model_inputs_metadata()
        model_inputs_version = metadata['output_path']
        if self.fh_subnationals:
            full_data_path = Path(model_inputs_version) / 'full_data_fh_subnationals.csv'
        else:
            full_data_path = Path(model_inputs_version) / 'full_data.csv'
        full_data = pd.read_csv(full_data_path)
        full_data['date'] = pd.to_datetime(full_data['Date'])
        full_data = full_data.drop(columns=['Date'])
        full_data['location_id'] = full_data['location_id'].astype(int)
        return full_data

    def load_population(self) -> pd.DataFrame:
        return self._get_regression_data_interface().load_population()

    def load_five_year_population(self) -> pd.DataFrame:
        return self._get_regression_data_interface().load_five_year_population()

    def load_total_deaths(self):
        """Load cumulative deaths by location."""
        location_ids = self.load_location_ids()
        full_data = self.load_full_data()
        total_deaths = full_data.groupby('location_id')['Deaths'].max().rename('deaths')
        return total_deaths.loc[location_ids]

    #####################
    # Forecast data I/O #
    #####################

    def save_specification(self, specification: ForecastSpecification) -> None:
        io.dump(specification.to_dict(), self.forecast_root.specification())

    def load_specification(self) -> ForecastSpecification:
        spec_dict = io.load(self.forecast_root.specification())
        return ForecastSpecification.from_dict(spec_dict)

    def save_raw_covariates(self, covariates: pd.DataFrame, scenario: str, draw_id: int) -> None:
        io.dump(covariates, self.forecast_root.raw_covariates(scenario=scenario, draw_id=draw_id))

    def load_raw_covariates(self, scenario: str, draw_id: int) -> pd.DataFrame:
        return io.load(self.forecast_root.raw_covariates(scenario=scenario, draw_id=draw_id))

    def save_ode_params(self, ode_params: pd.DataFrame, scenario: str, draw_id: int) -> None:
        io.dump(ode_params, self.forecast_root.ode_params(scenario=scenario, draw_id=draw_id))

    def load_ode_params(self, scenario: str, draw_id: int) -> pd.DataFrame:
        return io.load(self.forecast_root.ode_params(scenario=scenario, draw_id=draw_id))

    def save_components(self, forecasts: pd.DataFrame, scenario: str, draw_id: int):
        io.dump(forecasts, self.forecast_root.component_draws(scenario=scenario, draw_id=draw_id))

    def load_components(self, scenario: str, draw_id: int):
        components = io.load(self.forecast_root.component_draws(scenario=scenario, draw_id=draw_id))
        components['date'] = pd.to_datetime(components['date'])
        return components.set_index(['location_id', 'date'])

    def save_beta_scales(self, scales: pd.DataFrame, scenario: str, draw_id: int):
        io.dump(scales, self.forecast_root.beta_scaling(scenario=scenario, draw_id=draw_id))

    def load_beta_scales(self, scenario: str, draw_id: int):
        return io.load(self.forecast_root.beta_scaling(scenario=scenario, draw_id=draw_id))

    def save_raw_outputs(self, raw_outputs: pd.DataFrame, scenario: str, draw_id: int):
        io.dump(raw_outputs, self.forecast_root.raw_outputs(scenario=scenario, draw_id=draw_id))

    def load_raw_outputs(self, scenario: str, draw_id: int):
        return io.load(self.forecast_root.raw_outputs(scenario=scenario, draw_id=draw_id))

    #########################
    # Non-interface helpers #
    #########################

    @staticmethod
    def _format_covariate_data(dataset: pd.DataFrame, location_ids: List[int], with_observed: bool = False):
        index_columns = ['location_id']
        if with_observed:
            index_columns.append('observed')
        if 'date' in dataset.columns:
            dataset['date'] = pd.to_datetime(dataset['date'])
            index_columns.append('date')
        dataset = dataset.loc[dataset.location_id.isin(location_ids)]
        dataset = dataset.set_index(index_columns)
        return dataset

    def _get_regression_data_interface(self) -> RegressionDataInterface:
        regression_spec = RegressionSpecification.from_dict(io.load(self.regression_root.specification()))
        regression_di = RegressionDataInterface.from_specification(regression_spec)
        return regression_di