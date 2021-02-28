from dataclasses import dataclass
from typing import Dict, List, Union

import pandas as pd

from covid_model_seiir_pipeline.lib import (
    utilities,
)

# This is just exposing these containers from this namespace so we're not
# importing from the regression stage everywhere.
from covid_model_seiir_pipeline.pipeline.regression.model.containers import (
    RatioData,
    HospitalCensusData,
    HospitalMetrics,
    HospitalCorrectionFactors,
)


class Indices:
    """Abstraction for building square datasets."""

    def __init__(self,
                 past_start_dates: pd.Series,
                 forecast_start_dates: pd.Series,
                 forecast_end_dates: pd.Series):
        self._past_index = self._build_index(past_start_dates, forecast_start_dates, pd.Timedelta(days=1))
        self._future_index = self._build_index(forecast_start_dates, forecast_end_dates)
        self._initial_condition_index = (
            forecast_start_dates
            .reset_index()
            .set_index(['location_id', 'date'])
            .sort_index()
            .index
        )
        self._full_index = self._build_index(past_start_dates, forecast_end_dates)

    @property
    def past(self) -> pd.MultiIndex:
        """Location-date index for the past."""
        return self._past_index.copy()

    @property
    def future(self) -> pd.MultiIndex:
        """Location-date index for the future."""
        return self._future_index.copy()

    @property
    def initial_condition(self) -> pd.MultiIndex:
        """Location-date index for the initial condition.
        This index has one date per location.
        """
        return self._initial_condition_index.copy()

    @property
    def full(self) -> pd.MultiIndex:
        """Location-date index for the full time series, past and future."""
        return self._full_index.copy()

    @staticmethod
    def _build_index(start: pd.Series,
                     end: pd.Series,
                     end_offset: pd.Timedelta = pd.Timedelta(days=0)) -> pd.MultiIndex:
        index = (pd.concat([start.rename('start'), end.rename('end')], axis=1)
                 .groupby('location_id')
                 .apply(lambda x: pd.date_range(x.iloc[0, 0], x.iloc[0, 1] - end_offset))
                 .explode()
                 .rename('date')
                 .reset_index()
                 .set_index(['location_id', 'date'])
                 .sort_index()
                 .index)
        return index


@dataclass
class ModelParameters:
    # Core parameters
    alpha: pd.Series
    beta: pd.Series
    sigma: pd.Series
    gamma1: pd.Series
    gamma2: pd.Series

    # Theta parameters
    theta_plus: pd.Series
    theta_minus: pd.Series

    # Vaccine parameters
    unprotected_lr: pd.Series
    protected_wild_type_lr: pd.Series
    protected_all_types_lr: pd.Series
    immune_wild_type_lr: pd.Series
    immune_all_types_lr: pd.Series

    unprotected_hr: pd.Series
    protected_wild_type_hr: pd.Series
    protected_all_types_hr: pd.Series
    immune_wild_type_hr: pd.Series
    immune_all_types_hr: pd.Series

    # Variant parameters
    beta_wild: pd.Series
    beta_variant: pd.Series
    p_wild: pd.Series
    p_variant: pd.Series
    p_all_variant: pd.Series
    probability_cross_immune: pd.Series

    def with_index(self, index: pd.MultiIndex):
        return ModelParameters(**{
            parameter_name: parameter.loc[index] for parameter_name, parameter in self.to_dict().items()
        })

    def to_dict(self) -> Dict[str, pd.Series]:
        return utilities.asdict(self)

    def to_df(self) -> pd.DataFrame:
        return pd.concat([v.rename(k) for k, v in self.to_dict().items()], axis=1)


@dataclass
class OutputMetrics:
    components: pd.DataFrame
    infections: pd.Series
    cases: pd.Series
    admissions: pd.Series
    deaths: pd.DataFrame
    r_controlled: pd.Series
    r_effective: pd.Series
    herd_immunity: pd.Series
    total_susceptible: pd.Series
    total_immune: pd.Series

    def to_dict(self) -> Dict[str, pd.DataFrame]:
        return utilities.asdict(self)


@dataclass
class CompartmentInfo:
    compartments: List[str]
    group_suffixes: List[str]

    def to_dict(self) -> Dict[str, List[str]]:
        return utilities.asdict(self)


@dataclass
class VariantScalars:
    beta: pd.Series
    ifr: pd.Series

    def to_dict(self) -> Dict[str, pd.Series]:
        return utilities.asdict(self)
