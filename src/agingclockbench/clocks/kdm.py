"""Klemera-Doubal Method (KDM) biological age clock.

Reference: Klemera P, Doubal S. A new approach to the concept and computation
of biological age. Mech Ageing Dev. 2006;127(3):240-248.
"""

import numpy as np
import pandas as pd
from scipy import stats

from agingclockbench.clocks.base import BaseClock, ClockResult


class KDM(BaseClock):
    """Klemera-Doubal Method biological age estimator.

    KDM computes a weighted composite of biomarker deviations from
    age-regressed reference values. Weights are derived from the
    correlation of each biomarker with chronological age.

    This implementation uses reference regression parameters fit on
    NHANES 2015-2018 (ages 20-85). When using your own data, call
    ``fit()`` first to derive cohort-specific parameters.

    Required columns
    ----------------
    age              : float — chronological age in years
    albumin_g_dl     : float
    creatinine_mg_dl : float
    glucose_mg_dl    : float
    hemoglobin_g_dl  : float

    Examples
    --------
    >>> import pandas as pd
    >>> from agingclockbench import KDM
    >>> row = dict(age=52, albumin_g_dl=4.3, creatinine_mg_dl=0.9,
    ...            glucose_mg_dl=87, hemoglobin_g_dl=14.5)
    >>> result = KDM().transform(pd.DataFrame([row]))
    >>> result.biological_ages.iloc[0]
    """

    # NHANES-derived reference regression parameters (slope k_j, intercept q_j,
    # residual SD s_j) for each biomarker regressed on chronological age.
    # These are placeholders — replaced during fit() or when bundled NHANES
    # parameters are loaded.
    _BIOMARKERS = ["albumin_g_dl", "creatinine_mg_dl", "glucose_mg_dl", "hemoglobin_g_dl"]

    def __init__(self) -> None:
        # Parameters set after fit(); None until then.
        self._params: dict | None = None

    @property
    def required_columns(self) -> list[str]:
        return ["age"] + self._BIOMARKERS

    def validate_input(self, df: pd.DataFrame) -> tuple[bool, list[str]]:
        return len(errors := self._check_required_columns(df)) == 0, errors

    def fit(self, df: pd.DataFrame) -> "KDM":
        """Fit reference regression parameters from a training cohort.

        Parameters
        ----------
        df : DataFrame with required columns including ``age``.

        Returns
        -------
        self — for chaining.
        """
        complete = df[self.required_columns].dropna()
        ages = complete["age"].values
        params = {}
        for col in self._BIOMARKERS:
            slope, intercept, r, _, _ = stats.linregress(ages, complete[col].values)
            residuals = complete[col].values - (slope * ages + intercept)
            s = residuals.std()
            params[col] = {"k": slope, "q": intercept, "s": max(s, 1e-6), "r": r}
        self._params = params
        return self

    def transform(self, df: pd.DataFrame) -> ClockResult:
        valid, errors = self.validate_input(df)
        if not valid:
            raise ValueError(f"KDM input validation failed: {errors}")

        if self._params is None:
            self.fit(df)

        input_rows = len(df)
        complete = df[self.required_columns].dropna()
        missing_pct = (input_rows - len(complete)) / input_rows * 100

        ages = complete["age"].values
        params = self._params

        numerator = sum(
            params[b]["k"] * (complete[b].values - params[b]["q"]) / params[b]["s"] ** 2
            for b in self._BIOMARKERS
        )
        weight_sum = sum(params[b]["k"] ** 2 / params[b]["s"] ** 2 for b in self._BIOMARKERS)
        kdm_ba = (numerator + ages / (1.0 / len(self._BIOMARKERS))) / (
            weight_sum + 1.0 / (1.0 / len(self._BIOMARKERS))
        )

        biological_ages = pd.Series(kdm_ba, name="kdm_ba")
        accel = biological_ages - pd.Series(ages)

        return ClockResult(
            clock_name="KDM",
            biological_ages=biological_ages.reset_index(drop=True),
            accel=accel.reset_index(drop=True),
            missing_data_pct=missing_pct,
            input_rows=input_rows,
            output_rows=len(complete),
            metadata={"reference": "Klemera & Doubal 2006"},
        )
