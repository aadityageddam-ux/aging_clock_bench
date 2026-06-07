"""Klemera-Doubal Method (KDM) biological age clock.

Reference: Klemera P, Doubal S. A new approach to the concept and computation
of biological age. Mech Ageing Dev. 2006;127(3):240-248.

Algorithm
---------
1. For each of m biomarkers, fit a linear regression: x_j = q_j + k_j * age.
2. Compute preliminary biological age (BA1) as the weighted maximum-likelihood
   estimate of age given observed biomarker values.
3. Estimate s_BA: the standard deviation of (BA1 - chronological_age) in the
   reference cohort.
4. Compute final KDM biological age incorporating chronological age as an
   additional "measurement" anchored with precision 1/s_BA^2.

Default reference parameters are derived from NHANES 1999-2000 (N=4,086
complete cases). Provide your own cohort via ``fit()`` before ``transform()``.
"""

import numpy as np
import pandas as pd
from scipy import stats

from agingclockbench.clocks.base import BaseClock, ClockResult

# NHANES 1999-2000 reference regression parameters (slope k, intercept q, residual SD s).
# Biomarkers are in their standard NHANES clinical units (g/dL, mg/dL, %, fL, U/L).
_NHANES_PARAMS: dict[str, dict] = {
    "albumin_g_dl":     {"k": -0.002775, "q":  4.5571, "s": 0.3431},
    "creatinine_mg_dl": {"k":  0.005381, "q":  0.4866, "s": 0.5684},
    "glucose_mg_dl":    {"k":  0.467131, "q": 74.8338, "s": 36.0172},
    "rdw_pct":          {"k":  0.012186, "q": 12.1414, "s": 1.0829},
    "mcv_fl":           {"k":  0.049082, "q": 87.9137, "s": 5.1742},
    "wbc_k_ul":         {"k": -0.012827, "q":  7.9449, "s": 2.1350},
    "alp_u_l":          {"k":  0.190120, "q": 74.7400, "s": 32.1450},
    "lymphocyte_pct":   {"k": -0.020487, "q": 30.5987, "s": 8.5504},
}
_NHANES_S_BA: float = 40.4491  # SD of (BA1 - chronological_age) in NHANES reference


class KDM(BaseClock):
    """Klemera-Doubal Method biological age estimator.

    KDM is a maximum-likelihood estimator of biological age from a set of
    biomarkers, each linearly regressed on chronological age in a reference
    population. Chronological age itself is incorporated as a final "measurement"
    with precision 1/s_BA^2, where s_BA is the variability of the preliminary
    estimate in the reference cohort.

    Default reference parameters are from NHANES 1999-2000 (N=4,086).
    For your own cohort, call ``fit(df)`` before ``transform(df)``.

    Required columns (NHANES clinical units)
    -----------------------------------------
    age              : float — chronological age in years
    albumin_g_dl     : float — g/dL
    creatinine_mg_dl : float — mg/dL
    glucose_mg_dl    : float — mg/dL
    rdw_pct          : float — %
    mcv_fl           : float — fL
    wbc_k_ul         : float — 10³/μL
    alp_u_l          : float — U/L
    lymphocyte_pct   : float — %

    Examples
    --------
    >>> import pandas as pd
    >>> from agingclockbench import KDM
    >>> row = dict(age=53, albumin_g_dl=4.1, creatinine_mg_dl=0.5, glucose_mg_dl=94,
    ...            rdw_pct=12.7, mcv_fl=87.8, wbc_k_ul=7.4, alp_u_l=98,
    ...            lymphocyte_pct=35.8)
    >>> result = KDM().transform(pd.DataFrame([row]))
    >>> result.biological_ages.iloc[0]
    """

    _BIOMARKERS = list(_NHANES_PARAMS.keys())

    def __init__(self) -> None:
        self._params: dict = _NHANES_PARAMS.copy()
        self._s_ba: float = _NHANES_S_BA

    @property
    def required_columns(self) -> list[str]:
        return ["age"] + self._BIOMARKERS

    def validate_input(self, df: pd.DataFrame) -> tuple[bool, list[str]]:
        return len(errors := self._check_required_columns(df)) == 0, errors

    def fit(self, df: pd.DataFrame) -> "KDM":
        """Fit reference regression parameters from a training cohort.

        Derives slopes, intercepts, and residual SDs by regressing each
        biomarker on chronological age, then estimates s_BA.

        Parameters
        ----------
        df : DataFrame with all required columns.

        Returns
        -------
        self — for method chaining.
        """
        complete = df[self.required_columns].dropna()
        if len(complete) < 30:
            raise ValueError(f"Need at least 30 complete rows to fit KDM; got {len(complete)}.")

        ages = complete["age"].values
        params = {}
        for col in self._BIOMARKERS:
            slope, intercept, _, _, _ = stats.linregress(ages, complete[col].values)
            resid = complete[col].values - (slope * ages + intercept)
            s = max(resid.std(), 1e-6)
            params[col] = {"k": slope, "q": intercept, "s": s}
        self._params = params

        # Compute preliminary BA1 and estimate s_BA
        ba1 = self._preliminary_ba(complete)
        self._s_ba = max(float((ba1 - complete["age"]).std()), 1e-6)
        return self

    def transform(self, df: pd.DataFrame) -> ClockResult:
        valid, errors = self.validate_input(df)
        if not valid:
            raise ValueError(f"KDM input validation failed: {errors}")

        input_rows = len(df)
        complete = df[self.required_columns].dropna()
        missing_pct = (input_rows - len(complete)) / input_rows * 100

        if len(complete) == 0:
            raise ValueError("No complete rows after dropping NaN values.")

        ba1 = self._preliminary_ba(complete)

        # Full KDM: incorporate chronological age as an additional measurement
        # BA = [Σ(k_j*(x_j - q_j)/s_j²) + CA/s_BA²] / [Σ(k_j²/s_j²) + 1/s_BA²]
        numerator = (
            sum(
                self._params[b]["k"] * (complete[b] - self._params[b]["q"]) / self._params[b]["s"] ** 2
                for b in self._BIOMARKERS
            )
            + complete["age"] / self._s_ba ** 2
        )
        denominator = (
            sum(self._params[b]["k"] ** 2 / self._params[b]["s"] ** 2 for b in self._BIOMARKERS)
            + 1.0 / self._s_ba ** 2
        )
        biological_ages = (numerator / denominator).reset_index(drop=True)
        accel = (biological_ages - complete["age"].reset_index(drop=True)).rename("accel")

        return ClockResult(
            clock_name="KDM",
            biological_ages=biological_ages,
            accel=accel,
            missing_data_pct=missing_pct,
            input_rows=input_rows,
            output_rows=len(complete),
            original_index=complete.index,
            metadata={
                "reference": "Klemera & Doubal 2006; params from NHANES 1999-2000",
                "s_ba": self._s_ba,
                "n_biomarkers": len(self._BIOMARKERS),
            },
        )

    def _preliminary_ba(self, complete: pd.DataFrame) -> pd.Series:
        """Weighted ML estimate of age without the chronological age anchor."""
        numerator = sum(
            self._params[b]["k"] * (complete[b] - self._params[b]["q"]) / self._params[b]["s"] ** 2
            for b in self._BIOMARKERS
        )
        denominator = sum(
            self._params[b]["k"] ** 2 / self._params[b]["s"] ** 2 for b in self._BIOMARKERS
        )
        return numerator / denominator
