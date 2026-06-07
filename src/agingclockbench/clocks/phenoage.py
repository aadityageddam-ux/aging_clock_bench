"""PhenoAge clock — Levine et al. 2018 (Aging Cell).

Reference: Levine ME, et al. An epigenetic biomarker of aging for lifespan and
healthspan. Aging Cell. 2018;17(4):e12759.

Unit notes
----------
Inputs are accepted in standard NHANES clinical units (g/dL, mg/dL, mg/L).
The transform() method converts internally before applying the Levine 2018
coefficients, which were calibrated on:
  albumin  → g/L   (×10)
  creatinine → μmol/L (×88.4)
  glucose  → mmol/L (×0.0555)
  crp      → mg/L  (input already; natural-log applied after +0.001 epsilon)
"""

import numpy as np
import pandas as pd

from agingclockbench.clocks.base import BaseClock, ClockResult

# Coefficients applied to CONVERTED units (Levine 2018).
_COEFFICIENTS: dict[str, float] = {
    "intercept": -19.907,
    "age": 0.0804,
    "albumin_g_l": -0.0336,        # after ×10 from g/dL
    "creatinine_umol_l": 0.0095,   # after ×88.4 from mg/dL
    "glucose_mmol_l": 0.1953,      # after ×0.0555 from mg/dL
    "ln_crp_mg_l": 0.0954,         # ln(mg/L + 0.001)
    "lymphocyte_pct": -0.0120,
    "mcv_fl": 0.0268,
    "rdw_pct": 0.3306,
    "alp_u_l": 0.00188,
    "wbc_k_ul": 0.0554,
}

# Gompertz parameters (10-year mortality, t=120 months)
_GAMMA: float = 0.0076927
_T_MONTHS: int = 120
_PHENOAGE_INTERCEPT: float = 141.50
_PHENOAGE_SLOPE: float = 0.090165
_MORT_CONSTANT: float = -0.00553


class PhenoAge(BaseClock):
    """Biological age calculator implementing the Levine 2018 PhenoAge algorithm.

    All inputs are in standard NHANES clinical units. Unit conversions to the
    Levine 2018 coefficient scale are applied internally.

    Required columns
    ----------------
    age              : float — chronological age in years
    albumin_g_dl     : float — albumin in g/dL  (converted internally to g/L)
    creatinine_mg_dl : float — creatinine in mg/dL (converted to μmol/L)
    glucose_mg_dl    : float — glucose in mg/dL (converted to mmol/L)
    crp_mg_l         : float — C-reactive protein in mg/L (ln-transformed)
    lymphocyte_pct   : float — lymphocyte percentage (%)
    mcv_fl           : float — mean corpuscular volume in fL
    rdw_pct          : float — red cell distribution width (%)
    alp_u_l          : float — alkaline phosphatase in U/L
    wbc_k_ul         : float — white blood cell count in 10³/μL

    Examples
    --------
    >>> import pandas as pd
    >>> from agingclockbench import PhenoAge
    >>> row = dict(age=52, albumin_g_dl=4.3, creatinine_mg_dl=0.9,
    ...            glucose_mg_dl=87, crp_mg_l=0.3, lymphocyte_pct=28,
    ...            mcv_fl=90, rdw_pct=13.0, alp_u_l=65, wbc_k_ul=6.0)
    >>> result = PhenoAge().transform(pd.DataFrame([row]))
    >>> round(result.biological_ages.iloc[0], 1)
    44.9
    """

    @property
    def required_columns(self) -> list[str]:
        return [
            "age",
            "albumin_g_dl",
            "creatinine_mg_dl",
            "glucose_mg_dl",
            "crp_mg_l",
            "lymphocyte_pct",
            "mcv_fl",
            "rdw_pct",
            "alp_u_l",
            "wbc_k_ul",
        ]

    def validate_input(self, df: pd.DataFrame) -> tuple[bool, list[str]]:
        errors = self._check_required_columns(df)
        if not errors and (df["crp_mg_l"].dropna() < 0).any():
            errors.append("crp_mg_l contains negative values — check units.")
        return len(errors) == 0, errors

    def transform(self, df: pd.DataFrame) -> ClockResult:
        valid, errors = self.validate_input(df)
        if not valid:
            raise ValueError(f"PhenoAge input validation failed: {errors}")

        input_rows = len(df)
        complete = df[self.required_columns].dropna()
        dropped = input_rows - len(complete)
        missing_pct = dropped / input_rows * 100

        if len(complete) == 0:
            raise ValueError("No complete rows after dropping NaN values.")

        # --- Unit conversions (applied before coefficients) ---
        albumin_g_l = complete["albumin_g_dl"] * 10.0
        creatinine_umol_l = complete["creatinine_mg_dl"] * 88.4
        glucose_mmol_l = complete["glucose_mg_dl"] * 0.0555
        # CRP is already in mg/L; +0.001 epsilon prevents ln(0)
        ln_crp = np.log(complete["crp_mg_l"].clip(lower=0.001))

        # --- Linear predictor (xb) ---
        c = _COEFFICIENTS
        xb = (
            c["intercept"]
            + c["age"] * complete["age"]
            + c["albumin_g_l"] * albumin_g_l
            + c["creatinine_umol_l"] * creatinine_umol_l
            + c["glucose_mmol_l"] * glucose_mmol_l
            + c["ln_crp_mg_l"] * ln_crp
            + c["lymphocyte_pct"] * complete["lymphocyte_pct"]
            + c["mcv_fl"] * complete["mcv_fl"]
            + c["rdw_pct"] * complete["rdw_pct"]
            + c["alp_u_l"] * complete["alp_u_l"]
            + c["wbc_k_ul"] * complete["wbc_k_ul"]
        )

        # --- 10-year mortality score (Gompertz, t=120 months) ---
        mortality_score = 1 - np.exp(
            -np.exp(xb) * (np.exp(_GAMMA * _T_MONTHS) - 1) / _GAMMA
        )
        mortality_score = mortality_score.clip(upper=0.9999)

        # --- Phenotypic age (Levine 2018 Eq. 2) ---
        biological_ages = (
            _PHENOAGE_INTERCEPT
            + np.log(_MORT_CONSTANT * np.log(1 - mortality_score)) / _PHENOAGE_SLOPE
        )

        accel = biological_ages - complete["age"].values

        return ClockResult(
            clock_name="PhenoAge",
            biological_ages=biological_ages.reset_index(drop=True),
            accel=pd.Series(accel, name="accel"),
            missing_data_pct=missing_pct,
            input_rows=input_rows,
            output_rows=len(complete),
            original_index=complete.index,
            metadata={"reference": "Levine 2018 Aging Cell", "coefficients": _COEFFICIENTS},
        )
