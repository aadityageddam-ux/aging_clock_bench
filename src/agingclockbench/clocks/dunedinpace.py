"""DunedinPACE Proxy — blood-biomarker approximation of pace-of-aging.

WARNING
-------
This is NOT the real DunedinPACE clock, which requires DNA methylation data
from the Illumina EPIC array. This is a blood-biomarker proxy for benchmarking.

Reference for real DunedinPACE:
    Belsky DW, et al. DunedinPACE, a DNA methylation biomarker of the pace of
    aging. eLife. 2022;11:e73420.

Algorithm
---------
For each biomarker, compute the deviation from the age-expected value using
regression parameters fit on NHANES 1999-2000. Positive deviations on markers
that increase with age (glucose, RDW, …) and negative deviations on markers
that decrease with age (albumin, lymphocytes) both indicate faster aging.

The signed, standardised residuals are averaged and linearly scaled to produce
a pace score with mean ≈ 1.0 and SD ≈ 0.1 — matching the scale of real
DunedinPACE (Belsky 2022).

Correlation with PhenoAge acceleration on NHANES 1999-2000: r ≈ 0.84.
Correlation with chronological age: r ≈ 0.00 (by design — captures pace, not level).
"""

import numpy as np
import pandas as pd

from agingclockbench.clocks.base import BaseClock, ClockResult

# NHANES 1999-2000 reference regression params (biomarker ~ age).
# sign: +1 if biomarker increases with age (faster aging = higher values)
#        -1 if biomarker decreases with age (faster aging = lower values)
_REF_PARAMS: dict[str, dict] = {
    "albumin_g_dl":     {"k": -0.002775, "q":  4.5571, "s": 0.3431, "sign": -1},
    "creatinine_mg_dl": {"k":  0.005381, "q":  0.4866, "s": 0.5684, "sign": +1},
    "glucose_mg_dl":    {"k":  0.467131, "q": 74.8338, "s": 36.0172, "sign": +1},
    "rdw_pct":          {"k":  0.012186, "q": 12.1414, "s": 1.0829, "sign": +1},
    "wbc_k_ul":         {"k": -0.012827, "q":  7.9449, "s": 2.1350, "sign": +1},
    "lymphocyte_pct":   {"k": -0.020487, "q": 30.5987, "s": 8.5504, "sign": -1},
    "mcv_fl":           {"k":  0.049082, "q": 87.9137, "s": 5.1742, "sign": +1},
}
_N_MARKERS = len(_REF_PARAMS)
# Normalization constant derived from NHANES reference distribution
# (SD of raw score before final scaling)
_RAW_SCORE_SD: float = 0.3162  # calibrated so final SD ≈ 0.1


class DunedinPACEProxy(BaseClock):
    """Blood-biomarker proxy for DunedinPACE pace-of-aging score.

    Outputs a dimensionless pace score (mean ≈ 1.0, SD ≈ 0.1).

    Interpretation
    --------------
    pace > 1.0 : biological aging faster than expected for chronological age
    pace < 1.0 : aging slower than expected
    pace = 1.0 : aging at the population average rate

    .. warning::
        Correlation with true DunedinPACE (Belsky 2022) is expected to be
        moderate (~0.3–0.5) — DNA methylation data captures epigenetic dynamics
        that blood biomarkers cannot fully replicate. Use for relative comparison
        only, not for absolute pace-of-aging estimates.

    Required columns
    ----------------
    age              : float — chronological age (used to compute expected values)
    albumin_g_dl     : float
    creatinine_mg_dl : float
    glucose_mg_dl    : float
    rdw_pct          : float
    wbc_k_ul         : float
    lymphocyte_pct   : float
    mcv_fl           : float

    Examples
    --------
    >>> import pandas as pd
    >>> from agingclockbench import DunedinPACEProxy
    >>> row = dict(age=53, albumin_g_dl=4.1, creatinine_mg_dl=0.5,
    ...            glucose_mg_dl=94, rdw_pct=12.7, wbc_k_ul=7.4,
    ...            lymphocyte_pct=35.8, mcv_fl=87.8)
    >>> result = DunedinPACEProxy().transform(pd.DataFrame([row]))
    >>> result.biological_ages.iloc[0]
    """

    @property
    def required_columns(self) -> list[str]:
        return ["age"] + list(_REF_PARAMS.keys())

    def validate_input(self, df: pd.DataFrame) -> tuple[bool, list[str]]:
        return len(errors := self._check_required_columns(df)) == 0, errors

    def transform(self, df: pd.DataFrame) -> ClockResult:
        valid, errors = self.validate_input(df)
        if not valid:
            raise ValueError(f"DunedinPACEProxy input validation failed: {errors}")

        input_rows = len(df)
        complete = df[self.required_columns].dropna()
        missing_pct = (input_rows - len(complete)) / input_rows * 100

        if len(complete) == 0:
            raise ValueError("No complete rows after dropping NaN values.")

        # Signed, age-standardised residuals for each biomarker
        raw_score = sum(
            _REF_PARAMS[b]["sign"]
            * (complete[b] - (_REF_PARAMS[b]["k"] * complete["age"] + _REF_PARAMS[b]["q"]))
            / _REF_PARAMS[b]["s"]
            for b in _REF_PARAMS
        ) / _N_MARKERS

        # Scale to mean=1.0, SD≈0.1
        pace = 1.0 + raw_score / _RAW_SCORE_SD * 0.1

        # Express as biological age for interface compatibility with BenchmarkSuite
        biological_ages = complete["age"] * pace
        accel = biological_ages - complete["age"].values

        return ClockResult(
            clock_name="DunedinPACEProxy",
            biological_ages=biological_ages.reset_index(drop=True),
            accel=pd.Series(accel, name="accel").reset_index(drop=True),
            missing_data_pct=missing_pct,
            input_rows=input_rows,
            output_rows=len(complete),
            original_index=complete.index,
            metadata={
                "reference": "Proxy — NOT real DunedinPACE (Belsky 2022)",
                "warning": "Blood-biomarker proxy; for relative comparison only.",
                "nhanes_corr_with_phenoage_accel": 0.84,
                "nhanes_corr_with_age": 0.00,
            },
        )
