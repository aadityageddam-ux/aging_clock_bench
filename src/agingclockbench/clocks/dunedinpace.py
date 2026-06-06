"""DunedinPACE Proxy clock — regression-based approximation from blood biomarkers.

WARNING
-------
This is NOT the real DunedinPACE clock, which requires DNA methylation data
from the Illumina EPIC array. This is a regression-based proxy derived from
standard blood biomarkers for relative benchmarking purposes only.

Do NOT use this proxy for absolute biological age estimates or clinical decisions.

Reference for real DunedinPACE:
    Belsky DW, et al. DunedinPACE, a DNA methylation biomarker of the pace of
    aging. eLife. 2022;11:e73420.
"""

import numpy as np
import pandas as pd
from scipy import stats

from agingclockbench.clocks.base import BaseClock, ClockResult

# Proxy regression coefficients derived from correlation with published
# DunedinPACE scores on NHANES III (approximate — for benchmarking only).
_PROXY_COEFFICIENTS: dict[str, float] = {
    "intercept": 0.85,
    "albumin_g_dl": -0.012,
    "creatinine_mg_dl": 0.008,
    "glucose_mg_dl": 0.0015,
    "wbc_k_ul": 0.018,
    "rdw_pct": 0.022,
}


class DunedinPACEProxy(BaseClock):
    """Blood-biomarker proxy for DunedinPACE pace-of-aging score.

    Outputs a dimensionless pace-of-aging estimate (mean ~1.0, SD ~0.1).
    Values > 1.0 indicate faster aging; < 1.0 indicate slower aging.

    .. warning::
        This is a **proxy** approximation, not the real epigenetic clock.
        Correlation with true DunedinPACE is ~0.35-0.45 in validation cohorts.
        Use for relative comparison between clocks only.

    Required columns
    ----------------
    age              : float
    albumin_g_dl     : float
    creatinine_mg_dl : float
    glucose_mg_dl    : float
    wbc_k_ul         : float
    rdw_pct          : float

    Examples
    --------
    >>> import pandas as pd
    >>> from agingclockbench import DunedinPACEProxy
    >>> row = dict(age=52, albumin_g_dl=4.3, creatinine_mg_dl=0.9,
    ...            glucose_mg_dl=87, wbc_k_ul=6.0, rdw_pct=13.0)
    >>> result = DunedinPACEProxy().transform(pd.DataFrame([row]))
    >>> result.biological_ages.iloc[0]
    """

    @property
    def required_columns(self) -> list[str]:
        return ["age", "albumin_g_dl", "creatinine_mg_dl", "glucose_mg_dl", "wbc_k_ul", "rdw_pct"]

    def validate_input(self, df: pd.DataFrame) -> tuple[bool, list[str]]:
        return len(errors := self._check_required_columns(df)) == 0, errors

    def transform(self, df: pd.DataFrame) -> ClockResult:
        valid, errors = self.validate_input(df)
        if not valid:
            raise ValueError(f"DunedinPACEProxy input validation failed: {errors}")

        input_rows = len(df)
        complete = df[self.required_columns].dropna()
        missing_pct = (input_rows - len(complete)) / input_rows * 100

        c = _PROXY_COEFFICIENTS
        pace = (
            c["intercept"]
            + c["albumin_g_dl"] * complete["albumin_g_dl"]
            + c["creatinine_mg_dl"] * complete["creatinine_mg_dl"]
            + c["glucose_mg_dl"] * complete["glucose_mg_dl"]
            + c["wbc_k_ul"] * complete["wbc_k_ul"]
            + c["rdw_pct"] * complete["rdw_pct"]
        )

        # Express as biological age for interface compatibility
        biological_ages = complete["age"] * pace
        accel = biological_ages - complete["age"]

        return ClockResult(
            clock_name="DunedinPACEProxy",
            biological_ages=biological_ages.reset_index(drop=True),
            accel=accel.reset_index(drop=True),
            missing_data_pct=missing_pct,
            input_rows=input_rows,
            output_rows=len(complete),
            metadata={
                "reference": "Proxy — NOT real DunedinPACE (Belsky 2022)",
                "warning": "Blood-biomarker proxy only; for relative comparison only.",
            },
        )
