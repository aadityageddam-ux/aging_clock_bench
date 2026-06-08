"""Dataset loaders for bundled reference data."""

from pathlib import Path

import pandas as pd

_DATA_DIR = Path(__file__).parent


def load_nhanes_sample() -> pd.DataFrame:
    """Load the bundled NHANES 1999-2000 sample (N=4,086 complete cases).

    Source: CDC National Health and Nutrition Examination Survey 1999-2000,
    with mortality linkage from NCHS Public-Use Linked Mortality Files.
    All biomarkers are complete (no missing values). Participants are 20-85 years
    old.

    Columns
    -------
    seqn             : int   — NHANES participant sequence number
    age              : float — chronological age in years
    sex              : str   — 'male' or 'female'
    albumin_g_dl     : float — albumin (g/dL)
    creatinine_mg_dl : float — creatinine (mg/dL)
    glucose_mg_dl    : float — glucose (mg/dL)
    crp_mg_l         : float — C-reactive protein (mg/L)
    lymphocyte_pct   : float — lymphocyte percentage (%)
    mcv_fl           : float — mean corpuscular volume (fL)
    rdw_pct          : float — red cell distribution width (%)
    alp_u_l          : float — alkaline phosphatase (U/L)
    wbc_k_ul         : float — white blood cell count (10³/μL)
    mortstat         : int   — vital status at follow-up (1=deceased, 0=alive/censored)
    permth_exm       : float — months from examination to death or censoring

    Returns
    -------
    pd.DataFrame with 4,086 rows and 14 columns.

    Examples
    --------
    >>> df = load_nhanes_sample()
    >>> len(df)
    4086
    >>> list(df.columns[:4])
    ['seqn', 'age', 'sex', 'albumin_g_dl']
    """
    parquet_path = _DATA_DIR / "nhanes_sample.parquet"
    if not parquet_path.exists():
        raise FileNotFoundError(
            f"Bundled NHANES sample not found at {parquet_path}.\n"
            "This file should be included in the installed package. "
            "If you installed from source, run: python -m agingclockbench.datasets.build"
        )
    return pd.read_parquet(parquet_path)
