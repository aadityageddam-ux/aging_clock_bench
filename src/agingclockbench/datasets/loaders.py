"""Dataset loaders for bundled reference data."""

from pathlib import Path

import pandas as pd

_DATA_DIR = Path(__file__).parent


def load_nhanes_sample() -> pd.DataFrame:
    """Load the bundled NHANES 2015-2018 sample (preprocessed, ~5000 rows).

    Returns a DataFrame with all columns required by PhenoAge, KDM, and
    DunedinPACEProxy, plus mortality linkage columns from NCHS.

    Columns
    -------
    age, albumin_g_dl, creatinine_mg_dl, glucose_mg_dl, crp_mg_l,
    lymphocyte_pct, mcv_fl, rdw_pct, alp_u_l, wbc_k_ul,
    hemoglobin_g_dl, sex, mortstat, permth_exm

    Returns
    -------
    pd.DataFrame

    Raises
    ------
    FileNotFoundError : if the bundled parquet file is not present.
        Run ``agingclockbench datasets download`` to fetch it.
    """
    parquet_path = _DATA_DIR / "nhanes_sample.parquet"
    if not parquet_path.exists():
        raise FileNotFoundError(
            f"Bundled NHANES sample not found at {parquet_path}.\n"
            "Run: agingclockbench datasets download"
        )
    return pd.read_parquet(parquet_path)
