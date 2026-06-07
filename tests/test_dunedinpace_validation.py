"""Week 2 validation tests for DunedinPACEProxy.

DunedinPACEProxy is a blood-biomarker approximation, NOT the real DunedinPACE.
These tests verify the proxy's statistical properties and biological plausibility.
"""

import numpy as np
import pandas as pd
import pytest

from agingclockbench.clocks.dunedinpace import DunedinPACEProxy
from agingclockbench.clocks.base import ClockResult


@pytest.fixture(scope="module")
def nhanes_dunedin_df():
    parquet = "C:/Users/aadit/phenoage-engine/nhanes_cache/nhanes_phenoage_processed.parquet"
    try:
        ref = pd.read_parquet(parquet)
    except FileNotFoundError:
        pytest.skip("NHANES reference parquet not available")
    return pd.DataFrame({
        "age":             ref["age"],
        "albumin_g_dl":    ref["albumin"],
        "creatinine_mg_dl":ref["creatinine"],
        "glucose_mg_dl":   ref["glucose"],
        "rdw_pct":         ref["rdw"],
        "wbc_k_ul":        ref["wbc"],
        "lymphocyte_pct":  ref["lymphocyte_pct"],
        "mcv_fl":          ref["mcv"],
        "_phenoage_accel": ref["delta_age"].values,
    })


def test_dunedin_returns_clock_result(dunedin_df):
    result = DunedinPACEProxy().transform(dunedin_df)
    assert isinstance(result, ClockResult)
    assert result.clock_name == "DunedinPACEProxy"


def test_dunedin_pace_near_one(nhanes_dunedin_df):
    """Pace scores should be centered near 1.0 on NHANES."""
    df = nhanes_dunedin_df.drop(columns=["_phenoage_accel"])
    result = DunedinPACEProxy().transform(df)
    # Extract pace from biological_ages / chronological_ages
    pace = result.biological_ages / nhanes_dunedin_df["age"].reset_index(drop=True)
    assert abs(pace.mean() - 1.0) < 0.05, f"Mean pace={pace.mean():.3f} not near 1.0"


def test_dunedin_pace_std_near_01(nhanes_dunedin_df):
    """Pace SD should be close to 0.1 (like real DunedinPACE scale)."""
    df = nhanes_dunedin_df.drop(columns=["_phenoage_accel"])
    result = DunedinPACEProxy().transform(df)
    pace = result.biological_ages / nhanes_dunedin_df["age"].reset_index(drop=True)
    assert 0.05 < pace.std() < 0.20, f"Pace SD={pace.std():.3f} not in [0.05, 0.20]"


def test_dunedin_correlated_with_phenoage_accel(nhanes_dunedin_df):
    """Proxy should correlate meaningfully with PhenoAge acceleration (r >= 0.60)."""
    df = nhanes_dunedin_df.drop(columns=["_phenoage_accel"])
    result = DunedinPACEProxy().transform(df)
    ref_accel = pd.Series(nhanes_dunedin_df["_phenoage_accel"].values)
    r = result.accel.corr(ref_accel)
    assert r >= 0.60, f"Proxy~PhenoAge accel correlation r={r:.4f} below 0.60"


def test_dunedin_low_age_correlation(nhanes_dunedin_df):
    """Pace should have low correlation with age — it captures rate not level."""
    df = nhanes_dunedin_df.drop(columns=["_phenoage_accel"])
    result = DunedinPACEProxy().transform(df)
    pace = result.biological_ages / nhanes_dunedin_df["age"].reset_index(drop=True)
    r = abs(pace.corr(nhanes_dunedin_df["age"].reset_index(drop=True)))
    assert r < 0.20, f"Proxy has too-high age correlation r={r:.4f} (expected < 0.20)"


def test_dunedin_faster_aging_higher_rdw():
    """Higher RDW (aging marker) should increase pace score."""
    base = dict(age=55, albumin_g_dl=4.0, creatinine_mg_dl=1.0, glucose_mg_dl=95,
                wbc_k_ul=6.5, lymphocyte_pct=25, mcv_fl=90)
    rows = [dict(rdw_pct=r, **base) for r in [11.5, 13.0, 14.5, 16.0]]
    result = DunedinPACEProxy().transform(pd.DataFrame(rows))
    ages = result.biological_ages.values
    assert all(ages[i] < ages[i + 1] for i in range(len(ages) - 1)), \
        "Higher RDW should increase pace score"


def test_dunedin_missing_column_raises(dunedin_df):
    bad_df = dunedin_df.drop(columns=["rdw_pct"])
    with pytest.raises(ValueError, match="Missing required column"):
        DunedinPACEProxy().transform(bad_df)


def test_dunedin_original_index_populated(dunedin_df):
    result = DunedinPACEProxy().transform(dunedin_df)
    assert result.original_index is not None
    assert len(result.original_index) == result.output_rows
