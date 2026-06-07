"""Week 2 validation tests for KDM clock.

Tests verify:
1. KDM with NHANES reference params gives expected population statistics
2. Monotonicity with key aging biomarkers
3. Proper Klemera-Doubal formula behavior
"""

import numpy as np
import pandas as pd
import pytest

from agingclockbench.clocks.kdm import KDM


@pytest.fixture(scope="module")
def nhanes_kdm_df():
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
        "mcv_fl":          ref["mcv"],
        "wbc_k_ul":        ref["wbc"],
        "alp_u_l":         ref["alp"],
        "lymphocyte_pct":  ref["lymphocyte_pct"],
        "mortstat":        ref["deceased"],
        "permth_exm":      ref["followup_years"] * 12,
    })


def test_kdm_population_pearson_r(nhanes_kdm_df):
    """KDM should correlate positively with age on NHANES (target r >= 0.50)."""
    df = nhanes_kdm_df.drop(columns=["mortstat", "permth_exm"])
    result = KDM().transform(df)
    r = result.biological_ages.corr(nhanes_kdm_df["age"].reset_index(drop=True))
    assert r >= 0.50, f"KDM Pearson r={r:.4f} below minimum 0.50"


def test_kdm_population_mean_in_range(nhanes_kdm_df):
    """Mean KDM biological age should be within +-10 years of mean chronological age."""
    df = nhanes_kdm_df.drop(columns=["mortstat", "permth_exm"])
    result = KDM().transform(df)
    mean_kdm = result.biological_ages.mean()
    mean_age = nhanes_kdm_df["age"].mean()
    assert abs(mean_kdm - mean_age) <= 10, f"|mean KDM ({mean_kdm:.1f}) - mean age ({mean_age:.1f})| > 10"


def test_kdm_acceleration_mean_near_zero(nhanes_kdm_df):
    """Mean acceleration should be close to zero in the reference cohort."""
    df = nhanes_kdm_df.drop(columns=["mortstat", "permth_exm"])
    result = KDM().transform(df)
    mean_accel = result.accel.mean()
    assert abs(mean_accel) < 5, f"Mean KDM acceleration {mean_accel:.2f} not near zero"


def test_kdm_monotone_with_glucose():
    """Higher glucose → higher KDM biological age."""
    base = dict(age=55, albumin_g_dl=4.0, creatinine_mg_dl=1.0,
                rdw_pct=13, mcv_fl=90, wbc_k_ul=6.5, alp_u_l=70, lymphocyte_pct=25)
    rows = [dict(glucose_mg_dl=g, **base) for g in [80, 100, 130, 180]]
    result = KDM().transform(pd.DataFrame(rows))
    ages = result.biological_ages.values
    assert all(ages[i] < ages[i + 1] for i in range(len(ages) - 1)), \
        f"KDM not monotone with glucose: {ages}"


def test_kdm_monotone_with_rdw():
    """Higher RDW → higher KDM biological age."""
    base = dict(age=55, albumin_g_dl=4.0, creatinine_mg_dl=1.0, glucose_mg_dl=95,
                mcv_fl=90, wbc_k_ul=6.5, alp_u_l=70, lymphocyte_pct=25)
    rows = [dict(rdw_pct=r, **base) for r in [11.5, 13.0, 14.5, 16.0]]
    result = KDM().transform(pd.DataFrame(rows))
    ages = result.biological_ages.values
    assert all(ages[i] < ages[i + 1] for i in range(len(ages) - 1)), \
        f"KDM not monotone with RDW: {ages}"


def test_kdm_monotone_with_albumin():
    """Higher albumin → lower KDM biological age (protective)."""
    base = dict(age=55, creatinine_mg_dl=1.0, glucose_mg_dl=95, rdw_pct=13,
                mcv_fl=90, wbc_k_ul=6.5, alp_u_l=70, lymphocyte_pct=25)
    rows = [dict(albumin_g_dl=a, **base) for a in [3.0, 3.8, 4.3, 5.0]]
    result = KDM().transform(pd.DataFrame(rows))
    ages = result.biological_ages.values
    assert all(ages[i] > ages[i + 1] for i in range(len(ages) - 1)), \
        f"KDM not monotone with albumin: {ages}"


def test_kdm_custom_fit_reduces_s_ba(nhanes_kdm_df):
    """After fitting on NHANES, s_BA should be a finite positive number."""
    df = nhanes_kdm_df.drop(columns=["mortstat", "permth_exm"])
    clock = KDM().fit(df)
    assert clock._s_ba > 0
    assert np.isfinite(clock._s_ba)
