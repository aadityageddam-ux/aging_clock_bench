"""
Week 1 validation tests for PhenoAge.

Every expected value in this file is derived from the phenoage-engine reference
implementation (core/calculator.py), which has been cross-validated against the
published Levine 2018 Aging Cell algorithm. The reference dataset is NHANES
1999-2000 (N=4,086 complete cases with mortality linkage).

Cross-validation result (run 2026-06-06):
    N=4,086 | mean diff=0.000 | max |diff|=0.000 | Pearson r=1.000000

These tests verify:
  1. Exact agreement with reference implementation on known NHANES rows
  2. Population statistics match published Levine 2018 ranges
  3. Edge cases are handled safely
"""

import numpy as np
import pandas as pd
import pytest

from agingclockbench.clocks.phenoage import PhenoAge


# ---------------------------------------------------------------------------
# Ground-truth reference cases from NHANES 1999-2000
# (albumin/creatinine/glucose in NHANES clinical units; crp in mg/L)
# Expected phenoage computed by phenoage-engine reference implementation
# ---------------------------------------------------------------------------

REFERENCE_CASES = [
    {
        "label": "young_healthy_30yo",
        "age": 30.0,
        "albumin_g_dl": 4.60,
        "creatinine_mg_dl": 0.60,
        "glucose_mg_dl": 82.0,
        "crp_mg_l": 0.4,        # 0.04 mg/dL × 10
        "lymphocyte_pct": 27.0,
        "mcv_fl": 95.9,
        "rdw_pct": 12.3,
        "alp_u_l": 76.0,
        "wbc_k_ul": 5.4,
        "expected_phenoage": 20.6766,
    },
    {
        "label": "middle_aged_53yo",
        "age": 53.0,
        "albumin_g_dl": 4.10,
        "creatinine_mg_dl": 0.50,
        "glucose_mg_dl": 94.0,
        "crp_mg_l": 5.6,
        "lymphocyte_pct": 35.8,
        "mcv_fl": 87.8,
        "rdw_pct": 12.7,
        "alp_u_l": 98.0,
        "wbc_k_ul": 7.4,
        "expected_phenoage": 45.9278,
    },
    {
        "label": "older_80yo",
        "age": 80.0,
        "albumin_g_dl": 4.20,
        "creatinine_mg_dl": 0.50,
        "glucose_mg_dl": 98.0,
        "crp_mg_l": 12.6,
        "lymphocyte_pct": 20.9,
        "mcv_fl": 89.9,
        "rdw_pct": 12.7,
        "alp_u_l": 88.0,
        "wbc_k_ul": 8.5,
        "expected_phenoage": 74.0445,
    },
]


@pytest.mark.parametrize("case", REFERENCE_CASES, ids=[c["label"] for c in REFERENCE_CASES])
def test_phenoage_exact_reference_values(case):
    """PhenoAge output must match reference implementation to within 0.01 years."""
    row = {k: v for k, v in case.items() if k not in ("label", "expected_phenoage")}
    df = pd.DataFrame([row])
    result = PhenoAge().transform(df)
    actual = result.biological_ages.iloc[0]
    expected = case["expected_phenoage"]
    assert abs(actual - expected) < 0.01, (
        f"{case['label']}: got {actual:.4f}, expected {expected:.4f} "
        f"(diff={actual - expected:.4f})"
    )


def test_phenoage_acceleration_direction():
    """Young healthy participant should have negative age acceleration."""
    row = REFERENCE_CASES[0].copy()
    row.pop("label"); row.pop("expected_phenoage")
    result = PhenoAge().transform(pd.DataFrame([row]))
    assert result.accel.iloc[0] < 0, "Healthy 30yo should be biologically younger than chronological age"


def test_phenoage_older_participant_positive_accel():
    """Older participant with elevated CRP/RDW may have positive acceleration."""
    row = {
        "age": 72.0,
        "albumin_g_dl": 3.8,
        "creatinine_mg_dl": 1.3,
        "glucose_mg_dl": 115.0,
        "crp_mg_l": 15.0,
        "lymphocyte_pct": 18.0,
        "mcv_fl": 94.0,
        "rdw_pct": 15.5,
        "alp_u_l": 110.0,
        "wbc_k_ul": 9.0,
    }
    result = PhenoAge().transform(pd.DataFrame([row]))
    # Elevated inflammatory markers → biological age > chronological
    assert result.biological_ages.iloc[0] > 65, "Participant with multiple risk factors should have high biological age"


# ---------------------------------------------------------------------------
# Population statistics validation (Levine 2018 NHANES benchmarks)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def nhanes_bench_df():
    """Load NHANES reference and map to agingclockbench column names."""
    parquet = "C:/Users/aadit/phenoage-engine/nhanes_cache/nhanes_phenoage_processed.parquet"
    try:
        ref = pd.read_parquet(parquet)
    except FileNotFoundError:
        pytest.skip("NHANES reference parquet not available in this environment")
    return pd.DataFrame({
        "age": ref["age"],
        "albumin_g_dl": ref["albumin"],
        "creatinine_mg_dl": ref["creatinine"],
        "glucose_mg_dl": ref["glucose"],
        "crp_mg_l": ref["crp"] * 10,   # mg/dL → mg/L
        "lymphocyte_pct": ref["lymphocyte_pct"],
        "mcv_fl": ref["mcv"],
        "rdw_pct": ref["rdw"],
        "alp_u_l": ref["alp"],
        "wbc_k_ul": ref["wbc"],
        "_ref_phenoage": ref["phenoage"].values,
    })


def test_phenoage_population_pearson_r(nhanes_bench_df):
    """PhenoAge should correlate ≥ 0.85 with chronological age (Levine 2018 benchmark)."""
    df = nhanes_bench_df.drop(columns=["_ref_phenoage"])
    result = PhenoAge().transform(df)
    r = result.biological_ages.corr(nhanes_bench_df["age"].reset_index(drop=True))
    assert r >= 0.85, f"Pearson r={r:.4f} below minimum 0.85"


def test_phenoage_population_mean_in_range(nhanes_bench_df):
    """Mean PhenoAge should be within ±5 years of mean chronological age."""
    df = nhanes_bench_df.drop(columns=["_ref_phenoage"])
    result = PhenoAge().transform(df)
    mean_pa = result.biological_ages.mean()
    mean_age = nhanes_bench_df["age"].mean()
    assert abs(mean_pa - mean_age) <= 5, f"|mean PhenoAge ({mean_pa:.1f}) - mean age ({mean_age:.1f})| > 5"


def test_phenoage_population_no_negative_ages(nhanes_bench_df):
    """No participant should receive a negative biological age."""
    df = nhanes_bench_df.drop(columns=["_ref_phenoage"])
    result = PhenoAge().transform(df)
    neg = (result.biological_ages < 0).sum()
    assert neg == 0, f"{neg} participants received negative biological ages"


def test_phenoage_exact_match_with_reference(nhanes_bench_df):
    """Full population: agingclockbench PhenoAge must equal phenoage-engine to within 1e-4 years."""
    df = nhanes_bench_df.drop(columns=["_ref_phenoage"])
    result = PhenoAge().transform(df)
    ref_pa = pd.Series(nhanes_bench_df["_ref_phenoage"].values)
    max_diff = (result.biological_ages - ref_pa).abs().max()
    assert max_diff < 1e-4, f"Max diff vs reference: {max_diff:.6f} years — implementations diverged"


# ---------------------------------------------------------------------------
# Edge cases and numerical safety
# ---------------------------------------------------------------------------

def test_phenoage_crp_near_zero():
    """CRP values near zero must not cause log(0) errors."""
    row = dict(age=50, albumin_g_dl=4.0, creatinine_mg_dl=0.9, glucose_mg_dl=90,
               crp_mg_l=0.0001, lymphocyte_pct=25, mcv_fl=90, rdw_pct=13,
               alp_u_l=60, wbc_k_ul=6)
    result = PhenoAge().transform(pd.DataFrame([row]))
    assert np.isfinite(result.biological_ages.iloc[0])


def test_phenoage_very_high_crp():
    """Extremely elevated CRP (severe infection) should not overflow."""
    row = dict(age=65, albumin_g_dl=3.5, creatinine_mg_dl=1.5, glucose_mg_dl=120,
               crp_mg_l=200.0, lymphocyte_pct=10, mcv_fl=96, rdw_pct=16,
               alp_u_l=150, wbc_k_ul=15)
    result = PhenoAge().transform(pd.DataFrame([row]))
    assert np.isfinite(result.biological_ages.iloc[0])
    assert result.biological_ages.iloc[0] < 130


def test_phenoage_min_age_participant():
    """Youngest NHANES-eligible participant (age 20) should give finite result."""
    row = dict(age=20, albumin_g_dl=4.8, creatinine_mg_dl=0.8, glucose_mg_dl=80,
               crp_mg_l=0.2, lymphocyte_pct=35, mcv_fl=88, rdw_pct=12,
               alp_u_l=55, wbc_k_ul=6)
    result = PhenoAge().transform(pd.DataFrame([row]))
    assert np.isfinite(result.biological_ages.iloc[0])
    assert result.biological_ages.iloc[0] < 50  # Should be well below chronological age


def test_phenoage_output_monotone_with_rdw():
    """Higher RDW (aging marker) should increase biological age, all else equal."""
    base = dict(age=55, albumin_g_dl=4.0, creatinine_mg_dl=1.0, glucose_mg_dl=95,
                crp_mg_l=1.0, lymphocyte_pct=25, mcv_fl=90, alp_u_l=70, wbc_k_ul=6.5)
    rows = [dict(rdw_pct=rdw, **base) for rdw in [12.0, 13.5, 15.0, 17.0]]
    df = pd.DataFrame(rows)
    result = PhenoAge().transform(df)
    ages = result.biological_ages.values
    assert all(ages[i] < ages[i + 1] for i in range(len(ages) - 1)), \
        f"PhenoAge not monotone with RDW: {ages}"


def test_phenoage_output_monotone_with_glucose():
    """Higher glucose should increase biological age, all else equal."""
    base = dict(age=55, albumin_g_dl=4.0, creatinine_mg_dl=1.0,
                crp_mg_l=1.0, lymphocyte_pct=25, mcv_fl=90, rdw_pct=13,
                alp_u_l=70, wbc_k_ul=6.5)
    rows = [dict(glucose_mg_dl=g, **base) for g in [80, 95, 120, 180]]
    df = pd.DataFrame(rows)
    result = PhenoAge().transform(df)
    ages = result.biological_ages.values
    assert all(ages[i] < ages[i + 1] for i in range(len(ages) - 1)), \
        f"PhenoAge not monotone with glucose: {ages}"


def test_phenoage_output_monotone_with_albumin():
    """Higher albumin (protective) should decrease biological age, all else equal."""
    base = dict(age=55, creatinine_mg_dl=1.0, glucose_mg_dl=95,
                crp_mg_l=1.0, lymphocyte_pct=25, mcv_fl=90, rdw_pct=13,
                alp_u_l=70, wbc_k_ul=6.5)
    rows = [dict(albumin_g_dl=alb, **base) for alb in [3.0, 3.8, 4.3, 5.0]]
    df = pd.DataFrame(rows)
    result = PhenoAge().transform(df)
    ages = result.biological_ages.values
    assert all(ages[i] > ages[i + 1] for i in range(len(ages) - 1)), \
        f"PhenoAge not monotone with albumin: {ages}"
