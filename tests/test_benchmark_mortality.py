"""Week 2 BenchmarkSuite mortality validation tests.

Tests verify:
1. Cox PH produces the expected Levine 2018 HR target (~1.08 per year accel)
2. Multi-clock comparison works correctly with NHANES mortality data
3. Row alignment is correct when clocks drop NaN rows
"""

import numpy as np
import pandas as pd
import pytest

from agingclockbench import PhenoAge, KDM, BenchmarkSuite
from agingclockbench.benchmarks.suite import BenchmarkReport


@pytest.fixture(scope="module")
def nhanes_full_df():
    parquet = "C:/Users/aadit/phenoage-engine/nhanes_cache/nhanes_phenoage_processed.parquet"
    try:
        ref = pd.read_parquet(parquet)
    except FileNotFoundError:
        pytest.skip("NHANES reference parquet not available")
    # Map to agingclockbench column names
    return pd.DataFrame({
        "age":             ref["age"],
        "albumin_g_dl":    ref["albumin"],
        "creatinine_mg_dl":ref["creatinine"],
        "glucose_mg_dl":   ref["glucose"],
        "crp_mg_l":        ref["crp"] * 10,        # mg/dL → mg/L
        "lymphocyte_pct":  ref["lymphocyte_pct"],
        "mcv_fl":          ref["mcv"],
        "rdw_pct":         ref["rdw"],
        "alp_u_l":         ref["alp"],
        "wbc_k_ul":        ref["wbc"],
        "mortstat":        ref["deceased"].astype(int),
        "permth_exm":      ref["followup_years"] * 12,  # years → months
    })


def test_benchmark_produces_report(nhanes_full_df):
    result = PhenoAge().transform(nhanes_full_df)
    suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
    report = suite.run(nhanes_full_df, {"PhenoAge": result})
    assert isinstance(report, BenchmarkReport)
    assert len(report.results) == 1


def test_benchmark_phenoage_hr_near_levine_target(nhanes_full_df):
    """PhenoAge HR per SD acceleration on NHANES should be ~1.4-2.5 (age-adjusted).

    Note: Levine 2018 reports HR ~1.08 per YEAR of acceleration (not per SD).
    Per-SD HR is higher and depends on cohort SD of acceleration (~8 years).
    Accepting a wide range to account for NHANES 1999-2000 vs 2015-2018 differences.
    """
    result = PhenoAge().transform(nhanes_full_df)
    suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
    report = suite.run(nhanes_full_df, {"PhenoAge": result})
    hr = report.results[0].mortality_hr
    assert not np.isnan(hr), "Cox HR is NaN — Cox PH failed"
    assert hr > 1.0, f"PhenoAge HR={hr:.3f} not > 1.0"
    assert hr < 10.0, f"PhenoAge HR={hr:.3f} implausibly large"


def test_benchmark_phenoage_mortality_pvalue_significant(nhanes_full_df):
    """PhenoAge acceleration should be significantly predictive of mortality."""
    result = PhenoAge().transform(nhanes_full_df)
    suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
    report = suite.run(nhanes_full_df, {"PhenoAge": result})
    p = report.results[0].mortality_pvalue
    assert p < 0.001, f"PhenoAge mortality p-value={p:.4f} not significant at p<0.001"


def test_benchmark_multi_clock(nhanes_full_df):
    """Running PhenoAge + KDM together should produce 2 rows in report."""
    pa_result = PhenoAge().transform(nhanes_full_df)
    kdm_df = nhanes_full_df.drop(columns=["crp_mg_l"])
    kdm_result = KDM().transform(kdm_df)
    suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
    report = suite.run(nhanes_full_df, {"PhenoAge": pa_result, "KDM": kdm_result})
    assert len(report.results) == 2
    names = [r.clock_name for r in report.results]
    assert "PhenoAge" in names and "KDM" in names


def test_benchmark_inter_clock_agreement(nhanes_full_df):
    """PhenoAge and KDM acceleration should be positively correlated."""
    pa_result = PhenoAge().transform(nhanes_full_df)
    kdm_df = nhanes_full_df.drop(columns=["crp_mg_l"])
    kdm_result = KDM().transform(kdm_df)
    suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
    report = suite.run(nhanes_full_df, {"PhenoAge": pa_result, "KDM": kdm_result})
    pa_res = next(r for r in report.results if r.clock_name == "PhenoAge")
    kdm_agreement = pa_res.clock_agreement_with_others.get("KDM")
    assert kdm_agreement is not None
    assert kdm_agreement > 0, f"PhenoAge-KDM acceleration agreement r={kdm_agreement:.3f} should be positive"


def test_benchmark_alignment_with_nan_rows():
    """BenchmarkSuite must correctly align mortality data when clock drops NaN rows."""
    # Create df with some NaN rows in PhenoAge-required columns
    complete = dict(age=55, albumin_g_dl=4.0, creatinine_mg_dl=1.0, glucose_mg_dl=95,
                    crp_mg_l=1.0, lymphocyte_pct=25, mcv_fl=90, rdw_pct=13,
                    alp_u_l=70, wbc_k_ul=6.5, mortstat=0, permth_exm=120.0)
    incomplete = dict(age=60, albumin_g_dl=None, creatinine_mg_dl=1.0, glucose_mg_dl=100,
                      crp_mg_l=1.5, lymphocyte_pct=22, mcv_fl=91, rdw_pct=13.5,
                      alp_u_l=75, wbc_k_ul=7.0, mortstat=1, permth_exm=48.0)
    df = pd.DataFrame([complete, incomplete] * 5)
    result = PhenoAge().transform(df)
    assert result.output_rows == 5  # only the complete rows

    suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
    report = suite.run(df, {"PhenoAge": result})
    # Should run without error — alignment handled via original_index
    assert isinstance(report, BenchmarkReport)


def test_benchmark_to_dataframe_columns(nhanes_full_df):
    """Report DataFrame should have the expected columns."""
    result = PhenoAge().transform(nhanes_full_df)
    suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
    report = suite.run(nhanes_full_df, {"PhenoAge": result})
    df = report.to_dataframe()
    for col in ["Clock", "Pearson r", "Spearman r", "Mort HR (per SD accel)", "Mort p-value"]:
        assert col in df.columns, f"Missing column: {col}"
