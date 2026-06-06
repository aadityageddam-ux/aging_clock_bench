"""Tests for BenchmarkSuite."""

import numpy as np
import pandas as pd
import pytest

from agingclockbench import PhenoAge, BenchmarkSuite
from agingclockbench.benchmarks.suite import BenchmarkReport


def test_benchmark_suite_returns_report(phenoage_df):
    result = PhenoAge().transform(phenoage_df)
    suite = BenchmarkSuite()
    report = suite.run(phenoage_df, {"PhenoAge": result})
    assert isinstance(report, BenchmarkReport)


def test_benchmark_to_dataframe_shape(phenoage_df):
    result = PhenoAge().transform(phenoage_df)
    suite = BenchmarkSuite()
    report = suite.run(phenoage_df, {"PhenoAge": result})
    df = report.to_dataframe()
    assert len(df) == 1
    assert "Clock" in df.columns
    assert "Pearson r" in df.columns


def test_benchmark_pearson_r_valid_range():
    # Need varied ages to get a defined Pearson r
    rows = [
        dict(age=40, albumin_g_dl=4.5, creatinine_mg_dl=0.8, glucose_mg_dl=85,
             crp_mg_l=0.2, lymphocyte_pct=30, mcv_fl=88, rdw_pct=12.5, alp_u_l=55, wbc_k_ul=5.5),
        dict(age=55, albumin_g_dl=4.1, creatinine_mg_dl=1.0, glucose_mg_dl=95,
             crp_mg_l=0.5, lymphocyte_pct=25, mcv_fl=91, rdw_pct=13.5, alp_u_l=70, wbc_k_ul=6.5),
        dict(age=70, albumin_g_dl=3.8, creatinine_mg_dl=1.2, glucose_mg_dl=110,
             crp_mg_l=1.0, lymphocyte_pct=20, mcv_fl=93, rdw_pct=14.5, alp_u_l=90, wbc_k_ul=7.5),
    ] * 5
    df = pd.DataFrame(rows)
    result = PhenoAge().transform(df)
    suite = BenchmarkSuite()
    report = suite.run(df, {"PhenoAge": result})
    r = report.results[0].pearson_r
    assert -1.0 <= r <= 1.0


def test_benchmark_cv_positive(phenoage_df):
    result = PhenoAge().transform(phenoage_df)
    suite = BenchmarkSuite()
    report = suite.run(phenoage_df, {"PhenoAge": result})
    cv = report.results[0].cv
    assert cv >= 0 or np.isnan(cv)
