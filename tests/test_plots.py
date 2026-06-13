"""Tests for visualization functions (Week 3)."""

import pandas as pd
import pytest

from agingclockbench import PhenoAge, KDM, BenchmarkSuite
from agingclockbench.datasets import load_nhanes_sample


@pytest.fixture(scope="module")
def nhanes_with_results():
    df = load_nhanes_sample()
    pa_result = PhenoAge().transform(df)
    kdm_df = df.drop(columns=["crp_mg_l"])
    kdm_result = KDM().transform(kdm_df)
    suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
    report = suite.run(df, {"PhenoAge": pa_result, "KDM": kdm_result})
    return df, {"PhenoAge": pa_result, "KDM": kdm_result}, report


def test_plot_comparison_returns_figure(nhanes_with_results):
    df, results, report = nhanes_with_results
    import matplotlib
    matplotlib.use("Agg")
    fig = report.plot_comparison()
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close("all")


def test_plot_km_survival_returns_figure(nhanes_with_results):
    df, results, report = nhanes_with_results
    import matplotlib
    matplotlib.use("Agg")
    fig = report.plot_km_survival()
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close("all")


def test_plot_correlation_heatmap_returns_figure(nhanes_with_results):
    df, results, report = nhanes_with_results
    import matplotlib
    matplotlib.use("Agg")
    fig = report.plot_correlation_heatmap()
    assert fig is not None
    import matplotlib.pyplot as plt
    plt.close("all")


def test_to_html_creates_file(nhanes_with_results, tmp_path):
    df, results, report = nhanes_with_results
    html_path = str(tmp_path / "report.html")
    report.to_html(html_path)
    import os
    assert os.path.exists(html_path)
    with open(html_path) as f:
        content = f.read()
    assert "AgingClockBench" in content
    # Altair (vega-embed) is the primary renderer; Plotly CDN only loads on fallback
    assert "vega" in content.lower()


def test_to_html_contains_table(nhanes_with_results, tmp_path):
    df, results, report = nhanes_with_results
    html_path = str(tmp_path / "report2.html")
    report.to_html(html_path)
    with open(html_path) as f:
        content = f.read()
    assert "PhenoAge" in content
    assert "KDM" in content
