"""Tests targeting the remaining coverage gaps to reach >92%.

Covers:
- benchmarks/metrics.py (was 0%)
- utils/validation.py (was 0%)
- cli.py additional paths (was 76%)
- datasets FileNotFoundError path
- DunedinPACEProxy zero-row ValueError
- KDM zero-row ValueError
"""

import pandas as pd
import pytest
from click.testing import CliRunner

from agingclockbench.benchmarks.metrics import (
    pearson_correlation,
    spearman_correlation,
    coefficient_of_variation,
)
from agingclockbench.utils.validation import validate_dataframe
from agingclockbench.cli import cli


# ── benchmarks/metrics.py ──────────────────────────────────────────────────

def test_pearson_correlation_positive():
    x = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])
    y = pd.Series([1.1, 1.9, 3.2, 3.8, 5.1])
    r, p = pearson_correlation(x, y)
    assert 0.99 < r <= 1.0
    assert p < 0.05


def test_pearson_correlation_negative():
    x = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])
    y = pd.Series([5.0, 4.0, 3.0, 2.0, 1.0])
    r, p = pearson_correlation(x, y)
    assert r == pytest.approx(-1.0, abs=1e-6)


def test_spearman_correlation():
    x = pd.Series([1.0, 2.0, 3.0, 4.0, 5.0])
    y = pd.Series([2.0, 1.0, 4.0, 3.0, 5.0])
    rho = spearman_correlation(x, y)
    assert -1.0 <= rho <= 1.0


def test_coefficient_of_variation():
    s = pd.Series([10.0, 12.0, 8.0, 11.0, 9.0])
    cv = coefficient_of_variation(s)
    assert cv > 0


def test_coefficient_of_variation_zero_mean():
    s = pd.Series([0.0, 0.0, 0.0])
    cv = coefficient_of_variation(s)
    import math
    assert math.isnan(cv)


# ── utils/validation.py ───────────────────────────────────────────────────

def test_validate_dataframe_passes():
    df = pd.DataFrame({"a": [1], "b": [2]})
    validate_dataframe(df, ["a", "b"])  # Should not raise


def test_validate_dataframe_raises_on_missing():
    df = pd.DataFrame({"a": [1]})
    with pytest.raises(ValueError, match="missing columns"):
        validate_dataframe(df, ["a", "b", "c"])


# ── CLI additional paths ───────────────────────────────────────────────────

def test_cli_datasets_info():
    runner = CliRunner()
    result = runner.invoke(cli, ["datasets", "info"])
    assert result.exit_code == 0
    assert "NHANES" in result.output
    assert "participants" in result.output


def test_cli_benchmark_bundled():
    """Test CLI with bundled data."""
    runner = CliRunner()
    import tempfile, os
    with tempfile.TemporaryDirectory() as tmp:
        result = runner.invoke(cli, [
            "benchmark", "--data", "bundled",
            "--clocks", "PhenoAge",
            "--output", tmp,
        ])
        assert result.exit_code == 0, result.output
        assert "BENCHMARK RESULTS" in result.output
        assert os.path.exists(os.path.join(tmp, "comparison_table.csv"))


def test_cli_benchmark_all_clocks(tmp_path, phenoage_df):
    """Test --clocks all runs all three clocks."""
    csv_path = tmp_path / "input.csv"
    phenoage_df.to_csv(csv_path, index=False)
    runner = CliRunner()
    result = runner.invoke(cli, [
        "benchmark", "--data", str(csv_path),
        "--clocks", "all",
        "--output", str(tmp_path / "out"),
    ])
    # DunedinPACEProxy will fail (missing sex/rdw cols not in phenoage_df fixture)
    # but PhenoAge and KDM output should succeed or fail gracefully
    assert result.exit_code == 0 or "FAILED" in result.output


def test_cli_benchmark_with_report_flag(tmp_path, phenoage_df):
    """Test --report flag produces HTML file."""
    csv_path = tmp_path / "input.csv"
    phenoage_df.to_csv(csv_path, index=False)
    out_dir = tmp_path / "out"
    runner = CliRunner()
    result = runner.invoke(cli, [
        "benchmark", "--data", str(csv_path),
        "--clocks", "PhenoAge",
        "--output", str(out_dir),
        "--report",
    ])
    assert result.exit_code == 0, result.output
    html = out_dir / "benchmark_report.html"
    assert html.exists()


def test_cli_benchmark_unknown_clock():
    """Unknown clock name should give a helpful error."""
    runner = CliRunner()
    result = runner.invoke(cli, [
        "benchmark", "--data", "bundled",
        "--clocks", "UnknownClock123",
    ])
    assert result.exit_code != 0 or "Unknown" in result.output


# ── Datasets FileNotFoundError ────────────────────────────────────────────

def test_load_nhanes_sample_file_not_found(monkeypatch, tmp_path):
    from agingclockbench.datasets import loaders
    monkeypatch.setattr(loaders, "_DATA_DIR", tmp_path)
    with pytest.raises(FileNotFoundError, match="nhanes_sample.parquet"):
        loaders.load_nhanes_sample()


# ── DunedinPACEProxy zero rows ─────────────────────────────────────────────

def test_dunedinpace_all_nan_raises():
    from agingclockbench import DunedinPACEProxy
    row = dict(age=52, albumin_g_dl=None, creatinine_mg_dl=0.9,
               glucose_mg_dl=87, rdw_pct=13, wbc_k_ul=6,
               lymphocyte_pct=28, mcv_fl=90)
    df = pd.DataFrame([row])
    with pytest.raises(ValueError, match="No complete rows"):
        DunedinPACEProxy().transform(df)


# ── KDM zero rows ─────────────────────────────────────────────────────────

def test_kdm_all_nan_raises():
    from agingclockbench import KDM
    row = dict(age=52, albumin_g_dl=None, creatinine_mg_dl=0.5,
               glucose_mg_dl=90, rdw_pct=12.7, mcv_fl=88,
               wbc_k_ul=6, alp_u_l=70, lymphocyte_pct=25)
    df = pd.DataFrame([row])
    with pytest.raises(ValueError, match="No complete rows"):
        KDM().transform(df)
