"""Tests for Click CLI."""

from click.testing import CliRunner

from agingclockbench.cli import cli


def test_cli_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "benchmark" in result.output


def test_cli_version():
    runner = CliRunner()
    result = runner.invoke(cli, ["version"])
    assert result.exit_code == 0
    assert "agingclockbench" in result.output


def test_cli_datasets_list():
    runner = CliRunner()
    result = runner.invoke(cli, ["datasets", "list"])
    assert result.exit_code == 0
    assert "nhanes_sample" in result.output


def test_cli_benchmark_with_csv(tmp_path, phenoage_df):
    import pandas as pd

    csv_path = tmp_path / "test_input.csv"
    phenoage_df.to_csv(csv_path, index=False)
    out_dir = tmp_path / "out"

    runner = CliRunner()
    result = runner.invoke(cli, [
        "benchmark",
        "--data", str(csv_path),
        "--clocks", "PhenoAge",
        "--output", str(out_dir),
    ])
    assert result.exit_code == 0, result.output
    assert (out_dir / "comparison_table.csv").exists()
