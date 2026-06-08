"""Command-line interface for AgingClockBench."""

import sys
from pathlib import Path

import click

from agingclockbench.config import VERSION

_ALL_CLOCKS = ["PhenoAge", "KDM", "DunedinPACEProxy"]


def _load_clocks(names: tuple[str]) -> dict:
    from agingclockbench import PhenoAge, KDM, DunedinPACEProxy
    registry = {"PhenoAge": PhenoAge, "KDM": KDM, "DunedinPACEProxy": DunedinPACEProxy}
    if "all" in names:
        return {k: v() for k, v in registry.items()}
    unknown = [n for n in names if n not in registry]
    if unknown:
        raise click.BadParameter(
            f"Unknown clock(s): {unknown}. Choose from {_ALL_CLOCKS} or 'all'."
        )
    return {n: registry[n]() for n in names}


@click.group()
@click.version_option(version=VERSION, prog_name="agingclockbench")
def cli():
    """AgingClockBench — benchmark biological aging clocks on your data.

    \b
    Quick start:
        agingclockbench benchmark --data bundled --clocks all
    """


@cli.command()
@click.option("--data", required=True,
              help="Path to input CSV, or 'bundled' to use the NHANES 1999-2000 sample.")
@click.option("--clocks", multiple=True, default=["PhenoAge"], show_default=True,
              help="Clocks to run. Repeat for multiple, or pass 'all'. "
                   f"Choices: {_ALL_CLOCKS}")
@click.option("--mortality-col", default="mortstat", show_default=True,
              help="Column for vital status (1=deceased, 0=censored).")
@click.option("--followup-col", default="permth_exm", show_default=True,
              help="Column for follow-up time in months.")
@click.option("--output", default="./results", show_default=True,
              help="Directory for output files.")
@click.option("--report", is_flag=True, default=False,
              help="Generate an interactive HTML report (requires plotly).")
@click.option("--verbose", is_flag=True, default=False,
              help="Print detailed progress.")
def benchmark(data, clocks, mortality_col, followup_col, output, report, verbose):
    """Run a benchmark comparison of aging clocks on your data.

    \b
    Examples:
        agingclockbench benchmark --data bundled --clocks all
        agingclockbench benchmark --data my_data.csv --clocks PhenoAge KDM --report
        agingclockbench benchmark --data my_data.csv --clocks all \\
            --mortality-col vital_status --followup-col followup_months
    """
    import pandas as pd
    from agingclockbench import BenchmarkSuite
    from agingclockbench.datasets import load_nhanes_sample

    # --- Load data ---
    if data == "bundled":
        df = load_nhanes_sample()
        click.echo(f"Loaded bundled NHANES 1999-2000 sample: {len(df)} participants.")
    else:
        try:
            df = pd.read_csv(data)
        except Exception as e:
            raise click.ClickException(f"Could not read {data}: {e}")
        click.echo(f"Loaded {len(df):,} rows from {data}.")

    # --- Run clocks ---
    selected = _load_clocks(clocks)
    results = {}
    for name, clock in selected.items():
        try:
            result = clock.transform(df)
            results[name] = result
            click.echo(
                f"  {name}: {result.output_rows}/{result.input_rows} rows "
                f"(mean BA = {result.biological_ages.mean():.1f} yr, "
                f"mean accel = {result.accel.mean():.1f} yr)"
            )
        except Exception as e:
            click.echo(f"  {name}: FAILED — {e}", err=True)

    if not results:
        raise click.ClickException("No clocks produced results. Exiting.")

    # --- Benchmark ---
    suite = BenchmarkSuite(mortality_col=mortality_col, followup_col=followup_col)
    bench_report = suite.run(df, results)

    # --- Output ---
    out_dir = Path(output)
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = bench_report.to_dataframe()
    csv_path = out_dir / "comparison_table.csv"
    summary.to_csv(csv_path, index=False)

    click.echo(f"\n{'='*60}")
    click.echo("BENCHMARK RESULTS")
    click.echo("="*60)
    click.echo(summary.to_string(index=False))
    click.echo(f"\nComparison table saved to: {csv_path}")

    if report:
        html_path = out_dir / "benchmark_report.html"
        try:
            bench_report.to_html(str(html_path))
            click.echo(f"Interactive report saved to: {html_path}")
        except ImportError:
            click.echo("plotly not installed — skipping HTML report. Run: pip install plotly",
                       err=True)


@cli.group()
def datasets():
    """Manage bundled reference datasets."""


@datasets.command("list")
def datasets_list():
    """List available bundled datasets."""
    click.echo("Available bundled datasets:")
    click.echo("  nhanes_sample  — NHANES 1999-2000 (N=4,086, mortality-linked)")
    click.echo("                   Columns: age, sex, all 9 PhenoAge biomarkers,")
    click.echo("                            mortstat, permth_exm")


@datasets.command("info")
def datasets_info():
    """Print summary statistics for the bundled NHANES sample."""
    from agingclockbench.datasets import load_nhanes_sample
    df = load_nhanes_sample()
    click.echo(f"NHANES 1999-2000 sample — {len(df):,} participants")
    click.echo(f"  Age: {df.age.min():.0f}–{df.age.max():.0f} yr  "
               f"(mean {df.age.mean():.1f}, SD {df.age.std():.1f})")
    click.echo(f"  Deaths: {df.mortstat.sum():,} ({df.mortstat.mean()*100:.1f}%)")
    click.echo(f"  Median follow-up: {df.permth_exm.median()/12:.1f} yr")


@cli.command()
def version():
    """Print version and exit."""
    click.echo(f"agingclockbench {VERSION}")
