"""Command-line interface for AgingClockBench."""

import click

from agingclockbench.config import VERSION


@click.group()
@click.version_option(version=VERSION)
def cli():
    """AgingClockBench: benchmark biological aging clocks on your data."""


@cli.command()
@click.option("--data", required=True, help="Path to input CSV or 'bundled' for NHANES sample.")
@click.option("--clocks", multiple=True, default=["PhenoAge"], show_default=True,
              help="Clocks to run: PhenoAge, KDM, DunedinPACEProxy, or 'all'.")
@click.option("--mortality-col", default="mortstat", show_default=True,
              help="Column name for vital status (1=dead, 0=censored).")
@click.option("--followup-col", default="permth_exm", show_default=True,
              help="Column name for follow-up time in months.")
@click.option("--output", default="./results", show_default=True,
              help="Output directory for results.")
@click.option("--report", is_flag=True, default=False,
              help="Generate an HTML benchmark report.")
@click.option("--verbose", is_flag=True, default=False)
def benchmark(data, clocks, mortality_col, followup_col, output, report, verbose):
    """Run a benchmark comparison of aging clocks on your data."""
    import pandas as pd
    from pathlib import Path
    from agingclockbench import PhenoAge, KDM, DunedinPACEProxy, BenchmarkSuite
    from agingclockbench.datasets import load_nhanes_sample

    if data == "bundled":
        df = load_nhanes_sample()
        click.echo("Loaded bundled NHANES 2015-2018 sample.")
    else:
        df = pd.read_csv(data)
        click.echo(f"Loaded {len(df)} rows from {data}.")

    all_clocks = {"PhenoAge": PhenoAge(), "KDM": KDM(), "DunedinPACEProxy": DunedinPACEProxy()}
    selected = all_clocks if "all" in clocks else {k: v for k, v in all_clocks.items() if k in clocks}

    results = {}
    for name, clock in selected.items():
        click.echo(f"Running {name}...")
        results[name] = clock.transform(df)

    suite = BenchmarkSuite(mortality_col=mortality_col, followup_col=followup_col)
    bench_report = suite.run(df, results)

    out_dir = Path(output)
    out_dir.mkdir(parents=True, exist_ok=True)

    df_out = bench_report.to_dataframe()
    csv_path = out_dir / "comparison_table.csv"
    df_out.to_csv(csv_path, index=False)
    click.echo(f"\n{df_out.to_string(index=False)}")
    click.echo(f"\nResults saved to {csv_path}")


@cli.group()
def datasets():
    """Manage bundled reference datasets."""


@datasets.command("list")
def datasets_list():
    """List available bundled datasets."""
    click.echo("Available datasets:\n  nhanes_sample — NHANES 2015-2018 (5000 rows, mortality-linked)")


@cli.command()
def version():
    """Print version."""
    click.echo(f"agingclockbench {VERSION}")
