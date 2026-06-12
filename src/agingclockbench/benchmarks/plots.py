"""Visualization functions for BenchmarkReport.

All functions return matplotlib/plotly Figure objects so callers can
save, display, or embed them as needed.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from agingclockbench.benchmarks.suite import BenchmarkReport
    from agingclockbench.clocks.base import ClockResult


def plot_comparison(
    report: "BenchmarkReport",
    df: pd.DataFrame,
    results: dict[str, "ClockResult"],
):
    """Scatter plot of biological age vs chronological age for each clock.

    Parameters
    ----------
    report : BenchmarkReport from BenchmarkSuite.run()
    df : original input DataFrame (must contain 'age')
    results : dict mapping clock name -> ClockResult

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5), squeeze=False)
    palette = sns.color_palette("husl", n)

    for ax, (name, result), color in zip(axes[0], results.items(), palette):
        if result.original_index is not None:
            age = df.loc[result.original_index, "age"].values
        else:
            age = df["age"].iloc[: result.output_rows].values

        bio_age = result.biological_ages.values

        # Scatter with alpha for density
        ax.scatter(age, bio_age, alpha=0.3, s=8, color=color)

        # Identity line (biological age = chronological age)
        lo, hi = min(age.min(), bio_age.min()), max(age.max(), bio_age.max())
        ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="y = x")

        # Pearson r from benchmark results
        br = next((r for r in report.results if r.clock_name == name), None)
        r_str = f"r = {br.pearson_r:.3f}" if br and not np.isnan(br.pearson_r) else ""
        ax.set_title(f"{name}\n{r_str}", fontsize=12)
        ax.set_xlabel("Chronological Age (years)")
        ax.set_ylabel("Biological Age (years)")

    fig.suptitle("Biological Age vs Chronological Age", fontsize=14, y=1.02)
    plt.tight_layout()
    return fig


def plot_km_survival(
    df: pd.DataFrame,
    results: dict[str, "ClockResult"],
    mortality_col: str = "mortstat",
    followup_col: str = "permth_exm",
    n_quartiles: int = 4,
):
    """Kaplan-Meier survival curves stratified by age-acceleration quartile.

    Parameters
    ----------
    df : DataFrame with mortality columns.
    results : dict mapping clock name -> ClockResult.
    mortality_col : event indicator column (1=event, 0=censored).
    followup_col : time-to-event/censoring column (months).
    n_quartiles : number of strata (default 4).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter
    import seaborn as sns

    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5), squeeze=False)
    palette = sns.color_palette("RdYlGn_r", n_quartiles)

    for ax, (name, result) in zip(axes[0], results.items()):
        if result.original_index is not None:
            aligned = df.loc[result.original_index].reset_index(drop=True)
        else:
            aligned = df.iloc[: result.output_rows].reset_index(drop=True)

        if mortality_col not in aligned.columns or followup_col not in aligned.columns:
            ax.text(0.5, 0.5, "No mortality data", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(name)
            continue

        analysis = pd.DataFrame({
            "accel": result.accel.values,
            "event": aligned[mortality_col].values,
            "time": aligned[followup_col].values,
        }).dropna()

        quartile_labels = [f"Q{i+1}" for i in range(n_quartiles)]
        analysis["quartile"] = pd.qcut(analysis["accel"], n_quartiles,
                                       labels=quartile_labels)

        kmf = KaplanMeierFitter()
        for label, color in zip(quartile_labels, palette):
            mask = analysis["quartile"] == label
            kmf.fit(
                analysis.loc[mask, "time"] / 12,  # months → years
                analysis.loc[mask, "event"],
                label=label,
            )
            kmf.plot_survival_function(ax=ax, color=color, ci_show=False)

        ax.set_title(f"{name}\nKaplan-Meier by Accel Quartile", fontsize=11)
        ax.set_xlabel("Follow-up (years)")
        ax.set_ylabel("Survival Probability")
        ax.legend(title="Accel\nQuartile", fontsize=8)
        ax.set_ylim(0, 1)

    fig.suptitle("Survival by Biological Age Acceleration Quartile", fontsize=14, y=1.02)
    plt.tight_layout()
    return fig


def plot_correlation_heatmap(results: dict[str, "ClockResult"]):
    """Heatmap of Pearson correlations between clock accelerations.

    Parameters
    ----------
    results : dict mapping clock name -> ClockResult.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    names = list(results.keys())
    n = len(names)
    corr = np.eye(n)

    for i, n1 in enumerate(names):
        for j, n2 in enumerate(names):
            if i != j:
                a1 = results[n1].accel
                a2 = results[n2].accel
                min_len = min(len(a1), len(a2))
                if min_len > 2:
                    corr[i, j] = a1.iloc[:min_len].corr(a2.iloc[:min_len])

    corr_df = pd.DataFrame(corr, index=names, columns=names)
    fig, ax = plt.subplots(figsize=(max(4, n * 1.5), max(3, n * 1.5)))
    sns.heatmap(
        corr_df,
        annot=True,
        fmt=".3f",
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        ax=ax,
        square=True,
        cbar_kws={"shrink": 0.8},
    )
    ax.set_title("Inter-Clock Acceleration Correlations (Pearson r)", fontsize=12)
    plt.tight_layout()
    return fig


def to_html(
    report: "BenchmarkReport",
    df: pd.DataFrame,
    results: dict[str, "ClockResult"],
    filename: str,
    mortality_col: str = "mortstat",
    followup_col: str = "permth_exm",
) -> None:
    """Export an interactive Plotly HTML benchmark report.

    Parameters
    ----------
    report : BenchmarkReport from BenchmarkSuite.run()
    df : original input DataFrame
    results : dict mapping clock name -> ClockResult
    filename : output .html path
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.express as px

    # --- Benchmark table (built first; needed by Altair scatter) ---
    summary_df = report.to_dataframe()
    def _fmt(col):
        s = summary_df[col]
        if pd.api.types.is_numeric_dtype(s):
            return s.round(4).astype(str).tolist()
        return s.astype(str).tolist()

    fig_table = go.Figure(
        data=[go.Table(
            header=dict(
                values=list(summary_df.columns),
                fill_color="#2c3e50",
                font=dict(color="white", size=12),
                align="left",
            ),
            cells=dict(
                values=[_fmt(c) for c in summary_df.columns],
                fill_color="lavender",
                align="left",
            ),
        )]
    )
    # Calculate table height dynamically: ~40px per row + header + margin
    _table_height = max(300, 100 + 40 * (len(summary_df) + 1))
    fig_table.update_layout(title="Benchmark Summary", height=_table_height)

    # --- Scatter: Altair (preferred) with Plotly fallback ---
    altair_html_section: str | None = None
    try:
        from .altair_plots import generate_scatter_heatmap
        altair_html_section = generate_scatter_heatmap(
            df=df,
            summary_df=summary_df,
            results=results,
            top_n_clocks=2,
        )
    except Exception as _altair_err:
        import warnings
        warnings.warn(
            f"Altair scatter failed ({_altair_err}); falling back to Plotly.",
            stacklevel=2,
        )

    # Plotly fallback scatter (used if Altair is unavailable or errors)
    n = len(results)
    colors = px.colors.qualitative.Plotly
    fig_scatter = make_subplots(
        rows=1, cols=n,
        subplot_titles=list(results.keys()),
        shared_yaxes=False,
    )
    for col, (name, result) in enumerate(results.items(), start=1):
        if result.original_index is not None:
            age = df.loc[result.original_index, "age"].values
        else:
            age = df["age"].iloc[: result.output_rows].values
        bio_age = result.biological_ages.values
        br = next((r for r in report.results if r.clock_name == name), None)
        r_val = br.pearson_r if br else float("nan")
        fig_scatter.add_trace(
            go.Scatter(
                x=age, y=bio_age, mode="markers",
                marker=dict(size=6, color=colors[col - 1], opacity=0.65),
                name=f"{name} (r={r_val:.3f})",
            ),
            row=1, col=col,
        )
        lo = min(float(age.min()), float(bio_age.min()))
        hi = max(float(age.max()), float(bio_age.max()))
        fig_scatter.add_trace(
            go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines",
                       line=dict(color="black", dash="dash", width=1),
                       showlegend=False),
            row=1, col=col,
        )
    fig_scatter.update_layout(
        title="Biological Age vs Chronological Age",
        height=520, template="plotly_white",
    )

    # Combine into single HTML
    html_table = fig_table.to_html(full_html=False, include_plotlyjs=False)
    html_scatter_fallback = fig_scatter.to_html(full_html=False, include_plotlyjs=False)

    scatter_section = (
        f'<h2>Biological Age vs Chronological Age</h2>{altair_html_section}'
        if altair_html_section is not None
        else f'<h2>Biological Age vs Chronological Age</h2>{html_scatter_fallback}'
    )

    html = f"""<!DOCTYPE html>
<html>
<head>
  <title>AgingClockBench Report</title>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <style>
    body {{ font-family: Arial, sans-serif; max-width: 1200px; margin: auto; padding: 20px; }}
    h1 {{ color: #2c3e50; }}
    h2 {{ color: #34495e; border-bottom: 1px solid #bdc3c7; padding-bottom: 6px; }}
  </style>
</head>
<body>
  <h1>AgingClockBench Report</h1>
  <h2>Benchmark Summary</h2>
  {html_table}
  {scatter_section}
</body>
</html>"""

    with open(filename, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Report saved to {filename}")
