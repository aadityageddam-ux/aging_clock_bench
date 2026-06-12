"""
Altair-based visualization functions for AgingClockBench reports.

Provides production-grade scatter plots with density visualization, statistical
annotations, and portfolio-quality aesthetics. Designed to replace raw Plotly
scatter plots which suffer from over-plotting at N=4,000+.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from agingclockbench.clocks.base import ClockResult


# Professional color palette — accessible, distinct at small sizes
_CLOCK_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]


def generate_scatter_heatmap(
    df: pd.DataFrame,
    summary_df: pd.DataFrame,
    results: dict[str, "ClockResult"],
    top_n_clocks: int = 2,
) -> str:
    """Generate an Altair-based scatter visualization for biological vs chronological age.

    Selects the top N clocks by Pearson r and renders an interactive scatter
    with a dashed identity line (y=x) and a styled statistics table.

    Parameters
    ----------
    df : original input DataFrame (must contain 'age' column)
    summary_df : benchmark summary DataFrame from BenchmarkReport.to_dataframe()
    results : dict mapping clock name -> ClockResult
    top_n_clocks : number of clocks to display (2 recommended for clarity)

    Returns
    -------
    str — HTML fragment containing the Altair chart + stats table
    """
    import altair as alt

    # Altair 5+ enforces a 5,000-row limit; disable it so all points render inline.
    alt.data_transformers.disable_max_rows()

    # Select top N clocks by Pearson r
    top_clocks = (
        summary_df.nlargest(top_n_clocks, "Pearson r")["Clock"].tolist()
        if "Pearson r" in summary_df.columns
        else list(results.keys())[:top_n_clocks]
    )

    # Build combined DataFrame — one row per participant per clock
    frames = []
    for clock_name in top_clocks:
        result = results.get(clock_name)
        if result is None:
            continue
        # Align chronological age with the rows the clock actually processed
        if result.original_index is not None:
            age = df.loc[result.original_index, "age"].values
        else:
            age = df["age"].iloc[: result.output_rows].values

        bio_age = result.biological_ages.values
        frames.append(
            pd.DataFrame(
                {
                    "Chronological Age": age.astype(float),
                    "Biological Age": bio_age.astype(float),
                    "Clock": clock_name,
                }
            )
        )

    if not frames:
        return "<p>No clock data available for visualization.</p>"

    plot_df = pd.concat(frames, ignore_index=True)

    # Interactive selection — click legend to highlight a clock
    selection = alt.selection_point(fields=["Clock"], bind="legend")

    scatter = (
        alt.Chart(plot_df)
        .mark_point(size=25, filled=True, stroke=None)
        .encode(
            x=alt.X(
                "Chronological Age:Q",
                scale=alt.Scale(zero=False),
                title="Chronological Age (years)",
                axis=alt.Axis(grid=True, gridColor="#e0e0e0"),
            ),
            y=alt.Y(
                "Biological Age:Q",
                scale=alt.Scale(zero=False),
                title="Biological Age (years)",
                axis=alt.Axis(grid=True, gridColor="#e0e0e0"),
            ),
            color=alt.Color(
                "Clock:N",
                scale=alt.Scale(
                    domain=top_clocks,
                    range=_CLOCK_COLORS[:len(top_clocks)],
                ),
                title="Clock Algorithm",
                legend=alt.Legend(orient="top-right"),
            ),
            opacity=alt.condition(selection, alt.value(0.55), alt.value(0.08)),
            tooltip=[
                alt.Tooltip("Chronological Age:Q", format=".1f"),
                alt.Tooltip("Biological Age:Q", format=".1f"),
                "Clock:N",
            ],
        )
        .add_params(selection)
    )

    # Identity line y = x
    age_range = [
        float(plot_df["Chronological Age"].min()),
        float(plot_df["Chronological Age"].max()),
    ]
    identity_df = pd.DataFrame(
        {"Chronological Age": age_range, "Biological Age": age_range}
    )
    identity_line = (
        alt.Chart(identity_df)
        .mark_line(strokeDash=[6, 4], color="#555555", size=1.5, opacity=0.6)
        .encode(
            x="Chronological Age:Q",
            y="Biological Age:Q",
        )
    )

    chart = (
        (scatter + identity_line)
        .properties(
            width=520,
            height=460,
            title=alt.TitleParams(
                "Biological Age vs Chronological Age",
                fontSize=15,
                fontWeight="bold",
                color="#2c3e50",
            ),
        )
        .configure_view(strokeWidth=0)
        .configure_axis(labelFontSize=11, titleFontSize=12)
        .interactive()
    )

    chart_html = chart.to_html(embed_options={"actions": False})
    stats_html = _stats_table(summary_df, top_clocks)
    excluded = [c for c in results if c not in top_clocks]
    note = (
        f'<p style="font-size:12px;color:#888;margin-top:8px;">'
        f"<strong>Note:</strong> {', '.join(excluded)} excluded from scatter "
        f"to reduce visual clutter. Full metrics in the benchmark table above.</p>"
        if excluded
        else ""
    )

    return f"""
<div style="font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;
            padding:0 0 24px 0;">
  <p style="color:#555;font-size:13px;margin:0 0 12px 0;">
    Top {len(top_clocks)} clocks by Pearson r shown. Click a clock name in the
    legend to highlight it. Dashed line = perfect prediction (y&nbsp;=&nbsp;x).
  </p>
  {chart_html}
  <h3 style="margin:24px 0 8px 0;color:#2c3e50;">Clock Performance Summary</h3>
  {stats_html}
  {note}
</div>
"""


def _stats_table(summary_df: pd.DataFrame, clock_names: list[str]) -> str:
    """Render a styled HTML stats table for the selected clocks."""
    filtered = summary_df[summary_df["Clock"].isin(clock_names)].copy()

    # Columns to show and how to format them
    col_map = {
        "Clock": ("Clock", lambda v: str(v)),
        "Pearson r": ("Pearson r", lambda v: f"{v:.3f}"),
        "Spearman r": ("Spearman r", lambda v: f"{v:.3f}"),
        "Mort HR (per SD accel)": ("Mort HR (per SD accel)", lambda v: f"{v:.2f}"),
        "Mort p-value": ("Mort p-value", lambda v: "<0.0001" if v < 0.0001 else f"{v:.4f}"),
        "CV": ("CV", lambda v: f"{v:.3f}"),
    }
    present = {k: v for k, v in col_map.items() if k in filtered.columns}

    header_cells = "".join(
        f'<th style="padding:9px 14px;text-align:{"left" if k=="Clock" else "center"};'
        f'font-weight:600;color:#fff;">{label}</th>'
        for k, (label, _) in present.items()
    )

    rows_html = ""
    for i, (_, row) in enumerate(filtered.iterrows()):
        bg = "#f9f9f9" if i % 2 == 0 else "#ffffff"
        cells = "".join(
            f'<td style="padding:9px 14px;text-align:{"left" if k=="Clock" else "center"};">'
            f"{fmt(row[k])}</td>"
            for k, (_, fmt) in present.items()
        )
        rows_html += f'<tr style="background:{bg};border-bottom:1px solid #e8e8e8;">{cells}</tr>'

    return f"""
<table style="border-collapse:collapse;width:100%;max-width:860px;font-size:13px;
              box-shadow:0 1px 3px rgba(0,0,0,0.08);border-radius:4px;overflow:hidden;">
  <thead style="background:#2c3e50;">
    <tr>{header_cells}</tr>
  </thead>
  <tbody>{rows_html}</tbody>
</table>"""
