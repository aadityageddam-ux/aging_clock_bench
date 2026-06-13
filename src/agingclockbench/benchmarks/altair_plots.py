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

# Publication palette: vibrant teal (tight/confident) + bright orange (more variance)
_PUB_COLORS = {"PhenoAge": "#0277bd", "KDM": "#ff6f00", "DunedinPACEProxy": "#546e7a"}
_PUB_FALLBACK = ["#0277bd", "#ff6f00", "#546e7a", "#2c3e50"]

# Per-clock ellipse style: (hex_color, strokeDash, strokeWidth, opacity)
# Slot 0 = "better" clock (solid, thick); Slot 1 = comparison clock (dashed, slightly thinner)
_ELLIPSE_STYLES = [
    ("#0277bd", [],     3.5, 0.92),   # PhenoAge — teal, solid, prominent
    ("#ff6f00", [5, 5], 2.5, 0.88),   # KDM — orange, dashed, secondary
]


def compute_confidence_ellipse(
    x: np.ndarray,
    y: np.ndarray,
    confidence: float = 0.95,
    n_points: int = 120,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute a 2D confidence ellipse for a set of (x, y) points.

    Uses eigendecomposition of the covariance matrix scaled by the chi-squared
    critical value for the requested confidence level.

    Parameters
    ----------
    x, y : 1-D arrays of equal length
    confidence : probability mass enclosed (default 0.95 = 95 %)
    n_points : number of boundary points (higher = smoother curve)

    Returns
    -------
    (ellipse_x, ellipse_y) : boundary coordinates, each length n_points
    """
    from scipy import stats

    mean_x, mean_y = float(np.mean(x)), float(np.mean(y))
    cov = np.cov(x, y)

    # eigh is for symmetric matrices — returns REAL eigenvalues in ascending order
    eigenvalues, eigenvectors = np.linalg.eigh(cov)

    chi2_val = stats.chi2.ppf(confidence, df=2)
    # eigenvalues are [minor, major] after eigh (ascending order)
    b = np.sqrt(np.abs(eigenvalues[0]) * chi2_val)   # minor semi-axis
    a = np.sqrt(np.abs(eigenvalues[1]) * chi2_val)   # major semi-axis

    # Rotation angle from major eigenvector
    angle = np.arctan2(eigenvectors[1, 1], eigenvectors[0, 1])
    cos_a, sin_a = np.cos(angle), np.sin(angle)

    theta = np.linspace(0, 2 * np.pi, n_points)
    ex = a * np.cos(theta)
    ey = b * np.sin(theta)

    return (
        mean_x + ex * cos_a - ey * sin_a,
        mean_y + ex * sin_a + ey * cos_a,
    )


def generate_publication_chart(
    df: pd.DataFrame,
    summary_df: pd.DataFrame,
    results: dict[str, "ClockResult"],
    top_n_clocks: int = 2,
) -> str:
    """Publication-grade scatter + 95% confidence ellipses per clock.

    Designed to make the performance difference between clocks immediately visible:
    a tight ellipse = consistent predictions; a wide ellipse = high variance.
    Uses low-opacity scatter (density implied by stacking) + solid ellipse boundaries.

    Parameters
    ----------
    df : original DataFrame (must contain 'age' column)
    summary_df : BenchmarkReport.to_dataframe() output
    results : dict mapping clock name -> ClockResult
    top_n_clocks : number of clocks to show (default 2 = PhenoAge + KDM)

    Returns
    -------
    str — self-contained HTML fragment (Altair chart + stats table + caption)
    """
    import altair as alt

    alt.data_transformers.disable_max_rows()

    # Prefer PhenoAge + KDM (most scientifically distinct clocks)
    _PREFERRED_ORDER = ["PhenoAge", "KDM", "DunedinPACEProxy"]
    available = [c for c in _PREFERRED_ORDER if c in results]
    top_clocks = (
        available[:top_n_clocks] if len(available) >= top_n_clocks
        else list(results.keys())[:top_n_clocks]
    )

    # Assign colors: teal = tighter clock, rust = wider spread
    color_range = [
        _PUB_COLORS.get(name, _PUB_FALLBACK[i])
        for i, name in enumerate(top_clocks)
    ]

    # ── Build scatter data ──────────────────────────────────────────────────
    scatter_frames = []
    ellipse_frames = []

    for clock_name in top_clocks:
        result = results[clock_name]
        age = (
            df.loc[result.original_index, "age"].values
            if result.original_index is not None
            else df["age"].iloc[: result.output_rows].values
        )
        bio_age = result.biological_ages.values

        scatter_frames.append(pd.DataFrame({
            "Chronological Age": age.astype(float),
            "Biological Age": bio_age.astype(float),
            "Clock": clock_name,
        }))

        # Compute + store 95% confidence ellipse
        ex, ey = compute_confidence_ellipse(age, bio_age)
        ellipse_frames.append(pd.DataFrame({
            "Chronological Age": ex,
            "Biological Age": ey,
            "Clock": clock_name,
        }))

    scatter_df = pd.concat(scatter_frames, ignore_index=True)
    ellipse_df = pd.concat(ellipse_frames, ignore_index=True)

    color_enc = alt.Color(
        "Clock:N",
        scale=alt.Scale(domain=top_clocks, range=color_range),
        title="Clock Algorithm",
        legend=alt.Legend(orient="bottom", direction="horizontal"),
    )

    # ── Scatter layer (low opacity = density by stacking) ──────────────────
    scatter_layer = (
        alt.Chart(scatter_df)
        .mark_point(size=18, filled=True, stroke=None)
        .encode(
            x=alt.X("Chronological Age:Q", scale=alt.Scale(zero=False),
                    axis=alt.Axis(title="Chronological Age (years)", grid=True,
                                  gridColor="#ebebeb")),
            y=alt.Y("Biological Age:Q", scale=alt.Scale(zero=False),
                    axis=alt.Axis(title="Biological Age (years)", grid=True,
                                  gridColor="#ebebeb")),
            color=color_enc,
            opacity=alt.value(0.15),
            tooltip=[
                alt.Tooltip("Chronological Age:Q", format=".1f"),
                alt.Tooltip("Biological Age:Q", format=".1f"),
                "Clock:N",
            ],
        )
    )

    # ── Confidence ellipse layers (one per clock, explicit color + stroke style) ──
    # Splitting into separate layers guarantees each clock gets its own color and
    # stroke pattern — a shared color scale can blend or collapse in Altair.
    # Layer order: KDM first (bottom), PhenoAge on top (more prominent).
    ellipse_layers = []
    for i, clock_name in enumerate(reversed(top_clocks)):  # reversed = KDM first
        style_idx = top_clocks.index(clock_name)           # keep style aligned to position
        color, dash, width, opacity = _ELLIPSE_STYLES[min(style_idx, len(_ELLIPSE_STYLES) - 1)]
        clock_ellipse_df = ellipse_df[ellipse_df["Clock"] == clock_name]
        layer = (
            alt.Chart(clock_ellipse_df)
            .mark_line(
                strokeDash=dash if dash else alt.Undefined,
                size=width,
                opacity=opacity,
            )
            .encode(
                x=alt.X("Chronological Age:Q", scale=alt.Scale(zero=False)),
                y=alt.Y("Biological Age:Q", scale=alt.Scale(zero=False)),
                color=alt.value(color),
            )
        )
        ellipse_layers.append(layer)

    # ── Identity line (y = x) ───────────────────────────────────────────────
    age_min = float(scatter_df["Chronological Age"].min())
    age_max = float(scatter_df["Chronological Age"].max())
    identity_df = pd.DataFrame({
        "Chronological Age": [age_min, age_max],
        "Biological Age": [age_min, age_max],
    })
    identity_line = (
        alt.Chart(identity_df)
        .mark_line(strokeDash=[6, 4], color="#444444", size=1.8, opacity=0.7)
        .encode(x="Chronological Age:Q", y="Biological Age:Q")
    )

    # ── Compose layers: scatter → KDM ellipse → PhenoAge ellipse → identity ──
    base = scatter_layer
    for layer in ellipse_layers:
        base = base + layer
    chart = (
        (base + identity_line)
        .properties(
            width=530, height=480,
            title=alt.TitleParams(
                "Biological Age vs Chronological Age — 95% Confidence Ellipses",
                fontSize=14, fontWeight="bold", color="#2c3e50",
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
        f'<p style="font-size:12px;color:#888;margin-top:6px;">'
        f"<strong>Note:</strong> {', '.join(excluded)} excluded from scatter. "
        f"Full metrics in the table above.</p>"
    ) if excluded else ""

    return f"""
<div style="font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;
            padding:0 0 28px 0;">
  <p style="color:#555;font-size:13px;margin:0 0 10px 0;">
    Points at 15% opacity — darker regions = higher point density.
    Dashed ellipses enclose 95% of each clock&apos;s predictions.
    A <em>narrower</em> ellipse means more consistent biological age estimates.
    Dashed black line = perfect prediction (y&nbsp;=&nbsp;x).
  </p>
  {chart_html}
  <h3 style="margin:22px 0 8px 0;color:#2c3e50;">Clock Performance Summary</h3>
  {stats_html}
  {note}
</div>
"""


def generate_bland_altman_plots(
    df: pd.DataFrame,
    summary_df: pd.DataFrame,
    results: dict[str, "ClockResult"],
    top_n_clocks: int = 2,
) -> str:
    """Generate side-by-side Bland-Altman plots for clock validation.

    Bland-Altman is the standard method-comparison format in biomedical research.
    Shows systematic bias (mean error) and limits of agreement (±1.96 SD) for each
    clock, plotted against mean age to reveal age-dependent heteroscedasticity.

    Parameters
    ----------
    df : original DataFrame (must contain 'age' column)
    summary_df : BenchmarkReport.to_dataframe() output
    results : dict mapping clock name -> ClockResult
    top_n_clocks : number of clocks to show (default 2 = PhenoAge + KDM)

    Returns
    -------
    str — self-contained HTML fragment (hconcat Altair chart + stats table + legend)
    """
    import altair as alt

    alt.data_transformers.disable_max_rows()

    # Always prefer PhenoAge + KDM (most scientifically distinct)
    _PREFERRED_ORDER = ["PhenoAge", "KDM", "DunedinPACEProxy"]
    available = [c for c in _PREFERRED_ORDER if c in results]
    top_clocks = (
        available[:top_n_clocks] if len(available) >= top_n_clocks
        else list(results.keys())[:top_n_clocks]
    )

    # ── Per-clock computation (Python, not Altair transforms) ─────────────
    clock_charts = []
    interp_rows = []   # for the interpretation block below the chart

    for i, clock_name in enumerate(top_clocks):
        result = results[clock_name]
        age = (
            df.loc[result.original_index, "age"].values
            if result.original_index is not None
            else df["age"].iloc[: result.output_rows].values
        ).astype(float)
        bio = result.biological_ages.values.astype(float)

        mean_age = (age + bio) / 2.0
        diff = bio - age

        # Stats computed on FULL data (standard BA practice)
        mean_bias  = float(np.mean(diff))
        upper_loa  = mean_bias + 1.96 * float(np.std(diff, ddof=1))
        lower_loa  = mean_bias - 1.96 * float(np.std(diff, ddof=1))

        # X-axis domain: pin to 1st–99th percentile of chronological age
        # so extreme bio-age outliers (KDM can predict 200+ yrs) don't blow the scale
        x_lo = float(np.percentile(age, 1))
        x_hi = float(np.percentile(age, 99))

        color = _PUB_COLORS.get(clock_name, _PUB_FALLBACK[i])

        # Scatter data — clip to x-axis domain (stats already computed on full data)
        all_pts = pd.DataFrame({
            "Mean Age":         mean_age,
            "Prediction Error": diff,
        })
        pts_df = all_pts[(all_pts["Mean Age"] >= x_lo) & (all_pts["Mean Age"] <= x_hi)]

        # Reference lines — single-row DataFrames so Altair doesn't embed 4K rows per rule
        bias_df  = pd.DataFrame({"y": [mean_bias]})
        upper_df = pd.DataFrame({"y": [upper_loa]})
        lower_df = pd.DataFrame({"y": [lower_loa]})
        zero_df  = pd.DataFrame({"y": [0.0]})

        # Shared x/y encoding for scatter
        # x: fixed to chronological age range (clips extreme bio-age outliers)
        # y: fixed domain so both clocks are directly visually comparable
        x_enc = alt.X("Mean Age:Q", scale=alt.Scale(domain=[x_lo, x_hi]),
                      title="Mean Age (years)",
                      axis=alt.Axis(grid=True, gridColor="#ebebeb"))
        # Only show y-axis title on the leftmost plot to reduce clutter
        y_title = "Prediction Error: Bio − Chron (years)" if i == 0 else None
        y_enc = alt.Y("Prediction Error:Q",
                      scale=alt.Scale(domain=[-100, 100]),
                      title=y_title,
                      axis=alt.Axis(grid=True, gridColor="#ebebeb"))

        # ── Scatter points ──────────────────────────────────────────────
        pts = (
            alt.Chart(pts_df)
            .mark_point(size=28, filled=True, opacity=0.40)
            .encode(
                x=x_enc, y=y_enc,
                color=alt.value(color),
                tooltip=[
                    alt.Tooltip("Mean Age:Q", format=".1f", title="Mean Age"),
                    alt.Tooltip("Prediction Error:Q", format="+.1f", title="Error (yrs)"),
                ],
            )
        )

        # ── Mean bias line (solid, prominent) ──────────────────────────
        mean_line = (
            alt.Chart(bias_df)
            .mark_rule(size=2.5, opacity=0.92)
            .encode(y=alt.Y("y:Q"), color=alt.value(color))
        )

        # ── Upper LoA (dashed) ──────────────────────────────────────────
        upper_line = (
            alt.Chart(upper_df)
            .mark_rule(size=1.6, opacity=0.70, strokeDash=[6, 3])
            .encode(y=alt.Y("y:Q"), color=alt.value(color))
        )

        # ── Lower LoA (dashed) ──────────────────────────────────────────
        lower_line = (
            alt.Chart(lower_df)
            .mark_rule(size=1.6, opacity=0.70, strokeDash=[6, 3])
            .encode(y=alt.Y("y:Q"), color=alt.value(color))
        )

        # ── Zero reference (light grey) ─────────────────────────────────
        zero_line = (
            alt.Chart(zero_df)
            .mark_rule(size=1.0, opacity=0.30, strokeDash=[2, 2], color="#888888")
            .encode(y=alt.Y("y:Q"))
        )

        plot = (
            (pts + zero_line + lower_line + upper_line + mean_line)
            .properties(
                width=390, height=360,
                title=alt.TitleParams(
                    f"{clock_name} — Bland-Altman",
                    fontSize=13, fontWeight="bold", color="#2c3e50",
                ),
            )
        )
        clock_charts.append(plot)

        interp_rows.append({
            "name": clock_name,
            "color": color,
            "bias": mean_bias,
            "upper": upper_loa,
            "lower": lower_loa,
        })

    # ── hconcat: one Vega embed for both plots ─────────────────────────
    combined = (
        alt.hconcat(*clock_charts)
        .configure_view(strokeWidth=0)
        .configure_axis(labelFontSize=11, titleFontSize=12)
    )
    chart_html = combined.to_html(embed_options={"actions": False})

    # ── Interpretation block ────────────────────────────────────────────
    interp_html = ""
    for row in interp_rows:
        interp_html += f"""
<div style="text-align:center;min-width:180px;">
  <span style="display:inline-block;width:12px;height:12px;border-radius:2px;
               background:{row['color']};margin-right:6px;vertical-align:middle;"></span>
  <strong style="color:{row['color']};">{row['name']}</strong><br>
  <span style="font-size:12px;color:#555;line-height:1.8;">
    Mean bias: <strong>{row['bias']:+.2f} yrs</strong><br>
    95% LoA: <strong>{row['lower']:+.2f}</strong> to <strong>{row['upper']:+.2f} yrs</strong>
  </span>
</div>"""

    stats_html = _stats_table(summary_df, top_clocks)

    return f"""
<div style="font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;
            padding:0 0 28px 0;">
  <p style="color:#555;font-size:13px;margin:0 0 12px 0;">
    Solid line = mean bias. Dashed lines = 95% limits of agreement (&plusmn;1.96&nbsp;SD).
    Grey dotted line = zero error. Points should cluster around zero;
    systematic drift with age indicates heteroscedasticity.
    Independent y-axes: KDM has wider spread than PhenoAge by design.
  </p>
  {chart_html}
  <div style="display:flex;justify-content:center;gap:48px;
              flex-wrap:wrap;margin:14px 0 20px 0;">{interp_html}</div>
  <h3 style="margin:20px 0 8px 0;color:#2c3e50;">Clock Performance Summary</h3>
  {stats_html}
</div>
"""


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

    # Prefer PhenoAge + KDM for meaningful contrast (different algorithm families).
    # Fall back to Pearson r ranking only if those names aren't present.
    _PREFERRED_ORDER = ["PhenoAge", "KDM", "DunedinPACEProxy"]
    available = [c for c in _PREFERRED_ORDER if c in results]
    if len(available) >= top_n_clocks:
        top_clocks = available[:top_n_clocks]
    elif "Pearson r" in summary_df.columns:
        top_clocks = summary_df.nlargest(top_n_clocks, "Pearson r")["Clock"].tolist()
    else:
        top_clocks = list(results.keys())[:top_n_clocks]

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
                legend=alt.Legend(orient="bottom", direction="horizontal"),
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
        .mark_line(strokeDash=[6, 4], color="#555555", size=2, opacity=0.8)
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
