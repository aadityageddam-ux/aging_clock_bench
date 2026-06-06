"""BenchmarkSuite — runs validation metrics across multiple aging clocks."""

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from agingclockbench.clocks.base import ClockResult


@dataclass
class BenchmarkResult:
    """Validation metrics for a single clock."""

    clock_name: str
    pearson_r: float = float("nan")
    spearman_r: float = float("nan")
    pearson_pvalue: float = float("nan")
    mortality_hr: float = float("nan")
    mortality_hr_ci_lower: float = float("nan")
    mortality_hr_ci_upper: float = float("nan")
    mortality_pvalue: float = float("nan")
    cox_nobs: int = 0
    cv: float = float("nan")
    clock_agreement_with_others: dict = field(default_factory=dict)


class BenchmarkSuite:
    """Run a standardized validation benchmark on one or more aging clocks.

    Parameters
    ----------
    mortality_col : str
        Column name for vital status (1 = dead, 0 = censored).
    followup_col : str
        Column name for follow-up time in months.

    Examples
    --------
    >>> suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
    >>> report = suite.run(df, results={"PhenoAge": phenoage_result})
    >>> print(report.to_dataframe())
    """

    def __init__(self, mortality_col: str = "mortstat", followup_col: str = "permth_exm") -> None:
        self.mortality_col = mortality_col
        self.followup_col = followup_col

    def run(
        self,
        df: pd.DataFrame,
        results: dict[str, ClockResult],
    ) -> "BenchmarkReport":
        """Compute benchmark metrics for each clock result.

        Parameters
        ----------
        df : DataFrame — original input, must include ``age`` and optionally
             mortality columns.
        results : dict mapping clock name to ClockResult.

        Returns
        -------
        BenchmarkReport
        """
        benchmark_results: list[BenchmarkResult] = []
        accel_series: dict[str, pd.Series] = {}

        for name, result in results.items():
            br = BenchmarkResult(clock_name=name)

            # Align age to result rows (complete cases only)
            age = df["age"].iloc[: result.output_rows].reset_index(drop=True)

            # Pearson / Spearman correlation with chronological age
            r, p = stats.pearsonr(age, result.biological_ages)
            br.pearson_r = round(float(r), 4)
            br.pearson_pvalue = round(float(p), 6)
            br.spearman_r = round(float(stats.spearmanr(age, result.biological_ages).statistic), 4)

            # Coefficient of variation
            mean_ba = result.biological_ages.mean()
            std_ba = result.biological_ages.std()
            br.cv = round(float(std_ba / mean_ba), 4) if mean_ba != 0 else float("nan")

            # Cox PH mortality prediction (requires lifelines)
            if self.mortality_col in df.columns and self.followup_col in df.columns:
                br = self._run_cox(df, result, br)

            benchmark_results.append(br)
            accel_series[name] = result.accel

        # Inter-clock agreement
        for br in benchmark_results:
            others = {k: v for k, v in accel_series.items() if k != br.clock_name}
            for other_name, other_accel in others.items():
                min_len = min(len(accel_series[br.clock_name]), len(other_accel))
                r, _ = stats.pearsonr(
                    accel_series[br.clock_name].iloc[:min_len],
                    other_accel.iloc[:min_len],
                )
                br.clock_agreement_with_others[other_name] = round(float(r), 4)

        return BenchmarkReport(results=benchmark_results)

    def _run_cox(
        self,
        df: pd.DataFrame,
        result: ClockResult,
        br: BenchmarkResult,
    ) -> BenchmarkResult:
        try:
            from lifelines import CoxPHFitter

            analysis_df = df[[self.mortality_col, self.followup_col, "age"]].iloc[
                : result.output_rows
            ].copy().reset_index(drop=True)
            analysis_df["clock_acceleration"] = result.accel.values

            # Standardise acceleration for per-SD HR
            sd = analysis_df["clock_acceleration"].std()
            if sd > 0:
                analysis_df["clock_acceleration"] /= sd

            analysis_df = analysis_df.dropna()
            if len(analysis_df) < 10 or analysis_df[self.mortality_col].sum() == 0:
                return br

            cph = CoxPHFitter()
            cph.fit(
                analysis_df,
                duration_col=self.followup_col,
                event_col=self.mortality_col,
                formula="clock_acceleration + age",
            )
            summary = cph.summary
            row = summary.loc["clock_acceleration"]
            br.mortality_hr = round(float(np.exp(row["coef"])), 4)
            br.mortality_hr_ci_lower = round(float(np.exp(row["coef lower 95%"])), 4)
            br.mortality_hr_ci_upper = round(float(np.exp(row["coef upper 95%"])), 4)
            br.mortality_pvalue = round(float(row["p"]), 6)
            br.cox_nobs = int(cph.event_observed.sum())
        except Exception:
            pass  # mortality data unavailable or insufficient
        return br


class BenchmarkReport:
    """Container for all benchmark results with display and export methods."""

    def __init__(self, results: list[BenchmarkResult]) -> None:
        self.results = results

    def to_dataframe(self) -> pd.DataFrame:
        """Return a summary DataFrame — one row per clock."""
        rows = []
        for r in self.results:
            rows.append(
                {
                    "Clock": r.clock_name,
                    "Pearson r": r.pearson_r,
                    "Spearman r": r.spearman_r,
                    "Mort HR": r.mortality_hr,
                    "Mort p-value": r.mortality_pvalue,
                    "CV": r.cv,
                    "Cox N (events)": r.cox_nobs,
                }
            )
        return pd.DataFrame(rows)

    def plot_comparison(self):
        """Scatter plot of biological age vs chronological age per clock."""
        raise NotImplementedError("Plotting implemented in Week 3.")

    def plot_km_survival(self):
        """Kaplan-Meier survival curves stratified by acceleration quartile."""
        raise NotImplementedError("Plotting implemented in Week 3.")

    def to_html(self, filename: str) -> None:
        """Export interactive Plotly HTML report."""
        raise NotImplementedError("HTML export implemented in Week 3.")
