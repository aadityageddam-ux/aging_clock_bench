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
    """Run a standardised validation benchmark on one or more aging clocks.

    Parameters
    ----------
    mortality_col : str
        Column name for vital status (1 = dead, 0 = censored).
    followup_col : str
        Column name for follow-up time. Units must be consistent — the HR
        interpretation assumes months if using NHANES permth_exm.

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

        Uses ``ClockResult.original_index`` to align clock outputs with
        the correct rows in ``df`` (handles missing-data row drops).

        Parameters
        ----------
        df : Original input DataFrame.
        results : Mapping of clock name → ClockResult.

        Returns
        -------
        BenchmarkReport
        """
        benchmark_results: list[BenchmarkResult] = []
        accel_series: dict[str, pd.Series] = {}

        for name, result in results.items():
            br = BenchmarkResult(clock_name=name)

            # Align df rows to the rows the clock actually processed
            if result.original_index is not None:
                aligned_df = df.loc[result.original_index].reset_index(drop=True)
            else:
                aligned_df = df.iloc[: result.output_rows].reset_index(drop=True)

            age = aligned_df["age"]

            # Pearson / Spearman correlation with chronological age
            if age.nunique() > 1:
                r, p = stats.pearsonr(age, result.biological_ages)
                br.pearson_r = round(float(r), 4)
                br.pearson_pvalue = round(float(p), 6)
                spr = stats.spearmanr(age, result.biological_ages).statistic
                br.spearman_r = round(float(spr), 4)

            # Coefficient of variation
            mean_ba = result.biological_ages.mean()
            std_ba = result.biological_ages.std()
            br.cv = round(float(std_ba / mean_ba), 4) if mean_ba != 0 else float("nan")

            # Cox PH mortality prediction
            if self.mortality_col in aligned_df.columns and self.followup_col in aligned_df.columns:
                br = self._run_cox(aligned_df, result, br)

            benchmark_results.append(br)
            accel_series[name] = result.accel

        # Inter-clock agreement (Pearson r of accelerations)
        for br in benchmark_results:
            others = {k: v for k, v in accel_series.items() if k != br.clock_name}
            for other_name, other_accel in others.items():
                min_len = min(len(accel_series[br.clock_name]), len(other_accel))
                if min_len > 2 and accel_series[br.clock_name].iloc[:min_len].nunique() > 1:
                    r, _ = stats.pearsonr(
                        accel_series[br.clock_name].iloc[:min_len],
                        other_accel.iloc[:min_len],
                    )
                    br.clock_agreement_with_others[other_name] = round(float(r), 4)

        return BenchmarkReport(
            results=benchmark_results,
            df=df,
            clock_results=results,
            mortality_col=self.mortality_col,
            followup_col=self.followup_col,
        )

    def _run_cox(
        self,
        aligned_df: pd.DataFrame,
        result: ClockResult,
        br: BenchmarkResult,
    ) -> BenchmarkResult:
        """Fit a Cox PH model: mortality ~ clock_acceleration_sd + age."""
        try:
            from lifelines import CoxPHFitter

            analysis_df = aligned_df[[self.mortality_col, self.followup_col, "age"]].copy()
            analysis_df["clock_acceleration"] = result.accel.values

            # Standardise acceleration to per-SD hazard ratio
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
            br.cox_nobs = int(analysis_df[self.mortality_col].sum())
        except Exception:
            pass
        return br


class BenchmarkReport:
    """Container for all benchmark results with display and export methods.

    Attributes
    ----------
    results : list[BenchmarkResult]
    _df : the original input DataFrame (set by BenchmarkSuite.run)
    _clock_results : dict mapping clock name -> ClockResult
    """

    def __init__(
        self,
        results: list[BenchmarkResult],
        df: pd.DataFrame | None = None,
        clock_results: dict | None = None,
        mortality_col: str = "mortstat",
        followup_col: str = "permth_exm",
    ) -> None:
        self.results = results
        self._df = df
        self._clock_results = clock_results or {}
        self._mortality_col = mortality_col
        self._followup_col = followup_col

    def to_dataframe(self) -> pd.DataFrame:
        """Return a summary DataFrame — one row per clock."""
        rows = []
        for r in self.results:
            rows.append({
                "Clock": r.clock_name,
                "Pearson r": r.pearson_r,
                "Spearman r": r.spearman_r,
                "Mort HR (per SD accel)": r.mortality_hr,
                "HR 95% CI lower": r.mortality_hr_ci_lower,
                "HR 95% CI upper": r.mortality_hr_ci_upper,
                "Mort p-value": r.mortality_pvalue,
                "Cox N (events)": r.cox_nobs,
                "CV": r.cv,
            })
        return pd.DataFrame(rows)

    def plot_comparison(self, df: pd.DataFrame | None = None,
                        results: dict | None = None):
        """Scatter plot of biological age vs chronological age per clock.

        Returns matplotlib Figure. Pass ``df`` and ``results`` only if you
        did not run via BenchmarkSuite.run().
        """
        from agingclockbench.benchmarks.plots import plot_comparison
        return plot_comparison(
            self,
            df if df is not None else self._df,
            results if results is not None else self._clock_results,
        )

    def plot_km_survival(self, df: pd.DataFrame | None = None,
                         results: dict | None = None,
                         n_quartiles: int = 4):
        """Kaplan-Meier survival by age-acceleration quartile.

        Returns matplotlib Figure.
        """
        from agingclockbench.benchmarks.plots import plot_km_survival
        return plot_km_survival(
            df if df is not None else self._df,
            results if results is not None else self._clock_results,
            mortality_col=self._mortality_col,
            followup_col=self._followup_col,
            n_quartiles=n_quartiles,
        )

    def plot_correlation_heatmap(self, results: dict | None = None):
        """Heatmap of Pearson correlations between clock accelerations.

        Returns matplotlib Figure.
        """
        from agingclockbench.benchmarks.plots import plot_correlation_heatmap
        return plot_correlation_heatmap(
            results if results is not None else self._clock_results
        )

    def to_html(self, filename: str, df: pd.DataFrame | None = None,
                results: dict | None = None) -> None:
        """Export an interactive Plotly HTML benchmark report.

        Parameters
        ----------
        filename : output path (e.g. 'report.html')
        """
        from agingclockbench.benchmarks.plots import to_html
        to_html(
            self,
            df if df is not None else self._df,
            results if results is not None else self._clock_results,
            filename,
            mortality_col=self._mortality_col,
            followup_col=self._followup_col,
        )
