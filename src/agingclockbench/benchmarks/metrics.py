"""Individual metric functions used by BenchmarkSuite."""

import numpy as np
import pandas as pd
from scipy import stats


def pearson_correlation(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    """Return (r, p-value) Pearson correlation between x and y."""
    r, p = stats.pearsonr(x, y)
    return float(r), float(p)


def spearman_correlation(x: pd.Series, y: pd.Series) -> float:
    """Return Spearman rho between x and y."""
    return float(stats.spearmanr(x, y).statistic)


def coefficient_of_variation(series: pd.Series) -> float:
    """Return coefficient of variation (SD / mean)."""
    mean = series.mean()
    return float(series.std() / mean) if mean != 0 else float("nan")
