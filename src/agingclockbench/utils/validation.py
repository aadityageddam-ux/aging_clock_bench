"""Input validation helpers."""

import pandas as pd


def validate_dataframe(df: pd.DataFrame, required_cols: list[str]) -> None:
    """Raise ValueError if any required column is missing from df."""
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Input DataFrame is missing columns: {missing}")
