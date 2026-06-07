"""Abstract base class for all aging clocks."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd


@dataclass
class ClockResult:
    """Output from a clock's .transform() call.

    Attributes
    ----------
    original_index : pd.Index of the rows in the original input DataFrame that
        were actually processed (i.e., had complete data). Used by BenchmarkSuite
        to align mortality/survival data with clock outputs.
    """

    clock_name: str
    biological_ages: pd.Series
    accel: pd.Series
    missing_data_pct: float
    input_rows: int
    output_rows: int
    original_index: Optional[pd.Index] = None
    metadata: dict = field(default_factory=dict)


class BaseClock(ABC):
    """Abstract interface all aging clocks must implement."""

    @property
    @abstractmethod
    def required_columns(self) -> list[str]:
        """Column names required in the input DataFrame."""

    @abstractmethod
    def validate_input(self, df: pd.DataFrame) -> tuple[bool, list[str]]:
        """Return (is_valid, list_of_error_messages)."""

    @abstractmethod
    def transform(self, df: pd.DataFrame) -> ClockResult:
        """Compute biological ages. Returns a ClockResult."""

    def _check_required_columns(self, df: pd.DataFrame) -> list[str]:
        missing = [c for c in self.required_columns if c not in df.columns]
        return [f"Missing required column: '{c}'" for c in missing]
