"""Tests for KDM clock."""

import numpy as np
import pytest

from agingclockbench.clocks.kdm import KDM
from agingclockbench.clocks.base import ClockResult


def test_kdm_returns_clock_result(kdm_df):
    result = KDM().transform(kdm_df)
    assert isinstance(result, ClockResult)
    assert result.clock_name == "KDM"


def test_kdm_output_length(kdm_df):
    result = KDM().transform(kdm_df)
    assert len(result.biological_ages) == len(kdm_df)


def test_kdm_reproducible(kdm_df):
    # KDM fits on same data both times — results must be identical
    r1 = KDM().transform(kdm_df)
    r2 = KDM().transform(kdm_df)
    np.testing.assert_array_almost_equal(r1.biological_ages.values, r2.biological_ages.values)


def test_kdm_missing_column_raises(kdm_df):
    import pandas as pd
    bad_df = kdm_df.drop(columns=["hemoglobin_g_dl"])
    with pytest.raises(ValueError, match="Missing required column"):
        KDM().transform(bad_df)
