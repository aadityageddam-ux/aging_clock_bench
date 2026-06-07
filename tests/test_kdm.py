"""Tests for KDM clock."""

import numpy as np
import pandas as pd
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
    # Uses NHANES reference params — results must be identical across calls
    r1 = KDM().transform(kdm_df)
    r2 = KDM().transform(kdm_df)
    np.testing.assert_array_almost_equal(r1.biological_ages.values, r2.biological_ages.values)


def test_kdm_missing_column_raises(kdm_df):
    bad_df = kdm_df.drop(columns=["rdw_pct"])
    with pytest.raises(ValueError, match="Missing required column"):
        KDM().transform(bad_df)


def test_kdm_biological_ages_plausible(kdm_df):
    result = KDM().transform(kdm_df)
    assert (result.biological_ages > 0).all()
    assert (result.biological_ages < 150).all()


def test_kdm_custom_fit():
    """Fitting on own data (>=30 rows) and transforming should produce finite ages."""
    rows = [
        dict(age=30 + i * 1.5, albumin_g_dl=4.6 - i * 0.02, creatinine_mg_dl=0.6 + i * 0.01,
             glucose_mg_dl=82 + i * 0.8, rdw_pct=12.3 + i * 0.05, mcv_fl=95 + i * 0.1,
             wbc_k_ul=5.4 + i * 0.05, alp_u_l=76 + i * 0.5, lymphocyte_pct=27 - i * 0.1)
        for i in range(35)
    ]
    df = pd.DataFrame(rows)
    clock = KDM()
    clock.fit(df)
    result = clock.transform(df)
    assert np.isfinite(result.biological_ages.values).all()


def test_kdm_original_index_populated(kdm_df):
    result = KDM().transform(kdm_df)
    assert result.original_index is not None
    assert len(result.original_index) == result.output_rows


def test_kdm_nan_rows_dropped():
    row_complete = dict(age=53, albumin_g_dl=4.1, creatinine_mg_dl=0.5, glucose_mg_dl=94,
                        rdw_pct=12.7, mcv_fl=87.8, wbc_k_ul=7.4, alp_u_l=98, lymphocyte_pct=35.8)
    row_incomplete = dict(age=53, albumin_g_dl=None, creatinine_mg_dl=0.5, glucose_mg_dl=94,
                          rdw_pct=12.7, mcv_fl=87.8, wbc_k_ul=7.4, alp_u_l=98, lymphocyte_pct=35.8)
    df = pd.DataFrame([row_complete, row_incomplete])
    result = KDM().transform(df)
    assert result.output_rows == 1
    assert result.missing_data_pct == 50.0
