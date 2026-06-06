"""Tests for PhenoAge clock."""

import numpy as np
import pandas as pd
import pytest

from agingclockbench.clocks.phenoage import PhenoAge
from agingclockbench.clocks.base import ClockResult


def test_phenoage_returns_clock_result(phenoage_df):
    result = PhenoAge().transform(phenoage_df)
    assert isinstance(result, ClockResult)
    assert result.clock_name == "PhenoAge"


def test_phenoage_output_length(phenoage_df):
    result = PhenoAge().transform(phenoage_df)
    assert len(result.biological_ages) == len(phenoage_df)


def test_phenoage_biological_ages_plausible(phenoage_df):
    result = PhenoAge().transform(phenoage_df)
    # All biological ages should be within a plausible human lifespan range
    assert (result.biological_ages > 10).all()
    assert (result.biological_ages < 130).all()


def test_phenoage_acceleration_defined(phenoage_df):
    result = PhenoAge().transform(phenoage_df)
    assert len(result.accel) == len(result.biological_ages)
    # accel = biological_age - chronological_age
    expected_accel = result.biological_ages - phenoage_df["age"].values
    np.testing.assert_allclose(result.accel.values, expected_accel.values, atol=0.01)


def test_phenoage_missing_column_raises():
    bad_df = pd.DataFrame([{"age": 52, "albumin_g_dl": 4.3}])
    with pytest.raises(ValueError, match="Missing required column"):
        PhenoAge().transform(bad_df)


def test_phenoage_nan_rows_dropped():
    # A single row with a missing biomarker: all rows drop → should raise
    row = dict(age=52, albumin_g_dl=None, creatinine_mg_dl=0.9,
               glucose_mg_dl=87.0, crp_mg_l=0.3, lymphocyte_pct=28.0,
               mcv_fl=90.0, rdw_pct=13.0, alp_u_l=65.0, wbc_k_ul=6.0)
    df = pd.DataFrame([row])
    with pytest.raises(ValueError, match="No complete rows"):
        PhenoAge().transform(df)


def test_phenoage_partial_nan_rows_dropped():
    # Mix of complete and incomplete rows: incomplete rows should be dropped
    complete = dict(age=52, albumin_g_dl=4.3, creatinine_mg_dl=0.9,
                    glucose_mg_dl=87.0, crp_mg_l=0.3, lymphocyte_pct=28.0,
                    mcv_fl=90.0, rdw_pct=13.0, alp_u_l=65.0, wbc_k_ul=6.0)
    incomplete = dict(age=52, albumin_g_dl=None, creatinine_mg_dl=0.9,
                      glucose_mg_dl=87.0, crp_mg_l=0.3, lymphocyte_pct=28.0,
                      mcv_fl=90.0, rdw_pct=13.0, alp_u_l=65.0, wbc_k_ul=6.0)
    df = pd.DataFrame([complete, incomplete])
    result = PhenoAge().transform(df)
    assert result.output_rows == 1
    assert result.missing_data_pct == 50.0


def test_phenoage_reproducible(phenoage_df):
    r1 = PhenoAge().transform(phenoage_df)
    r2 = PhenoAge().transform(phenoage_df)
    np.testing.assert_array_equal(r1.biological_ages.values, r2.biological_ages.values)


def test_phenoage_missing_data_pct_zero_when_complete(phenoage_df):
    result = PhenoAge().transform(phenoage_df)
    assert result.missing_data_pct == 0.0


def test_phenoage_metadata_has_reference(phenoage_df):
    result = PhenoAge().transform(phenoage_df)
    assert "reference" in result.metadata
    assert "Levine" in result.metadata["reference"]
