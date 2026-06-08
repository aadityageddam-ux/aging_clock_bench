"""Tests for bundled NHANES dataset loader."""

import pandas as pd
import pytest

from agingclockbench.datasets import load_nhanes_sample


def test_load_nhanes_sample_shape():
    df = load_nhanes_sample()
    assert len(df) == 4086
    assert len(df.columns) == 14


def test_load_nhanes_sample_required_columns():
    df = load_nhanes_sample()
    required = [
        "age", "albumin_g_dl", "creatinine_mg_dl", "glucose_mg_dl",
        "crp_mg_l", "lymphocyte_pct", "mcv_fl", "rdw_pct", "alp_u_l",
        "wbc_k_ul", "mortstat", "permth_exm",
    ]
    for col in required:
        assert col in df.columns, f"Missing column: {col}"


def test_load_nhanes_sample_no_missing():
    df = load_nhanes_sample()
    assert df.isnull().sum().sum() == 0, "Bundled sample contains missing values"


def test_load_nhanes_sample_age_range():
    df = load_nhanes_sample()
    assert df.age.min() >= 20
    assert df.age.max() <= 90


def test_load_nhanes_sample_mortality_binary():
    df = load_nhanes_sample()
    assert set(df.mortstat.unique()).issubset({0, 1})
    assert df.mortstat.sum() > 0, "No deaths in mortality column"


def test_load_nhanes_sample_followup_positive():
    df = load_nhanes_sample()
    assert (df.permth_exm > 0).all()


def test_load_nhanes_sample_biomarker_ranges():
    df = load_nhanes_sample()
    assert df.albumin_g_dl.between(2.0, 6.0).all(), "Albumin out of range"
    assert df.crp_mg_l.ge(0).all(), "Negative CRP values"
    # Allow up to 700 mg/dL — NHANES includes uncontrolled diabetics (max ~561)
    assert df.glucose_mg_dl.ge(30).all(), "Glucose values below 30 mg/dL"
    assert df.glucose_mg_dl.le(700).all(), "Glucose values above 700 mg/dL"


def test_nhanes_works_with_phenoage():
    """Bundled NHANES must work seamlessly with PhenoAge.transform()."""
    from agingclockbench import PhenoAge
    df = load_nhanes_sample()
    result = PhenoAge().transform(df)
    assert result.output_rows == len(df)
    assert result.missing_data_pct == 0.0
