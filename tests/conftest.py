"""Shared pytest fixtures."""

import pandas as pd
import pytest


@pytest.fixture
def minimal_phenoage_row() -> dict:
    return dict(
        age=52,
        albumin_g_dl=4.3,
        creatinine_mg_dl=0.9,
        glucose_mg_dl=87.0,
        crp_mg_l=0.3,
        lymphocyte_pct=28.0,
        mcv_fl=90.0,
        rdw_pct=13.0,
        alp_u_l=65.0,
        wbc_k_ul=6.0,
    )


@pytest.fixture
def phenoage_df(minimal_phenoage_row) -> pd.DataFrame:
    return pd.DataFrame([minimal_phenoage_row] * 10)


@pytest.fixture
def kdm_df() -> pd.DataFrame:
    rows = [
        dict(age=45, albumin_g_dl=4.5, creatinine_mg_dl=0.85, glucose_mg_dl=90.0, hemoglobin_g_dl=14.2),
        dict(age=60, albumin_g_dl=4.0, creatinine_mg_dl=1.0, glucose_mg_dl=100.0, hemoglobin_g_dl=13.8),
        dict(age=75, albumin_g_dl=3.7, creatinine_mg_dl=1.2, glucose_mg_dl=115.0, hemoglobin_g_dl=13.0),
    ]
    return pd.DataFrame(rows * 5)


@pytest.fixture
def dunedin_df() -> pd.DataFrame:
    rows = [
        dict(age=45, albumin_g_dl=4.5, creatinine_mg_dl=0.85, glucose_mg_dl=90.0, wbc_k_ul=5.5, rdw_pct=12.8),
        dict(age=60, albumin_g_dl=4.0, creatinine_mg_dl=1.0, glucose_mg_dl=100.0, wbc_k_ul=6.5, rdw_pct=13.5),
    ]
    return pd.DataFrame(rows * 5)
