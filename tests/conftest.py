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
        dict(age=30, albumin_g_dl=4.6, creatinine_mg_dl=0.6, glucose_mg_dl=82,
             rdw_pct=12.3, mcv_fl=95.9, wbc_k_ul=5.4, alp_u_l=76, lymphocyte_pct=27),
        dict(age=53, albumin_g_dl=4.1, creatinine_mg_dl=0.5, glucose_mg_dl=94,
             rdw_pct=12.7, mcv_fl=87.8, wbc_k_ul=7.4, alp_u_l=98, lymphocyte_pct=35.8),
        dict(age=80, albumin_g_dl=4.2, creatinine_mg_dl=0.5, glucose_mg_dl=98,
             rdw_pct=12.7, mcv_fl=89.9, wbc_k_ul=8.5, alp_u_l=88, lymphocyte_pct=20.9),
    ]
    return pd.DataFrame(rows * 5)


@pytest.fixture
def dunedin_df() -> pd.DataFrame:
    rows = [
        dict(age=45, albumin_g_dl=4.5, creatinine_mg_dl=0.85, glucose_mg_dl=90,
             rdw_pct=12.8, wbc_k_ul=5.5, lymphocyte_pct=30, mcv_fl=89),
        dict(age=60, albumin_g_dl=4.0, creatinine_mg_dl=1.0, glucose_mg_dl=100,
             rdw_pct=13.5, wbc_k_ul=6.5, lymphocyte_pct=25, mcv_fl=91),
        dict(age=72, albumin_g_dl=3.8, creatinine_mg_dl=1.1, glucose_mg_dl=115,
             rdw_pct=14.2, wbc_k_ul=7.5, lymphocyte_pct=20, mcv_fl=93),
    ]
    return pd.DataFrame(rows * 5)
