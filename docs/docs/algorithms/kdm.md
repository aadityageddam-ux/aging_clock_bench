# KDM — Klemera & Doubal 2006

## Reference

> Klemera P, Doubal S. *A new approach to the concept and computation of biological age.*
> Mech Ageing Dev. 2006;127(3):240-248. PMID: 16318874

## Overview

The Klemera-Doubal Method (KDM) estimates biological age as a **maximum-likelihood weighted average** of how far each biomarker deviates from its age-expected value. It is more statistically principled than simple weighted scoring and naturally incorporates chronological age as a final "anchor" measurement.

## Algorithm

### Step 1 — Reference regression

For each biomarker $j$, regress on chronological age in a reference cohort:

$$x_j = q_j + k_j \cdot t + \varepsilon_j, \quad \varepsilon_j \sim \mathcal{N}(0, s_j^2)$$

Default parameters are from NHANES 1999-2000 (N=4,086). Override with `KDM().fit(your_df)`.

### Step 2 — Preliminary biological age (BA₁)

$$
BA_1 = \frac{\displaystyle\sum_j \frac{k_j (x_j - q_j)}{s_j^2}}{\displaystyle\sum_j \frac{k_j^2}{s_j^2}}
$$

### Step 3 — Estimate s_BA

$$s_{BA} = \text{SD}(BA_1 - t_{\text{reference}})$$

### Step 4 — Final KDM biological age

Incorporate chronological age as an additional measurement with precision $1/s_{BA}^2$:

$$
BA_{\text{KDM}} = \frac{\displaystyle\sum_j \frac{k_j (x_j - q_j)}{s_j^2} + \frac{t}{s_{BA}^2}}{\displaystyle\sum_j \frac{k_j^2}{s_j^2} + \frac{1}{s_{BA}^2}}
$$

## NHANES reference parameters

| Biomarker | k (slope) | q (intercept) | s (residual SD) |
|-----------|-----------|---------------|-----------------|
| albumin_g_dl | −0.002775 | 4.557 | 0.343 |
| creatinine_mg_dl | +0.005381 | 0.487 | 0.568 |
| glucose_mg_dl | +0.467131 | 74.83 | 36.02 |
| rdw_pct | +0.012186 | 12.14 | 1.083 |
| mcv_fl | +0.049082 | 87.91 | 5.174 |
| wbc_k_ul | −0.012827 | 7.945 | 2.135 |
| alp_u_l | +0.190120 | 74.74 | 32.15 |
| lymphocyte_pct | −0.020487 | 30.60 | 8.550 |

s_BA = 40.45 years (NHANES 1999-2000 reference)

## Usage

```python
from agingclockbench import KDM
from agingclockbench.datasets import load_nhanes_sample

df = load_nhanes_sample()
result = KDM().transform(df)   # uses NHANES reference params

# Or fit on your own cohort first:
clock = KDM().fit(df)
result = clock.transform(new_df)
```
