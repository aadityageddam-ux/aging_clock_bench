# PhenoAge — Levine 2018

## Reference

> Levine ME, et al. *An epigenetic biomarker of aging for lifespan and healthspan.*
> Aging Cell. 2018;17(4):e12759. PMID: 29987380

## Overview

PhenoAge is a biological age clock derived from 9 routine blood biomarkers measured in NHANES. It predicts **10-year all-cause mortality** using a Gompertz proportional hazards model, then converts that mortality probability into a biological age estimate on the chronological age scale.

## Mathematical formula

### Step 1 — Linear predictor (xb)

$$
xb = -19.907
  + (-0.0336 \times \text{albumin}_{g/L})
  + (0.0095 \times \text{creatinine}_{\mu mol/L})
  + (0.1953 \times \text{glucose}_{mmol/L})
  + (0.0954 \times \ln(\text{CRP}_{mg/L}))
  + (-0.0120 \times \text{lymphocyte\%})
  + (0.0268 \times \text{MCV}_{fL})
  + (0.3306 \times \text{RDW\%})
  + (0.00188 \times \text{ALP}_{U/L})
  + (0.0554 \times \text{WBC}_{10^3/\mu L})
  + (0.0804 \times \text{age})
$$

!!! note "Unit conversions (applied internally)"
    Inputs are accepted in standard NHANES clinical units. AgingClockBench converts:

    - Albumin: g/dL → g/L (×10)
    - Creatinine: mg/dL → μmol/L (×88.4)
    - Glucose: mg/dL → mmol/L (×0.0555)
    - CRP: mg/L input, apply ln(crp + 0.001)

### Step 2 — 10-year mortality score

$$
M = 1 - \exp\!\left(\frac{-e^{xb} \cdot (e^{\gamma t} - 1)}{\gamma}\right)
$$

where $\gamma = 0.0076927$ and $t = 120$ months (10 years).

### Step 3 — Phenotypic age

$$
\text{PhenoAge} = 141.50 + \frac{\ln(-0.00553 \cdot \ln(1 - M))}{0.090165}
$$

## Validation (NHANES 1999-2000, N=4,086)

| Metric | Value |
|--------|-------|
| Pearson r (PhenoAge vs age) | 0.93 |
| Mean PhenoAge | 45.1 yr |
| Mean age | 49.8 yr |
| Mortality HR per SD acceleration | 1.83 (p < 0.001) |

## Usage

```python
from agingclockbench import PhenoAge
from agingclockbench.datasets import load_nhanes_sample

df = load_nhanes_sample()
result = PhenoAge().transform(df)

print(f"Mean PhenoAge: {result.biological_ages.mean():.1f} yr")
print(f"Mean acceleration: {result.accel.mean():.1f} yr")
```
