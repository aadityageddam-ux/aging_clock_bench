# DunedinPACEProxy

!!! danger "This is NOT the real DunedinPACE"
    The real DunedinPACE (Belsky 2022) requires **DNA methylation data** from the
    Illumina EPIC array. This is a blood-biomarker **proxy** for benchmarking.
    
    Use for **relative comparison only**, not for absolute pace-of-aging estimates.

## Reference (real DunedinPACE)

> Belsky DW, et al. *DunedinPACE, a DNA methylation biomarker of the pace of aging.*
> eLife. 2022;11:e73420. PMID: 35580151

## Overview

DunedinPACEProxy approximates the **pace-of-aging** concept from DunedinPACE using standard blood biomarkers. Rather than estimating a biological *age level*, it estimates a biological *aging rate*:

- **Pace = 1.0** → aging at the population average rate for your age
- **Pace > 1.0** → aging faster than expected
- **Pace < 1.0** → aging slower than expected

## Algorithm

For each biomarker $j$ with reference regression parameters $(k_j, q_j, s_j)$:

1. Compute **age-standardised residual**: how far the biomarker deviates from the age-expected value
2. Apply a **sign** (+1 if higher = faster aging, −1 if lower = faster aging)
3. Average across all biomarkers
4. Normalise to mean = 1.0, SD ≈ 0.1

$$
\text{raw} = \frac{1}{m} \sum_j \text{sign}_j \cdot \frac{x_j - (k_j \cdot t + q_j)}{s_j}
$$

$$
\text{pace} = 1.0 + \frac{\text{raw}}{s_{\text{raw}}} \times 0.1
$$

## Validation

| Property | Value |
|----------|-------|
| Correlation with PhenoAge acceleration | r = 0.84 |
| Correlation with chronological age | r ≈ 0.00 (by design) |
| Mean pace score | 1.000 |
| SD of pace score | 0.100 |

The near-zero age correlation is by design — the proxy captures *rate* of aging, not *level*.

## Limitations

- Correlation with true DunedinPACE (methylation) is expected to be ~0.3–0.5
- Blood biomarkers cannot capture the epigenetic dynamics measured by EPIC array
- Do not use for clinical decisions or publication claims about DunedinPACE

## Usage

```python
from agingclockbench import DunedinPACEProxy

result = DunedinPACEProxy().transform(df)
# result.biological_ages = age × pace_score
# result.accel = age × (pace_score - 1) = "excess aging"
```
