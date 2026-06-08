# FAQ

## Which clock should I use?

It depends on your data and research question:

| Situation | Recommended clock |
|-----------|------------------|
| You have all 9 standard blood biomarkers | **PhenoAge** — best validated, highest mortality prediction |
| You want a relative pace score (not absolute age) | **DunedinPACEProxy** — captures rate, not level |
| You want to fit reference parameters to your own cohort | **KDM** with `clock.fit(your_df)` |
| You want a simple comparison baseline | **KDM** with default NHANES params |

For most research purposes, **start with PhenoAge** — it has the strongest mortality prediction (HR per SD ≈ 1.83 on NHANES) and is the most widely cited blood-based aging clock.

## What units should my biomarkers be in?

Always use **standard clinical / NHANES units**. AgingClockBench converts internally before applying coefficients.

| Biomarker | Input unit | Notes |
|-----------|-----------|-------|
| Albumin | g/dL | Converted to g/L internally |
| Creatinine | mg/dL | Converted to μmol/L internally |
| Glucose | mg/dL | Converted to mmol/L internally |
| CRP | **mg/L** | NOT mg/dL — already in mg/L |
| Lymphocytes | % | No conversion |
| MCV | fL | No conversion |
| RDW | % | No conversion |
| ALP | U/L | No conversion |
| WBC | 10³/μL | No conversion |

!!! warning "CRP units"
    CRP is **mg/L**, not mg/dL. NHANES reports CRP in mg/dL — multiply by 10 before passing to AgingClockBench.

## Why is my PhenoAge negative?

Negative biological ages are physiologically implausible but can occur with extreme biomarker values. Check:

1. CRP units — if CRP is in mg/dL instead of mg/L, PhenoAge will be severely underestimated
2. Very low glucose or RDW values

## Why doesn't KDM include CRP?

The Klemera-Doubal method works best with biomarkers that have a roughly linear relationship with age and consistent residual variance. CRP has high within-person variability and a skewed distribution, making it a poor fit for the KDM regression framework. PhenoAge handles CRP via a log transformation.

## Can I use AgingClockBench for clinical decisions?

**No.** AgingClockBench is a research tool for comparative benchmarking of aging algorithms. It is not validated for clinical use, medical diagnosis, or treatment decisions.

## How do I add a new clock?

Implement the `BaseClock` interface:

```python
from agingclockbench.clocks.base import BaseClock, ClockResult
import pandas as pd

class MyClock(BaseClock):
    @property
    def required_columns(self) -> list[str]:
        return ["age", "my_biomarker"]

    def validate_input(self, df: pd.DataFrame) -> tuple[bool, list[str]]:
        errors = self._check_required_columns(df)
        return len(errors) == 0, errors

    def transform(self, df: pd.DataFrame) -> ClockResult:
        complete = df[self.required_columns].dropna()
        biological_ages = complete["age"] + 0.5 * complete["my_biomarker"]
        accel = biological_ages - complete["age"]
        return ClockResult(
            clock_name="MyClock",
            biological_ages=biological_ages.reset_index(drop=True),
            accel=accel.reset_index(drop=True),
            missing_data_pct=(len(df) - len(complete)) / len(df) * 100,
            input_rows=len(df),
            output_rows=len(complete),
            original_index=complete.index,
        )
```

## How do I cite AgingClockBench?

```bibtex
@software{geddam2026agingclockbench,
  author = {Geddam, Aaditya},
  title  = {AgingClockBench: Benchmarking biological aging clocks},
  url    = {https://github.com/aadityageddam-ux/aging_clock_bench},
  year   = {2026}
}
```
