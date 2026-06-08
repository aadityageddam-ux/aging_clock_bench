# AgingClockBench

**Benchmark biological aging clocks on your data in minutes.**

[![PyPI version](https://badge.fury.io/py/agingclockbench.svg)](https://badge.fury.io/py/agingclockbench)
[![Tests](https://github.com/aadityageddam-ux/aging_clock_bench/actions/workflows/test.yml/badge.svg)](https://github.com/aadityageddam-ux/aging_clock_bench/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## The problem

Multiple biological aging clocks exist — **PhenoAge** (Levine 2018), **KDM** (Klemera & Doubal 2006), **DunedinPACE** (Belsky 2022) — but no standard tool lets researchers compare them side-by-side on their own data.

Researchers either:

- Implement each algorithm from scratch (error-prone, weeks of work)
- Use one clock in isolation (miss comparative validation)
- Pay for proprietary services (expensive, not reproducible)

## The solution

AgingClockBench is the **first open-source Python package** that implements multiple aging clocks with a unified interface and reproducible benchmarking against mortality.

```bash
pip install agingclockbench
```

```python
from agingclockbench import PhenoAge, KDM, BenchmarkSuite
from agingclockbench.datasets import load_nhanes_sample

df = load_nhanes_sample()   # bundled NHANES 1999-2000 (N=4,086)

results = {
    "PhenoAge": PhenoAge().transform(df),
    "KDM":      KDM().transform(df),
}

suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
report = suite.run(df, results)

print(report.to_dataframe())
# Clock     Pearson r  Mort HR  Mort p-value
# PhenoAge  0.93       1.83     < 0.001
# KDM       0.68       1.41     < 0.001

report.plot_km_survival()   # Kaplan-Meier survival curves
report.to_html("report.html")  # interactive Plotly report
```

## CLI

```bash
agingclockbench benchmark --data bundled --clocks all --report
agingclockbench benchmark --data my_data.csv --clocks PhenoAge KDM \
    --mortality-col vital_status --followup-col followup_months
```

## Implemented clocks

| Clock | Reference | Biomarkers | Notes |
|-------|-----------|-----------|-------|
| **PhenoAge** | Levine et al. 2018 *Aging Cell* | 9 blood biomarkers | Validated against NHANES reference |
| **KDM** | Klemera & Doubal 2006 *Mech Ageing Dev* | 8 blood biomarkers | NHANES reference params built-in; `fit()` for custom cohorts |
| **DunedinPACEProxy** | Proxy (NOT Belsky 2022) | 7 blood biomarkers | Blood-biomarker approximation; pace score mean≈1.0 |

!!! warning "DunedinPACEProxy is not the real DunedinPACE"
    The real DunedinPACE requires DNA methylation data (Illumina EPIC array).
    This is a blood-biomarker proxy for benchmarking purposes only.

## NHANES benchmark results

Run on NHANES 1999-2000 (N=4,086, 20-year mortality follow-up):

| Clock | Pearson r (vs age) | Mortality HR (per SD accel) | Mort p-value |
|-------|-------------------|---------------------------|-------------|
| PhenoAge | 0.93 | 1.83 | < 0.001 |
| KDM | 0.68 | 1.41 | < 0.001 |

## Citation

```bibtex
@software{geddam2026agingclockbench,
  author = {Geddam, Aaditya},
  title  = {AgingClockBench: Benchmarking biological aging clocks},
  url    = {https://github.com/aadityageddam-ux/aging_clock_bench},
  year   = {2026}
}
```
