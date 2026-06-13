# AgingClockBench

[![PyPI version](https://badge.fury.io/py/agingclockbench.svg)](https://badge.fury.io/py/agingclockbench)
[![Tests](https://github.com/aadityageddam-ux/aging_clock_bench/actions/workflows/test.yml/badge.svg)](https://github.com/aadityageddam-ux/aging_clock_bench/actions/workflows/test.yml)
[![Coverage](https://img.shields.io/badge/coverage-89%25-brightgreen)](https://github.com/aadityageddam-ux/aging_clock_bench)
[![Docs](https://img.shields.io/badge/docs-README-blue)](https://github.com/aadityageddam-ux/aging_clock_bench#quick-start)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/downloads/)

**Benchmark biological aging clocks on your data in minutes.**

Multiple biological aging clocks exist — PhenoAge, KDM, DunedinPACE — but no standard tool lets researchers compare them side-by-side. AgingClockBench is the **first open-source Python package** implementing multiple clocks with a unified interface and reproducible mortality-validated benchmarking.

📖 **[Documentation & Examples](#quick-start)**

---

## Install

```bash
pip install agingclockbench
```

Requires Python 3.11+.

---

## Quick Start

```python
from agingclockbench import PhenoAge, KDM, BenchmarkSuite
from agingclockbench.datasets import load_nhanes_sample

# Load bundled NHANES 1999-2000 (N=4,086, 20-year mortality follow-up)
df = load_nhanes_sample()

# Compute biological ages
results = {
    "PhenoAge": PhenoAge().transform(df),
    "KDM":      KDM().transform(df),
}

# Benchmark against mortality
suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
report = suite.run(df, results)

print(report.to_dataframe())
```

```
Clock     Pearson r  Mort HR (per SD accel)  Mort p-value
PhenoAge      0.930                    1.83      0.000001
KDM           0.677                    1.41      0.000001
```

```python
report.plot_km_survival()      # Kaplan-Meier by acceleration quartile
report.plot_comparison()       # biological age vs chronological age
report.to_html("report.html")  # interactive Plotly report
```

---

## CLI

```bash
# Benchmark on bundled NHANES with HTML report
agingclockbench benchmark --data bundled --clocks all --report

# Your own data
agingclockbench benchmark \
  --data my_cohort.csv \
  --clocks PhenoAge KDM \
  --mortality-col vital_status \
  --followup-col followup_months \
  --output ./results/
```

---

## Implemented Clocks

| Clock | Reference | Biomarkers | Key metric (NHANES) |
|-------|-----------|-----------|---------------------|
| **PhenoAge** | Levine et al. 2018 *Aging Cell* | 9 blood | Pearson r=0.93, HR=1.83 |
| **KDM** | Klemera & Doubal 2006 *Mech Ageing Dev* | 8 blood | Pearson r=0.68, HR=1.41 |
| **DunedinPACEProxy** | Proxy (NOT Belsky 2022) | 7 blood | pace corr w/ PhenoAge accel r=0.84 |

> **Note:** DunedinPACEProxy is a blood-biomarker approximation. The real DunedinPACE requires DNA methylation data (Illumina EPIC array).

---

## Features

- ✅ **Unified interface** — all clocks share the same `transform()` API
- ✅ **Validated** — PhenoAge implementation cross-validated against reference; zero numerical difference on N=4,086 NHANES participants  
- ✅ **Mortality benchmarking** — Cox PH hazard ratios, Kaplan-Meier curves
- ✅ **Bundled data** — NHANES 1999-2000 with 20-year mortality follow-up, ready to use
- ✅ **Interactive reports** — Plotly HTML with comparison plots and benchmark table
- ✅ **CLI tool** — `agingclockbench benchmark` works out of the box
- ✅ **89% test coverage** — 76+ tests, CI/CD on every push

---

## Documentation

| Resource | Link |
|----------|------|
| Full docs | [aadityageddam-ux.github.io/aging_clock_bench](https://aadityageddam-ux.github.io/aging_clock_bench/) |
| Quickstart | [docs/quickstart](https://aadityageddam-ux.github.io/aging_clock_bench/quickstart/) |
| PhenoAge algorithm | [docs/algorithms/phenoage](https://aadityageddam-ux.github.io/aging_clock_bench/algorithms/phenoage/) |
| FAQ | [docs/faq](https://aadityageddam-ux.github.io/aging_clock_bench/faq/) |
| Example notebooks | [examples/](examples/) |

---

## Citation

If you use AgingClockBench in your research, please cite:

```bibtex
@software{geddam2026agingclockbench,
  author = {Geddam, Aaditya},
  title  = {AgingClockBench: Benchmarking biological aging clocks},
  url    = {https://github.com/aadityageddam-ux/aging_clock_bench},
  year   = {2026}
}
```

**Also cite the underlying clock papers:**

- Levine ME, et al. *Aging Cell.* 2018. (PhenoAge)
- Klemera P, Doubal S. *Mech Ageing Dev.* 2006. (KDM)

---

## Contributing

Contributions welcome! See [CONTRIBUTING.md](CONTRIBUTING.md).

To add a new clock, implement the `BaseClock` interface — see the [FAQ](https://aadityageddam-ux.github.io/aging_clock_bench/faq/#how-do-i-add-a-new-clock).

---

## License

MIT © Aaditya Geddam
