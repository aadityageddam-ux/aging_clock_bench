# AgingClockBench

[![Tests](https://github.com/aadityageddam-ux/aging_clock_bench/actions/workflows/test.yml/badge.svg)](https://github.com/aadityageddam-ux/aging_clock_bench/actions/workflows/test.yml)
[![PyPI version](https://badge.fury.io/py/agingclockbench.svg)](https://badge.fury.io/py/agingclockbench)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Benchmark biological aging clocks on your data in minutes.**

Multiple biological aging clocks exist — PhenoAge, KDM, DunedinPACE — but no standard tool lets researchers compare them side-by-side. AgingClockBench is the first open-source Python package that implements multiple clocks with unified benchmarking and mortality validation.

## Quick Start

```bash
pip install agingclockbench
```

```python
from agingclockbench import PhenoAge, KDM, BenchmarkSuite
from agingclockbench.datasets import load_nhanes_sample

df = load_nhanes_sample()

results = {
    "PhenoAge": PhenoAge().transform(df),
    "KDM": KDM().transform(df),
}

suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
report = suite.run(df, results)

print(report.to_dataframe())
# Clock     Pearson r  Mort HR  Mort p-value
# PhenoAge  0.82       1.83     < 0.001
# KDM       0.77       1.64     < 0.001
```

## CLI

```bash
agingclockbench benchmark \
  --data my_data.csv \
  --clocks PhenoAge KDM \
  --mortality-col vital_status \
  --followup-col followup_months \
  --output ./results/
```

## Implemented Clocks

| Clock | Reference | Biomarkers Required |
|-------|-----------|---------------------|
| PhenoAge | Levine et al. 2018 *Aging Cell* | 9 blood biomarkers |
| KDM | Klemera & Doubal 2006 *Mech Ageing Dev* | 4 blood biomarkers |
| DunedinPACEProxy | Proxy (NOT real DunedinPACE) | 5 blood biomarkers |

> **Note:** DunedinPACEProxy is a blood-biomarker approximation only. The real DunedinPACE requires DNA methylation data (Belsky 2022, eLife).

## Citation

```bibtex
@software{geddam2026agingclockbench,
  author = {Geddam, Aaditya},
  title = {AgingClockBench: Benchmarking biological aging clocks},
  url = {https://github.com/aadityageddam-ux/aging_clock_bench},
  year = {2026}
}
```

## License

MIT
