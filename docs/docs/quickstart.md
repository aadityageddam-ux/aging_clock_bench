# Quickstart

## Installation

```bash
pip install agingclockbench
```

Requires Python 3.11+.

---

## 5-minute tutorial

### Load the bundled NHANES data

```python
from agingclockbench.datasets import load_nhanes_sample

df = load_nhanes_sample()
print(f"{len(df)} participants, age {df.age.min():.0f}–{df.age.max():.0f} yr")
# 4086 participants, age 20–85 yr
```

### Compute biological ages

```python
from agingclockbench import PhenoAge, KDM

results = {
    "PhenoAge": PhenoAge().transform(df),
    "KDM":      KDM().transform(df),
}

for name, res in results.items():
    print(f"{name}: mean BA = {res.biological_ages.mean():.1f} yr  "
          f"mean accel = {res.accel.mean():.1f} yr")
```

### Run the benchmark

```python
from agingclockbench import BenchmarkSuite

suite = BenchmarkSuite(mortality_col="mortstat", followup_col="permth_exm")
report = suite.run(df, results)

print(report.to_dataframe())
```

Output:

```
Clock     Pearson r  Mort HR (per SD accel)  Mort p-value
PhenoAge      0.930                    1.83      < 0.001
KDM           0.677                    1.41      < 0.001
```

### Visualize

```python
report.plot_comparison()     # biological age vs chronological age
report.plot_km_survival()    # Kaplan-Meier by acceleration quartile
report.to_html("report.html")  # interactive Plotly report
```

---

## Using your own data

Your CSV needs the columns listed below. Mortality columns are optional.

```python
import pandas as pd
from agingclockbench import PhenoAge, BenchmarkSuite

df = pd.read_csv("my_cohort.csv")
result = PhenoAge().transform(df)

suite = BenchmarkSuite(
    mortality_col="vital_status",    # your column names
    followup_col="followup_months",
)
report = suite.run(df, {"PhenoAge": result})
print(report.to_dataframe())
```

### Required columns per clock

| Column | Unit | PhenoAge | KDM | DunedinPACEProxy |
|--------|------|:--------:|:---:|:----------------:|
| `age` | years | ✓ | ✓ | ✓ |
| `albumin_g_dl` | g/dL | ✓ | ✓ | ✓ |
| `creatinine_mg_dl` | mg/dL | ✓ | ✓ | ✓ |
| `glucose_mg_dl` | mg/dL | ✓ | ✓ | ✓ |
| `crp_mg_l` | mg/L | ✓ | | |
| `lymphocyte_pct` | % | ✓ | ✓ | ✓ |
| `mcv_fl` | fL | ✓ | ✓ | ✓ |
| `rdw_pct` | % | ✓ | ✓ | ✓ |
| `alp_u_l` | U/L | ✓ | ✓ | |
| `wbc_k_ul` | 10³/μL | ✓ | ✓ | ✓ |

---

## CLI

```bash
# Use bundled NHANES sample
agingclockbench benchmark --data bundled --clocks all

# Your own CSV, with HTML report
agingclockbench benchmark \
  --data my_data.csv \
  --clocks PhenoAge KDM \
  --mortality-col vital_status \
  --followup-col followup_months \
  --report

# List bundled datasets
agingclockbench datasets list
agingclockbench datasets info
```
