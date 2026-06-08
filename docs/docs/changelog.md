# Changelog

## 0.1.0 (2026-06-08)

Initial release.

### Clocks

- **PhenoAge** (Levine 2018) — 9-biomarker Gompertz mortality model; validated against NHANES 1999-2000 reference implementation (zero numerical difference)
- **KDM** (Klemera & Doubal 2006) — full two-stage formula with chronological age anchor; NHANES reference parameters built-in; `fit()` for custom cohorts
- **DunedinPACEProxy** — age-standardised blood-biomarker proxy; mean=1.0, SD=0.1; corr with PhenoAge acceleration r=0.84

### Benchmarking

- `BenchmarkSuite` with Cox PH (age-adjusted), Pearson/Spearman correlation, coefficient of variation, inter-clock agreement
- Correct row alignment via `ClockResult.original_index` (handles NaN-dropped rows)

### Visualization

- `plot_comparison()` — scatter: biological age vs chronological age
- `plot_km_survival()` — Kaplan-Meier by acceleration quartile
- `plot_correlation_heatmap()` — inter-clock acceleration correlations
- `to_html()` — interactive Plotly HTML report

### Data

- Bundled NHANES 1999-2000 sample (N=4,086, mortality-linked)
- `load_nhanes_sample()` function

### CLI

- `agingclockbench benchmark` — full benchmark with optional HTML report
- `agingclockbench datasets list/info` — dataset management

### Tests

- 89 tests; 89% code coverage
