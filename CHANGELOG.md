# Changelog

## [0.1.0] - 2026-06-11

### Added
- **PhenoAge** — Levine 2018 9-biomarker biological age clock (Gompertz mortality model)
  - Cross-validated to zero difference on 4,086 NHANES participants
  - Full unit conversions (albumin, creatinine, glucose, CRP)
  
- **KDM** — Klemera-Doubal 2006 two-stage maximum-likelihood estimator
  - NHANES 1999-2000 reference parameters included
  - `fit()` method to derive custom parameters on user data
  
- **DunedinPACEProxy** — Age-standardized biomarker proxy for pace-of-aging
  - NOT the real epigenetic DunedinPACE (requires DNA methylation)
  - Mean=1.0, SD=0.1 scale; r=0.84 with PhenoAge acceleration
  
- **BenchmarkSuite** — Mortality-validated benchmarking
  - Cox proportional hazards (Hazard Ratio + CI + p-value)
  - Kaplan-Meier survival curves (by acceleration quartile)
  - Pearson/Spearman correlation, coefficient of variation
  
- **Bundled NHANES 1999-2000** — N=4,086 complete cases with 20-year mortality follow-up
  - Ready to use: `from agingclockbench.datasets import load_nhanes_sample`
  
- **CLI tool** — `agingclockbench benchmark --data bundled --clocks all --report`
  - Pretty-printed output with benchmark table and per-clock summary
  
- **Visualization** — Matplotlib scatter, KM curves, heatmaps + Plotly interactive HTML reports

- **Full documentation** — MkDocs with Material theme
  - Quickstart, algorithm explainers (with LaTeX), FAQ, contributing guide
  - Auto-generated API reference via mkdocstrings
  - Auto-deployed to GitHub Pages on every push

### Test Coverage
- **96% overall** — 91/91 tests passing
- Unit tests for all 3 clocks with validation against published examples
- Mortality benchmark tests with real NHANES data
- Integration tests for CLI, plotting, datasets

### Infrastructure
- **GitHub Actions CI/CD**
  - Tests run on every push
  - Docs auto-deploy to GitHub Pages
  - Package auto-publishes to PyPI on version tags
- **PyPI packaging** — wheel (86 KB) with bundled NHANES parquet
- **MIT License**

---

**References:**
- Levine ME, et al. (2018). Aging Cell. PhenoAge algorithm.
- Klemera P, Doubal S. (2006). Mech Ageing Dev. KDM algorithm.
