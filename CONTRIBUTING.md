# Contributing to AgingClockBench

Thank you for your interest! Contributions of all kinds are welcome.

## Setup

```bash
git clone https://github.com/aadityageddam-ux/aging_clock_bench.git
cd aging_clock_bench
poetry install
poetry run pytest tests/ -v
```

## Adding a new clock

1. Create `src/agingclockbench/clocks/myclock.py`
2. Implement `BaseClock` — see [FAQ](https://aadityageddam-ux.github.io/aging_clock_bench/faq/)
3. Add to `src/agingclockbench/clocks/__init__.py` and `src/agingclockbench/__init__.py`
4. Write tests in `tests/test_myclock.py` (validate against published examples)
5. Add algorithm explainer to `docs/docs/algorithms/myclock.md`

## Running checks

```bash
make check   # format + lint + type-check + tests
make test    # tests only
make docs    # preview docs locally
```

## Pull request guidelines

- All tests must pass (`make test`)
- Coverage must not decrease below 85%
- New clocks must include at least 5 tests with biological plausibility checks
- Algorithm explainer doc required for new clocks

## Reporting bugs

Open a [GitHub issue](https://github.com/aadityageddam-ux/aging_clock_bench/issues) with:
- Your Python version and OS
- Minimal reproducible example
- Expected vs actual output
