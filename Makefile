.PHONY: install test lint format mypy check docs clean publish

install:
	poetry install

test:
	poetry run pytest tests/ -v --cov=src/agingclockbench --cov-report=term-missing

lint:
	poetry run ruff check src/ tests/

format:
	poetry run black src/ tests/

mypy:
	poetry run mypy src/

check: format lint mypy test

docs:
	poetry run mkdocs serve

docs-build:
	poetry run mkdocs build

clean:
	find . -type d -name __pycache__ -exec rm -r {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	rm -rf dist/ build/ *.egg-info .coverage htmlcov/

publish:
	poetry publish
