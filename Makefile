.PHONY: install test lint format typecheck clean

install:
	pip install -e packages/dyadic-core
	pip install -e packages/dyadic-math --no-deps

test:
	pytest packages/ -v --tb=short

lint:
	ruff check packages/

format:
	ruff format packages/

typecheck:
	mypy packages/

clean:
	rm -rf .pytest_cache .ruff_cache .hypothesis
	find packages -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null; true
	find packages -type d -name '*.egg-info' -exec rm -rf {} + 2>/dev/null; true
