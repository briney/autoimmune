# Agent Instructions for autoimmune

## Architecture

This is a Python package using a `src/` layout with `hatchling` as the build backend.
Source code lives in `src/autoimmune/`, tests in `tests/`.

## Code Style

- Python 3.11+. Use modern typing syntax (`X | Y`, `type` aliases, `Self`).
- All tools configured in `pyproject.toml` — no separate config files.
- Ruff for linting and formatting. Run `ruff check` and `ruff format` before committing.
- Type hints required on all public function signatures.
- Google-style docstrings on public classes and functions.

## Testing

- Use pytest. Test files mirror source layout: `src/autoimmune/foo.py` -> `tests/test_foo.py`.
- Run `pytest` to verify changes. All tests must pass before committing.
- Mark slow tests with `@pytest.mark.slow`.

## What Not To Do

- Do not add `setup.py`, `setup.cfg`, or `requirements.txt`.
- Do not use `print()` in library code — use `logging`.
- Do not commit generated files, model weights, or large data.
- Do not add `# type: ignore` without a specific error code and justification.
- Do not use bare `except:` clauses.
