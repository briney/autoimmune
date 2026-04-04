# autoimmune

## Quick Reference

```bash
# Install (editable, with dev extras)
pip install -e ".[dev]"

# Run tests
pytest

# Lint & format
ruff check src/ tests/
ruff format src/ tests/

# Type check
mypy src/
```

## Project Structure

```
autoimmune/
├── src/
│   └── autoimmune/       # importable package
├── tests/                # pytest tests mirror src/ layout
├── pyproject.toml        # single source of truth for config
└── CLAUDE.md             # this file
```

## Conventions

- **Python 3.11+** with modern syntax (`match`, `X | Y` unions, `type` aliases).
- **src layout** — always import from `autoimmune`, never relative to the repo root.
- **Ruff** for linting and formatting. All config in `pyproject.toml`.
- **pytest** for testing. Mark slow tests with `@pytest.mark.slow`.
- Type hints on all public function signatures.
- Google-style docstrings on public APIs.
- Max line length: 100 characters.
