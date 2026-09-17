# Install

[中文](Install_zh)

## Runtime

Python **3.6+**. The only runtime dependency is NumPy. On CPython 3.6 the NumPy pin is **1.19.x** (last line that fully supports 3.6).

```bash
pip install pymeteo-kit
```

The PyPI project is `pymeteo-kit` because `pymeteo` was already taken and `py-meteo` was rejected as too similar. The import is still:

```python
import pymeteo as pm
```

## From a git checkout

Current `[build-system]` requires `hatchling>=1.18`, which needs **Python 3.8+**. Do not build the sdist with 3.6/3.7; those interpreters should install the `py3-none-any` wheel from PyPI.

```bash
pip install -e .
pip install -e ".[dev]"   # pytest; ruff on 3.8+
pytest
ruff check src tests
```

Lint in CI runs on Python 3.12 only. ruff’s oldest `target-version` in this tree is `py37`.

## CI note (3.6)

GitHub-hosted `ubuntu-20.04` and `actions/setup-python` 3.6 were removed in 2025. CI still runs 3.6 tests in the `python:3.6.15-buster` container against a wheel built on 3.12. Wheels set `core-metadata-version = "2.1"` so pip 21.3 (last pip on 3.6) can read the metadata.

## Publishing

A GitHub Release, or **Actions → Publish → Run workflow**, builds the sdist/wheel and uploads to PyPI from the `pypi` environment using the `PYPI_API_TOKEN` secret. Workflow: `.github/workflows/publish.yml`. See [Publishing](Publishing).
