# Publishing and contributing

[中文](Publishing_zh)

## Install for users

```bash
pip install pymeteo-kit
```

Import remains `import pymeteo`. Python 3.6+. NumPy only at runtime.

## PyPI

Distribution name: [`pymeteo-kit`](https://pypi.org/project/pymeteo-kit/). Publishing is GitHub Actions:

* Trigger: a GitHub **Release**, or **Actions → Publish → Run workflow**
* Workflow: `.github/workflows/publish.yml`
* Environment: `pypi` (URL https://pypi.org/p/pymeteo-kit)
* Secret: `PYPI_API_TOKEN`
* Build: Python 3.12, `python -m build` (hatchling ≥ 1.18)

Wheels/sdists set `core-metadata-version = "2.1"` so pip 21.3 on CPython 3.6 can still read metadata.

## CI

`.github/workflows/ci.yml`: ruff on 3.12; pytest on 3.8 / 3.10 / 3.12 / 3.13; a 3.6 job installs the 3.12-built wheel into `python:3.6.15-buster`.

```bash
pip install -e ".[dev]"
pytest
ruff check src tests
```

## Contributing (brief)

* Open issues and pull requests against `master` on https://github.com/IncubatorShokuhou/pyMeteo
* Keep the public API English snake_case at the package root. NCL names belong in `pymeteo.ncl` only; do not invent NCL names that NCL does not have (heat index, wind chill, Showalter, K, SWEAT, …).
* Units are strings. Do not add Pint, MetPy, or extra runtime dependencies.
* Pointwise diagnostics only: no I/O, plotting, map projections, FAO56 suite, grid advection, or a full CAPE/CIN sounding kit unless that is an explicit, reviewed expansion.
* Formulas should be reimplemented from public literature, not pasted from MetPy or NCL sources.
* Tests should cover the function you touch. NCL documentation example numbers belong in `tests/test_ncl_official_examples.py`.

## License

MIT. Copyright IncubatorShokuhou, 2019–2026. See `LICENSE` in the source tree.
