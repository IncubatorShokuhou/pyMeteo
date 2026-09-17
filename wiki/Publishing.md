# Publishing and contributing

[中文](Publishing_zh)

```bash
pip install pymeteo-kit
```

Import remains `import pymeteo`. Python 3.7+. NumPy only at runtime.

Distribution name: [`pymeteo-kit`](https://pypi.org/project/pymeteo-kit/). A GitHub Release (or **Actions → Publish → Run workflow**) uploads to PyPI. Workflow: `.github/workflows/publish.yml`.

```bash
pip install -e ".[dev]"
pytest
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
