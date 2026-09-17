# pymeteo

Meteorological diagnostic functions with explicit string unit parameters. Runtime dependency: NumPy only. Import: `pymeteo`. PyPI name: `pymeteo-kit`.

Public functions take string unit kwargs. There is no Pint, no Sounding/Wind object facade, no MetPy, and no NCL runtime. It is not an I/O, plotting, or full CAPE/CIN sounding package.

Parameter tables, unit aliases, and the NCL name map are in [README_zh.md](README_zh.md).

## Install

Python 3.6+ and NumPy. On 3.6 the pin is NumPy 1.19.x (last line that fully supports 3.6).

```bash
pip install pymeteo-kit
```

The PyPI project is `pymeteo-kit` because `pymeteo` was already taken and `py-meteo` was rejected as too similar. The import is still `import pymeteo`.

A GitHub Release (or **Actions → Publish → Run workflow**) builds the sdist/wheel and uploads to PyPI from the `pypi` environment using the `PYPI_API_TOKEN` secret. Workflow: [`.github/workflows/publish.yml`](.github/workflows/publish.yml).

From a git checkout, current hatchling (`>=1.18` in `[build-system]`) needs Python 3.8+. On 3.6/3.7 install the wheel from PyPI; do not build the sdist with those interpreters.

```bash
pip install -e .
pip install -e ".[dev]"   # pytest, ruff
```

## Quick start

```python
import pymeteo as pm

dewpoint = pm.dewpoint_from_relative_humidity(
    18.0, 46.5, temperature_unit="C", humidity_unit="%"
)
si = pm.showalter_index(16.6, 0.6, -15.9, temperature_unit="C")
km = pm.earth_distance(39.9, 116.4, 31.2, 121.5, output_distance_unit="km")
```

Default units: temperature `C`, pressure `hPa`, wind `m/s`, relative humidity `%`, mixing ratio `kg/kg`, distance `km` (thickness defaults to `m`), angles in degrees. Aliases such as `K`/`kelvin`, `hPa`/`mb`, `kt`, `fraction`, `g/kg`, `m`/`nmi` are documented in [README_zh.md](README_zh.md).

Top-level `pymeteo` uses English snake_case names. NCL builtin names live in `pymeteo.ncl` and are not re-exported from the package root.

## Versions

- **v2.2.0** adds pointwise thermodynamic, stability, kinematic and NWS comfort algorithms (potential/equivalent potential temperature, LCL, wet-bulb, lifted index, Coriolis, hypsometric thickness, ω↔w, heat index, wind chill) plus NCL name shims. See the “新增” section in [README_zh.md](README_zh.md). No MetPy, NCL, Pint, or extra runtime dependencies — NumPy only, string unit kwargs.
- **v2.2.1** adds `tests/test_ncl_official_examples.py`, which pins expected numbers taken from NCL documentation examples.
- **v2.2.2** attempted to publish as `py-meteo` (`pymeteo` was taken); PyPI rejected that name as too similar to existing `pymeteo`.
- **v2.2.3** publishes the package on PyPI as `pymeteo-kit`. The import remains `import pymeteo`.
- **v2.3.0** lowers supported Python to 3.6+.

## NCL names (`pymeteo.ncl`)

Optional thin wrappers with NCL builtin names and fixed NCL units (K, Pa/hPa, %, kg/kg). They only translate arguments and call the modern functions. These names are not re-exported from top-level `pymeteo`.

```python
import pymeteo as pm
from pymeteo.ncl import dewtemp_trh, relhum_ttd, mixhum_ptrh, wind_speed

td_k = dewtemp_trh(18.0 + 273.15, 46.5)   # tk in K, rh in % → dewpoint K
td_k = pm.ncl.dewtemp_trh(18.0 + 273.15, 46.5)
rh = pm.ncl.relhum_ttd(291.15, 279.45, 0)  # opt=0 → %, opt=1 → fraction
```

Full mapping (NCL name → modern function and units) is in [README_zh.md](README_zh.md).

## Known limits

- GitHub-hosted `ubuntu-20.04` and `actions/setup-python` 3.6 were removed in 2025. CI still runs 3.6 tests in the `python:3.6.15-buster` container against a wheel built on 3.12. That job can break if the image disappears or if a future hatchling default raises core metadata above 2.1 (pip 21.3 is the last pip on 3.6). Wheels in this repo set `core-metadata-version = "2.1"` for that reason.
- Building this tree (editable install or `python -m build`) needs Python 3.8+ because `[build-system]` requires `hatchling>=1.18`. 3.6/3.7 should install the `py3-none-any` wheel.
- 3.6 is capped at NumPy 1.19.x. Newer NumPy dropped 3.6.
- ruff’s oldest `target-version` is `py37`; lint runs on Python 3.12 only.
- No Pint. Units are strings. Scope limits (no FAO56 suite, no grid advection, no full CAPE/CIN) are listed in [README_zh.md](README_zh.md).

## License

GPL-3.0. See [LICENSE](LICENSE). Bug reports: [GitHub Issues](https://github.com/IncubatorShokuhou/pyMeteo/issues).
