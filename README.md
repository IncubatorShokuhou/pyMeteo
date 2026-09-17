# pymeteo

Lightweight meteorological diagnostic functions for Python 3.10+.

This is a **breaking rewrite** of the old single-file `pyMeteo.py`. The import name is now lowercase `pymeteo`. Public functions take **explicit string unit parameters** (no Pint, no Sounding/Wind object facade). Runtime dependency: **NumPy only**.

**Full documentation is in Chinese:** [README_zh.md](README_zh.md)

## Install

```bash
pip install pymeteo-kit
```

The PyPI project name is `pymeteo-kit` because `pymeteo` was already taken and `py-meteo` was rejected as too similar. The import is still `import pymeteo`.

For local development:

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

**v2.2.0** adds pointwise thermodynamic, stability, kinematic and NWS comfort algorithms (potential/equivalent potential temperature, LCL, wet-bulb, lifted index, Coriolis, hypsometric thickness, ω↔w, heat index, wind chill) plus NCL name shims. See the “新增” section in [README_zh.md](README_zh.md). No MetPy, NCL, Pint, or extra runtime dependencies — NumPy only, string unit kwargs.

**v2.2.1** adds NCL official-documentation example regression tests (`tests/test_ncl_official_examples.py`).

**v2.2.2** attempted to publish as `py-meteo` (`pymeteo` was taken); PyPI rejected that name as too similar to existing `pymeteo`.

**v2.2.3** publishes the package on PyPI as `pymeteo-kit`. The import remains `import pymeteo`.

The old names (`showalter`, `E_WATER`, `SWEAT_calculate`, …) are **gone** from the top-level modern API—no compatibility aliases there.

## NCL names (`pymeteo.ncl`)

Optional **thin wrappers** with NCL builtin names and **fixed NCL units** (K, Pa/hPa, %, kg/kg). They only translate arguments and call the modern functions. These names are **not** re-exported from top-level `pymeteo`.

```python
import pymeteo as pm
from pymeteo.ncl import dewtemp_trh, relhum_ttd, mixhum_ptrh, wind_speed

td_k = dewtemp_trh(18.0 + 273.15, 46.5)   # tk in K, rh in % → dewpoint K
td_k = pm.ncl.dewtemp_trh(18.0 + 273.15, 46.5)
rh = pm.ncl.relhum_ttd(291.15, 279.45, 0)  # opt=0 → %, opt=1 → fraction
```

Full mapping (NCL name → modern function and units) is in [README_zh.md](README_zh.md).

## License

GPL-3.0. See [LICENSE](LICENSE).
