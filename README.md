# pymeteo

Lightweight meteorological diagnostic functions for Python 3.10+.

This is a **breaking rewrite** of the old single-file `pyMeteo.py`. The import name is now lowercase `pymeteo`. Public functions take **explicit string unit parameters** (no Pint, no Sounding/Wind object facade). Runtime dependency: **NumPy only**.

**Full documentation is in Chinese:** [README_zh.md](README_zh.md)

## Install

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

Default units: temperature `C`, pressure `hPa`, wind `m/s`, relative humidity `%`, mixing ratio `kg/kg`, distance `km`, angles in degrees. Aliases such as `K`/`kelvin`, `hPa`/`mb`, `kt`, `fraction`, `g/kg`, `m`/`nmi` are documented in [README_zh.md](README_zh.md).

The old names (`showalter`, `E_WATER`, `SWEAT_calculate`, …) are **gone**—no compatibility aliases.

## License

GPL-3.0. See [LICENSE](LICENSE).
