# pymeteo

NumPy routines for meteorological diagnostics: dewpoint, relative humidity, potential temperature, stability indices, wind, great-circle distance, and related point quantities.
Meant for sounding points and gridded pointwise fields. No file I/O, no plotting, no sounding object model.

[中文说明](README_zh.md)

## Install

Python 3.7+ (PyPI package name: `pymeteo-kit`):

```bash
pip install pymeteo-kit
```

```python
import pymeteo as pm
```

## Examples

```python
import pymeteo as pm

td = pm.dewpoint_from_relative_humidity(
    18.0, 46.5, temperature_unit="C", humidity_unit="%"
)
si = pm.showalter_index(16.6, 0.6, -15.9, temperature_unit="C")
ws = pm.wind_speed(3.0, 4.0)
km = pm.earth_distance(39.9, 116.4, 31.2, 121.5)
```

## Units

Pass units as strings (`temperature_unit="C"`); conversion happens inside the function. Defaults are in each function’s docstring. No Pint.

## Main API

Public names are re-exported at the package root. See source docstrings for arguments.

| Module | Contents |
|--------|----------|
| `pymeteo.thermo` | dewpoint, relative humidity, mixing ratio, saturation vapor pressure, potential / equivalent potential temperature, LCL, wet-bulb, virtual temperature |
| `pymeteo.indices` | Showalter, K / A / TT, SWEAT, lifted index |
| `pymeteo.wind` | speed, direction, u/v, bulk shear |
| `pymeteo.geo` | great-circle distance, gravity, sea-level pressure, layer thickness |
| `pymeteo.dynamics` | Coriolis parameter, ω ↔ w |
| `pymeteo.comfort` | heat index, wind chill |

## NCL names

A subset of NCL builtins lives in `pymeteo.ncl`, with NCL’s unit conventions (temperature in K, humidity in %, …):

```python
from pymeteo.ncl import dewtemp_trh

td_k = dewtemp_trh(18.0 + 273.15, 46.5)
```

## License

MIT License. See [LICENSE](LICENSE).
