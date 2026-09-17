# pymeteo

Pointwise meteorological diagnostics with explicit string unit parameters (Python 3.7+).

Functions convert units internally and return plain `float` values or NumPy arrays. The only runtime dependency is NumPy. This is not an I/O, plotting, map-projection, or full CAPE/CIN sounding package, and it does not use Pint, MetPy, or an NCL runtime.

[中文说明](README_zh.md)

## Install

```bash
pip install pymeteo-kit
```

The PyPI project is `pymeteo-kit` because `pymeteo` is already taken; the import name is still `pymeteo`.

## Quick start

```python
import pymeteo as pm

dewpoint = pm.dewpoint_from_relative_humidity(
    18.0, 46.5, temperature_unit="C", humidity_unit="%"
)
si = pm.showalter_index(16.6, 0.6, -15.9, temperature_unit="C")
speed = pm.wind_speed(3.0, 4.0, speed_unit="m/s")
km = pm.earth_distance(39.9, 116.4, 31.2, 121.5, output_distance_unit="km")
```

## Units

Dimensional arguments take string kwargs (`temperature_unit`, `pressure_unit`, `output_*_unit`, …). Unrecognized strings raise `pymeteo.UnitError`. Defaults:

| Quantity | Default | Common aliases |
|----------|---------|----------------|
| temperature | `C` | `K`, `F` |
| pressure | `hPa` | `mb`, `Pa` |
| wind | `m/s` | `kt`, `km/h` |
| relative humidity | `%` | `fraction` |
| mixing ratio | `kg/kg` | `g/kg` |
| distance | `km` | `m`, `nmi` (`height_thickness` defaults to `m`) |
| angles | degrees | `rad` |

## Modules

`pymeteo.thermo`, `indices`, `wind`, `geo`, `dynamics`, and `comfort` are re-exported at the package root. Parameter details are in the source docstrings.

`pymeteo.ncl` is a thin shim: NCL builtin names and NCL’s fixed units (K, Pa/hPa, %, kg/kg). Those names are not re-exported from `pymeteo`.

```python
from pymeteo.ncl import dewtemp_trh

td_k = dewtemp_trh(18.0 + 273.15, 46.5)
```

## Development

From a checkout on Python 3.8+:

```bash
pip install -e ".[dev]"
pytest
ruff check src tests
```

## License

GPL-3.0. See [LICENSE](LICENSE).
