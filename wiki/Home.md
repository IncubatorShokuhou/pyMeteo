# pymeteo

pymeteo is a set of pointwise meteorological calculations on NumPy arrays.

It does not read files, plot maps, or wrap a full sounding object. You pass scalars or arrays in and get scalars or arrays back.

## Install

```bash
pip install pymeteo-kit
```

```python
import pymeteo as pm
```

Python 3.7+. The import name is `pymeteo`.

## Example

```python
import pymeteo as pm

td = pm.dewpoint_from_relative_humidity(18.0, 46.5)
# about 6.3 °C
```

Units are strings on the call (`temperature_unit="C"`). Defaults live in each function’s docstring.

## Pages

- [Thermo](Thermo.md) — vapour, θ, LCL, wet-bulb
- [Indices](Indices.md) — Showalter, K, TT, SWEAT, lifted index
- [Wind](Wind.md) — speed, from-direction, u/v, bulk shear
- [Geo](Geo.md) — great-circle distance, gravity, sea-level pressure
- [Comfort](Comfort.md) — heat index, wind chill
- [Dynamics](Dynamics.md) — Coriolis, ω ↔ w
- [Units](Units.md) — aliases and `UnitError`
- [NCL](NCL.md) — NCL names with NCL’s unit conventions
- [Agent](Agent.md) — skill / MCP / `MeteoEngine`
- [Install](Install.md)

中文：[首页](Home_zh.md)

MIT License. [GitHub](https://github.com/IncubatorShokuhou/pyMeteo) · [PyPI](https://pypi.org/project/pymeteo-kit/)
