# pymeteo

Dewpoint, vapour mixing, potential and equivalent potential temperature, LCL, wet-bulb, visibility, Showalter / K / A / TT / SWEAT / lifted index, wind speed and from-direction, geodetic distance, station-to-sea-level pressure, NWS heat index and wind chill, Coriolis parameter, and hydrostatic ω ↔ w.

Python **3.7+**. NumPy only. Units are strings on the call (`temperature_unit="C"`, `pressure_unit="hPa"`, …). Returns `float` or NumPy arrays.

```python
import pymeteo as pm

td = pm.dewpoint_from_relative_humidity(
    18.0, 46.5, temperature_unit="C", humidity_unit="%"
)
si = pm.showalter_index(16.6, 0.6, -15.9, temperature_unit="C")
u, v = pm.uv_from_speed_direction(10.0, 270.0, speed_unit="kt")
km = pm.earth_distance(39.9, 116.4, 31.2, 121.5, output_distance_unit="km")
```

Defaults: temperature `C`, pressure `hPa`, wind `m/s`, relative humidity `%`, mixing ratio `kg/kg`, distance `km` (thickness defaults to `m`), angles in degrees. Unknown strings raise `pymeteo.UnitError`. Details: [Units](Units.md).

```bash
pip install pymeteo-kit
```

`import pymeteo`. NCL builtin names (fixed NCL units) are in [`pymeteo.ncl`](NCL.md), not at the package root.

中文：[首页](Home_zh.md).

| Page | Contents |
|------|----------|
| [Thermo](Thermo.md) | Vapour, θ, θe, LCL, wet-bulb, visibility |
| [Indices](Indices.md) | Showalter, K, A, TT, SWEAT, LI |
| [Wind](Wind.md) | Speed, direction, components, bulk shear |
| [Geo](Geo.md) | Distance, gravity, SLP, thickness |
| [Comfort](Comfort.md) | NWS heat index, wind chill |
| [Dynamics](Dynamics.md) | Coriolis, ω ↔ w |
| [Units](Units.md) | Aliases and `UnitError` |
| [NCL](NCL.md) | NCL name shims |
| [Install](Install.md) | `pip install pymeteo-kit` |

MIT License. [GitHub](https://github.com/IncubatorShokuhou/pyMeteo) · [PyPI](https://pypi.org/project/pymeteo-kit/)
