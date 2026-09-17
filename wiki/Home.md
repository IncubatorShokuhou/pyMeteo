# pymeteo

Pointwise meteorological diagnostics for **Python 3.7+**, with **explicit string unit parameters**. Runtime dependency: **NumPy only**. Import name: `pymeteo`. PyPI distribution: [`pymeteo-kit`](https://pypi.org/project/pymeteo-kit/).

There is no Pint unit object, no Sounding/Wind class facade, no MetPy, and no NCL runtime. The library is not an I/O, plotting, map-projection, or full CAPE/CIN sounding package.

```python
import pymeteo as pm

td = pm.dewpoint_from_relative_humidity(18.0, 46.5, temperature_unit="C", humidity_unit="%")
si = pm.showalter_index(16.6, 0.6, -15.9, temperature_unit="C")
km = pm.earth_distance(39.9, 116.4, 31.2, 121.5, output_distance_unit="km")
```

Default public units: temperature `C`, pressure `hPa`, wind `m/s`, relative humidity `%`, mixing ratio `kg/kg`, distance `km` (layer thickness defaults to `m`), angles in degrees. Unrecognised unit strings raise `pymeteo.UnitError`.

Top-level names are English snake_case. NCL builtin names live in [`pymeteo.ncl`](NCL) and are **not** re-exported from the package root.

中文入口：[首页](Home_zh)。

## Pages

| Page | Module | Contents |
|------|--------|----------|
| [Install](Install) | — | `pip install pymeteo-kit`, Python 3.7+, development |
| [Units](Units) | `pymeteo.units` | Aliases, internals, `UnitError` |
| [Thermo](Thermo) | `pymeteo.thermo` | Vapour, θ, θe, LCL, wet-bulb, visibility |
| [Indices](Indices) | `pymeteo.indices` | Showalter, K, A, TT, SWEAT, LI |
| [Wind](Wind) | `pymeteo.wind` | Speed, direction, components, bulk shear |
| [Geo](Geo) | `pymeteo.geo` | Distance, gravity, SLP, thickness |
| [Comfort](Comfort) | `pymeteo.comfort` | NWS heat index, wind chill |
| [Dynamics](Dynamics) | `pymeteo.dynamics` | Coriolis, ω ↔ w |
| [NCL](NCL) | `pymeteo.ncl` | Fixed-unit name shims |
| [Publishing](Publishing) | — | PyPI release, tests, contributing |

## Public top-level names

Re-exported from `import pymeteo as pm` (see `pymeteo.__all__`):

`UnitError`, `__version__`, `ncl`, and the modern functions listed on the module pages above.

## Scope (what this is not)

No FAO56 radiation/ET suite, no grid spherical-harmonic advection, no cross-section, no complete CAPE/CIN sounding kit. Heat index and wind chill have **no** NCL-name shims; the ncl module does not invent those names.

## Links

* Repository: https://github.com/IncubatorShokuhou/pyMeteo
* Site: https://incubatorshokuhou.github.io/pyMeteo/
* PyPI: https://pypi.org/project/pymeteo-kit/
* License: MIT (see the `LICENSE` file in the source tree)
