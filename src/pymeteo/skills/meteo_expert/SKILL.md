---
name: meteo-expert
description: Meteorological diagnostics helper for pymeteo. Activate when the user asks which pymeteo function to call, or about dewpoint, relative humidity, mixing ratio, potential temperature, LCL, wet-bulb, virtual temperature, visibility, Showalter/K/TT/SWEAT/lifted index, wind speed/direction/u-v/shear, great-circle distance, gravity, sea-level pressure, layer thickness, Coriolis, omega versus w, heat index, wind chill, NCL names, or pymeteo unit strings. Use to pick a function, check units, and avoid documented pitfalls. The classic import pymeteo as pm API is unchanged.
---

You are helping with **pymeteo**, a small NumPy library of pointwise meteorological diagnostics. Import name is `pymeteo`; PyPI name is `pymeteo-kit`. There is no file I/O, no plotting, and no sounding object.

Classic call, unchanged:

```python
import pymeteo as pm

td = pm.dewpoint_from_relative_humidity(
    18.0, 46.5, temperature_unit="C", humidity_unit="%"
)
```

Do not invent MetPy, Pint, or NCL runtime dependencies. Do not wrap results in a session object unless the user asked for `MeteoEngine`.

## When to activate

- “Which pymeteo function computes …?”
- Dewpoint, RH, mixing ratio, θ, θe, LCL, wet-bulb, virtual T, visibility
- Stability indices (Showalter, K, A, TT, SWEAT, lifted index)
- Wind speed / from-direction / u-v / bulk shear
- Great-circle distance, gravity, sea-level pressure, layer thickness
- Coriolis, ω ↔ w
- Heat index, wind chill
- NCL names (`dewtemp_trh`, `relhum`, `mixhum_ptrh`, …)
- Unit strings, `UnitError`, “is this hPa or Pa?”

## What to load

| File | Open when |
|------|-----------|
| `references/workflow.md` | Choosing a module or function |
| `references/units.md` | Any call with a unit kwarg |
| `references/pitfalls.md` | Before reporting a number |
| `references/ncl.md` | User said NCL, or a name like `dewtemp_trh` |
| `references/api.md` | Need a compact public-name table |

## How to pick a function

1. Decide the domain (thermo / indices / wind / geo / dynamics / comfort) using `references/workflow.md`.
2. Prefer the modern root API (`pm.dewpoint_from_relative_humidity`) over `pymeteo.ncl` unless the user is porting NCL and wants NCL’s fixed units.
3. Pass units as strings. Defaults are in each docstring; they are **not** all SI. Potential temperature output defaults to `K`. `earth_distance` output defaults to `km`; `height_thickness` defaults to `m`.
4. If the user mentioned an NCL builtin, resolve it with `references/ncl.md` or `MeteoEngine.ncl_lookup`. Some NCL pressure units are Pa, some hPa; `omega_to_w` argument order differs.
5. For a one-off scalar check you may use `MeteoEngine.run_calc` (JSON scalars/lists only). For real work, emit ordinary Python that calls `pymeteo`.

Planning from Python (no MCP required):

```python
from pymeteo.engine import MeteoEngine

engine = MeteoEngine()
engine.recommend("dewpoint from RH")
engine.explain("showalter_index")
engine.ncl_lookup("dewtemp_trh")
engine.unit_help("temperature")
```

MCP (optional extra): `python -m pymeteo.mcp_server` exposes `list_functions`, `explain_function`, `recommend`, `ncl_lookup`, `unit_help`, `run_calc`.

## Hard rules

- Do not call NCL names from `import pymeteo as pm` — they are not in `__all__`. Use `from pymeteo.ncl import dewtemp_trh`.
- Do not assume ice-phase saturation; `saturation_vapor_pressure` is liquid water (李社宏 1994).
- Do not silently convert °C as if it were K. If units are unspecified, use the function default and say so.
- `UnitError` is a `ValueError` subclass raised for unknown unit strings.
- This library does not compute CAPE, CIN, or a full sounding.

Read `references/pitfalls.md` before you quote a numeric result.
