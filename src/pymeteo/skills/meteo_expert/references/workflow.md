# Workflow: which module?

pymeteo is a bag of functions, not a sounding pipeline. Pick the module from the *quantity the user wants back*, then pick the function from the *inputs they actually have*.

## Decision

```
Want a water-vapour or temperature-derived state variable
(dewpoint, RH, mixing ratio, e, es, θ, θe, LCL, Tw, Tv, visibility)?
  → pymeteo.thermo  (re-exported at the package root)

Want a named stability index (Showalter, K, A, TT, SWEAT, LI)?
  → pymeteo.indices
  Showalter needs T850, Td850, T500.
  K needs 850/700/500 temperatures and 850/700 dewpoints.
  SWEAT also needs 850 and 500 u, v.
  lifted_index_from_surface lifts a near-surface parcel to 500 hPa.
  lifted_index is just Tenv(500) − Tparcel(500); you already did the lift.

Want speed, meteorological from-direction, u/v, or two-level vector shear?
  → pymeteo.wind
  wind_components is an alias of uv_from_speed_direction.

Want great-circle distance, g(φ), station→MSLP, or isobaric thickness?
  → pymeteo.geo
  height_thickness wants a layer-mean temperature (use virtual temperature if moist).

Want f = 2Ω sinφ or ω ↔ w under hydrostatic ideal-gas?
  → pymeteo.dynamics
  Modern omega_to_w / w_to_omega take (omega_or_w, temperature, pressure).
  NCL wrappers take (omega_or_w, p, t) with p in Pa.

Want NWS heat index or wind chill?
  → pymeteo.comfort
  No NCL names for these.
```

NCL builtin **names** live only in `pymeteo.ncl`. They translate arguments and then call the modern functions. Use them when porting NCL; use the modern API when you care about unit kwargs.

## Inputs you must not invent

- If the user has T and RH, you can get dewpoint (`dewpoint_from_relative_humidity`). You do not need mixing ratio first.
- If the user has T and Td, use `relative_humidity_from_dewpoint`, not the mixing-ratio path.
- Mixing ratio from RH needs **pressure** (`mixing_ratio_from_relative_humidity`). Mixing ratio from dewpoint also needs pressure.
- `relative_humidity_from_mixing_ratio` needs T, r, and p.
- SWEAT and TT: pass **either** `dewpoint_850` **or** `relative_humidity_850`, not both.
- `lifted_index_from_surface` needs near-surface p, T, Td and T500. The 500 hPa target is fixed; `pressure_unit` does not move it.
- Wind direction is meteorological **from**-direction, north 0°, clockwise. Calm → 0 in the modern API.

## Planning helpers

```python
from pymeteo.engine import MeteoEngine
engine = MeteoEngine()
engine.list_functions(module="thermo")
engine.recommend("wet-bulb from T and RH")
engine.explain("wet_bulb_temperature")
```

`recommend` is keyword overlap, not a physical solver. Still read the docstring of the function you emit.

## Out of scope

No GRIB/NetCDF readers, no maps, no hodographs, no CAPE/CIN, no hydrometeor classification, no ice-bulb. If the user needs those, say so and stay on the point diagnostics pymeteo actually implements.
