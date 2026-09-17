# Pitfalls

Grounded in current `src/pymeteo` behaviour. Check these before you quote a number.

## Water vapour

- **Liquid only.** `saturation_vapor_pressure` uses 李社宏 (1994) *water* formula. Docstring: not ice. 0 °C → 6.1078 hPa. Do not treat this as frost point.
- **Several saturation formulae.** Dewpoint ↔ RH uses Dutton latent heat. Mixing ratio from RH uses Tetens (es0 = 6.11, NCL `mixhum_ptrh`). RH from mixing ratio uses the NCL `relhum` lookup table (173.16–375.16 K). They will not match to many decimals. Say which function you called.
- **Td > T.** `relative_humidity_from_dewpoint` will return RH > 100%; it does not clip. `dewpoint_from_relative_humidity` with RH > 100% yields Td > T. Bolton LCL clamps Td to T (already saturated). `condensation_temperature` (李社宏 iterator) may fail to converge on Td > T or negative mixing ratio and then returns the last iterate.
- **RH = 0.** `dewpoint_from_relative_humidity` returns `nan` for RH ≤ 0.
- **RH from mixing ratio.** Values **above 100% are kept**. Values below 0 become 0.0001%. Temperature is clipped to the lookup table range.
- **Specific humidity ≥ 1.** `convert_humidity` to mixing ratio returns `nan` (impossible).
- **Wet-bulb.** `wet_bulb_temperature` is Stull (2011) **sea-level** (~1013.25 hPa). Valid roughly −20–50 °C and RH 5–99%. Not a pressure-dependent psychrometric wet-bulb.
- **Visibility.** `method` is `RUC` or `FSL` only. RUC is an exponential decay from 60 km, not the old implementation that multiplied by 1000 twice.

## Indices

- **SWEAT.** Td850 is in °C; wind terms are in **knots**. Negative terms are zeroed. The shear term is kept only when 850° is 130–250, 500° is 210–310, 500°−850° > 0, and both speeds ≥ 15 kt. Do not pass Kelvin dewpoint as if it were Celsius.
- **Showalter.** Mixing ratio in the iterator uses **T850 saturation**, not the dewpoint (李社宏 format). Grid points with Td > T or nonsense mixing ratios can stall the iteration.
- **TT / SWEAT humidity.** Provide exactly one of `dewpoint_850` or `relative_humidity_850`.
- **`lifted_index_from_surface`.** Target is always 500 hPa. `pressure_unit` only interprets the *surface* pressure you passed. This is not CAPE.

## Wind

- Meteorological **from**-direction. `u = −speed sin(dir)`, `v = −speed cos(dir)`.
- Modern `wind_direction`: calm → **0**. NCL `wind_direction` `opt=1` → `nan`, other scalar → fill.
- `bulk_wind_shear` is the vector difference modulus at two points, not a height-integrated 0–6 km shear.

## Geo / dynamics

- **`earth_distance`** is WGS84 Vincenty; antipodes may fall back to mean-radius haversine. Default output **km**.
- **`gravity`**. Latitude is converted from `deg` (default) to radians. Do not pass degrees into a raw `sin`.
- **`height_thickness`**. Default output **m**. Pass virtual temperature as `mean_temperature` if the layer is moist.
- **`omega_to_w`**. Hydrostatic, *dry* ideal gas (Rd = 287.058). Modern argument order is (omega, temperature, pressure) with pressure default **hPa**. NCL is (omega, p, t) with p in **Pa**. Mixing these up is a 100× pressure error.

## Comfort

- Heat index: Rothfusz/NWS. Unreliable outside the Steadman table. Extreme T/RH should be flagged.
- Wind chill: NWS 2001. Intended for T ≤ 50 °F and V ≥ 3 mph; the formula is still evaluated outside that range.

## Packaging / names

- Top-level `pm.relhum` / `pm.dewtemp_trh` **do not exist**. Those names are only on `pymeteo.ncl`.
- Do not add MetPy or Pint “to be safe”. pymeteo will not read their unit objects.
