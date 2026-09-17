# Units

[中文](Units_zh.md)

Each dimensional argument has a string unit keyword (or a pair of input/output units). Values are converted to an internal SI-ish scale, then converted back. Results are ordinary `float` or NumPy arrays (`restore_shape` returns a Python float when every original input was a scalar).

Unrecognised or inapplicable strings raise **`pymeteo.UnitError`** (a `ValueError` subclass).

## Defaults and internals

| Quantity | Internal | Default public unit | Aliases (not exhaustive) |
|----------|----------|---------------------|--------------------------|
| Temperature | K | `C` | `degC` / `celsius` / `centigrade`, `K` / `kelvin`, `F` / `degF` / `fahrenheit` |
| Pressure | Pa | `hPa` | `mb` / `mbar` / `millibar`, `Pa` / `pascal`, `kPa`, `atm` |
| Wind speed | m/s | `m/s` | `ms-1` / `mps`, `kt` / `knot` / `knots` / `kn`, `km/h` / `kmh` / `kph`, `mph` |
| Relative humidity | fraction 0–1 | `%` | `percent` / `percentage`, `fraction` / `1` / `ratio` |
| Mixing ratio / specific humidity | kg/kg | `kg/kg` | `g/kg`, `mg/kg`, `kgkg-1`, `g/g` |
| Distance | m | `km` | `m` / `meter`, `cm`, `mi` / `mile`, `ft` / `foot`, `nmi` / `nm` / `nauticalmile` |
| Angle (lat/lon, wind from-direction) | converted to rad before trig | `deg` | `degree` / `degrees`, `rad` / `radian` |

Layer **thickness** (`height_thickness`) defaults to **`m`**, not `km`.

Stability indices (Showalter, K, A, TT, SWEAT) are temperature differences or dimensionless combinations. TT / K / A / Showalter have the **same numeric value** in Celsius and Kelvin.

Output unit kwargs follow the pattern `output_temperature_unit`, `output_pressure_unit`, `output_speed_unit`, `output_distance_unit`, `output_humidity_unit`, `output_omega_unit`.

## Conversion helpers (`pymeteo.units`)

These are used internally. They are not re-exported from the package root (except `UnitError`).

Canonicalisers: `canonical_temperature_unit`, `canonical_pressure_unit`, `canonical_speed_unit`, `canonical_rh_unit`, `canonical_mass_humidity_unit`, `canonical_distance_unit`, `canonical_angle_unit`.

Converters: `to_kelvin` / `from_kelvin`, `to_pascal` / `from_pascal`, `to_mps` / `from_mps`, `to_rh_fraction` / `from_rh_fraction`, `to_kgkg` / `from_kgkg`, `to_meters` / `from_meters`, `to_radians`.

Knots use `1852/3600` m/s. Statute mile is `1609.344` m. Nautical mile is `1852` m. Atmosphere is `101325` Pa.

```python
from pymeteo.units import UnitError, to_kelvin, from_kelvin

k = to_kelvin(18.0, "C")          # 291.15
c = from_kelvin(k, "celsius")     # 18.0
```
