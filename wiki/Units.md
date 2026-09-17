# Units

[中文](Units_zh.md)

Each dimensional argument takes a string unit (sometimes a pair of input/output units). Values are converted to an internal SI-ish scale, then converted back. Scalars stay Python `float`; arrays stay arrays.

Unknown or inapplicable strings raise `pymeteo.UnitError`, a `ValueError` subclass.

There is no Pint.

## Defaults

| Quantity | Internal | Usual default | Common aliases |
|----------|----------|---------------|----------------|
| Temperature | K | `C` | `degC`, `celsius`, `K`, `kelvin`, `F`, `fahrenheit` |
| Pressure | Pa | `hPa` | `mb`, `mbar`, `Pa`, `kPa`, `atm` |
| Wind speed | m/s | `m/s` | `kt`, `knots`, `km/h`, `mph` |
| Relative humidity | 0–1 | `%` | `percent`, `fraction` |
| Mixing ratio / specific humidity | kg/kg | `kg/kg` | `g/kg`, `mg/kg` |
| Distance | m | `km` | `m`, `mi`, `ft`, `nmi` |
| Angle | rad before trig | `deg` | `degree`, `rad` |

Layer thickness (`height_thickness`) defaults to **`m`**, not `km`. Potential temperature defaults to **K** on the way out.

Showalter, K, A, TT, and lifted index are temperature differences. `C` and `K` give the same number.

Output kwargs look like `output_temperature_unit`, `output_pressure_unit`, `output_speed_unit`, `output_distance_unit`, `output_humidity_unit`, `output_omega_unit`.

## Helpers in `pymeteo.units`

Used internally. Only `UnitError` is re-exported at the package root.

Canonicalisers: `canonical_temperature_unit`, `canonical_pressure_unit`, `canonical_speed_unit`, `canonical_rh_unit`, `canonical_mass_humidity_unit`, `canonical_distance_unit`, `canonical_angle_unit`.

Converters: `to_kelvin` / `from_kelvin`, `to_pascal` / `from_pascal`, `to_mps` / `from_mps`, `to_rh_fraction` / `from_rh_fraction`, `to_kgkg` / `from_kgkg`, `to_meters` / `from_meters`, `to_radians`.

Knots are `1852/3600` m/s. A statute mile is 1609.344 m. A nautical mile is 1852 m. One atmosphere is 101325 Pa.

```python
from pymeteo.units import to_kelvin, from_kelvin

k = to_kelvin(18.0, "C")       # 291.15
c = from_kelvin(k, "celsius")  # 18.0
```
