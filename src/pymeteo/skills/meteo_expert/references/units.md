# Units

Every dimensional argument is a Python number plus a **string** unit kwarg. Conversion happens inside the function. There is no Pint.

Unknown or inapplicable strings raise `UnitError`, a `ValueError` subclass from `pymeteo.units` (re-exported as `pymeteo.UnitError`).

Converters in `pymeteo.units` are internal. Agents should pass the public kwargs, not call those helpers, unless the user is writing unit tests.

## Internal scale

| Quantity | Internal | Usual public default | Common aliases |
|----------|----------|----------------------|----------------|
| Temperature | K | `C` | `degC`, `celsius`, `K`, `kelvin`, `F`, `fahrenheit` |
| Pressure | Pa | `hPa` | `mb`, `mbar`, `Pa`, `kPa`, `atm` |
| Wind speed | m/s | `m/s` | `kt`, `knots`, `kn`, `km/h`, `mph` |
| Relative humidity | 0–1 | `%` | `percent`, `fraction`, `1`, `ratio` |
| Mixing ratio / specific humidity | kg/kg | `kg/kg` | `g/kg`, `mg/kg`, `kgkg-1` |
| Distance | m | `km` for `earth_distance`; **`m` for `height_thickness`** | `meter`, `mi`, `ft`, `nmi` / `nm` |
| Angle | rad before trig | `deg` | `degree`, `rad` |
| ω | Pa/s | `Pa/s` | `hPa/s`, `mb/s` (in `omega_to_w` / `w_to_omega` only) |

Output kwargs: `output_temperature_unit`, `output_pressure_unit`, `output_speed_unit`, `output_distance_unit`, `output_humidity_unit`, `output_omega_unit`.

If `output_*` is omitted, wind speed/heat index/wind chill/virtual temperature usually echo the input unit. **Potential temperature output defaults to `K`**, not `C`.

## Defaults that surprise people

- Modern thermo pressure kwargs default to **`hPa`**. Several NCL shims use **Pa** (`pot_temp`, `mixhum_ptd`, `relhum`, `omega_to_w`). `mixhum_ptrh` and `lclvl` use **hPa**. See `ncl.md`.
- Index functions (Showalter, K, A, TT, LI) are temperature *differences*. `C` and `K` give the same number. `F` for a depression multiplies by 9/5.
- `sea_level_pressure` `lapse_rate` default `0.005` is **degrees per metre in `temperature_unit`**. Leave it if temperatures are °C. Do not reuse 0.005 with `temperature_unit="F"`.
- Knots are `1852/3600` m/s. Statute mile 1609.344 m. Nautical mile 1852 m. 1 atm = 101325 Pa.

## Lookup from Python

```python
from pymeteo.engine import MeteoEngine
MeteoEngine().unit_help("pressure")
```

`mb` is an alias of `hPa`. `%` is a real key in the RH table; stripping spaces does not remove it.
