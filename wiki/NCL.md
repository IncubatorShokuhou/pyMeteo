# NCL compatibility (`pymeteo.ncl`)

[中文](NCL_zh)

Thin wrappers with NCL builtin **names** and **fixed NCL units**. They only translate arguments and call the modern functions. They are **not** re-exported from top-level `pymeteo` (`from pymeteo import dewtemp_trh` fails). Need flexible units? Call the modern API.

```python
import pymeteo as pm
from pymeteo.ncl import dewtemp_trh, relhum_ttd, mixhum_ptrh, wind_speed

td_k = dewtemp_trh(18.0 + 273.15, 46.5)        # tk K, rh % → dewpoint K (~279.45)
td_k = pm.ncl.dewtemp_trh(18.0 + 273.15, 46.5)
rh = pm.ncl.relhum_ttd(291.15, 279.45, 0)      # opt=0 → %, opt=1 → fraction
```

Showalter / K / SWEAT / heat index / wind chill have no matching NCL builtins here; those names are not invented.

Regression tests that pin numbers from NCL documentation examples live in `tests/test_ncl_official_examples.py`. Wiring and unit flags are in `tests/test_ncl.py`.

## Map

| NCL name | Fixed units / flags | Calls |
|----------|---------------------|--------|
| `dewtemp_trh(tk, rh)` | `tk` K, `rh` %, return dewpoint K | `dewpoint_from_relative_humidity` |
| `relhum_ttd(t, td, opt)` | `t`/`td` K; `opt=0` → %, `opt=1` → fraction | `relative_humidity_from_dewpoint` |
| `relhum(t, w, p)` | `t` K, `w` kg/kg, `p` Pa, return % | `relative_humidity_from_mixing_ratio` |
| `mixhum_ptrh(p, tk, rh, iswit)` | `p` **hPa**, `tk` K, `rh` %; `\|iswit\|=1` mixing ratio, `2` specific humidity; sign − → g/kg, + → kg/kg | `mixing_ratio_from_relative_humidity` / `specific_humidity_from_relative_humidity` |
| `mixhum_ptd(p, tdk, iswit)` | `p` **Pa**, `tdk` K; `iswit` as above | `mixing_ratio_from_dewpoint` (+ convert) |
| `mixhum_convert(wq, wqType, iounit)` | `wqType` `"w"` mixing→specific, `"q"` reverse; `iounit=(in,out)` 0=kg/kg, 1=g/kg | `convert_humidity` |
| `vapor_pres_rh(rh, es)` | `rh` %; `es` and return share a unit | `RH/100 · e_s` |
| `pot_temp(p, t, dim=-1, opt=False)` | `p` Pa, `t` K, return K. `dim`/`opt` ignored | `potential_temperature` |
| `pot_temp_equiv(p, t, w, dim=-1, humVarType="r")` | `p` Pa, `t` K; `humVarType` `"r"`/`"w"` mixing kg/kg, `"q"` specific humidity, `"rh"` RH %. Internally Bolton with LCL (closer to `pot_temp_equiv_tlcl` than NCL 6.4’s no-LCL approx). `dim` ignored | `equivalent_potential_temperature` |
| `temp_virtual(t, w, iounit)` | `iounit` length 3: T in 0=°C/1=K/2=°F, mixing 0=kg/kg or 1=g/kg, output T. Uses `T(1+r/ε)/(1+r)`, not NCL’s `T(1+0.61 r)` approx | `virtual_temperature` |
| `wetbulb_stull(t, rh, iounit, opt=False)` | `rh` %; `iounit` length 2 (0=°C, 1=K, 2=°F). `opt` unused | `wet_bulb_temperature` |
| `lclvl(p, tk, tdk)` | `p` hPa, temperatures K, return LCL pressure hPa only | `lifting_condensation_level` |
| `wind_speed(u, v)` | m/s → m/s | `wind_speed` |
| `wind_direction(u, v, opt=0)` | from-direction degrees; calm: `opt=0` → 0, `opt=1` → nan, other scalar → fill | `wind_direction` |
| `wind_component(wspd, wdir, opt=0)` | m/s and from-direction → `(u, v)` tuple. NCL `opt` unused | `uv_from_speed_direction` |
| `coriolis_param(lat)` | latitude degrees → s⁻¹ | `coriolis_parameter` |
| `omega_to_w(omega, p, t)` | Pa/s, **Pa**, K → m/s. Order `(omega, p, t)` | `omega_to_w` (modern order is temperature before pressure) |
| `w_to_omega(w, p, t)` | inverse of the above | `w_to_omega` |

## Examples that match the test suite

```python
# mixhum_ptrh, Wallace & Hobbs-style 1000 hPa / 18 °C / 46.5 %
mix_kg = pm.ncl.mixhum_ptrh(1000.0, 18.0 + 273.15, 46.5, 1)   # ~0.006018 kg/kg
mix_g  = pm.ncl.mixhum_ptrh(1000.0, 18.0 + 273.15, 46.5, -1)  # ~6.018 g/kg

# pot_temp: 100000 Pa, 301.25 K → 301.25 K
pm.ncl.pot_temp(100000.0, 301.25)

# wetbulb_stull: 20 °C, 50 %, iounit (0,0) → ~13.70 °C
pm.ncl.wetbulb_stull(20.0, 50.0, (0, 0), False)

# lclvl: 1000 hPa, 15 °C, Td 4 °C → ~847–849 hPa (Bolton vs Stipanuk)
pm.ncl.lclvl(1000.0, 15.0 + 273.15, 4.0 + 273.15)
```
