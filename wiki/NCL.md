# NCL

[中文](NCL_zh.md)

Thin wrappers that keep NCL builtin **names** and NCL’s **fixed units**. They translate arguments, then call the modern functions. They are not re-exported from top-level `pymeteo`.

Need flexible units? Use the modern API.

```python
from pymeteo.ncl import dewtemp_trh

td_k = dewtemp_trh(18.0 + 273.15, 46.5)  # tk in K, rh in %, dewpoint in K
```

Showalter, K, SWEAT, heat index, and wind chill have no NCL builtins here. Those names are not invented.

Numbers from NCL documentation examples are pinned in `tests/test_ncl_official_examples.py`.

## Map

| NCL name | Units / flags | Calls |
|----------|---------------|--------|
| `dewtemp_trh(tk, rh)` | `tk` K, `rh` %, dewpoint K | `dewpoint_from_relative_humidity` |
| `relhum_ttd(t, td, opt)` | `t`/`td` K; `opt=0` → %, `opt=1` → fraction | `relative_humidity_from_dewpoint` |
| `relhum(t, w, p)` | `t` K, `w` kg/kg, `p` Pa, result % | `relative_humidity_from_mixing_ratio` |
| `mixhum_ptrh(p, tk, rh, iswit)` | `p` **hPa**, `tk` K, `rh` %; `\|iswit\|=1` mixing ratio, `2` specific humidity; minus → g/kg | mixing / specific humidity from RH |
| `mixhum_ptd(p, tdk, iswit)` | `p` **Pa**, `tdk` K; `iswit` as above | `mixing_ratio_from_dewpoint` |
| `mixhum_convert(wq, wqType, iounit)` | `"w"` mixing→specific, `"q"` reverse; `iounit=(in,out)` 0=kg/kg, 1=g/kg | `convert_humidity` |
| `vapor_pres_rh(rh, es)` | `rh` %; `es` and the result share a unit | `RH/100 · e_s` |
| `pot_temp(p, t, dim=-1, opt=False)` | `p` Pa, `t` K → K. `dim`/`opt` ignored | `potential_temperature` |
| `pot_temp_equiv(p, t, w, dim=-1, humVarType="r")` | `p` Pa, `t` K; `"r"`/`"w"` mixing kg/kg, `"q"` specific humidity, `"rh"` RH %. Bolton with LCL (closer to NCL `pot_temp_equiv_tlcl` than the no-LCL approx). `dim` ignored | `equivalent_potential_temperature` |
| `temp_virtual(t, w, iounit)` | `iounit` length 3: T 0=°C/1=K/2=°F, mixing 0=kg/kg or 1=g/kg, output T. Uses `T(1+r/ε)/(1+r)`, not `T(1+0.61 r)` | `virtual_temperature` |
| `wetbulb_stull(t, rh, iounit, opt=False)` | `rh` %; `iounit` length 2 (0=°C, 1=K, 2=°F). `opt` unused | `wet_bulb_temperature` |
| `lclvl(p, tk, tdk)` | `p` hPa, temperatures K, returns LCL pressure only | `lifting_condensation_level` |
| `wind_speed(u, v)` | m/s | `wind_speed` |
| `wind_direction(u, v, opt=0)` | from-direction degrees; calm: `opt=0` → 0, `opt=1` → nan, other scalar → fill | `wind_direction` |
| `wind_component(wspd, wdir, opt=0)` | m/s, from-direction → `(u, v)` tuple. NCL `opt` unused | `uv_from_speed_direction` |
| `coriolis_param(lat)` | latitude degrees → s⁻¹ | `coriolis_parameter` |
| `omega_to_w(omega, p, t)` | Pa/s, **Pa**, K → m/s. Order `(omega, p, t)` | modern `omega_to_w` is `(omega, temperature, pressure)` |
| `w_to_omega(w, p, t)` | inverse | `w_to_omega` |

```python
import pymeteo as pm

pm.ncl.mixhum_ptrh(1000.0, 18.0 + 273.15, 46.5, 1)   # ~0.00602 kg/kg
pm.ncl.pot_temp(100000.0, 301.25)                     # 301.25 K
pm.ncl.wetbulb_stull(20.0, 50.0, (0, 0), False)       # ~13.7 °C
pm.ncl.lclvl(1000.0, 15.0 + 273.15, 4.0 + 273.15)     # ~847 hPa
```
