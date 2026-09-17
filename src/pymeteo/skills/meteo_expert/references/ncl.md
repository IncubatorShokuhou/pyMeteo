# NCL names

`pymeteo.ncl` is a thin compatibility layer: **NCL builtin names** and **NCL’s fixed units**. It does not re-implement the physics. Functions are **not** re-exported at `import pymeteo as pm`.

```python
from pymeteo.ncl import dewtemp_trh

td_k = dewtemp_trh(18.0 + 273.15, 46.5)  # tk in K, rh in %, dewpoint in K
```

Need flexible units? Call the modern name instead.

There are **no** NCL shims for Showalter, K, A, TT, SWEAT, lifted index, heat index, or wind chill.

## Map

| NCL name | Units / flags | Modern callable |
|----------|---------------|-----------------|
| `dewtemp_trh(tk, rh)` | `tk` K, `rh` %, dewpoint K | `dewpoint_from_relative_humidity` |
| `relhum_ttd(t, td, opt)` | `t`/`td` K; `opt=0` → %, `opt=1` → fraction | `relative_humidity_from_dewpoint` |
| `relhum(t, w, p)` | `t` K, `w` kg/kg, `p` **Pa**, result % | `relative_humidity_from_mixing_ratio` |
| `mixhum_ptrh(p, tk, rh, iswit)` | `p` **hPa**, `tk` K, `rh` %; `\|iswit\|=1` mixing ratio, `2` specific humidity; minus → g/kg | `mixing_ratio_from_relative_humidity` / `specific_humidity_from_relative_humidity` |
| `mixhum_ptd(p, tdk, iswit)` | `p` **Pa**, `tdk` K; `iswit` as above | `mixing_ratio_from_dewpoint` |
| `mixhum_convert(wq, wqType, iounit)` | `"w"` mixing→specific, `"q"` reverse; `iounit=(in,out)` 0=kg/kg, 1=g/kg | `convert_humidity` |
| `vapor_pres_rh(rh, es)` | `rh` %; `es` and result share a unit; no temperature | `e = (RH%/100)·e_s` (not the T-based modern helper) |
| `pot_temp(p, t, dim=-1, opt=False)` | `p` **Pa**, `t` K → K. `dim`/`opt` ignored | `potential_temperature` |
| `pot_temp_equiv(p, t, w, dim=-1, humVarType="r")` | `p` Pa, `t` K; `"r"`/`"w"` mixing kg/kg, `"q"` specific humidity, `"rh"` RH % | `equivalent_potential_temperature` (Bolton with LCL) |
| `temp_virtual(t, w, iounit)` | `iounit` length 3: T 0=°C/1=K/2=°F, mixing 0=kg/kg or 1=g/kg, output T | `virtual_temperature` — uses `T(1+r/ε)/(1+r)`, not `T(1+0.61 r)` |
| `wetbulb_stull(t, rh, iounit, opt=False)` | `rh` %; `iounit` length 2 (0=°C, 1=K, 2=°F) | `wet_bulb_temperature` |
| `lclvl(p, tk, tdk)` | `p` hPa, temperatures K, **pressure only** | `lifting_condensation_level` (modern also returns T_LCL) |
| `wind_speed(u, v)` | m/s | `wind_speed` |
| `wind_direction(u, v, opt=0)` | from-direction degrees; calm: `opt=0` → 0, `opt=1` → nan | `wind_direction` |
| `wind_component(wspd, wdir, opt=0)` | m/s, from-direction → `(u, v)` tuple | `uv_from_speed_direction` |
| `coriolis_param(lat)` | latitude degrees → s⁻¹ | `coriolis_parameter` |
| `omega_to_w(omega, p, t)` | Pa/s, **Pa**, K → m/s. Order `(omega, p, t)` | modern `omega_to_w(omega, temperature, pressure)` |
| `w_to_omega(w, p, t)` | inverse, same order | `w_to_omega` |

## Pressure unit traps

- `mixhum_ptrh`: p in **hPa**
- `mixhum_ptd`, `relhum`, `pot_temp`, `pot_temp_equiv`, NCL `omega_to_w`: p in **Pa**

```python
from pymeteo.engine import MeteoEngine
MeteoEngine().ncl_lookup("mixhum_ptrh")
```
