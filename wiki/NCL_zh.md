# NCL

[English](NCL.md)

一层薄封装：保住 NCL 内建**函数名**和 NCL 的**固定单位**。它们只翻译参数，再调用现代接口。不会从顶层 `pymeteo` 再导出。

需要灵活单位时，直接用现代函数。

```python
from pymeteo.ncl import dewtemp_trh

td_k = dewtemp_trh(18.0 + 273.15, 46.5)  # tk 为 K，rh 为 %，露点为 K
```

沙氏、K、SWEAT、热指数、风寒这里没有对应的 NCL 内建名，也不会去编。

NCL 文档例题的数字钉在 `tests/test_ncl_official_examples.py`。

## 对照

| NCL 名 | 单位 / 开关 | 转调 |
|--------|-------------|------|
| `dewtemp_trh(tk, rh)` | `tk` K，`rh` %，露点 K | `dewpoint_from_relative_humidity` |
| `relhum_ttd(t, td, opt)` | `t`/`td` K；`opt=0` → %，`opt=1` → 小数 | `relative_humidity_from_dewpoint` |
| `relhum(t, w, p)` | `t` K，`w` kg/kg，`p` Pa，结果 % | `relative_humidity_from_mixing_ratio` |
| `mixhum_ptrh(p, tk, rh, iswit)` | `p` 为 **hPa**，`tk` K，`rh` %；`\|iswit\|=1` 混合比，`2` 比湿；负号 → g/kg | 由 RH 求混合比 / 比湿 |
| `mixhum_ptd(p, tdk, iswit)` | `p` 为 **Pa**，`tdk` K；`iswit` 同上 | `mixing_ratio_from_dewpoint` |
| `mixhum_convert(wq, wqType, iounit)` | `"w"` 混合比→比湿，`"q"` 反过来；`iounit=(in,out)` 0=kg/kg，1=g/kg | `convert_humidity` |
| `vapor_pres_rh(rh, es)` | `rh` %；`es` 与返回值单位相同 | `RH/100 · e_s` |
| `pot_temp(p, t, dim=-1, opt=False)` | `p` Pa，`t` K → K。`dim`/`opt` 忽略 | `potential_temperature` |
| `pot_temp_equiv(p, t, w, dim=-1, humVarType="r")` | `p` Pa，`t` K；`"r"`/`"w"` 混合比 kg/kg，`"q"` 比湿，`"rh"` RH %。内部仍用带 LCL 的 Bolton（比 NCL 6.4 无 LCL 近似更接近 `pot_temp_equiv_tlcl`）。`dim` 忽略 | `equivalent_potential_temperature` |
| `temp_virtual(t, w, iounit)` | `iounit` 长度 3：温度 0=°C/1=K/2=°F，混合比 0=kg/kg 或 1=g/kg，输出温度。用 `T(1+r/ε)/(1+r)`，不用 `T(1+0.61 r)` | `virtual_temperature` |
| `wetbulb_stull(t, rh, iounit, opt=False)` | `rh` %；`iounit` 长度 2（0=°C，1=K，2=°F）。`opt` 未使用 | `wet_bulb_temperature` |
| `lclvl(p, tk, tdk)` | `p` hPa，温度为 K，只返回 LCL 气压 | `lifting_condensation_level` |
| `wind_speed(u, v)` | m/s | `wind_speed` |
| `wind_direction(u, v, opt=0)` | 来向，度；静风：`opt=0` → 0，`opt=1` → nan，其它标量当填充值 | `wind_direction` |
| `wind_component(wspd, wdir, opt=0)` | m/s、来向 → `(u, v)` 元组。NCL 的 `opt` 未使用 | `uv_from_speed_direction` |
| `coriolis_param(lat)` | 纬度（度）→ s⁻¹ | `coriolis_parameter` |
| `omega_to_w(omega, p, t)` | Pa/s，**Pa**，K → m/s。顺序 `(omega, p, t)` | 现代 `omega_to_w` 是 `(omega, temperature, pressure)` |
| `w_to_omega(w, p, t)` | 逆变换 | `w_to_omega` |

```python
import pymeteo as pm

pm.ncl.mixhum_ptrh(1000.0, 18.0 + 273.15, 46.5, 1)   # 约 0.00602 kg/kg
pm.ncl.pot_temp(100000.0, 301.25)                     # 301.25 K
pm.ncl.wetbulb_stull(20.0, 50.0, (0, 0), False)       # 约 13.7 °C
pm.ncl.lclvl(1000.0, 15.0 + 273.15, 4.0 + 273.15)     # 约 847 hPa
```
