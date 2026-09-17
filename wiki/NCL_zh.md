# NCL 兼容（`pymeteo.ncl`）

[English](NCL)

一层薄封装：沿用 NCL 内建**函数名**和**固定单位**，只翻译参数再调用现代函数。这些名字**不会**导出到 `pymeteo` 顶层（`from pymeteo import dewtemp_trh` 会失败）。需要灵活单位时请直接调现代 API。

```python
import pymeteo as pm
from pymeteo.ncl import dewtemp_trh, relhum_ttd, mixhum_ptrh, wind_speed

td_k = dewtemp_trh(18.0 + 273.15, 46.5)        # tk 为 K，rh 为 % → 露点 K（约 279.45）
td_k = pm.ncl.dewtemp_trh(18.0 + 273.15, 46.5)
rh = pm.ncl.relhum_ttd(291.15, 279.45, 0)      # opt=0 → %，opt=1 → 小数
```

沙氏 / K / SWEAT / 热指数 / 风寒在这里没有对应的 NCL 内建名，本模块不伪造。

对照 NCL 文档例题数值的回归在 `tests/test_ncl_official_examples.py`；接线与单位开关在 `tests/test_ncl.py`。

## 对照表

| NCL 名 | 固定单位 / 开关 | 调用的现代函数 |
|--------|-----------------|----------------|
| `dewtemp_trh(tk, rh)` | `tk` 为 K，`rh` 为 %，返回露点 K | `dewpoint_from_relative_humidity` |
| `relhum_ttd(t, td, opt)` | `t`/`td` 为 K；`opt=0` → %，`opt=1` → 小数 | `relative_humidity_from_dewpoint` |
| `relhum(t, w, p)` | `t` 为 K，`w` 为 kg/kg，`p` 为 Pa，返回 % | `relative_humidity_from_mixing_ratio` |
| `mixhum_ptrh(p, tk, rh, iswit)` | `p` 为 **hPa**，`tk` 为 K，`rh` 为 %；`\|iswit\|=1` 混合比、`2` 比湿；负号 → g/kg，正号 → kg/kg | `mixing_ratio_from_relative_humidity` / `specific_humidity_from_relative_humidity` |
| `mixhum_ptd(p, tdk, iswit)` | `p` 为 **Pa**，`tdk` 为 K；`iswit` 同上 | `mixing_ratio_from_dewpoint`（再转换） |
| `mixhum_convert(wq, wqType, iounit)` | `wqType` 为 `"w"` 混合比→比湿、`"q"` 相反；`iounit=(in,out)` 中 0=kg/kg、1=g/kg | `convert_humidity` |
| `vapor_pres_rh(rh, es)` | `rh` 为 %；`es` 与返回同单位 | `RH/100 · e_s` |
| `pot_temp(p, t, dim=-1, opt=False)` | `p` 为 Pa，`t` 为 K，返回 K。`dim`/`opt` 忽略 | `potential_temperature` |
| `pot_temp_equiv(p, t, w, dim=-1, humVarType="r")` | `p` 为 Pa，`t` 为 K；`humVarType`：`r`/`w` 混合比 kg/kg、`q` 比湿、`rh` 相对湿度 %。内部 Bolton 含 LCL（比 NCL 6.4 无 LCL 近似更接近 `pot_temp_equiv_tlcl`）。`dim` 忽略 | `equivalent_potential_temperature` |
| `temp_virtual(t, w, iounit)` | `iounit` 长度 3：温度 0=°C/1=K/2=°F，混合比 0=kg/kg 或 1=g/kg，输出温度。用 `T(1+r/ε)/(1+r)`，不是 NCL 文档里的 `T(1+0.61 r)` 近似 | `virtual_temperature` |
| `wetbulb_stull(t, rh, iounit, opt=False)` | `rh` 为 %；`iounit` 长度 2（0=°C、1=K、2=°F）。`opt` 未使用 | `wet_bulb_temperature` |
| `lclvl(p, tk, tdk)` | `p` 为 hPa，温度为 K，只返回 LCL 气压 hPa | `lifting_condensation_level` |
| `wind_speed(u, v)` | m/s → m/s | `wind_speed` |
| `wind_direction(u, v, opt=0)` | 来向（度）；静风：`opt=0` 为 0，`opt=1` 为 nan，其它标量为填充值 | `wind_direction` |
| `wind_component(wspd, wdir, opt=0)` | m/s + 来向 → `(u, v)` 元组。NCL 的 `opt` 未使用 | `uv_from_speed_direction` |
| `coriolis_param(lat)` | 纬度度 → s⁻¹ | `coriolis_parameter` |
| `omega_to_w(omega, p, t)` | Pa/s、**Pa**、K → m/s。顺序 `(omega, p, t)` | `omega_to_w`（现代 API 温度在气压前） |
| `w_to_omega(w, p, t)` | 上式之逆 | `w_to_omega` |

## 与测试套件一致的例子

```python
# mixhum_ptrh，1000 hPa / 18 °C / 46.5 %
mix_kg = pm.ncl.mixhum_ptrh(1000.0, 18.0 + 273.15, 46.5, 1)   # 约 0.006018 kg/kg
mix_g  = pm.ncl.mixhum_ptrh(1000.0, 18.0 + 273.15, 46.5, -1)  # 约 6.018 g/kg

# pot_temp：100000 Pa、301.25 K → 301.25 K
pm.ncl.pot_temp(100000.0, 301.25)

# wetbulb_stull：20 °C、50 %、iounit (0,0) → 约 13.70 °C
pm.ncl.wetbulb_stull(20.0, 50.0, (0, 0), False)

# lclvl：1000 hPa、15 °C、露点 4 °C → 约 847–849 hPa（Bolton 与 Stipanuk 之差）
pm.ncl.lclvl(1000.0, 15.0 + 273.15, 4.0 + 273.15)
```
