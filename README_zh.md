# pymeteo

轻量气象诊断函数库（Python 3.10+）。由原先单文件 `pyMeteo.py` **不兼容重写** 而来：包名改为 `pymeteo`，公开函数使用明确的英文蛇形命名，并用**字符串单位参数**在函数内部完成换算（不依赖 Pint，也不提供 Sounding / Wind 面向对象封装）。

更短的英文说明见 [README.md](README.md)。

## 安装

需要 Python ≥ 3.10 与 NumPy。运行时**只有** `numpy` 依赖，不再需要 `geopy`。

```bash
pip install -e .
```

开发（测试与 lint）：

```bash
pip install -e ".[dev]"
pytest
ruff check src tests
```

`tests/test_ncl_official_examples.py` **含 NCL 官网例题回归**（黄金值硬编码自公开文档，不依赖 MetPy / NCL / Pint）；`tests/test_ncl.py` 只测兼容层接线与单位。

## 快速开始

```python
import pymeteo as pm

# 旧 __main__ 示例：18 °C、相对湿度 46.5% → 露点约 6.30 °C
dewpoint = pm.dewpoint_from_relative_humidity(
    18.0, 46.5, temperature_unit="C", humidity_unit="%"
)

# 沙氏指数：850 hPa 温度 / 露点、500 hPa 温度
si = pm.showalter_index(16.6, 0.6, -15.9, temperature_unit="C")

# 风速风向（气象学来向，度）
speed = pm.wind_speed(3.0, 4.0, speed_unit="m/s")
u, v = pm.uv_from_speed_direction(10.0, 270.0, speed_unit="kt", output_speed_unit="m/s")

# 两点大地线距离（WGS84 Vincenty，无 geopy）
km = pm.earth_distance(39.9, 116.4, 31.2, 121.5, output_distance_unit="km")
```

## 单位约定

每个带量纲的量都有对应的单位关键字参数（或成对的输入 / 输出单位）。未写出时使用下表默认值。函数先换算到内部量纲再计算，返回普通 `float` 或 NumPy 数组。

| 物理量 | 内部量纲 | 默认公开单位 | 常用别名 |
|--------|----------|--------------|----------|
| 温度 | K | `C` | `degC` / `celsius`、`K` / `kelvin`、`F` / `fahrenheit` |
| 气压 | Pa | `hPa` | `mb` / `mbar`、`Pa`、`kPa`、`atm` |
| 风速 | m/s | `m/s` | `kt` / `knot` / `knots`、`km/h`、`mph` |
| 相对湿度 | 小数 0–1 | `%` | `percent`、`fraction` |
| 混合比 / 比湿 | kg/kg | `kg/kg` | `g/kg`、`mg/kg` |
| 距离 | m | `km` | `m` / `meter`、`mi`、`nmi` / `nm`、`ft` |
| 经纬度 / 风向 | 计算前转弧度 | `deg` | `degree`、`rad` |

无法识别的单位会抛出 `pymeteo.UnitError`。

返回值单位写在各函数的中文文档字符串中；需要换单位时使用 `output_*_unit`（如 `output_temperature_unit`、`output_pressure_unit`、`output_speed_unit`、`output_distance_unit`、`output_humidity_unit`）。

稳定度指数（沙氏、K、A、TT、SWEAT）是温差或无量纲组合：在摄氏度与开尔文下 **TT / K / A / 沙氏的数值相同**。

## 模块与 API

导入名一律为小写：

```python
from pymeteo import showalter_index, wind_speed
```

`from pymeteo import showalter` 之类的旧名**已从顶层删除**。若需要 NCL 原名与固定单位，请使用 ``pymeteo.ncl``（见下一节），不要从 ``pymeteo`` 顶层导入。

## 新增（2.2.0）

在 2.1.0 水汽 / 指数 / 风 / 几何之外，补充一批**点上或一维廓线**算法。全部按公开文献公式用 NumPy 自行实现（Bolton 1980、Tetens、Stull 2011、压高方程、NWS Rothfusz / 风寒 2001 等），**不拷贝** MetPy（BSD-3）或 NCL 源码，也**不引入** Pint 单位对象。运行时依赖仍只有 `numpy`。

刻意不做：I/O、绘图、地图投影、网格球谐平流、FAO56 辐射/蒸散全套、剖面 cross-section、完整 CAPE / CIN 探空套件。

### 现代 API（`import pymeteo as pm`）

| 函数 | 模块 | 作用 |
|------|------|------|
| `saturation_mixing_ratio` | thermo | 饱和混合比 \(w_s(p,T)\) |
| `mixing_ratio_from_dewpoint` | thermo | 气压 + 露点 → 混合比 |
| `vapor_pressure_from_mixing_ratio` | thermo | 气压 + 混合比 → 水汽压 |
| `vapor_pressure_from_relative_humidity` | thermo | 温度 + 相对湿度 → 水汽压 |
| `potential_temperature` | thermo | 位温 θ（Poisson，κ=0.286） |
| `equivalent_potential_temperature` | thermo | 相当位温 θ_e（Bolton 1980 式 43） |
| `virtual_temperature` | thermo | 虚温 \(T_v=T(1+r/\varepsilon)/(1+r)\) |
| `wet_bulb_temperature` | thermo | 湿球温度（Stull 2011，仅近海平面） |
| `lifting_condensation_level` | thermo | LCL 气压与温度（Bolton） |
| `parcel_temperature_at_pressure` | thermo | 气块干+湿绝热抬到目标气压 |
| `lifted_index` | indices | 已知 \(T_{500}\) 与气块 \(T_{500}\) 的抬升指数 |
| `lifted_index_from_surface` | indices | 由地面 \(p,T,T_d\) 抬到 500 hPa 再算 LI |
| `bulk_wind_shear` | wind | 两层风矢量差的模 |
| `coriolis_parameter` | dynamics | \(f=2\Omega\sin\phi\) |
| `height_thickness` | geo | 压高方程气层厚度 |
| `omega_to_w` / `w_to_omega` | dynamics | 静力近似 ω ↔ w（现代参数顺序为温度在气压前） |
| `heat_index` | comfort | NWS 热指数（无对应 NCL 名） |
| `wind_chill` | comfort | NWS 2001 风寒（无对应 NCL 名） |

### 新增 NCL 名（`pymeteo.ncl`，固定 NCL 单位）

| NCL 名 | NCL 单位 | 调用的现代函数 |
|--------|----------|----------------|
| `mixhum_ptd(p, tdk, iswit)` | `p` 为 Pa，`tdk` 为 K；`iswit` 同 `mixhum_ptrh` | `mixing_ratio_from_dewpoint` |
| `vapor_pres_rh(rh, es)` | `rh` 为 %，`es` 与返回值同单位 | `RH/100 · e_s` |
| `pot_temp(p, t)` | `p` 为 Pa，`t` 为 K，返回 K | `potential_temperature` |
| `pot_temp_equiv(p, t, w, dim=-1, humVarType="r")` | `p` 为 Pa，`t` 为 K；`humVarType` 为 `"r"` 混合比 kg/kg、`"q"` 比湿、`"rh"` 相对湿度 %。内部 Bolton 含 LCL | `equivalent_potential_temperature` |
| `temp_virtual(t, w, iounit)` | `iounit` 长度 3：温度 C/K/F、混合比 kg/kg 或 g/kg、输出温度 | `virtual_temperature` |
| `wetbulb_stull(t, rh, iounit, opt=False)` | `rh` 为 %；`iounit` 长度 2 指定输入/输出温度（0=°C、1=K、2=°F） | `wet_bulb_temperature` |
| `lclvl(p, tk, tdk)` | `p` 为 hPa，温度 K，返回 LCL 气压 hPa | `lifting_condensation_level` |
| `coriolis_param(lat)` | 纬度度，返回 s⁻¹ | `coriolis_parameter` |
| `omega_to_w(omega, p, t)` | ω 为 Pa/s，`p` 为 Pa，`t` 为 K，返回 m/s | `omega_to_w`（注意参数顺序与现代 API 不同） |
| `w_to_omega(w, p, t)` | 上式之逆 | `w_to_omega` |

热指数、风寒、抬升指数在 NCL 中没有与本库一一对应的同名内建函数，**不伪造** NCL 名字。

### `pymeteo.thermo` 水汽

| 函数 | 作用 |
|------|------|
| `saturation_vapor_pressure` | 水面饱和水汽压（李社宏 1994） |
| `condensation_temperature` | 抬升凝结温度 |
| `relative_humidity_from_dewpoint` | 温度 + 露点 → 相对湿度（Dutton） |
| `dewpoint_from_relative_humidity` | 温度 + 相对湿度 → 露点 |
| `relative_humidity_from_mixing_ratio` | 温度 + 混合比 + 气压 → 相对湿度（NCL 查表） |
| `mixing_ratio_from_relative_humidity` | 气压 + 温度 + 相对湿度 → 混合比（Tetens） |
| `specific_humidity_from_relative_humidity` | 同上，返回比湿 |
| `convert_humidity` | 混合比 ↔ 比湿，并换算 `kg/kg` / `g/kg` |
| `visibility` | RUC / FSL 能见度估算 |
| `saturation_mixing_ratio` | 饱和混合比（2.2.0） |
| `mixing_ratio_from_dewpoint` | 气压 + 露点 → 混合比（2.2.0） |
| `vapor_pressure_from_mixing_ratio` | 气压 + 混合比 → 水汽压（2.2.0） |
| `vapor_pressure_from_relative_humidity` | 温度 + 相对湿度 → 水汽压（2.2.0） |
| `potential_temperature` | 位温 θ（2.2.0） |
| `equivalent_potential_temperature` | 相当位温 θ_e，Bolton 1980（2.2.0） |
| `virtual_temperature` | 虚温（2.2.0） |
| `wet_bulb_temperature` | Stull 湿球温度（2.2.0） |
| `lifting_condensation_level` | LCL 气压与温度（2.2.0） |
| `parcel_temperature_at_pressure` | 气块抬升到目标气压（2.2.0） |

### `pymeteo.indices` 指数

| 函数 | 作用 |
|------|------|
| `showalter_index` | 沙氏指数 |
| `k_index` | K 指数 |
| `a_index` | A 指数 |
| `total_totals_index` | 全总指数 TT（可给 850 hPa 露点或相对湿度） |
| `sweat_index` | SWEAT（露点用 °C，风速用节；见下方订正说明） |
| `temperature_dewpoint_depression` | 温度露点差（替代旧 `ttd850` 等） |
| `layer_temperature_difference` | 两层温度差（替代旧 `tt500`） |
| `lifted_index` / `lifted_index_from_surface` | 抬升指数 LI（2.2.0） |

### `pymeteo.wind` 风

| 函数 | 作用 |
|------|------|
| `wind_speed` | u、v → 风速 |
| `wind_direction` | u、v → 气象风向（来向，度） |
| `uv_from_speed_direction` / `wind_components` | 风速 + 风向 → u、v |
| `bulk_wind_shear` | 两层风矢量差模（2.2.0） |

### `pymeteo.dynamics` 轻量动力学（2.2.0）

| 函数 | 作用 |
|------|------|
| `coriolis_parameter` | 科里奥利参数 f(φ) |
| `omega_to_w` / `w_to_omega` | 静力近似垂直速度换算 |

### `pymeteo.comfort` 体感（2.2.0）

| 函数 | 作用 |
|------|------|
| `heat_index` | NWS 热指数 |
| `wind_chill` | NWS 2001 风寒 |

### `pymeteo.geo` 几何与站点

| 函数 | 作用 |
|------|------|
| `earth_distance` | WGS84 Vincenty 大地线距离 |
| `gravity` | 重力加速度 g(φ)，纬度默认按**度** |
| `sea_level_pressure` | 本站气压订正到海平面 |
| `height_thickness` | 压高方程气层厚度（2.2.0） |

完整参数说明见源码中的中文文档字符串。

## NCL 兼容层（`pymeteo.ncl`）

部分用户习惯 NCL 内建函数名与**固定单位**。`pymeteo.ncl` 提供一层薄封装：只翻译参数和单位，再调用上表中的现代函数，**不重复实现物理公式**。NCL 名字**不会**再导出到 `pymeteo` 顶层，以免污染现代 API。

```python
import pymeteo as pm
from pymeteo.ncl import dewtemp_trh, relhum_ttd, mixhum_ptrh, wind_speed

# 温度开尔文、相对湿度百分数 → 露点开尔文
td_k = dewtemp_trh(18.0 + 273.15, 46.5)
td_k = pm.ncl.dewtemp_trh(18.0 + 273.15, 46.5)

# opt=0 → 百分数；opt=1 → 0–1 小数
rh = pm.ncl.relhum_ttd(18.0 + 273.15, 6.3 + 273.15, 0)
```

| NCL 名 | NCL 单位 / 开关 | 调用的现代函数 |
|--------|-----------------|----------------|
| `dewtemp_trh(tk, rh)` | `tk` 为 K，`rh` 为 %，返回露点 K | `dewpoint_from_relative_humidity` |
| `relhum_ttd(t, td, opt)` | `t`/`td` 为 K；`opt=0` → %，`opt=1` → 小数 | `relative_humidity_from_dewpoint` |
| `relhum(t, w, p)` | `t` 为 K，`w` 为 kg/kg，`p` 为 Pa，返回 % | `relative_humidity_from_mixing_ratio` |
| `mixhum_ptrh(p, tk, rh, iswit)` | `p` 为 hPa，`tk` 为 K，`rh` 为 %；`iswit` ±1 混合比、±2 比湿；负号 → g/kg，正号 → kg/kg | `mixing_ratio_from_relative_humidity` / `specific_humidity_from_relative_humidity` |
| `mixhum_convert(wq, wqType, iounit)` | `wqType` 为 `"w"` 混合比→比湿、`"q"` 相反；`iounit=(in,out)` 中 0=kg/kg、1=g/kg | `convert_humidity` |
| `wind_speed(u, v)` | m/s → m/s | `wind_speed` |
| `wind_direction(u, v, opt=0)` | 气象学来向（度）；静风时 `opt=0` 为 0，`opt=1` 为 nan | `wind_direction` |
| `wind_component(wspd, wdir, opt=0)` | 风速 + 来向 → `(u, v)`（Python 返回元组；NCL 的 `opt` 未使用） | `uv_from_speed_direction` |
| `mixhum_ptd(p, tdk, iswit)` | `p` 为 Pa，`tdk` 为 K；`iswit` 同 `mixhum_ptrh` | `mixing_ratio_from_dewpoint` |
| `vapor_pres_rh(rh, es)` | `rh` 为 %，`es` 与返回同单位 | `RH/100 · e_s` |
| `pot_temp(p, t)` | `p` 为 Pa，`t` 为 K，返回 K | `potential_temperature` |
| `pot_temp_equiv(p, t, w, dim=-1, humVarType="r")` | `p` 为 Pa，`t` 为 K；`humVarType`：`r` 混合比、`q` 比湿、`rh` 相对湿度 % | `equivalent_potential_temperature` |
| `temp_virtual(t, w, iounit)` | `iounit` 长度 3，见上文新增节 | `virtual_temperature` |
| `wetbulb_stull(t, rh, iounit, opt=False)` | `rh` 为 %；`iounit` 长度 2（0=°C、1=K、2=°F） | `wet_bulb_temperature` |
| `lclvl(p, tk, tdk)` | `p` 为 hPa，温度 K，返回 LCL 气压 | `lifting_condensation_level` |
| `coriolis_param(lat)` | 纬度度 → s⁻¹ | `coriolis_parameter` |
| `omega_to_w(omega, p, t)` | Pa/s、Pa、K → m/s | `omega_to_w` |
| `w_to_omega(w, p, t)` | m/s、Pa、K → Pa/s | `w_to_omega` |

沙氏 / K / SWEAT 等指数在 NCL 中没有与本库一一对应的同名内建函数，因此**只保留现代 API**，本模块不伪造 NCL 名字。需要灵活单位时请直接调用现代函数。

数值回归**含 NCL 官网例题回归**：`tests/test_ncl_official_examples.py` 对照 dewtemp_trh、mixhum_ptrh、pot_temp、wetbulb_stull、lclvl 文档打印值；接线测试仍在 `tests/test_ncl.py`。

## 科学来源与相对旧代码的订正

公式意图仍来自：

- NCL 相关例程（`relhum`、`mixhum_ptrh`、`dewtemp_trh` 等）
- 李社宏. 用 C 语言开发的气象常用参数和物理量计算函数库（一）[J]. 陕西气象, 1994(03):42-45.
- Bolton, D., 1980: The computation of equivalent potential temperature. *Mon. Wea. Rev.*, 108, 1046–1053.
- Stull, R., 2011: Wet-bulb temperature from relative humidity and air temperature. *J. Appl. Meteor. Climatol.*, 50, 2267–2269.
- Rothfusz, L. P., 1990: The heat index equation. NWS Technical Attachment SR 90-23.
- NWS / Environment Canada, 2001: 风寒公式。

2.2.0 新增算法均为按上述文献**重新实现**，未粘贴 MetPy 或 NCL 源码。2.2.1 增加 NCL 官网例题回归测试，公式未改。

在保持上述公式意图的前提下，重写时修正了若干会误导结果的问题：

1. **SWEAT**：标准形式为 `12·Td850(°C) + 20·(TT−49) + 2·f850(kt) + f500(kt) + 125·(S+0.2)`。旧代码把开尔文露点直接乘 12，并用 m/s 风速去套接近“两倍系数”的 4 与 2，结果可偏大约一个量级。新实现按 NWS / Miller (1972) 使用 °C 与节；负项置零；切变项还要求风向差为正且两层风速 ≥ 15 kt。
2. **`earth_distance`**：去掉 `geopy` 与 `eval("a."+unit)`，改为 Vincenty。
3. **`gravity`**：纬度按度输入并先化为弧度。旧代码把度数直接交给 `sin`。
4. **能见度 RUC**：按 60 km × 指数衰减给出千米量级，再换算到 `output_distance_unit`，不再把千米结果误标成再乘 1000。
5. **数组**：闭合公式用 NumPy 广播；迭代型沙氏指数 / 凝结温度按元素迭代并设最大步数，避免无限循环。
6. **混合比 / 比湿**：用 `from_quantity` / `to_quantity` 与单位字符串，取代含义与文档互相矛盾的整型 `ISWIT` / `wqType`。

## 从旧单文件 API 迁移

| 旧（`pyMeteo.py`） | 新 |
|--------------------|----|
| `showalter` | `showalter_index(..., temperature_unit="C")` |
| `E_WATER` | `saturation_vapor_pressure` |
| `Tc` | `condensation_temperature` |
| `K` / `A` / `TT` | `k_index` / `a_index` / `total_totals_index` |
| `dewtemp_trh` / `relhum_ttd` | `dewpoint_from_relative_humidity` / `relative_humidity_from_dewpoint` |
| `relhum` | `relative_humidity_from_mixing_ratio` |
| `mixhum_ptrh` | `mixing_ratio_from_relative_humidity` 或 `specific_humidity_from_relative_humidity` |
| `mixhum_convert` | `convert_humidity` |
| `ws` / `wd` / `u` / `v` | `wind_speed` / `wind_direction` / `uv_from_speed_direction` |
| `SWEAT_calculate` | `sweat_index`（单位语义已修正，数值不可与旧结果逐点对比） |
| `earth_distance` | 同名，但无 geopy；用 `output_distance_unit` |
| `g` | `gravity(..., latitude_unit="deg")` |
| `ttd850` 等 | `temperature_dewpoint_depression` |

旧文件 `pyMeteo.py` 已删除。`import pyMeteo` 不再可用。

## 许可

GNU General Public License v3.0，见 [LICENSE](LICENSE)。
