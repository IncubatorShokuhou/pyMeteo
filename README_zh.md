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

### `pymeteo.wind` 风

| 函数 | 作用 |
|------|------|
| `wind_speed` | u、v → 风速 |
| `wind_direction` | u、v → 气象风向（来向，度） |
| `uv_from_speed_direction` / `wind_components` | 风速 + 风向 → u、v |

### `pymeteo.geo` 几何与站点

| 函数 | 作用 |
|------|------|
| `earth_distance` | WGS84 Vincenty 大地线距离 |
| `gravity` | 重力加速度 g(φ)，纬度默认按**度** |
| `sea_level_pressure` | 本站气压订正到海平面 |

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

沙氏 / K / SWEAT 等指数在 NCL 中没有与本库一一对应的同名内建函数，因此**只保留现代 API**，本模块不伪造 NCL 名字。需要灵活单位时请直接调用现代函数。

## 科学来源与相对旧代码的订正

公式意图仍来自：

- NCL 相关例程（`relhum`、`mixhum_ptrh`、`dewtemp_trh` 等）
- 李社宏. 用 C 语言开发的气象常用参数和物理量计算函数库（一）[J]. 陕西气象, 1994(03):42-45.

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
