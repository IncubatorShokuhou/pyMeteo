# 单位

[English](Units.md)

每个带量纲的量都有对应的字符串单位关键字（或成对的输入 / 输出单位）。函数先换算到内部量纲再计算，返回普通 `float` 或 NumPy 数组（全部原始输入都是标量时，`restore_shape` 会还原成 Python `float`）。

无法识别或不适用于当前物理量的单位会抛出 **`pymeteo.UnitError`**（`ValueError` 子类）。

## 默认值与内部量纲

| 物理量 | 内部量纲 | 默认公开单位 | 常用别名（不完全） |
|--------|----------|--------------|--------------------|
| 温度 | K | `C` | `degC` / `celsius` / `centigrade`、`K` / `kelvin`、`F` / `degF` / `fahrenheit` |
| 气压 | Pa | `hPa` | `mb` / `mbar` / `millibar`、`Pa` / `pascal`、`kPa`、`atm` |
| 风速 | m/s | `m/s` | `ms-1` / `mps`、`kt` / `knot` / `knots` / `kn`、`km/h` / `kmh` / `kph`、`mph` |
| 相对湿度 | 小数 0–1 | `%` | `percent` / `percentage`、`fraction` / `1` / `ratio` |
| 混合比 / 比湿 | kg/kg | `kg/kg` | `g/kg`、`mg/kg`、`kgkg-1`、`g/g` |
| 距离 | m | `km` | `m` / `meter`、`cm`、`mi` / `mile`、`ft` / `foot`、`nmi` / `nm` / `nauticalmile` |
| 角度（经纬度、风向来向） | 三角函数前转弧度 | `deg` | `degree` / `degrees`、`rad` / `radian` |

气层**厚度**（`height_thickness`）默认是 **`m`**，不是 `km`。

稳定度指数（沙氏、K、A、TT、SWEAT）是温差或无量纲组合：在摄氏度与开尔文下 **TT / K / A / 沙氏的数值相同**。

输出单位关键字形如 `output_temperature_unit`、`output_pressure_unit`、`output_speed_unit`、`output_distance_unit`、`output_humidity_unit`、`output_omega_unit`。

## 换算函数（`pymeteo.units`）

这些函数给内部用。除 `UnitError` 外，不从包根再导出。

规范化：`canonical_temperature_unit`、`canonical_pressure_unit`、`canonical_speed_unit`、`canonical_rh_unit`、`canonical_mass_humidity_unit`、`canonical_distance_unit`、`canonical_angle_unit`。

换算：`to_kelvin` / `from_kelvin`、`to_pascal` / `from_pascal`、`to_mps` / `from_mps`、`to_rh_fraction` / `from_rh_fraction`、`to_kgkg` / `from_kgkg`、`to_meters` / `from_meters`、`to_radians`。

节用 `1852/3600` m/s。法定英里 `1609.344` m。海里 `1852` m。标准大气压 `101325` Pa。

```python
from pymeteo.units import UnitError, to_kelvin, from_kelvin

k = to_kelvin(18.0, "C")          # 291.15
c = from_kelvin(k, "celsius")     # 18.0
```
