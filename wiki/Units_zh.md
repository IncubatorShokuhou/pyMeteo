# 单位

[English](Units.md)

带量纲的参数用字符串声明单位（有时是一对输入/输出单位）。先换到内部接近 SI 的量，再换回来。标量仍是 Python `float`，数组仍是数组。

无法识别或不适用于该物理量的字符串会抛 `pymeteo.UnitError`，它是 `ValueError` 的子类。

不用 Pint。

## 默认值

| 物理量 | 内部 | 常见默认 | 常用别名 |
|--------|------|----------|----------|
| 温度 | K | `C` | `degC`、`celsius`、`K`、`kelvin`、`F`、`fahrenheit` |
| 气压 | Pa | `hPa` | `mb`、`mbar`、`Pa`、`kPa`、`atm` |
| 风速 | m/s | `m/s` | `kt`、`knots`、`km/h`、`mph` |
| 相对湿度 | 0–1 | `%` | `percent`、`fraction` |
| 混合比 / 比湿 | kg/kg | `kg/kg` | `g/kg`、`mg/kg` |
| 距离 | m | `km` | `m`、`mi`、`ft`、`nmi` |
| 角度 | 三角函数前转弧度 | `deg` | `degree`、`rad` |

气层厚度（`height_thickness`）默认是 **`m`**，不是 `km`。位温输出默认是 **K**。

沙氏、K、A、TT、抬升指数是温度差。`C` 和 `K` 得到同一个数。

输出参数长这样：`output_temperature_unit`、`output_pressure_unit`、`output_speed_unit`、`output_distance_unit`、`output_humidity_unit`、`output_omega_unit`。

## `pymeteo.units` 里的辅助函数

给内部用的。包根只再导出 `UnitError`。

规范化：`canonical_temperature_unit`、`canonical_pressure_unit`、`canonical_speed_unit`、`canonical_rh_unit`、`canonical_mass_humidity_unit`、`canonical_distance_unit`、`canonical_angle_unit`。

换算：`to_kelvin` / `from_kelvin`、`to_pascal` / `from_pascal`、`to_mps` / `from_mps`、`to_rh_fraction` / `from_rh_fraction`、`to_kgkg` / `from_kgkg`、`to_meters` / `from_meters`、`to_radians`。

节是 `1852/3600` m/s。法定英里 1609.344 m。海里 1852 m。标准大气压 101325 Pa。

```python
from pymeteo.units import to_kelvin, from_kelvin

k = to_kelvin(18.0, "C")       # 291.15
c = from_kelvin(k, "celsius")  # 18.0
```
