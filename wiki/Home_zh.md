# pymeteo

露点、水汽混合比、位温与相当位温、LCL、湿球、能见度、沙氏 / K / A / TT / SWEAT / 抬升指数、风速与来向、大地线距离、本站气压订正海平面、NWS 热指数与风寒、科里奥利参数、静力近似 ω ↔ w。

**Python 3.7+**。运行时只有 NumPy。单位写在调用参数里（`temperature_unit="C"`、`pressure_unit="hPa"` 等），返回 `float` 或 NumPy 数组。

```python
import pymeteo as pm

td = pm.dewpoint_from_relative_humidity(
    18.0, 46.5, temperature_unit="C", humidity_unit="%"
)
si = pm.showalter_index(16.6, 0.6, -15.9, temperature_unit="C")
u, v = pm.uv_from_speed_direction(10.0, 270.0, speed_unit="kt")
km = pm.earth_distance(39.9, 116.4, 31.2, 121.5, output_distance_unit="km")
```

默认：温度 `C`、气压 `hPa`、风 `m/s`、相对湿度 `%`、混合比 `kg/kg`、距离 `km`（气层厚度默认 `m`）、角度为度。无法识别的单位抛出 `pymeteo.UnitError`。详见 [单位](Units_zh)。

```bash
pip install pymeteo-kit
```

`import pymeteo`。NCL 内建名（固定 NCL 单位）在 [`pymeteo.ncl`](NCL_zh)，不在包根。

English: [Home](Home).

| 页面 | 内容 |
|------|------|
| [热力学](Thermo_zh) | 水汽、θ、θe、LCL、湿球、能见度 |
| [指数](Indices_zh) | 沙氏、K、A、TT、SWEAT、LI |
| [风](Wind_zh) | 风速、风向、分量、体风切变 |
| [几何与气压](Geo_zh) | 距离、重力、海平面气压、厚度 |
| [体感](Comfort_zh) | NWS 热指数、风寒 |
| [动力学](Dynamics_zh) | 科里奥利、ω ↔ w |
| [单位](Units_zh) | 别名与 `UnitError` |
| [NCL 兼容](NCL_zh) | NCL 名封装 |
| [安装](Install_zh) | `pip install pymeteo-kit` |

MIT License。[GitHub](https://github.com/IncubatorShokuhou/pyMeteo) · [PyPI](https://pypi.org/project/pymeteo-kit/)
