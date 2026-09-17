# pymeteo

点上气象诊断函数（**Python 3.7+**），使用**字符串单位参数**。运行时依赖只有 **NumPy**。导入名：`pymeteo`。PyPI 发行名：[`pymeteo-kit`](https://pypi.org/project/pymeteo-kit/)。

不提供 Pint 单位对象、Sounding / Wind 类封装、MetPy 或 NCL 运行时。本库不是 I/O、绘图、地图投影或完整 CAPE/CIN 探空套件。

```python
import pymeteo as pm

td = pm.dewpoint_from_relative_humidity(18.0, 46.5, temperature_unit="C", humidity_unit="%")
si = pm.showalter_index(16.6, 0.6, -15.9, temperature_unit="C")
km = pm.earth_distance(39.9, 116.4, 31.2, 121.5, output_distance_unit="km")
```

默认公开单位：温度 `C`、气压 `hPa`、风 `m/s`、相对湿度 `%`、混合比 `kg/kg`、距离 `km`（气层厚度默认 `m`）、角度为度。无法识别的单位字符串会抛出 `pymeteo.UnitError`。

顶层是英文蛇形名。NCL 内建名在 [`pymeteo.ncl`](NCL_zh)，**不会**再导出到包根。

English: [Home](Home).

## 页面

| 页面 | 模块 | 内容 |
|------|------|------|
| [安装](Install_zh) | — | `pip install pymeteo-kit`，Python 3.7+，开发 |
| [单位](Units_zh) | `pymeteo.units` | 别名、内部量纲、`UnitError` |
| [热力学](Thermo_zh) | `pymeteo.thermo` | 水汽、θ、θe、LCL、湿球、能见度 |
| [指数](Indices_zh) | `pymeteo.indices` | 沙氏、K、A、TT、SWEAT、LI |
| [风](Wind_zh) | `pymeteo.wind` | 风速、风向、分量、体风切变 |
| [几何与气压](Geo_zh) | `pymeteo.geo` | 距离、重力、海平面气压、厚度 |
| [体感](Comfort_zh) | `pymeteo.comfort` | NWS 热指数、风寒 |
| [动力学](Dynamics_zh) | `pymeteo.dynamics` | 科里奥利、ω ↔ w |
| [NCL 兼容](NCL_zh) | `pymeteo.ncl` | 固定单位的名字封装 |
| [发布与贡献](Publishing_zh) | — | PyPI 发布、测试、贡献 |

## 顶层公开名

`import pymeteo as pm` 再导出（见 `pymeteo.__all__`）：

`UnitError`、`__version__`、`ncl`，以及上表各模块页面列出的现代函数。

## 范围（不做的事）

不做 FAO56 辐射/蒸散全套、网格球谐平流、剖面 cross-section、完整 CAPE / CIN 探空套件。热指数与风寒**没有**对应 NCL 名；`ncl` 模块不伪造这些名字。

## 链接

* 仓库：https://github.com/IncubatorShokuhou/pyMeteo
* 站点：https://incubatorshokuhou.github.io/pyMeteo/zh/
* PyPI：https://pypi.org/project/pymeteo-kit/
* 许可：MIT（见源码树中的 `LICENSE`）
