# pymeteo

点上气象诊断函数库，单位用字符串关键字传入（Python 3.7+）。

函数内部完成单位换算，返回普通 `float` 或 NumPy 数组。运行时只依赖 NumPy。不做 I/O、绘图、地图投影或完整 CAPE/CIN 探空套件，也不引入 Pint、MetPy 或 NCL 运行时。

[English](README.md)

## 安装

```bash
pip install pymeteo-kit
```

PyPI 发行名是 `pymeteo-kit`（`pymeteo` 已被占用）；导入仍为 `import pymeteo`。

## 快速开始

```python
import pymeteo as pm

dewpoint = pm.dewpoint_from_relative_humidity(
    18.0, 46.5, temperature_unit="C", humidity_unit="%"
)
si = pm.showalter_index(16.6, 0.6, -15.9, temperature_unit="C")
speed = pm.wind_speed(3.0, 4.0, speed_unit="m/s")
km = pm.earth_distance(39.9, 116.4, 31.2, 121.5, output_distance_unit="km")
```

## 单位

带量纲的量用字符串关键字（`temperature_unit`、`pressure_unit`、`output_*_unit` 等）。无法识别时抛出 `pymeteo.UnitError`。默认值：

| 物理量 | 默认 | 常用别名 |
|--------|------|----------|
| 温度 | `C` | `K`、`F` |
| 气压 | `hPa` | `mb`、`Pa` |
| 风 | `m/s` | `kt`、`km/h` |
| 相对湿度 | `%` | `fraction` |
| 混合比 | `kg/kg` | `g/kg` |
| 距离 | `km` | `m`、`nmi`（气层厚度默认 `m`） |
| 角度 | 度 | `rad` |

## 模块

`pymeteo.thermo`、`indices`、`wind`、`geo`、`dynamics`、`comfort` 的公开函数从包顶层再导出。参数说明见源码中的中文文档字符串。

`pymeteo.ncl` 用 NCL 内建函数名和固定单位（K、Pa/hPa、%、kg/kg）做薄封装，这些名字不会出现在 `pymeteo` 顶层。

```python
from pymeteo.ncl import dewtemp_trh

td_k = dewtemp_trh(18.0 + 273.15, 46.5)
```

## 开发

在 Python 3.8+ 的源码目录中：

```bash
pip install -e ".[dev]"
pytest
ruff check src tests
```

## 许可

GPL-3.0，见 [LICENSE](LICENSE)。
