# pymeteo

一组用 NumPy 算的气象诊断量：露点、相对湿度、位温、稳定度指数、风速风向、大圆距离等。
适合探空点算、批量网格点诊断；不做读资料、不做绘图、不做完整探空对象模型。

[English](README.md)

## 安装

Python 3.7+（PyPI 包名是 `pymeteo-kit`）：

```bash
pip install pymeteo-kit
```

```python
import pymeteo as pm
```

## 例子

```python
import pymeteo as pm

td = pm.dewpoint_from_relative_humidity(
    18.0, 46.5, temperature_unit="C", humidity_unit="%"
)
si = pm.showalter_index(16.6, 0.6, -15.9, temperature_unit="C")
ws = pm.wind_speed(3.0, 4.0)
km = pm.earth_distance(39.9, 116.4, 31.2, 121.5)
```

## 给 agent 用

可以用随包装的 `meteo-expert` skill，或可选的 MCP 服务，来选函数、核对单位、对照 NCL 名字。原来的 `import pymeteo as pm` 调用方式不变。

```bash
pip install pymeteo-kit
pymeteo install-skill              # Claude Code：~/.claude/skills/meteo-expert
pymeteo install-skill --project    # Codex：./skills/meteo-expert
pip install "pymeteo-kit[mcp]"
pymeteo mcp serve                  # 或：python -m pymeteo.mcp_server
pymeteo info
```

## 单位

函数用字符串参数声明单位（如 `temperature_unit="C"`），在函数内换算。默认单位写在各函数的文档字符串里。不引入 Pint。

## 主要接口

公开函数从包顶层导入。参数说明见源码文档字符串。

| 模块 | 内容 |
|------|------|
| `pymeteo.thermo` | 露点、相对湿度、混合比、饱和水汽压、位温、相当位温、LCL、湿球、虚温 |
| `pymeteo.indices` | 沙氏指数、K / A / TT、SWEAT、抬升指数 |
| `pymeteo.wind` | 风速、风向、uv 分量、风切变 |
| `pymeteo.geo` | 大圆距离、重力、海平面气压、气层厚度 |
| `pymeteo.dynamics` | 科氏参数、ω ↔ w |
| `pymeteo.comfort` | 热指数、风寒 |

## NCL 名字

一部分 NCL 内建名在 `pymeteo.ncl` 里，单位按 NCL 约定（温度 K、湿度 % 等）：

```python
from pymeteo.ncl import dewtemp_trh

td_k = dewtemp_trh(18.0 + 273.15, 46.5)
```

## 许可

MIT License，见 [LICENSE](LICENSE)。
