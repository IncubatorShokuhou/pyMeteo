# pymeteo

pymeteo 是一组在 NumPy 上做点诊断的气象计算。

不读资料、不绘图，也不包一层完整的探空对象。传入标量或数组，返回标量或数组。

## 安装

```bash
pip install pymeteo-kit
```

```python
import pymeteo as pm
```

Python 3.7+。导入名是 `pymeteo`。

## 例子

```python
import pymeteo as pm

td = pm.dewpoint_from_relative_humidity(18.0, 46.5)
# 大约 6.3 °C
```

单位写在调用参数里（`temperature_unit="C"`）。默认值在各函数的文档字符串中。

## 页面

- [热力学](Thermo_zh.md) — 水汽、位温、LCL、湿球
- [指数](Indices_zh.md) — 沙氏、K、TT、SWEAT、抬升指数
- [风](Wind_zh.md) — 风速、来向、uv、体风切变
- [几何与气压](Geo_zh.md) — 大圆距离、重力、海平面气压
- [体感](Comfort_zh.md) — 热指数、风寒
- [动力学](Dynamics_zh.md) — 科氏参数、ω ↔ w
- [单位](Units_zh.md) — 别名与 `UnitError`
- [NCL](NCL_zh.md) — NCL 函数名和它的单位约定
- [安装](Install_zh.md)

English: [Home](Home.md)

MIT License。[GitHub](https://github.com/IncubatorShokuhou/pyMeteo) · [PyPI](https://pypi.org/project/pymeteo-kit/)
