# 指数

[English](Indices.md)

用几层规定高度上的温度（SWEAT 还要风）算稳定度指数。

沙氏、K、A、TT、抬升指数用 `C` 或 `K` 传入，数值一样。SWEAT 内部按摄氏度和节来定义。这些指数没有 NCL 名封装——NCL 本身就没有对应内建函数，`pymeteo.ncl` 也不会去编。

## temperature_dewpoint_depression

`T − Td`。输出用华氏度时，差值乘 9/5。

```python
import pymeteo as pm
pm.temperature_dewpoint_depression(7.0, -2.0)  # 9.0
```

## layer_temperature_difference

下层减上层。常见用法是 850 hPa 减 500 hPa。

```python
pm.layer_temperature_difference(16.6, -15.9)  # 32.5
```

## k_index

`K = T850 − T500 + Td850 − (T700 − Td700)`，按摄氏度组合。

```python
pm.k_index(16.6, 0.6, 7.0, -2.0, -15.9)  # 24.1
```

文档字符串里的经验（不是预报结论）：K < 20 基本无雷暴；20–25 零星；25–30 分散；30–35 成片。

## a_index

`A = (T850 − T500) − (T850 − Td850) − (T700 − Td700) − (T500 − Td500)`。还需要 Td500。

```python
pm.a_index(16.6, 0.6, 7.0, -2.0, -15.9, -20.0)  # 3.4
```

## total_totals_index

`TT = T850 + Td850 − 2·T500`。`dewpoint_850` 和 `relative_humidity_850` 只能给一个。湿度走 Dutton 公式反推露点。

```python
pm.total_totals_index(18.0, -15.9, dewpoint_850=6.3)
```

## showalter_index

从 850 hPa 抬到 500 hPa 的气块。SI = T500 − T_parcel(500)。湿绝热段是李社宏（1994）。

```python
pm.showalter_index(16.6, 0.6, -15.9)  # 大约 1.1
```

网格插值出负混合比时迭代会难受。那是资料问题，不是单位问题。

## sweat_index

Miller（1972）/ NWS：

`12·Td850(°C) + 20·(TT−49) + 2·f850(kt) + f500(kt) + 125·(S+0.2)`

负的项置零。切变项只在这些条件同时成立时保留：850 风向 130–250°、500 风向 210–310°、风向差为正、两层风速都 ≥ 15 kt。

Td850 **或** RH850 二选一。风分量默认 `m/s`，内部换成节。文档里的经验：>300 有强对流潜势，>400 有龙卷潜势。

## lifted_index

已经有 500 hPa 气块温度时：`LI = T500 − T_parcel(500)`。负值表示气块比环境暖。

## lifted_index_from_surface

近地层气压、温度、露点，再加上 T500。气块先按 Bolton 求 LCL，再按李社宏湿熵抬到 500 hPa。目标层固定是 500 hPa；`pressure_unit` 只说明你传入的地面气压用什么单位。

这是常用的地面抬升指数。CAPE、最不稳定气块、混合层气块都不在这个库里。
