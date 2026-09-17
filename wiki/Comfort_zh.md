# 体感

[English](Comfort.md)

两个 NWS 业务回归：热指数和风寒。

内部先换成华氏度和英里每小时，再换回你要的单位。NCL 没有对应内建名，所以 `pymeteo.ncl` 里也没有这两项。

## heat_index

Rothfusz（1990）/ NWS SR 90-23，源自 Steadman（1979）。

先用一个和气温平均的简化式。若该均值 ≥ 80 °F，改用 Rothfusz 全式，并加上常见的低湿、高湿订正。

超出 Steadman 原表范围（又热又湿的极端）时，公式不可靠。输出单位默认跟 `temperature_unit` 相同。

```python
import pymeteo as pm
pm.heat_index(90.0, 60.0, temperature_unit="F")  # 大约 100 °F
```

## wind_chill

NWS / Environment Canada 2001，针对大约 10 m 高度的风：

`WC = 35.74 + 0.6215 T − 35.75 V^0.16 + 0.4275 T V^0.16`

（`T` 为 °F，`V` 为 mph。）大致适用于 `T ≤ 50 °F` 且 `V ≥ 3 mph`。超出范围仍会返回公式值，只是物理意义弱一些。

```python
pm.wind_chill(0.0, 10.0, temperature_unit="F", speed_unit="mph")  # 大约 −16 °F
```
