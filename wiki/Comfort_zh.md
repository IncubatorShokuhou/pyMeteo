# 体感（`pymeteo.comfort`）

[English](Comfort.md)

NWS 业务回归式，自行实现（不拷贝 MetPy）。内部先换到华氏度与英里每小时再算，再换回调用方单位。**没有 NCL 名封装**——NCL 无对应内建函数，`pymeteo.ncl` 也不伪造。

---

### `heat_index(temperature, relative_humidity, *, temperature_unit="C", humidity_unit="%", output_temperature_unit=None)`

Rothfusz（1990）/ NWS SR 90-23，基于 Steadman（1979）：

1. 先算简化式 `HI = 0.5 {T + 61 + (T−68)·1.2 + RH·0.094}`（`T` 为 °F，`RH` 为百分数），再与气温平均；
2. 若该均值 ≥ 80 °F，改用 Rothfusz 全式，并在低湿（RH<13%、80–112 °F）或高湿（RH>85%、80–87 °F）时加减订正项。

Steadman 原表范围外不可靠。90 °F、60% 时约 100 °F。

```python
pm.heat_index(90.0, 60.0, temperature_unit="F", output_temperature_unit="F")
```

输出单位默认与 `temperature_unit` 相同。

### `wind_chill(temperature, wind_speed, *, temperature_unit="C", speed_unit="m/s", output_temperature_unit=None)`

NWS / Environment Canada 2001：

`WC = 35.74 + 0.6215 T − 35.75 V^0.16 + 0.4275 T V^0.16`

（`T` 为 °F，`V` 为 mph；业务公式针对约 10 m 高度的风）。有效范围约 `T ≤ 50 °F` 且 `V ≥ 3 mph`；超出范围仍返回公式值。0 °F、10 mph 时约 −16 °F。

```python
pm.wind_chill(0.0, 10.0, temperature_unit="F", speed_unit="mph", output_temperature_unit="F")
```
