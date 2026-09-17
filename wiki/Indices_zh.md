# 指数（`pymeteo.indices`）

[English](Indices)

稳定度指数。沙氏、K、A、TT 无论温度用 `C` 还是 `K` 传入，数值相同。SWEAT 在内部按 °C 与节定义。

这些指数在 `pymeteo.ncl` 里**没有**对应的 NCL 内建名。

---

### `temperature_dewpoint_depression(temperature, dewpoint, *, temperature_unit="C", output_temperature_unit="C")`

`T − Td`。选用 `F` 时按华氏度差（乘 9/5）返回。

```python
pm.temperature_dewpoint_depression(7.0, -2.0)  # 9.0
```

### `layer_temperature_difference(temperature_lower, temperature_upper, *, temperature_unit="C", output_temperature_unit="C")`

下层减上层。

```python
pm.layer_temperature_difference(16.6, -15.9)  # 32.5
```

### `k_index(temperature_850, dewpoint_850, temperature_700, dewpoint_700, temperature_500, *, temperature_unit="C")`

`K = T850 − T500 + Td850 − (T700 − Td700)`，按摄氏度计算。

```python
pm.k_index(16.6, 0.6, 7.0, -2.0, -15.9)  # 24.1
```

文档字符串中的经验：K<20 无雷暴；20–25 零星；25–30 分散；30–35 成片。

### `a_index(temperature_850, dewpoint_850, temperature_700, dewpoint_700, temperature_500, dewpoint_500, *, temperature_unit="C")`

`A = (T850 − T500) − (T850 − Td850) − (T700 − Td700) − (T500 − Td500)`。

```python
pm.a_index(16.6, 0.6, 7.0, -2.0, -15.9, -20.0)  # 3.4
```

### `total_totals_index(temperature_850, temperature_500, *, dewpoint_850=None, relative_humidity_850=None, temperature_unit="C", humidity_unit="%")`

`TT = T850 + Td850 − 2·T500`。`dewpoint_850` 与 `relative_humidity_850` **必须提供其一**。相对湿度按 Dutton 公式反推露点。

```python
pm.total_totals_index(18.0, -15.9, relative_humidity_850=46.5)
```

### `showalter_index(temperature_850, dewpoint_850, temperature_500, *, temperature_unit="C")`

850 hPa 气块抬到 500 hPa；SI = T500 − T_parcel(500)。按元素迭代并设最大步数。

```python
pm.showalter_index(16.6, 0.6, -15.9)  # 约 1.1
```

### `sweat_index(temperature_850, temperature_500, u_850, v_850, u_500, v_500, *, dewpoint_850=None, relative_humidity_850=None, temperature_unit="C", humidity_unit="%", speed_unit="m/s")`

NWS / Miller（1972）：

`12·Td850(°C) + 20·(TT−49) + 2·f850(kt) + f500(kt) + 125·(S+0.2)`

负项置零。切变项还要求风向差为正且两层风速 ≥ 15 kt。提供 Td850 **或** RH850。风分量默认 `m/s`，内部换成节。

### `lifted_index(temperature_500, parcel_temperature_500, *, temperature_unit="C")`

已知 500 hPa 气块温度时：`LI = T500 − T_parcel(500)`。

### `lifted_index_from_surface(pressure, temperature, dewpoint, temperature_500, *, pressure_unit="hPa", temperature_unit="C")`

把近地层气块抬到 500 hPa（目标层固定 500 hPa，与 `pressure_unit` 无关），再调用 `lifted_index`。
