# 热力学（`pymeteo.thermo`）

[English](Thermo)

水汽、位温、LCL、湿球、能见度。水面饱和水汽压沿用李社宏（1994）；露点与相对湿度互换用 Dutton 经验潜热；由混合比求相对湿度用 NCL `relhum` 查表；由相对湿度求混合比用 Tetens（同 NCL `mixhum_ptrh`）。位温为 Poisson，κ = 0.286。相当位温与 LCL 用 Bolton（1980）。湿球为 Stull（2011）海平面经验式。均为按文献自行实现，不拷贝 MetPy / NCL 源码。

可从包根（`import pymeteo as pm`）或 `pymeteo.thermo` 导入。

---

### `saturation_vapor_pressure(temperature, *, temperature_unit="C", output_pressure_unit="hPa")`

水面饱和水汽压，未区分冰面。0 °C 时结果为 6.1078 hPa。

```python
import pymeteo as pm
pm.saturation_vapor_pressure(0.0)                 # 6.1078 hPa
pm.saturation_vapor_pressure(0.0, output_pressure_unit="Pa")  # 610.78
```

### `condensation_temperature(pressure, temperature, dewpoint, *, pressure_unit="hPa", temperature_unit="C", output_temperature_unit="C")`

抬升凝结高度上的温度（李社宏 1994 迭代）。网格插值导致露点高于温度或出现负混合比时可能不收敛，此时返回最后一次迭代值。

```python
pm.condensation_temperature(850.0, 16.6, 0.6)
```

### `relative_humidity_from_dewpoint(temperature, dewpoint, *, temperature_unit="C", output_humidity_unit="%")`

温度 + 露点 → 相对湿度（Dutton）。

### `dewpoint_from_relative_humidity(temperature, relative_humidity, *, temperature_unit="C", humidity_unit="%", output_temperature_unit="C")`

温度 + 相对湿度 → 露点。相对湿度为 0 时返回 `nan`。

```python
pm.dewpoint_from_relative_humidity(18.0, 46.5)   # 约 6.30 °C
```

### `relative_humidity_from_mixing_ratio(temperature, mixing_ratio, pressure, *, temperature_unit="C", mixing_ratio_unit="kg/kg", pressure_unit="hPa", output_humidity_unit="%")`

NCL `relhum` 饱和水汽压表（自 173.16 K 起每隔 1 K）。

### `mixing_ratio_from_relative_humidity(pressure, temperature, relative_humidity, *, pressure_unit="hPa", temperature_unit="C", humidity_unit="%", output_humidity_unit="kg/kg")`

Tetens 混合比。

```python
pm.mixing_ratio_from_relative_humidity(1000.0, 18.0, 46.5)  # 约 0.006018 kg/kg
```

### `specific_humidity_from_relative_humidity(...)`

参数同由相对湿度求混合比。`q = w / (1 + w)`。1000 hPa、18 °C、46.5%、`output_humidity_unit="g/kg"` 时约 5.982 g/kg。

### `convert_humidity(value, *, from_quantity="mixing_ratio", to_quantity="specific_humidity", humidity_unit="kg/kg", output_humidity_unit="kg/kg")`

`from_quantity` / `to_quantity` 为 `"mixing_ratio"` 或 `"specific_humidity"`。`q = w/(1+w)`，`w = q/(1-q)`。比湿 ≥ 1 时对应元素为 `nan`。

### `visibility(relative_humidity, temperature, method="RUC", *, humidity_unit="%", temperature_unit="C", output_distance_unit="km")`

`method` 为 `"RUC"`（60 km × 指数衰减）或 `"FSL"`（露点差公式）。例：相对湿度 80%、18 °C、RUC 约 11.8 km。

### `saturation_mixing_ratio(pressure, temperature, *, pressure_unit="hPa", temperature_unit="C", output_humidity_unit="kg/kg")`

饱和混合比 `w_s(p, T)`。

### `mixing_ratio_from_dewpoint(pressure, dewpoint, *, pressure_unit="hPa", temperature_unit="C", output_humidity_unit="kg/kg")`

与饱和混合比相同，温度实参为露点。

### `vapor_pressure_from_mixing_ratio(pressure, mixing_ratio, *, pressure_unit="hPa", mixing_ratio_unit="kg/kg", output_pressure_unit="hPa")`

### `vapor_pressure_from_relative_humidity(temperature, relative_humidity, *, temperature_unit="C", humidity_unit="%", output_pressure_unit="hPa")`

`e = RH · e_s(T)`。

### `potential_temperature(pressure, temperature, *, pressure_unit="hPa", temperature_unit="C", output_temperature_unit="K")`

Poisson 位温，κ = 0.286，p0 = 1000 hPa。默认输出为 **K**。1000 hPa、301.25 K → 301.25 K。

### `equivalent_potential_temperature(pressure, temperature, dewpoint, *, pressure_unit="hPa", temperature_unit="C", output_temperature_unit="K")`

Bolton（1980）式 (43)，含 LCL。默认输出 **K**。

### `virtual_temperature(temperature, mixing_ratio, *, temperature_unit="C", mixing_ratio_unit="kg/kg", output_temperature_unit=None)`

`T_v = T (1 + r/ε) / (1 + r)`。输出单位默认与 `temperature_unit` 相同。

### `wet_bulb_temperature(temperature, relative_humidity, *, temperature_unit="C", humidity_unit="%", output_temperature_unit="C")`

Stull 2011 海平面经验式。20 °C、50% → 约 13.70 °C。

### `lifting_condensation_level(pressure, temperature, dewpoint, *, pressure_unit="hPa", temperature_unit="C", output_pressure_unit="hPa", output_temperature_unit="C")`

返回 `(p_lcl, T_lcl)`。Bolton LCL 温度再沿干绝热求气压。1000 hPa、15 °C、露点 4 °C 时 p_LCL 约 847–849 hPa（Bolton 与 Stipanuk 差约 1.5 hPa）。

### `parcel_temperature_at_pressure(pressure, temperature, dewpoint, pressure_target, *, pressure_unit="hPa", temperature_unit="C", output_temperature_unit="C")`

气块干绝热抬到 LCL，再湿绝热抬到 `pressure_target`。`lifted_index_from_surface` 用它。
