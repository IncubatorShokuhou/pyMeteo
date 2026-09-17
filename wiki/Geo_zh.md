# 几何与气压（`pymeteo.geo`）

[English](Geo.md)

大地测量、重力、本站气压订正到海平面、压高厚度。

---

### `earth_distance(latitude_1, longitude_1, latitude_2, longitude_2, *, angle_unit="deg", output_distance_unit="km")`

WGS84 Vincenty 反解。极少数对跖点迭代失败时回退到平均地球半径的 haversine。相同点 → 0。

```python
pm.earth_distance(0.0, 0.0, 0.0, 1.0, output_distance_unit="m")  # 约 111319.5 m
pm.earth_distance(39.9, 116.4, 31.2, 121.5, output_distance_unit="km")
```

### `gravity(latitude, *, latitude_unit="deg")`

`g = 9.7803 · (1 + 0.0053024 sin²φ − 0.000005 sin²2φ)`，单位 m/s²。纬度先化为弧度。赤道：9.7803。

```python
pm.gravity(0.0)
pm.gravity(45.0, latitude_unit="deg")
```

### `sea_level_pressure(station_pressure, station_height, temperature, temperature_12h_ago, *, lapse_rate=0.005, pressure_unit="hPa", height_unit="m", temperature_unit="C", output_pressure_unit="hPa")`

气柱平均温度 `tm = (t + t12)/2 + lapse_rate·h/2`，再按 `p0 = ph · 10^(h / (18400 · (1 + tm/273)))` 订正。273 与 18400 沿用原经验常数，温度按摄氏度代入。

`lapse_rate` 的单位是 **`temperature_unit` 每米**。默认 0.005 即每 100 m 降低 0.5 °C。若温度用华氏度，递减率也要改。

```python
pm.sea_level_pressure(1000.0, 100.0, 20.0, 18.0)  # 约 1011.76 hPa
```

### `height_thickness(pressure_bottom, pressure_top, mean_temperature, *, pressure_unit="hPa", temperature_unit="C", output_distance_unit="m")`

压高方程：`Δz = (Rd T̄ / g0) ln(p_bottom / p_top)`，`Rd = 287.058` J K⁻¹ kg⁻¹，`g0 = 9.80665` m s⁻²。默认输出为**米**。

```python
pm.height_thickness(1000.0, 500.0, 0.0, temperature_unit="C")
```
