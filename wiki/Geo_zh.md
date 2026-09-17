# 几何与气压

[English](Geo.md)

椭球上的距离、按纬度的重力、本站气压订正海平面，以及压高公式给的气层厚度。

## earth_distance

两点经纬度之间的大地线距离。WGS84 Vincenty 反解，纯 NumPy，不依赖 geopy。

`angle_unit` 默认 `deg`（也认 `rad`）。输出默认 `km`。重合点是 0。极少数对跖点迭代失败时，退回平均地球半径的 haversine。

```python
import pymeteo as pm
pm.earth_distance(0.0, 0.0, 0.0, 1.0, output_distance_unit="m")  # 约 111319 m
pm.earth_distance(39.9, 116.4, 31.2, 121.5)  # 约 1070 km
```

## gravity

`g = 9.7803 · (1 + 0.0053024 sin²φ − 0.000005 sin²2φ)`，单位 m/s²。纬度会先换成弧度。赤道是 9.7803。

```python
pm.gravity(0.0)
pm.gravity(45.0)
```

## sea_level_pressure

用一个很简单的经验气柱，把本站气压订正到海平面：

`tm = (t + t12) / 2 + lapse_rate · h / 2`

`p0 = ph · 10^(h / (18400 · (1 + tm/273)))`

273 和 18400 是原来的经验常数；温度按摄氏度代入。

`lapse_rate` 的量纲是 **`temperature_unit` 每米**。默认 0.005 就是每 100 m 降 0.5 °C。如果温度用华氏度，递减率也要跟着改。

```python
pm.sea_level_pressure(1000.0, 100.0, 20.0, 18.0)  # 大约 1011.8 hPa
```

## height_thickness

两层等压面之间的厚度（压高公式）：

`ΔZ = (Rd T̄ / g0) ln(p_bottom / p_top)`

`Rd = 287.058` J K⁻¹ kg⁻¹，`g0 = 9.80665` m s⁻²。输出默认是 **米**，不是千米。气层有水汽时，`mean_temperature` 应传入虚温。

```python
pm.height_thickness(1000.0, 500.0, 0.0)  # 大约 5542 m
```
