# 风

[English](Wind.md)

点上的水平风：风速、气象来向、uv 分量。

风向是风从哪来。0° 为正北，顺时针增加，范围 `[0, 360)`。静风（风速为 0）在现代接口里风向返回 0。

`u = −speed · sin(direction)`，`v = −speed · cos(direction)`，三角函数用弧度。

## wind_speed

`sqrt(u² + v²)`。`speed_unit` 是输入分量的单位；`output_speed_unit` 默认跟它相同。别名包括 `kt` / `knots`、`km/h`、`mph`。

```python
import pymeteo as pm
pm.wind_speed(3.0, 4.0)  # 5.0 m/s
pm.wind_speed(10.0, 0.0, speed_unit="m/s", output_speed_unit="kt")
```

## wind_direction

来向，单位度。`speed_unit` 只是为了和其它风函数签名一致，风向跟单位无关。

```python
pm.wind_direction(3.0, 4.0)  # 大约 216.9°
pm.wind_direction(0.0, 0.0)  # 0.0
```

NCL 封装（`pymeteo.ncl.wind_direction`）多一个 `opt`，用来填静风：0、`nan` 或自定义标量。见 [NCL](NCL_zh.md)。

## uv_from_speed_direction

风速和来向（度）→ `(u, v)`。90° 是东风：u = −speed，v = 0。0° 是北风：u = 0，v = −speed。

```python
u, v = pm.uv_from_speed_direction(10.0, 270.0, speed_unit="kt")
```

## wind_components

和 `uv_from_speed_direction` 是同一个函数。用哪个名字都行。

## bulk_wind_shear

两层水平风矢量差的模：

`S = sqrt((u_top − u_bottom)² + (v_top − v_bottom)²)`

有两层风就可以当 0–6 km 体切变用。它不沿高度积分。
