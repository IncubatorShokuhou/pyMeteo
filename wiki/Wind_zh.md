# 风（`pymeteo.wind`）

[English](Wind.md)

气象学**来向**，单位度，范围 `[0, 360)`。静风（风速为 0）时现代 API 返回风向 0。

---

### `wind_speed(u, v, *, speed_unit="m/s", output_speed_unit=None)`

`sqrt(u²+v²)`。输出单位默认与 `speed_unit` 相同。

```python
pm.wind_speed(3.0, 4.0)  # 5.0 m/s
pm.wind_speed(10.0, 0.0, speed_unit="m/s", output_speed_unit="kt")
```

### `wind_direction(u, v, *, speed_unit="m/s")`

来向，度。`speed_unit` 仅为接口一致；风向与单位无关。

```python
pm.wind_direction(3.0, 4.0)  # 约 216.87°
pm.wind_direction(0.0, 0.0)  # 0.0
```

[`pymeteo.ncl`](NCL_zh.md) 里的 `wind_direction` 另有 `opt`，用于静风填充（0 / nan / 自定义）。

### `uv_from_speed_direction(speed, direction, *, speed_unit="m/s", output_speed_unit=None)`

风向为**度**（来向）。气象学约定：90° → u = −speed，v = 0；0° → u = 0，v = −speed。

```python
u, v = pm.uv_from_speed_direction(10.0, 270.0, speed_unit="kt", output_speed_unit="m/s")
```

### `wind_components(speed, direction, *, speed_unit="m/s", output_speed_unit=None)`

与 `uv_from_speed_direction` 相同。

### `bulk_wind_shear(u_bottom, v_bottom, u_top, v_top, *, speed_unit="m/s", output_speed_unit=None)`

两层水平风矢量差的模。
