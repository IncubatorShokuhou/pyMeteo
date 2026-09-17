# Wind (`pymeteo.wind`)

[中文](Wind_zh)

Meteorological **from**-direction in degrees, range `[0, 360)`. Calm (speed 0) → direction 0 in the modern API.

---

### `wind_speed(u, v, *, speed_unit="m/s", output_speed_unit=None)`

`sqrt(u²+v²)`. Output unit defaults to `speed_unit`.

```python
pm.wind_speed(3.0, 4.0)  # 5.0 m/s
pm.wind_speed(10.0, 0.0, speed_unit="m/s", output_speed_unit="kt")
```

### `wind_direction(u, v, *, speed_unit="m/s")`

From-direction in degrees. `speed_unit` is kept for a consistent interface; direction does not depend on it.

```python
pm.wind_direction(3.0, 4.0)  # about 216.87°
pm.wind_direction(0.0, 0.0)  # 0.0
```

NCL `wind_direction` in [`pymeteo.ncl`](NCL) adds `opt` for the calm fill (0 / nan / custom).

### `uv_from_speed_direction(speed, direction, *, speed_unit="m/s", output_speed_unit=None)`

Direction in **degrees** (from). Meteorological convention: 90° → u = −speed, v = 0; 0° → u = 0, v = −speed.

```python
u, v = pm.uv_from_speed_direction(10.0, 270.0, speed_unit="kt", output_speed_unit="m/s")
```

### `wind_components(speed, direction, *, speed_unit="m/s", output_speed_unit=None)`

Alias of `uv_from_speed_direction`.

### `bulk_wind_shear(u_bottom, v_bottom, u_top, v_top, *, speed_unit="m/s", output_speed_unit=None)`

Magnitude of the vector difference between two layers.
