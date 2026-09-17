# Wind

[中文](Wind_zh.md)

Horizontal wind at a point: speed, meteorological from-direction, and u/v.

Direction is the direction the wind is coming from. 0° is north, increasing clockwise, range `[0, 360)`. Calm (speed 0) returns direction 0 in the modern API.

`u = −speed · sin(direction)`, `v = −speed · cos(direction)`, with direction in radians for the trig.

## wind_speed

`sqrt(u² + v²)`. `speed_unit` describes the incoming components; `output_speed_unit` defaults to the same. Aliases include `kt` / `knots`, `km/h`, `mph`.

```python
import pymeteo as pm
pm.wind_speed(3.0, 4.0)  # 5.0 m/s
pm.wind_speed(10.0, 0.0, speed_unit="m/s", output_speed_unit="kt")
```

## wind_direction

From-direction in degrees. `speed_unit` is there so the signature matches the other wind functions; the angle does not depend on it.

```python
pm.wind_direction(3.0, 4.0)  # about 216.9°
pm.wind_direction(0.0, 0.0)  # 0.0
```

The NCL shim (`pymeteo.ncl.wind_direction`) adds `opt` for the calm fill: 0, `nan`, or a custom scalar. See [NCL](NCL.md).

## uv_from_speed_direction

Speed and from-direction (degrees) → `(u, v)`. 90° is easterly: u = −speed, v = 0. 0° is northerly: u = 0, v = −speed.

```python
u, v = pm.uv_from_speed_direction(10.0, 270.0, speed_unit="kt")
```

## wind_components

Same function as `uv_from_speed_direction`. Keep whichever name you already type.

## bulk_wind_shear

Magnitude of the vector difference between two layers:

`S = sqrt((u_top − u_bottom)² + (v_top − v_bottom)²)`

Useful as a 0–6 km bulk shear once you have the two winds. It does not integrate a profile.
