# Geo

[中文](Geo_zh.md)

Distance on the ellipsoid, gravity from latitude, station pressure reduced to sea level, and layer thickness from the hypsometric equation.

## earth_distance

Geodesic distance between two lat/lon points. WGS84 Vincenty inverse, in NumPy — no geopy.

`angle_unit` default `deg` (`rad` accepted). Output default `km`. Identical points are 0. A handful of antipodal pairs that fail to converge fall back to a mean-Earth-radius haversine.

```python
import pymeteo as pm
pm.earth_distance(0.0, 0.0, 0.0, 1.0, output_distance_unit="m")  # ~111319 m
pm.earth_distance(39.9, 116.4, 31.2, 121.5)  # ~1070 km
```

## gravity

`g = 9.7803 · (1 + 0.0053024 sin²φ − 0.000005 sin²2φ)` m/s². Latitude is converted to radians first. Equator: 9.7803.

```python
pm.gravity(0.0)
pm.gravity(45.0)
```

## sea_level_pressure

Reduce station pressure to sea level with a simple empirical column:

`tm = (t + t12) / 2 + lapse_rate · h / 2`

`p0 = ph · 10^(h / (18400 · (1 + tm/273)))`

The 273 and 18400 are the original empirical constants; temperature is substituted in Celsius.

`lapse_rate` is **per metre in `temperature_unit`**. The default 0.005 is 0.5 °C / 100 m. If you pass Fahrenheit temperatures, change the lapse rate too.

```python
pm.sea_level_pressure(1000.0, 100.0, 20.0, 18.0)  # about 1011.8 hPa
```

## height_thickness

Hypsometric thickness between two isobaric surfaces:

`ΔZ = (Rd T̄ / g0) ln(p_bottom / p_top)`

`Rd = 287.058` J K⁻¹ kg⁻¹, `g0 = 9.80665` m s⁻². Output default is **metres**, not km. Pass virtual temperature as `mean_temperature` if the layer is moist.

```python
pm.height_thickness(1000.0, 500.0, 0.0)  # about 5542 m
```
