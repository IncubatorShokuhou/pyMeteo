# Geo (`pymeteo.geo`)

[中文](Geo_zh)

Geodesy, gravity, station-to-sea-level pressure, hypsometric thickness.

---

### `earth_distance(latitude_1, longitude_1, latitude_2, longitude_2, *, angle_unit="deg", output_distance_unit="km")`

WGS84 Vincenty inverse. A handful of antipodal pairs that fail iteration fall back to a mean-Earth-radius haversine. Identical points → 0.

```python
pm.earth_distance(0.0, 0.0, 0.0, 1.0, output_distance_unit="m")  # about 111319.5 m
pm.earth_distance(39.9, 116.4, 31.2, 121.5, output_distance_unit="km")
```

### `gravity(latitude, *, latitude_unit="deg")`

`g = 9.7803 · (1 + 0.0053024 sin²φ − 0.000005 sin²2φ)` m/s². Latitude is converted to radians first. Equator: 9.7803.

```python
pm.gravity(0.0)
pm.gravity(45.0, latitude_unit="deg")
```

### `sea_level_pressure(station_pressure, station_height, temperature, temperature_12h_ago, *, lapse_rate=0.005, pressure_unit="hPa", height_unit="m", temperature_unit="C", output_pressure_unit="hPa")`

Column mean temperature `tm = (t + t12)/2 + lapse_rate·h/2`, then `p0 = ph · 10^(h / (18400 · (1 + tm/273)))`. The 273 and 18400 are the original empirical constants; temperature is substituted in Celsius.

`lapse_rate` is **per metre in `temperature_unit`**. Default 0.005 is 0.5 °C / 100 m. If you pass Fahrenheit temperatures, change the lapse rate too.

```python
pm.sea_level_pressure(1000.0, 100.0, 20.0, 18.0)  # about 1011.76 hPa
```

### `height_thickness(pressure_bottom, pressure_top, mean_temperature, *, pressure_unit="hPa", temperature_unit="C", output_distance_unit="m")`

Hypsometric: `Δz = (Rd T̄ / g0) ln(p_bottom / p_top)` with `Rd = 287.058` J K⁻¹ kg⁻¹, `g0 = 9.80665` m s⁻². Default output is **metres**.

```python
pm.height_thickness(1000.0, 500.0, 0.0, temperature_unit="C")
```
