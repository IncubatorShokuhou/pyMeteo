# Comfort (`pymeteo.comfort`)

[中文](Comfort_zh.md)

NWS operational regressions, reimplemented (not copied from MetPy). Internally the code converts to Fahrenheit and miles per hour, then back. **No NCL-name wrappers** — NCL has no matching builtins, and `pymeteo.ncl` does not invent them.

---

### `heat_index(temperature, relative_humidity, *, temperature_unit="C", humidity_unit="%", output_temperature_unit=None)`

Rothfusz (1990) / NWS SR 90-23, after Steadman (1979):

1. Simple form `HI = 0.5 {T + 61 + (T−68)·1.2 + RH·0.094}` (`T` in °F, `RH` in percent), averaged with air temperature.
2. If that average ≥ 80 °F, switch to the full Rothfusz polynomial, with the low-humidity (RH < 13%, 80–112 °F) and high-humidity (RH > 85%, 80–87 °F) adjustments.

Outside Steadman’s original table the formula is not reliable. 90 °F, 60% → about 100 °F.

```python
pm.heat_index(90.0, 60.0, temperature_unit="F", output_temperature_unit="F")
```

Output unit defaults to `temperature_unit`.

### `wind_chill(temperature, wind_speed, *, temperature_unit="C", speed_unit="m/s", output_temperature_unit=None)`

NWS / Environment Canada 2001:

`WC = 35.74 + 0.6215 T − 35.75 V^0.16 + 0.4275 T V^0.16`

(`T` in °F, `V` in mph; wind intended near 10 m). Valid roughly for `T ≤ 50 °F` and `V ≥ 3 mph`; values are still returned outside that range. 0 °F, 10 mph → about −16 °F.

```python
pm.wind_chill(0.0, 10.0, temperature_unit="F", speed_unit="mph", output_temperature_unit="F")
```
