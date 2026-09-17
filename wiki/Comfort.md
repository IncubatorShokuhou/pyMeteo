# Comfort

[中文](Comfort_zh.md)

Two NWS operational regressions: heat index and wind chill.

Internally the code converts to Fahrenheit and miles per hour, then converts back to whatever you asked for. NCL has no matching builtins, so there is nothing under `pymeteo.ncl` for these.

## heat_index

Rothfusz (1990) / NWS SR 90-23, after Steadman (1979).

A simple average with air temperature is used first. If that average is at least 80 °F, the full Rothfusz polynomial is used, with the usual low-humidity and high-humidity adjustments.

Outside Steadman’s original table (very hot and humid) the formula is not trustworthy. Output unit defaults to `temperature_unit`.

```python
import pymeteo as pm
pm.heat_index(90.0, 60.0, temperature_unit="F")  # about 100 °F
```

## wind_chill

NWS / Environment Canada 2001, intended for wind near 10 m:

`WC = 35.74 + 0.6215 T − 35.75 V^0.16 + 0.4275 T V^0.16`

(`T` in °F, `V` in mph.) Roughly valid for `T ≤ 50 °F` and `V ≥ 3 mph`. Values are still returned outside that range; they just mean less.

```python
pm.wind_chill(0.0, 10.0, temperature_unit="F", speed_unit="mph")  # about −16 °F
```
