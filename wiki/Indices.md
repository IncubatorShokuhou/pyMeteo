# Indices (`pymeteo.indices`)

[中文](Indices_zh.md)

Stability indices. Showalter, K, A, and TT have the same numeric value whether temperatures are passed in `C` or `K`. SWEAT is defined in °C and knots internally.

There are **no** NCL builtin names for these indices in `pymeteo.ncl`.

---

### `temperature_dewpoint_depression(temperature, dewpoint, *, temperature_unit="C", output_temperature_unit="C")`

`T − Td`. In `F`, the difference is scaled by 9/5.

```python
pm.temperature_dewpoint_depression(7.0, -2.0)  # 9.0
```

### `layer_temperature_difference(temperature_lower, temperature_upper, *, temperature_unit="C", output_temperature_unit="C")`

Lower minus upper.

```python
pm.layer_temperature_difference(16.6, -15.9)  # 32.5
```

### `k_index(temperature_850, dewpoint_850, temperature_700, dewpoint_700, temperature_500, *, temperature_unit="C")`

`K = T850 − T500 + Td850 − (T700 − Td700)`, evaluated in °C.

```python
pm.k_index(16.6, 0.6, 7.0, -2.0, -15.9)  # 24.1
```

Docstring rule of thumb: K < 20 none; 20–25 isolated; 25–30 scattered; 30–35 numerous thunderstorms.

### `a_index(temperature_850, dewpoint_850, temperature_700, dewpoint_700, temperature_500, dewpoint_500, *, temperature_unit="C")`

`A = (T850 − T500) − (T850 − Td850) − (T700 − Td700) − (T500 − Td500)`.

```python
pm.a_index(16.6, 0.6, 7.0, -2.0, -15.9, -20.0)  # 3.4
```

### `total_totals_index(temperature_850, temperature_500, *, dewpoint_850=None, relative_humidity_850=None, temperature_unit="C", humidity_unit="%")`

`TT = T850 + Td850 − 2·T500`. Provide **exactly one** of `dewpoint_850` or `relative_humidity_850`. Humidity is converted to dewpoint with Dutton’s formula.

```python
pm.total_totals_index(18.0, -15.9, relative_humidity_850=46.5)
```

### `showalter_index(temperature_850, dewpoint_850, temperature_500, *, temperature_unit="C")`

Parcel from 850 hPa lifted to 500 hPa; SI = T500 − T_parcel(500). Iteration is per-element with a max step count.

```python
pm.showalter_index(16.6, 0.6, -15.9)  # about 1.1
```

### `sweat_index(temperature_850, temperature_500, u_850, v_850, u_500, v_500, *, dewpoint_850=None, relative_humidity_850=None, temperature_unit="C", humidity_unit="%", speed_unit="m/s")`

NWS / Miller (1972):

`12·Td850(°C) + 20·(TT−49) + 2·f850(kt) + f500(kt) + 125·(S+0.2)`

Negative terms are set to zero. The shear term requires a positive direction difference and both layer speeds ≥ 15 kt. Provide Td850 **or** RH850. Wind components default to `m/s` and are converted to knots inside.

### `lifted_index(temperature_500, parcel_temperature_500, *, temperature_unit="C")`

`LI = T500 − T_parcel(500)` when the parcel temperature at 500 hPa is already known.

### `lifted_index_from_surface(pressure, temperature, dewpoint, temperature_500, *, pressure_unit="hPa", temperature_unit="C")`

Lifts a near-surface parcel to 500 hPa (`pressure_target` is fixed at 500 hPa; independent of `pressure_unit`) then calls `lifted_index`.
