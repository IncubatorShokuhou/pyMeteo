# Thermo (`pymeteo.thermo`)

[中文](Thermo_zh)

Vapour, potential temperature, LCL, wet-bulb, visibility. Surface saturation vapour pressure follows Li Shehong (1994); dewpoint ↔ RH uses Dutton’s latent-heat relation; mixing-ratio RH uses the NCL `relhum` table; mixing ratio from RH uses Tetens as in NCL `mixhum_ptrh`. Potential temperature is Poisson with κ = 0.286. Equivalent potential temperature and LCL use Bolton (1980). Wet-bulb is Stull (2011), sea-level empirical. All of these are reimplemented; they do not copy MetPy or NCL source.

Import from the package root (`import pymeteo as pm`) or from `pymeteo.thermo`.

---

### `saturation_vapor_pressure(temperature, *, temperature_unit="C", output_pressure_unit="hPa")`

Water-surface saturation vapour pressure. Ice surface is not distinguished. At 0 °C the result is 6.1078 hPa.

```python
import pymeteo as pm
pm.saturation_vapor_pressure(0.0)                 # 6.1078 hPa
pm.saturation_vapor_pressure(0.0, output_pressure_unit="Pa")  # 610.78
```

### `condensation_temperature(pressure, temperature, dewpoint, *, pressure_unit="hPa", temperature_unit="C", output_temperature_unit="C")`

Temperature at the lifting condensation level (Li Shehong 1994 iteration). If interpolation yields Td > T or a negative mixing ratio, the loop may not converge and the last iterate is returned.

```python
pm.condensation_temperature(850.0, 16.6, 0.6)
```

### `relative_humidity_from_dewpoint(temperature, dewpoint, *, temperature_unit="C", output_humidity_unit="%")`

T and Td → RH (Dutton).

### `dewpoint_from_relative_humidity(temperature, relative_humidity, *, temperature_unit="C", humidity_unit="%", output_temperature_unit="C")`

T and RH → dewpoint. RH = 0 yields `nan`.

```python
pm.dewpoint_from_relative_humidity(18.0, 46.5)   # about 6.30 °C
```

### `relative_humidity_from_mixing_ratio(temperature, mixing_ratio, pressure, *, temperature_unit="C", mixing_ratio_unit="kg/kg", pressure_unit="hPa", output_humidity_unit="%")`

Uses the NCL `relhum` saturation table (173.16 K onward, 1 K steps).

### `mixing_ratio_from_relative_humidity(pressure, temperature, relative_humidity, *, pressure_unit="hPa", temperature_unit="C", humidity_unit="%", output_humidity_unit="kg/kg")`

Tetens mixing ratio.

```python
pm.mixing_ratio_from_relative_humidity(1000.0, 18.0, 46.5)  # about 0.006018 kg/kg
```

### `specific_humidity_from_relative_humidity(...)`

Same arguments as mixing ratio from RH. `q = w / (1 + w)`. At 1000 hPa, 18 °C, 46.5% the result is about 5.982 g/kg if `output_humidity_unit="g/kg"`.

### `convert_humidity(value, *, from_quantity="mixing_ratio", to_quantity="specific_humidity", humidity_unit="kg/kg", output_humidity_unit="kg/kg")`

`from_quantity` / `to_quantity`: `"mixing_ratio"` or `"specific_humidity"`. `q = w/(1+w)`, `w = q/(1-q)`. Specific humidity ≥ 1 → `nan`.

### `visibility(relative_humidity, temperature, method="RUC", *, humidity_unit="%", temperature_unit="C", output_distance_unit="km")`

`method` is `"RUC"` (60 km × exponential decay) or `"FSL"` (dewpoint-depression formula). Example: 80% RH, 18 °C, RUC ≈ 11.8 km.

### `saturation_mixing_ratio(pressure, temperature, *, pressure_unit="hPa", temperature_unit="C", output_humidity_unit="kg/kg")`

`w_s(p, T)`.

### `mixing_ratio_from_dewpoint(pressure, dewpoint, *, pressure_unit="hPa", temperature_unit="C", output_humidity_unit="kg/kg")`

Same as saturation mixing ratio with temperature = dewpoint.

### `vapor_pressure_from_mixing_ratio(pressure, mixing_ratio, *, pressure_unit="hPa", mixing_ratio_unit="kg/kg", output_pressure_unit="hPa")`

### `vapor_pressure_from_relative_humidity(temperature, relative_humidity, *, temperature_unit="C", humidity_unit="%", output_pressure_unit="hPa")`

`e = RH · e_s(T)`.

### `potential_temperature(pressure, temperature, *, pressure_unit="hPa", temperature_unit="C", output_temperature_unit="K")`

Poisson θ, κ = 0.286, p0 = 1000 hPa. Default output is **K**. At 1000 hPa, 301.25 K → 301.25 K.

### `equivalent_potential_temperature(pressure, temperature, dewpoint, *, pressure_unit="hPa", temperature_unit="C", output_temperature_unit="K")`

Bolton (1980) eq. (43), including LCL. Default output **K**.

### `virtual_temperature(temperature, mixing_ratio, *, temperature_unit="C", mixing_ratio_unit="kg/kg", output_temperature_unit=None)`

`T_v = T (1 + r/ε) / (1 + r)`. Output unit defaults to `temperature_unit`.

### `wet_bulb_temperature(temperature, relative_humidity, *, temperature_unit="C", humidity_unit="%", output_temperature_unit="C")`

Stull 2011 sea-level formula. 20 °C, 50% → about 13.70 °C.

### `lifting_condensation_level(pressure, temperature, dewpoint, *, pressure_unit="hPa", temperature_unit="C", output_pressure_unit="hPa", output_temperature_unit="C")`

Returns `(p_lcl, T_lcl)`. Bolton LCL temperature then dry-adiabatic pressure. 1000 hPa, 15 °C, Td 4 °C → p_LCL around 847–849 hPa (Bolton vs Stipanuk differs by ~1.5 hPa).

### `parcel_temperature_at_pressure(pressure, temperature, dewpoint, pressure_target, *, pressure_unit="hPa", temperature_unit="C", output_temperature_unit="C")`

Dry adiabatic to LCL, then moist adiabatic to `pressure_target`. Used by `lifted_index_from_surface`.
