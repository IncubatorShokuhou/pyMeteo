# Thermo

[中文](Thermo_zh.md)

Humidity and parcel temperature. Most of the library’s thermodynamics lives here.

Functions are re-exported at the package root (`import pymeteo as pm`). Defaults are °C, hPa, %, kg/kg unless noted. Formula citations are on the individual functions.

## saturation_vapor_pressure

Saturation vapour pressure over liquid water. Ice is not treated separately.

`temperature` default unit `C`; result default `hPa`. At 0 °C this is 6.1078 hPa.

```python
import pymeteo as pm
pm.saturation_vapor_pressure(0.0)
pm.saturation_vapor_pressure(0.0, output_pressure_unit="Pa")  # 610.78
```

## condensation_temperature

Temperature at the lifting condensation level, via Li Shehong’s 1994 iteration.

Takes starting `pressure`, `temperature`, and `dewpoint`. If Td > T or the mixing ratio goes negative (common with sloppy grid interpolation), the loop may not settle; the last iterate is returned.

```python
pm.condensation_temperature(850.0, 16.6, 0.6)
```

## relative_humidity_from_dewpoint

T and Td → RH. Output default `%`; pass `output_humidity_unit="fraction"` for 0–1.

## dewpoint_from_relative_humidity

The inverse. RH ≤ 0 yields `nan` on that element.

```python
pm.dewpoint_from_relative_humidity(18.0, 46.5)  # about 6.30 °C
```

## relative_humidity_from_mixing_ratio

T, mixing ratio, and pressure → RH. Saturation vapour pressure is interpolated from a table spanning 173.16–375.16 K (same layout as NCL `relhum`). Values above 100% are left alone; negatives are clipped to 0.0001%.

`mixing_ratio_unit` default `kg/kg` (`g/kg` is accepted).

## mixing_ratio_from_relative_humidity

Pressure, T, RH → mixing ratio. Tetens, matching NCL `mixhum_ptrh`.

```python
pm.mixing_ratio_from_relative_humidity(1000.0, 18.0, 46.5)  # ~0.00602 kg/kg
```

## specific_humidity_from_relative_humidity

Same arguments. Then `q = w / (1 + w)`. At 1000 hPa, 18 °C, 46.5% you get about 5.98 g/kg if you ask for `output_humidity_unit="g/kg"`.

## convert_humidity

Mixes mixing ratio and specific humidity, and can change the mass unit at the same time.

`from_quantity` / `to_quantity` are `"mixing_ratio"` or `"specific_humidity"`. Specific humidity ≥ 1 is `nan` when converting to mixing ratio.

## visibility

A rough visibility estimate from RH and T. `method="RUC"` is a 60 km exponential decay; `"FSL"` uses dewpoint depression. Default output is km.

80% RH at 18 °C, RUC, is about 11.8 km.

## saturation_mixing_ratio

`w_s = ε e_s(T) / (p − e_s)`, ε = 0.622. `e_s` is the Li Shehong water-surface formula. If `e_s ≥ p`, that element is `nan`.

## mixing_ratio_from_dewpoint

Same formula with dewpoint in place of temperature. Equivalent to NCL `mixhum_ptd` in spirit.

## vapor_pressure_from_mixing_ratio

`e = w p / (ε + w)`. Algebraic inverse of the mixing-ratio definition.

## vapor_pressure_from_relative_humidity

`e = RH · e_s(T)`, again with the Li Shehong `e_s`.

## potential_temperature

Poisson θ: `θ = T (1000 hPa / p)^0.286`. Output default is **K**, not °C.

```python
pm.potential_temperature(1000.0, 28.1)  # 301.25 K
```

## equivalent_potential_temperature

Bolton (1980) eq. (43), including the LCL temperature from his eq. (22). Output default K.

## virtual_temperature

`T_v = T (1 + r/ε) / (1 + r)`. Output unit defaults to `temperature_unit`. This is the exact moist equation of state, not NCL’s `T (1 + 0.61 r)` approximation.

## wet_bulb_temperature

Stull (2011) empirical wet-bulb. Meant for ~1013 hPa. Useful roughly −20–50 °C and 5–99% RH; very cold and dry is where it drifts.

```python
pm.wet_bulb_temperature(20.0, 50.0)  # about 13.7 °C
```

## lifting_condensation_level

Returns `(p_lcl, T_lcl)`. `T_L` from Bolton (1980) eq. (22), then `p_L = p (T_L / T)^{1/κ}` with κ = 0.286.

Wallace & Hobbs’ 1000 hPa / 15 °C / Td 4 °C example comes out near 847 hPa here. NCL `lclvl` uses Stipanuk (1973); a couple of hPa of disagreement is expected.

## parcel_temperature_at_pressure

Lift a parcel dry-adiabatically to the LCL, then moist-adiabatically to `pressure_target`. If the target is still below the LCL (`p_target ≥ p_LCL`), the whole path is dry. Used by `lifted_index_from_surface`.
