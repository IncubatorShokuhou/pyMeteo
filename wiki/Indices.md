# Indices

[中文](Indices_zh.md)

Stability indices from a few mandatory-level temperatures (and, for SWEAT, winds).

Showalter, K, A, TT, and lifted index have the same numeric value whether you pass `C` or `K`. SWEAT is defined in °C and knots internally. There are no NCL-name wrappers for these — NCL does not ship matching builtins, and `pymeteo.ncl` does not invent them.

## temperature_dewpoint_depression

`T − Td`. In Fahrenheit the difference is scaled by 9/5.

```python
import pymeteo as pm
pm.temperature_dewpoint_depression(7.0, -2.0)  # 9.0
```

## layer_temperature_difference

Lower minus upper. Typical use is 850 hPa minus 500 hPa.

```python
pm.layer_temperature_difference(16.6, -15.9)  # 32.5
```

## k_index

`K = T850 − T500 + Td850 − (T700 − Td700)`, evaluated in °C.

```python
pm.k_index(16.6, 0.6, 7.0, -2.0, -15.9)  # 24.1
```

A common rule of thumb (from the docstring, not a forecast): K < 20 little thunderstorm activity; 20–25 isolated; 25–30 scattered; 30–35 numerous.

## a_index

`A = (T850 − T500) − (T850 − Td850) − (T700 − Td700) − (T500 − Td500)`. Needs Td500 as well.

```python
pm.a_index(16.6, 0.6, 7.0, -2.0, -15.9, -20.0)  # 3.4
```

## total_totals_index

`TT = T850 + Td850 − 2·T500`. Give exactly one of `dewpoint_850` or `relative_humidity_850`. Humidity is turned into dewpoint with the Dutton formula.

```python
pm.total_totals_index(18.0, -15.9, dewpoint_850=6.3)
```

## showalter_index

Parcel from 850 hPa lifted to 500 hPa. SI = T500 − T_parcel(500). The moist-adiabatic piece is Li Shehong (1994).

```python
pm.showalter_index(16.6, 0.6, -15.9)  # about 1.1
```

Negative mixing ratios from interpolated grids can make the iteration unhappy. That is a data problem, not a unit problem.

## sweat_index

Miller (1972) / NWS:

`12·Td850(°C) + 20·(TT−49) + 2·f850(kt) + f500(kt) + 125·(S+0.2)`

Negative terms are zeroed. The shear term is kept only when 850° is 130–250, 500° is 210–310, the direction difference is positive, and both speeds are at least 15 kt.

Pass Td850 **or** RH850. Wind components default to `m/s` and are converted to knots inside. Rule of thumb in the docstring: >300 severe-convection potential, >400 tornado potential.

## lifted_index

`LI = T500 − T_parcel(500)` when you already have the parcel temperature at 500 hPa. Negative means the parcel is warmer than the environment.

## lifted_index_from_surface

Near-surface p, T, Td plus T500. The parcel is taken to LCL (Bolton), then to 500 hPa (Li Shehong moist entropy). The 500 hPa target is fixed; `pressure_unit` only describes the surface pressure you passed.

This is the usual surface-based LI. CAPE, most-unstable parcels, and mixed-layer parcels are out of scope.
