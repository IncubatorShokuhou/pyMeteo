# Dynamics

[中文](Dynamics_zh.md)

Pointwise only: Coriolis parameter and hydrostatic ω ↔ w. No grid advection, spherical harmonics, or finite differences.

## coriolis_parameter

`f = 2 Ω sin φ` with `Ω = 7.292×10⁻⁵` rad s⁻¹. Result in s⁻¹. Equator is 0; 45° is about `1.031×10⁻⁴`.

```python
import pymeteo as pm
pm.coriolis_parameter(35.0)
pm.coriolis_parameter(45.0)
```

NCL name: `pymeteo.ncl.coriolis_param(lat)` — degrees in, s⁻¹ out.

## omega_to_w

Pressure-coordinate vertical velocity to geometric `w`:

`w = −ω / (ρ g)`, `ρ = p / (Rd T)`

`Rd = 287.058` J K⁻¹ kg⁻¹, `g = 9.80665` m s⁻². Positive ω (subsidence) gives negative w.

Modern argument order is `(omega, temperature, pressure)`. The NCL shim is `(omega, p, t)` — see [NCL](NCL.md).

`omega_unit` default `Pa/s`; `hPa/s` is accepted. Temperature default `C`, pressure default `hPa`.

```python
pm.omega_to_w(0.1, 0.0, 850.0)
```

## w_to_omega

The inverse: `ω = −ρ g w`. Same modern order `(w, temperature, pressure)`.
