# Dynamics (`pymeteo.dynamics`)

[中文](Dynamics_zh.md)

Pointwise only: no grid advection, spherical harmonics, or finite differences. ω ↔ w uses the hydrostatic ideal-gas relation also used by NCL `omega_to_w` / `w_to_omega`, reimplemented here.

---

### `coriolis_parameter(latitude, *, latitude_unit="deg")`

`f = 2 Ω sin φ` with `Ω = 7.292×10⁻⁵` rad s⁻¹. Result in s⁻¹. Equator 0; 35° ≈ 8.365×10⁻⁵; 45° ≈ 1.031×10⁻⁴; north pole ≈ 1.458×10⁻⁴.

```python
pm.coriolis_parameter(35.0)
pm.coriolis_parameter(45.0, latitude_unit="deg")
```

NCL name: `pymeteo.ncl.coriolis_param(lat)` (degrees in, s⁻¹ out).

### `omega_to_w(omega, temperature, pressure, *, omega_unit="Pa/s", temperature_unit="C", pressure_unit="hPa", output_speed_unit="m/s")`

`w = −ω / (ρ g)` with `ρ = p / (Rd T)`, `Rd = 287.058` J K⁻¹ kg⁻¹, `g = 9.80665` m s⁻². Positive ω (subsidence) → negative w.

**Argument order** is `(omega, temperature, pressure)`. The NCL shim is `(omega, p, t)` — see [NCL](NCL.md).

`omega_unit` default `Pa/s`; `hPa/s` is accepted.

```python
w = pm.omega_to_w(0.1, 0.0, 850.0)
```

### `w_to_omega(w, temperature, pressure, *, speed_unit="m/s", temperature_unit="C", pressure_unit="hPa", output_omega_unit="Pa/s")`

Algebraic inverse: `ω = −ρ g w`. Same modern argument order `(w, temperature, pressure)`.
