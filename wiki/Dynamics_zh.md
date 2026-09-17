# 动力学（`pymeteo.dynamics`）

[English](Dynamics.md)

只做点上诊断：不做网格平流、球谐或有限差分。ω ↔ w 采用静力、理想气体近似，与 NCL `omega_to_w` / `w_to_omega` 同一关系式，按文献自行实现。

---

### `coriolis_parameter(latitude, *, latitude_unit="deg")`

`f = 2 Ω sin φ`，`Ω = 7.292×10⁻⁵` rad s⁻¹，结果单位 s⁻¹。赤道为 0；35° 约 8.365×10⁻⁵；45° 约 1.031×10⁻⁴；北极约 1.458×10⁻⁴。

```python
pm.coriolis_parameter(35.0)
pm.coriolis_parameter(45.0, latitude_unit="deg")
```

NCL 名：`pymeteo.ncl.coriolis_param(lat)`（纬度度，返回 s⁻¹）。

### `omega_to_w(omega, temperature, pressure, *, omega_unit="Pa/s", temperature_unit="C", pressure_unit="hPa", output_speed_unit="m/s")`

`w = −ω / (ρ g)`，`ρ = p / (Rd T)`，`Rd = 287.058` J K⁻¹ kg⁻¹，`g = 9.80665` m s⁻²。下沉 ω>0 → w<0。

**现代参数顺序**是 `(omega, temperature, pressure)`。NCL 封装是 `(omega, p, t)`，见 [NCL 兼容](NCL_zh.md)。

`omega_unit` 默认 `Pa/s`，也接受 `hPa/s`。

```python
w = pm.omega_to_w(0.1, 0.0, 850.0)
```

### `w_to_omega(w, temperature, pressure, *, speed_unit="m/s", temperature_unit="C", pressure_unit="hPa", output_omega_unit="Pa/s")`

代数逆：`ω = −ρ g w`。现代参数顺序同样是 `(w, temperature, pressure)`。
