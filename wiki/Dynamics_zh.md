# 动力学

[English](Dynamics.md)

只做点上的量：科氏参数，以及静力近似下的 ω ↔ w。没有网格平流、球谐或有限差分。

## coriolis_parameter

`f = 2 Ω sin φ`，`Ω = 7.292×10⁻⁵` rad s⁻¹。结果单位 s⁻¹。赤道是 0；45° 大约 `1.031×10⁻⁴`。

```python
import pymeteo as pm
pm.coriolis_parameter(35.0)
pm.coriolis_parameter(45.0)
```

NCL 名：`pymeteo.ncl.coriolis_param(lat)`，纬度用度，结果 s⁻¹。

## omega_to_w

气压坐标垂直速度换成几何 `w`：

`w = −ω / (ρ g)`，`ρ = p / (Rd T)`

`Rd = 287.058` J K⁻¹ kg⁻¹，`g = 9.80665` m s⁻²。ω 为正（下沉）时 w 为负。

现代接口参数顺序是 `(omega, temperature, pressure)`。NCL 封装是 `(omega, p, t)`，见 [NCL](NCL_zh.md)。

`omega_unit` 默认 `Pa/s`，也认 `hPa/s`。温度默认 `C`，气压默认 `hPa`。

```python
pm.omega_to_w(0.1, 0.0, 850.0)
```

## w_to_omega

代数逆：`ω = −ρ g w`。现代顺序同样是 `(w, temperature, pressure)`。
