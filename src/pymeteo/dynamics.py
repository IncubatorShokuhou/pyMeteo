"""点上的轻量动力学诊断：科里奥利参数与垂直速度换算。

不做网格平流、球谐或有限差分。ω ↔ w 采用静力、理想气体近似，与 NCL
``omega_to_w`` / ``w_to_omega`` 同一关系式，按文献自行实现。
"""

import numpy as np

from pymeteo.units import (
    ArrayLike,
    ArrayOrScalar,
    UnitError,
    as_float_array,
    from_mps,
    restore_shape,
    to_kelvin,
    to_mps,
    to_pascal,
    to_radians,
)

# NCL coriolis_param / omega_to_w 所用常数
_EARTH_OMEGA = 7.292e-5  # rad·s⁻¹
_RD = 287.058  # J·K⁻¹·kg⁻¹
_G0 = 9.80665  # m·s⁻²


def coriolis_parameter(
    latitude: ArrayLike,
    *,
    latitude_unit: str = "deg",
) -> ArrayOrScalar:
    """计算科里奥利参数 ``f = 2 Ω sin φ``。

    参数
    ----
    latitude:
        地理纬度。
    latitude_unit:
        纬度单位，默认 ``deg``。亦可为 ``rad``。

    返回
    ----
    float 或 ndarray
        科里奥利参数，单位 s⁻¹。``Ω = 7.292×10⁻⁵ rad·s⁻¹``（与 NCL
        ``coriolis_param`` 相同）。赤道为 0，45° 约为 ``1.031×10⁻⁴``，
        北极约为 ``1.458×10⁻⁴``。
    """

    phi = to_radians(latitude, latitude_unit)
    return restore_shape(2.0 * _EARTH_OMEGA * np.sin(phi), latitude)


def _to_pascal_per_second(value: ArrayLike, unit: str) -> np.ndarray:
    """把垂直速度 ω 换算为 Pa/s。"""

    key = unit.strip().lower().replace(" ", "")
    arr = as_float_array(value)
    if key in {"pa/s", "pas-1", "pa s-1"}:
        return arr
    if key in {"hpa/s", "mb/s"}:
        return arr * 100.0
    raise UnitError(f"无法识别的ω单位 {unit!r}。支持：Pa/s、hPa/s")


def _from_pascal_per_second(value: ArrayLike, unit: str) -> np.ndarray:
    """把 Pa/s 换算为指定 ω 单位。"""

    key = unit.strip().lower().replace(" ", "")
    arr = as_float_array(value)
    if key in {"pa/s", "pas-1", "pa s-1"}:
        return arr
    if key in {"hpa/s", "mb/s"}:
        return arr / 100.0
    raise UnitError(f"无法识别的ω单位 {unit!r}。支持：Pa/s、hPa/s")


def omega_to_w(
    omega: ArrayLike,
    temperature: ArrayLike,
    pressure: ArrayLike,
    *,
    omega_unit: str = "Pa/s",
    temperature_unit: str = "C",
    pressure_unit: str = "hPa",
    output_speed_unit: str = "m/s",
) -> ArrayOrScalar:
    """把气压坐标垂直速度 ω 转为几何垂直速度 w。

    参数
    ----
    omega:
        垂直速度 ω（下沉为正）。
    temperature:
        空气温度。
    pressure:
        气压。
    omega_unit:
        ω 单位，默认 ``Pa/s``。支持 ``hPa/s``。
    temperature_unit:
        温度单位，默认 ``C``。
    pressure_unit:
        气压单位，默认 ``hPa``。
    output_speed_unit:
        输出 w 单位，默认 ``m/s``。

    返回
    ----
    float 或 ndarray
        ``w = -ω / (ρ g)``，``ρ = p / (R_d T)``，``R_d = 287.058``，
        ``g = 9.80665``。静力、干空气理想气体近似；对应 NCL
        ``omega_to_w``。现代 API 参数顺序为 ``(omega, temperature, pressure)``，
        与 NCL ``(omega, p, t)`` 不同。
    """

    omega_pas = _to_pascal_per_second(omega, omega_unit)
    temperature_k = to_kelvin(temperature, temperature_unit)
    pressure_pa = to_pascal(pressure, pressure_unit)
    density = pressure_pa / (_RD * temperature_k)
    w_mps = -omega_pas / (density * _G0)
    return restore_shape(from_mps(w_mps, output_speed_unit), omega, temperature, pressure)


def w_to_omega(
    w: ArrayLike,
    temperature: ArrayLike,
    pressure: ArrayLike,
    *,
    speed_unit: str = "m/s",
    temperature_unit: str = "C",
    pressure_unit: str = "hPa",
    output_omega_unit: str = "Pa/s",
) -> ArrayOrScalar:
    """把几何垂直速度 w 转为气压坐标垂直速度 ω。

    为 :func:`omega_to_w` 的代数逆：``ω = -ρ g w``。参数顺序为
    ``(w, temperature, pressure)``。
    """

    w_mps = to_mps(w, speed_unit)
    temperature_k = to_kelvin(temperature, temperature_unit)
    pressure_pa = to_pascal(pressure, pressure_unit)
    density = pressure_pa / (_RD * temperature_k)
    omega_pas = -density * _G0 * w_mps
    return restore_shape(
        _from_pascal_per_second(omega_pas, output_omega_unit), w, temperature, pressure
    )
