"""地球几何与站点气压订正。

大圆距离采用 WGS84 Vincenty 反解，不依赖 geopy；重力加速度按纬度用国际
重力公式近似；海平面气压沿用原经验公式。
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from pymeteo.units import (
    ArrayOrScalar,
    as_float_array,
    from_kelvin,
    from_meters,
    from_pascal,
    restore_shape,
    to_kelvin,
    to_meters,
    to_pascal,
    to_radians,
)

# WGS84 椭球
_WGS84_A = 6378137.0
_WGS84_F = 1.0 / 298.257223563
_WGS84_B = _WGS84_A * (1.0 - _WGS84_F)
_MEAN_EARTH_RADIUS = 6371008.8


def earth_distance(
    latitude_1: ArrayLike,
    longitude_1: ArrayLike,
    latitude_2: ArrayLike,
    longitude_2: ArrayLike,
    *,
    angle_unit: str = "deg",
    output_distance_unit: str = "km",
) -> ArrayOrScalar:
    """计算两点间的大地线距离（WGS84 Vincenty 反解）。

    参数
    ----
    latitude_1, longitude_1:
        第一点纬度和经度。
    latitude_2, longitude_2:
        第二点纬度和经度。
    angle_unit:
        经纬度单位，默认 ``deg``。亦可为 ``rad``。
    output_distance_unit:
        输出距离单位，默认 ``km``。支持 ``m`` / ``meter``、``mi``、``nmi`` / ``nm``、
        ``ft`` 等。

    返回
    ----
    float 或 ndarray
        两点距离。

    备注
    ----
    旧实现通过 ``eval("a." + unit)`` 读取 geopy 属性，既不安全也引入硬依赖。
    此处纯 NumPy 实现 Vincenty；对极少数对跖点迭代失败时回退到平均地球半径
    的 haversine 公式。
    """

    lat1 = to_radians(latitude_1, angle_unit)
    lon1 = to_radians(longitude_1, angle_unit)
    lat2 = to_radians(latitude_2, angle_unit)
    lon2 = to_radians(longitude_2, angle_unit)
    meters = _vincenty_inverse(lat1, lon1, lat2, lon2)
    return restore_shape(
        from_meters(meters, output_distance_unit),
        latitude_1,
        longitude_1,
        latitude_2,
        longitude_2,
    )


def _haversine(
    lat1: np.ndarray, lon1: np.ndarray, lat2: np.ndarray, lon2: np.ndarray
) -> np.ndarray:
    """平均半径球面 haversine 距离（米），用作 Vincenty 失败时的回退。"""

    dlat = lat2 - lat1
    dlon = lon2 - lon1
    hav = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    return 2.0 * _MEAN_EARTH_RADIUS * np.arcsin(np.sqrt(np.clip(hav, 0.0, 1.0)))


def _vincenty_inverse(
    lat1: ArrayLike,
    lon1: ArrayLike,
    lat2: ArrayLike,
    lon2: ArrayLike,
) -> np.ndarray:
    """向量化 Vincenty 反解，返回米。"""

    lat1, lon1, lat2, lon2 = np.broadcast_arrays(
        as_float_array(lat1),
        as_float_array(lon1),
        as_float_array(lat2),
        as_float_array(lon2),
    )
    flattened = np.column_stack([lat1.ravel(), lon1.ravel(), lat2.ravel(), lon2.ravel()])
    # 逐点迭代：Vincenty λ 更新难以整阵收敛；测站对数在数千以内足够快。
    out = np.empty(flattened.shape[0], dtype=float)
    for i, (phi1, lambda1, phi2, lambda2) in enumerate(flattened):
        out[i] = _vincenty_scalar(float(phi1), float(lambda1), float(phi2), float(lambda2))
    return out.reshape(lat1.shape)


def _vincenty_scalar(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """单点对 Vincenty 反解。输入为弧度，输出为米。"""

    if lat1 == lat2 and lon1 == lon2:
        return 0.0

    a = _WGS84_A
    f = _WGS84_F
    b = _WGS84_B
    l_diff = lon2 - lon1
    u1 = np.arctan((1.0 - f) * np.tan(lat1))
    u2 = np.arctan((1.0 - f) * np.tan(lat2))
    sin_u1, cos_u1 = np.sin(u1), np.cos(u1)
    sin_u2, cos_u2 = np.sin(u2), np.cos(u2)

    lam = l_diff
    cos2_alpha = 0.0
    sin_sigma = 0.0
    cos_sigma = 0.0
    sigma = 0.0
    cos2_sigma_m = 0.0
    converged = False
    for _ in range(200):
        sin_lam = np.sin(lam)
        cos_lam = np.cos(lam)
        sin_sigma = np.sqrt(
            (cos_u2 * sin_lam) ** 2 + (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lam) ** 2
        )
        if sin_sigma == 0.0:
            return 0.0
        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lam
        sigma = np.arctan2(sin_sigma, cos_sigma)
        sin_alpha = cos_u1 * cos_u2 * sin_lam / sin_sigma
        cos2_alpha = 1.0 - sin_alpha**2
        if cos2_alpha != 0.0:
            cos2_sigma_m = cos_sigma - 2.0 * sin_u1 * sin_u2 / cos2_alpha
        else:
            cos2_sigma_m = 0.0
        c_val = f / 16.0 * cos2_alpha * (2.0 + f * (4.0 - 3.0 * cos2_alpha))
        lam_prev = lam
        lam = l_diff + (1.0 - c_val) * f * sin_alpha * (
            sigma
            + c_val
            * sin_sigma
            * (cos2_sigma_m + c_val * cos_sigma * (-1.0 + 2.0 * cos2_sigma_m**2))
        )
        if abs(lam - lam_prev) < 1.0e-12:
            converged = True
            break

    if not converged:
        return float(_haversine(np.array(lat1), np.array(lon1), np.array(lat2), np.array(lon2)))

    u2_val = cos2_alpha * (a**2 - b**2) / b**2
    a_series = 1.0 + u2_val / 16384.0 * (
        4096.0 + u2_val * (-768.0 + u2_val * (320.0 - 175.0 * u2_val))
    )
    b_series = u2_val / 1024.0 * (256.0 + u2_val * (-128.0 + u2_val * (74.0 - 47.0 * u2_val)))
    delta_sigma = (
        b_series
        * sin_sigma
        * (
            cos2_sigma_m
            + b_series
            / 4.0
            * (
                cos_sigma * (-1.0 + 2.0 * cos2_sigma_m**2)
                - b_series
                / 6.0
                * cos2_sigma_m
                * (-3.0 + 4.0 * sin_sigma**2)
                * (-3.0 + 4.0 * cos2_sigma_m**2)
            )
        )
    )
    return float(b * a_series * (sigma - delta_sigma))


def gravity(
    latitude: ArrayLike,
    *,
    latitude_unit: str = "deg",
) -> ArrayOrScalar:
    """按纬度计算重力加速度。

    参数
    ----
    latitude:
        地理纬度。
    latitude_unit:
        纬度单位，默认 ``deg``。亦可为 ``rad``。

    返回
    ----
    float 或 ndarray
        重力加速度，单位 m/s²。

    备注
    ----
    采用 ``g = 9.7803 · (1 + 0.0053024 sin²φ - 0.000005 sin²2φ)``。
    旧实现把纬度直接送入 ``math.sin``，相当于把度数当弧度，45° 处会算错。
    此处按地理纬度（度）先化为弧度。
    """

    phi = to_radians(latitude, latitude_unit)
    g_val = 9.7803 * (1.0 + 0.0053024 * np.sin(phi) ** 2 - 0.000005 * np.sin(2.0 * phi) ** 2)
    return restore_shape(g_val, latitude)


def sea_level_pressure(
    station_pressure: ArrayLike,
    station_height: ArrayLike,
    temperature: ArrayLike,
    temperature_12h_ago: ArrayLike,
    *,
    lapse_rate: float = 0.005,
    pressure_unit: str = "hPa",
    height_unit: str = "m",
    temperature_unit: str = "C",
    output_pressure_unit: str = "hPa",
) -> ArrayOrScalar:
    """由本站气压订正到海平面气压。

    参数
    ----
    station_pressure:
        本站气压。
    station_height:
        测站海拔。
    temperature:
        当前气温。
    temperature_12h_ago:
        12 小时前气温。
    lapse_rate:
        气柱订正用递减率，单位为 ``temperature_unit`` **每米**，默认 0.005。
        在默认摄氏度下即每 100 m 降低 0.5 °C。若 ``temperature_unit="F"``，
        请把递减率也改成华氏度每米，否则不要沿用 0.005。
    pressure_unit:
        输入气压单位，默认 ``hPa``。
    height_unit:
        海拔单位，默认 ``m``。
    temperature_unit:
        温度单位，默认 ``C``。
    output_pressure_unit:
        输出气压单位，默认 ``hPa``。

    返回
    ----
    float 或 ndarray
        海平面气压。

    备注
    ----
    气柱平均温度 ``tm = (t + t12) / 2 + lapse_rate · h / 2``，再按
    ``p0 = ph · 10^(h / (18400 · (1 + tm/273)))`` 订正。公式中的 273 与 18400
    沿用原经验常数，温度按摄氏度代入。
    """

    pressure_hpa = from_pascal(to_pascal(station_pressure, pressure_unit), "hPa")
    height_m = to_meters(station_height, height_unit)
    t_c = from_kelvin(to_kelvin(temperature, temperature_unit), "C")
    t12_c = from_kelvin(to_kelvin(temperature_12h_ago, temperature_unit), "C")
    lapse_c_per_m = float(
        from_kelvin(to_kelvin(lapse_rate, temperature_unit), "C")
        - from_kelvin(to_kelvin(0.0, temperature_unit), "C")
    )
    tm = (t_c + t12_c) / 2.0 + lapse_c_per_m * height_m / 2.0
    slp_hpa = pressure_hpa * 10.0 ** (height_m / (18400.0 * (1.0 + tm / 273.0)))
    return restore_shape(
        from_pascal(to_pascal(slp_hpa, "hPa"), output_pressure_unit),
        station_pressure,
        station_height,
        temperature,
        temperature_12h_ago,
    )


# 干空气气体常数与标准重力，用于压高公式（SI）
_RD = 287.058  # J·K⁻¹·kg⁻¹
_G0 = 9.80665  # m·s⁻²


def height_thickness(
    pressure_bottom: ArrayLike,
    pressure_top: ArrayLike,
    mean_temperature: ArrayLike,
    *,
    pressure_unit: str = "hPa",
    temperature_unit: str = "C",
    output_distance_unit: str = "m",
) -> ArrayOrScalar:
    """由压高方程计算两等压面之间的厚度。

    参数
    ----
    pressure_bottom:
        下层气压（较大）。
    pressure_top:
        上层气压（较小）。
    mean_temperature:
        气层平均温度；有水汽时应传入**虚温**。
    pressure_unit:
        气压单位，默认 ``hPa``。
    temperature_unit:
        温度单位，默认 ``C``。
    output_distance_unit:
        输出厚度单位，默认 ``m``。

    返回
    ----
    float 或 ndarray
        气层厚度 ``ΔZ = (R_d T̄ / g) ln(p_bottom / p_top)``，
        ``R_d = 287.058 J·K⁻¹·kg⁻¹``，``g = 9.80665 m·s⁻²``。
        静力平衡下的经典压高公式（hypsometric equation）。
    """

    p_bottom = to_pascal(pressure_bottom, pressure_unit)
    p_top = to_pascal(pressure_top, pressure_unit)
    t_mean = to_kelvin(mean_temperature, temperature_unit)
    thickness_m = (_RD * t_mean / _G0) * np.log(p_bottom / p_top)
    return restore_shape(
        from_meters(thickness_m, output_distance_unit),
        pressure_bottom,
        pressure_top,
        mean_temperature,
    )
