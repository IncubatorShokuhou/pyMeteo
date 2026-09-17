"""风场诊断：风速、风向及 u/v 分量互换。

风向采用气象学惯例：来向，正北为 0°，顺时针增加。分量满足

* ``u = -speed · sin(direction)``
* ``v = -speed · cos(direction)``
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from pymeteo.units import (
    ArrayOrScalar,
    as_float_array,
    from_mps,
    restore_shape,
    to_mps,
)


def wind_speed(
    u: ArrayLike,
    v: ArrayLike,
    *,
    speed_unit: str = "m/s",
    output_speed_unit: str | None = None,
) -> ArrayOrScalar:
    """由 u、v 分量计算风速。

    参数
    ----
    u, v:
        纬向与经向风分量。
    speed_unit:
        输入分量单位，默认 ``m/s``。支持 ``kt`` / ``knots``、``km/h``、``mph``。
    output_speed_unit:
        输出风速单位；默认与 ``speed_unit`` 相同。

    返回
    ----
    float 或 ndarray
        水平风速 ``sqrt(u² + v²)``。
    """

    if output_speed_unit is None:
        output_speed_unit = speed_unit
    u_mps = to_mps(u, speed_unit)
    v_mps = to_mps(v, speed_unit)
    speed = np.hypot(u_mps, v_mps)
    return restore_shape(from_mps(speed, output_speed_unit), u, v)


def wind_direction(
    u: ArrayLike,
    v: ArrayLike,
    *,
    speed_unit: str = "m/s",
) -> ArrayOrScalar:
    """由 u、v 分量计算气象风向（度，来向）。

    参数
    ----
    u, v:
        风分量。
    speed_unit:
        分量单位，仅用于接口一致；风向与单位无关。

    返回
    ----
    float 或 ndarray
        风向，范围 ``[0, 360)``。静风（风速为 0）时返回 0。
    """

    u_mps = to_mps(u, speed_unit)
    v_mps = to_mps(v, speed_unit)
    speed = np.hypot(u_mps, v_mps)
    direction = np.degrees(np.arctan2(u_mps, v_mps)) + 180.0
    direction = np.mod(direction, 360.0)
    direction = np.where(speed == 0.0, 0.0, direction)
    return restore_shape(direction, u, v)


def uv_from_speed_direction(
    speed: ArrayLike,
    direction: ArrayLike,
    *,
    speed_unit: str = "m/s",
    output_speed_unit: str | None = None,
) -> tuple[ArrayOrScalar, ArrayOrScalar]:
    """由风速和气象风向计算 u、v 分量。

    参数
    ----
    speed:
        水平风速。
    direction:
        风向，单位为度（来向）。
    speed_unit:
        输入风速单位，默认 ``m/s``。
    output_speed_unit:
        输出分量单位；默认与 ``speed_unit`` 相同。

    返回
    ----
    (u, v)
        纬向与经向分量。
    """

    if output_speed_unit is None:
        output_speed_unit = speed_unit
    speed_mps = to_mps(speed, speed_unit)
    direction_rad = np.deg2rad(as_float_array(direction))
    u_mps = -speed_mps * np.sin(direction_rad)
    v_mps = -speed_mps * np.cos(direction_rad)
    u_out = from_mps(u_mps, output_speed_unit)
    v_out = from_mps(v_mps, output_speed_unit)
    return restore_shape(u_out, speed, direction), restore_shape(v_out, speed, direction)


def wind_components(
    speed: ArrayLike,
    direction: ArrayLike,
    *,
    speed_unit: str = "m/s",
    output_speed_unit: str | None = None,
) -> tuple[ArrayOrScalar, ArrayOrScalar]:
    """由风速风向得到 u、v 分量，与 :func:`uv_from_speed_direction` 相同。"""

    return uv_from_speed_direction(
        speed,
        direction,
        speed_unit=speed_unit,
        output_speed_unit=output_speed_unit,
    )
