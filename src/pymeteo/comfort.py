"""体感温度：热指数与风寒（NWS 业务公式）。

按 NOAA / NWS 公布的回归式自行实现，不依赖 MetPy。内部先换到华氏度与
英里每小时再计算，再换回调用方单位。NCL 无对应内建名，故不提供 ncl 封装。
"""

from typing import Optional

import numpy as np

from pymeteo.units import (
    ArrayLike,
    ArrayOrScalar,
    as_float_array,
    from_kelvin,
    from_mps,
    from_rh_fraction,
    restore_shape,
    to_kelvin,
    to_mps,
    to_rh_fraction,
)


def heat_index(
    temperature: ArrayLike,
    relative_humidity: ArrayLike,
    *,
    temperature_unit: str = "C",
    humidity_unit: str = "%",
    output_temperature_unit: Optional[str] = None,
) -> ArrayOrScalar:
    """计算热指数（Heat Index，体感炎热程度）。

    参数
    ----
    temperature:
        气温。
    relative_humidity:
        相对湿度。
    temperature_unit:
        温度单位，默认 ``C``。
    humidity_unit:
        相对湿度单位，默认 ``%``。
    output_temperature_unit:
        输出单位；默认与 ``temperature_unit`` 相同。

    返回
    ----
    float 或 ndarray
        热指数。算法来自 Rothfusz（1990）对 Steadman（1979）的多元回归，
        以及 NWS 技术附件 SR 90-23 的订正：

        1. 先算简化式 ``HI = 0.5 {T + 61 + (T-68)·1.2 + RH·0.094}``（``T`` 为
           °F，``RH`` 为百分数），再与气温平均；
        2. 若该均值 ≥ 80 °F，改用 Rothfusz 全式，并在低湿（RH<13%、
           80–112 °F）或高湿（RH>85%、80–87 °F）时加减 NWS 订正项。

        公式在 Steadman 原表范围外（极端高温高湿）不可靠。90 °F、60% 时
        全式约 100 °F。
    """

    if output_temperature_unit is None:
        output_temperature_unit = temperature_unit
    t_f = from_kelvin(to_kelvin(temperature, temperature_unit), "F")
    rh = from_rh_fraction(to_rh_fraction(relative_humidity, humidity_unit), "%")
    t_f, rh = np.broadcast_arrays(as_float_array(t_f), as_float_array(rh))

    simple = 0.5 * (t_f + 61.0 + ((t_f - 68.0) * 1.2) + (rh * 0.094))
    averaged = 0.5 * (simple + t_f)

    roth = (
        -42.379
        + 2.04901523 * t_f
        + 10.14333127 * rh
        - 0.22475541 * t_f * rh
        - 0.00683783 * t_f**2
        - 0.05481717 * rh**2
        + 0.00122874 * t_f**2 * rh
        + 0.00085282 * t_f * rh**2
        - 0.00000199 * t_f**2 * rh**2
    )
    adj_low = ((13.0 - rh) / 4.0) * np.sqrt(np.maximum(17.0 - np.abs(t_f - 95.0), 0.0) / 17.0)
    adj_high = ((rh - 85.0) / 10.0) * ((87.0 - t_f) / 5.0)
    use_low = (rh < 13.0) & (t_f >= 80.0) & (t_f <= 112.0)
    use_high = (rh > 85.0) & (t_f >= 80.0) & (t_f <= 87.0)
    roth = roth - np.where(use_low, adj_low, 0.0) + np.where(use_high, adj_high, 0.0)
    hi_f = np.where(averaged >= 80.0, roth, averaged)
    return restore_shape(
        from_kelvin(to_kelvin(hi_f, "F"), output_temperature_unit),
        temperature,
        relative_humidity,
    )


def wind_chill(
    temperature: ArrayLike,
    wind_speed: ArrayLike,
    *,
    temperature_unit: str = "C",
    speed_unit: str = "m/s",
    output_temperature_unit: Optional[str] = None,
) -> ArrayOrScalar:
    """计算风寒温度（NWS / Environment Canada 2001 公式）。

    参数
    ----
    temperature:
        气温。
    wind_speed:
        风速（业务公式针对约 10 m 高度的风）。
    temperature_unit:
        温度单位，默认 ``C``。
    speed_unit:
        风速单位，默认 ``m/s``。
    output_temperature_unit:
        输出单位；默认与 ``temperature_unit`` 相同。

    返回
    ----
    float 或 ndarray
        风寒温度。华氏度 / 英里每小时形式为

        ``WC = 35.74 + 0.6215 T - 35.75 V^{0.16} + 0.4275 T V^{0.16}``

        有效范围约 ``T ≤ 50 °F`` 且 ``V ≥ 3 mph``；超出范围仍返回公式值，
        但物理意义减弱。0 °F、10 mph 时约 -16 °F。
    """

    if output_temperature_unit is None:
        output_temperature_unit = temperature_unit
    t_f = from_kelvin(to_kelvin(temperature, temperature_unit), "F")
    v_mph = from_mps(to_mps(wind_speed, speed_unit), "mph")
    t_f, v_mph = np.broadcast_arrays(as_float_array(t_f), as_float_array(v_mph))
    speed_term = np.power(np.maximum(v_mph, 0.0), 0.16)
    wc_f = 35.74 + 0.6215 * t_f - 35.75 * speed_term + 0.4275 * t_f * speed_term
    return restore_shape(
        from_kelvin(to_kelvin(wc_f, "F"), output_temperature_unit),
        temperature,
        wind_speed,
    )
