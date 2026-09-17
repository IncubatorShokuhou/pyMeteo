"""对流稳定度指数。

K 指数、A 指数、全总指数为温度与露点的线性组合；沙氏指数沿用李社宏（1994）
湿绝热迭代；SWEAT 采用 Miller（1972）/ NWS 标准形式，并修正旧代码把露点当
开尔文、把风速当米每秒却套用节系数的问题。
"""

from typing import Optional, Union

import numpy as np

from pymeteo.thermo import (
    _condensation_temperature_c,
    _saturation_vapor_pressure_hpa,
    dewpoint_from_relative_humidity,
    parcel_temperature_at_pressure,
)
from pymeteo.units import (
    ArrayLike,
    ArrayOrScalar,
    as_float_array,
    canonical_temperature_unit,
    from_kelvin,
    from_pascal,
    restore_shape,
    to_kelvin,
    to_mps,
    to_pascal,
)
from pymeteo.wind import wind_direction, wind_speed

_T0 = 273.15


def temperature_dewpoint_depression(
    temperature: ArrayLike,
    dewpoint: ArrayLike,
    *,
    temperature_unit: str = "C",
    output_temperature_unit: str = "C",
) -> ArrayOrScalar:
    """计算温度露点差（T - Td）。

    参数
    ----
    temperature:
        空气温度。
    dewpoint:
        露点温度。
    temperature_unit:
        二者单位，默认 ``C``。
    output_temperature_unit:
        输出温差单位。温差在 ``C`` 与 ``K`` 下数值相同；选用 ``F`` 时按华氏度差
        （乘 9/5）返回。

    返回
    ----
    float 或 ndarray
        温度露点差。可替代旧的 ``ttd850`` / ``ttd700`` / ``ttd925``。
    """

    temperature_k = to_kelvin(temperature, temperature_unit)
    dewpoint_k = to_kelvin(dewpoint, temperature_unit)
    delta_k = temperature_k - dewpoint_k
    delta = _delta_temperature(delta_k, output_temperature_unit)
    return restore_shape(delta, temperature, dewpoint)


def layer_temperature_difference(
    temperature_lower: ArrayLike,
    temperature_upper: ArrayLike,
    *,
    temperature_unit: str = "C",
    output_temperature_unit: str = "C",
) -> ArrayOrScalar:
    """计算两层温度差（下层减上层）。

    典型用法是 850 hPa 与 500 hPa 的温度差，对应旧接口 ``tt500``。
    """

    lower_k = to_kelvin(temperature_lower, temperature_unit)
    upper_k = to_kelvin(temperature_upper, temperature_unit)
    delta = _delta_temperature(lower_k - upper_k, output_temperature_unit)
    return restore_shape(delta, temperature_lower, temperature_upper)


def _delta_temperature(delta_k: ArrayLike, output_unit: str) -> np.ndarray:
    """把开尔文温差换算为输出单位下的温差。"""

    kind = canonical_temperature_unit(output_unit)
    arr = as_float_array(delta_k)
    if kind in {"C", "K"}:
        return arr
    return arr * (9.0 / 5.0)


def k_index(
    temperature_850: ArrayLike,
    dewpoint_850: ArrayLike,
    temperature_700: ArrayLike,
    dewpoint_700: ArrayLike,
    temperature_500: ArrayLike,
    *,
    temperature_unit: str = "C",
) -> ArrayOrScalar:
    """计算 K 指数。

    K = T850 - T500 + Td850 - (T700 - Td700)。

    参数
    ----
    temperature_850, dewpoint_850:
        850 hPa 温度与露点。
    temperature_700, dewpoint_700:
        700 hPa 温度与露点。
    temperature_500:
        500 hPa 温度。
    temperature_unit:
        全部温度单位，默认 ``C``。

    返回
    ----
    float 或 ndarray
        K 指数。公式按摄氏度定义（Td850 以绝对温度加入），输入会先换算到
        °C，因此 ``temperature_unit="K"`` 与 ``"C"`` 得到同一指数。经验上：
        K<20 无雷暴；20–25 零星雷暴；25–30 分散雷暴；30–35 成片雷暴。
    """

    # K 指数按摄氏度定义（Td850 以绝对温度加入），先统一换到 °C 再组合。
    t850 = from_kelvin(to_kelvin(temperature_850, temperature_unit), "C")
    td850 = from_kelvin(to_kelvin(dewpoint_850, temperature_unit), "C")
    t700 = from_kelvin(to_kelvin(temperature_700, temperature_unit), "C")
    td700 = from_kelvin(to_kelvin(dewpoint_700, temperature_unit), "C")
    t500 = from_kelvin(to_kelvin(temperature_500, temperature_unit), "C")
    result = t850 - t500 + td850 - (t700 - td700)
    return restore_shape(
        result, temperature_850, dewpoint_850, temperature_700, dewpoint_700, temperature_500
    )


def a_index(
    temperature_850: ArrayLike,
    dewpoint_850: ArrayLike,
    temperature_700: ArrayLike,
    dewpoint_700: ArrayLike,
    temperature_500: ArrayLike,
    dewpoint_500: ArrayLike,
    *,
    temperature_unit: str = "C",
) -> ArrayOrScalar:
    """计算 A 指数。

    A = (T850 - T500) - (T850 - Td850) - (T700 - Td700) - (T500 - Td500)。
    """

    t850 = to_kelvin(temperature_850, temperature_unit)
    td850 = to_kelvin(dewpoint_850, temperature_unit)
    t700 = to_kelvin(temperature_700, temperature_unit)
    td700 = to_kelvin(dewpoint_700, temperature_unit)
    t500 = to_kelvin(temperature_500, temperature_unit)
    td500 = to_kelvin(dewpoint_500, temperature_unit)
    result = (t850 - t500) - (t850 - td850) - (t700 - td700) - (t500 - td500)
    return restore_shape(
        result,
        temperature_850,
        dewpoint_850,
        temperature_700,
        dewpoint_700,
        temperature_500,
        dewpoint_500,
    )


def total_totals_index(
    temperature_850: ArrayLike,
    temperature_500: ArrayLike,
    *,
    dewpoint_850: Optional[ArrayLike] = None,
    relative_humidity_850: Optional[ArrayLike] = None,
    temperature_unit: str = "C",
    humidity_unit: str = "%",
) -> ArrayOrScalar:
    """计算全总指数 TT = T850 + Td850 - 2·T500。

    参数
    ----
    temperature_850, temperature_500:
        850 hPa 与 500 hPa 温度。
    dewpoint_850:
        850 hPa 露点。与 ``relative_humidity_850`` 必须提供其一。
    relative_humidity_850:
        850 hPa 相对湿度；若给出则按 Dutton 公式反推露点。
    temperature_unit:
        温度单位，默认 ``C``。
    humidity_unit:
        相对湿度单位，默认 ``%``。

    返回
    ----
    float 或 ndarray
        全总指数。温度无论用摄氏度还是开尔文，TT 数值相同。
    """

    t850 = to_kelvin(temperature_850, temperature_unit)
    t500 = to_kelvin(temperature_500, temperature_unit)
    td850 = _dewpoint_850_k(
        t850, dewpoint_850, relative_humidity_850, temperature_unit, humidity_unit
    )
    result = t850 + td850 - 2.0 * t500
    extras = [item for item in (dewpoint_850, relative_humidity_850) if item is not None]
    return restore_shape(result, temperature_850, temperature_500, *extras)


def _dewpoint_850_k(
    temperature_850_k: ArrayLike,
    dewpoint_850: Optional[ArrayLike],
    relative_humidity_850: Optional[ArrayLike],
    temperature_unit: str,
    humidity_unit: str,
) -> np.ndarray:
    """把 850 hPa 露点或相对湿度统一为开尔文露点。"""

    if dewpoint_850 is not None and relative_humidity_850 is not None:
        raise ValueError("dewpoint_850 与 relative_humidity_850 只能提供其中一个")
    if dewpoint_850 is not None:
        return to_kelvin(dewpoint_850, temperature_unit)
    if relative_humidity_850 is None:
        raise ValueError("必须提供 dewpoint_850 或 relative_humidity_850")
    return to_kelvin(
        dewpoint_from_relative_humidity(
            temperature_850_k,
            relative_humidity_850,
            temperature_unit="K",
            humidity_unit=humidity_unit,
            output_temperature_unit="K",
        ),
        "K",
    )


def showalter_index(
    temperature_850: ArrayLike,
    dewpoint_850: ArrayLike,
    temperature_500: ArrayLike,
    *,
    temperature_unit: str = "C",
) -> ArrayOrScalar:
    """计算沙氏指数。

    参数
    ----
    temperature_850, dewpoint_850:
        850 hPa 温度与露点。
    temperature_500:
        500 hPa 环境温度。
    temperature_unit:
        温度单位，默认 ``C``。

    返回
    ----
    float 或 ndarray
        沙氏指数，等于 500 hPa 环境温度减去从 850 hPa 沿湿绝热抬升到 500 hPa
        的气块温度，单位为 K（与摄氏度温差相同）。

    备注
    ----
    算法来自李社宏（1994）。旧说明指出：网格插值若给出负混合比，迭代会失败。
    示例：T850=16.6 °C、Td850=0.6 °C、T500=-15.9 °C 时约为 1.1。
    """

    t8_c = from_kelvin(to_kelvin(temperature_850, temperature_unit), "C")
    td8_c = from_kelvin(to_kelvin(dewpoint_850, temperature_unit), "C")
    t5_c = from_kelvin(to_kelvin(temperature_500, temperature_unit), "C")
    result = _showalter_c(t8_c, td8_c, t5_c)
    return restore_shape(result, temperature_850, dewpoint_850, temperature_500)


def _showalter_c(t8_c: ArrayLike, td8_c: ArrayLike, t5_c: ArrayLike) -> np.ndarray:
    """沙氏指数内部实现，输入输出均为摄氏度温差。"""

    t8, td8, t5 = np.broadcast_arrays(
        as_float_array(t8_c), as_float_array(td8_c), as_float_array(t5_c)
    )
    p5 = 500.0
    p8 = 850.0
    cpd = 0.2403
    cpv = 0.445
    rd = 6.85578 * 0.01
    rw = 11.017874 * 0.01
    t0 = _T0
    l0 = 597.40
    c_const = 1.002
    c1 = 0.57
    t8_k = t0 + t8
    # 沿用旧实现：混合比用 850 hPa 气温的水面饱和水汽压，而非露点（与李社宏迭代格式一致）。
    etd8 = _saturation_vapor_pressure_hpa(t8)
    mixing = (rd / rw) * etd8 / (p8 - etd8)
    ml = (cpd * (1.0 + cpv * mixing / cpd)) / (rd * (1.0 + mixing / (rd / rw)))
    ta = _condensation_temperature_c(p8, t8, td8) + t0
    pa = p8 * np.power(ta / t8_k, ml)
    m2 = (cpd / rd) * (1.0 + c_const * mixing / cpd)

    def moist_entropy(pressure: Union[np.ndarray, float], temperature_k: np.ndarray) -> np.ndarray:
        e_val = _saturation_vapor_pressure_hpa(temperature_k - t0)
        return np.log((pressure - e_val) / np.power(temperature_k, m2)) - (0.622 / rd) * (
            ((l0 + c1 * (t0 - temperature_k)) / temperature_k) * (e_val / (pressure - e_val))
        )

    qc = moist_entropy(pa, ta)
    out = t8_k.copy()
    step = np.full(out.shape, 10.0, dtype=float)
    for _ in range(10_000):
        q5 = moist_entropy(p5, out)
        residual = np.abs(qc - q5)
        done = residual <= 0.0001
        if bool(np.all(done)):
            break
        too_low = (~done) & (qc > q5)
        too_high = (~done) & (qc <= q5)
        out = np.where(too_low, out - step, out)
        out = np.where(too_high, out + step - step / 5.0, out)
        step = np.where(too_high, step / 5.0, step)
    return t5 - (out - t0)


def sweat_index(
    temperature_850: ArrayLike,
    temperature_500: ArrayLike,
    u_850: ArrayLike,
    v_850: ArrayLike,
    u_500: ArrayLike,
    v_500: ArrayLike,
    *,
    dewpoint_850: Optional[ArrayLike] = None,
    relative_humidity_850: Optional[ArrayLike] = None,
    temperature_unit: str = "C",
    humidity_unit: str = "%",
    speed_unit: str = "m/s",
) -> ArrayOrScalar:
    """计算强天气威胁指数 SWEAT。

    SWEAT = 12·Td850 + 20·(TT-49) + 2·f850 + f500 + 125·(S+0.2)

    参数
    ----
    temperature_850, temperature_500:
        850 hPa 与 500 hPa 温度。
    u_850, v_850, u_500, v_500:
        对应层次的风分量。
    dewpoint_850, relative_humidity_850:
        850 hPa 湿度信息，二者择一。
    temperature_unit:
        温度单位，默认 ``C``。
    humidity_unit:
        相对湿度单位，默认 ``%``。
    speed_unit:
        风分量单位，默认 ``m/s``。支持 ``kt`` / ``knots``。

    返回
    ----
    float 或 ndarray
        SWEAT 指数（无量纲）。经验上：>300 有强对流潜势，>400 有龙卷潜势。

    备注
    ----
    * Td850 使用摄氏度；风速项使用节。旧实现把开尔文露点直接乘 12、把 m/s
      风速套 4 与 2 的系数，结果会偏大约一个量级，此处已按 NWS / Miller 1972
      修正。
    * 任一项为负则置零。切变项仅在以下条件同时满足时保留：850 hPa 风向
      130°–250°、500 hPa 风向 210°–310°、500 hPa 风向减 850 hPa 风向为正、
      两层风速均 ≥ 15 kt。
    """

    t850 = to_kelvin(temperature_850, temperature_unit)
    t500 = to_kelvin(temperature_500, temperature_unit)
    td850_k = _dewpoint_850_k(
        t850, dewpoint_850, relative_humidity_850, temperature_unit, humidity_unit
    )
    td850_c = from_kelvin(td850_k, "C")
    tt_val = t850 + td850_k - 2.0 * t500

    u8 = to_mps(u_850, speed_unit)
    v8 = to_mps(v_850, speed_unit)
    u5 = to_mps(u_500, speed_unit)
    v5 = to_mps(v_500, speed_unit)
    f850 = wind_speed(u8, v8, speed_unit="m/s", output_speed_unit="kt")
    f500 = wind_speed(u5, v5, speed_unit="m/s", output_speed_unit="kt")
    d850 = wind_direction(u8, v8, speed_unit="m/s")
    d500 = wind_direction(u5, v5, speed_unit="m/s")

    term_td = 12.0 * td850_c
    term_tt = 20.0 * (tt_val - 49.0)
    term_f8 = 2.0 * as_float_array(f850)
    term_f5 = as_float_array(f500)
    shear_angle = np.deg2rad(as_float_array(d500) - as_float_array(d850))
    term_shear = 125.0 * (np.sin(shear_angle) + 0.2)

    d850_a = as_float_array(d850)
    d500_a = as_float_array(d500)
    shear_ok = (
        (d850_a >= 130.0)
        & (d850_a <= 250.0)
        & (d500_a >= 210.0)
        & (d500_a <= 310.0)
        & ((d500_a - d850_a) > 0.0)
        & (term_f8 / 2.0 >= 15.0)
        & (term_f5 >= 15.0)
    )
    term_shear = np.where(shear_ok, term_shear, 0.0)

    parts = [term_td, term_tt, term_f8, term_f5, term_shear]
    parts = [np.where(part < 0.0, 0.0, part) for part in parts]
    result = parts[0] + parts[1] + parts[2] + parts[3] + parts[4]
    extras = [item for item in (dewpoint_850, relative_humidity_850) if item is not None]
    return restore_shape(
        result, temperature_850, temperature_500, u_850, v_850, u_500, v_500, *extras
    )


def lifted_index(
    temperature_500: ArrayLike,
    parcel_temperature_500: ArrayLike,
    *,
    temperature_unit: str = "C",
) -> ArrayOrScalar:
    """由 500 hPa 环境温度与气块温度计算抬升指数。

    ``LI = T_env(500) - T_parcel(500)``。负值表示气块暖于环境、层结不稳定。
    温度无论用摄氏度还是开尔文，LI 数值相同（都是开尔文温差）。

    参数
    ----
    temperature_500:
        500 hPa 环境温度。
    parcel_temperature_500:
        气块到达 500 hPa 时的温度（需由调用方完成抬升，或使用
        :func:`lifted_index_from_surface`）。
    temperature_unit:
        二者单位，默认 ``C``。
    """

    env_k = to_kelvin(temperature_500, temperature_unit)
    parcel_k = to_kelvin(parcel_temperature_500, temperature_unit)
    return restore_shape(env_k - parcel_k, temperature_500, parcel_temperature_500)


def lifted_index_from_surface(
    pressure: ArrayLike,
    temperature: ArrayLike,
    dewpoint: ArrayLike,
    temperature_500: ArrayLike,
    *,
    pressure_unit: str = "hPa",
    temperature_unit: str = "C",
) -> ArrayOrScalar:
    """由近地层气压、温度、露点与 500 hPa 环境温度计算抬升指数。

    气块先按 Bolton（1980）求 LCL，未饱和段干绝热、饱和段按李社宏（1994）
    湿熵抬到 500 hPa，再 ``LI = T_500 - T_parcel(500)``。这是最常用的地面
    抬升指数定义；完整 CAPE / 最不稳定气块等探空套件不在本库范围。
    目标层固定为 500 hPa，与 ``pressure_unit`` 无关。
    """

    target_500 = from_pascal(to_pascal(500.0, "hPa"), pressure_unit)
    parcel_500 = parcel_temperature_at_pressure(
        pressure,
        temperature,
        dewpoint,
        target_500,
        pressure_unit=pressure_unit,
        temperature_unit=temperature_unit,
        output_temperature_unit="K",
    )
    return lifted_index(
        to_kelvin(temperature_500, temperature_unit),
        parcel_500,
        temperature_unit="K",
    )
