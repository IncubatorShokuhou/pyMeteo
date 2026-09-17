"""水汽热力学诊断。

饱和水汽压、抬升凝结温度沿用李社宏（1994）水面公式；相对湿度与露点之间的
互换沿用 Dutton 经验潜热公式；由混合比求相对湿度沿用 NCL ``relhum`` 查表；
由相对湿度求混合比 / 比湿沿用 NCL ``mixhum_ptrh`` 的 Tetens 公式。
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike

from pymeteo.units import (
    ArrayOrScalar,
    as_float_array,
    from_kelvin,
    from_kgkg,
    from_meters,
    from_pascal,
    from_rh_fraction,
    restore_shape,
    to_kelvin,
    to_kgkg,
    to_pascal,
    to_rh_fraction,
)

# 李社宏（1994）水面饱和水汽压常数（气压结果为 hPa）
_E0 = 6.1078
_T0 = 273.15
_CL = 0.57
_RW_CAL = 0.1101787372
_L0 = 597.4

# Dutton 水汽气体常数与单位换算
_GC_J = 461.5  # J K^-1 kg^-1
_GC_CAL = _GC_J / (1000.0 * 4.186)  # cal g^-1 K^-1

# NCL mixhum_ptrh / Tetens
_EP = 0.622
_ONEMEP = 0.378
_ES0 = 6.11
_TETENS_A = 17.269
_TETENS_B = 35.86

# NCL relhum 饱和水汽压表，对应 173.16 K 起每隔 1 K；表值乘 0.1 后为 Pa
_RELHUM_ES_TABLE = np.array(
    [
        0.01403,
        0.01719,
        0.02101,
        0.02561,
        0.03117,
        0.03784,
        0.04584,
        0.05542,
        0.06685,
        0.08049,
        0.09672,
        0.1160,
        0.1388,
        0.1658,
        0.1977,
        0.2353,
        0.2796,
        0.3316,
        0.3925,
        0.4638,
        0.5472,
        0.6444,
        0.7577,
        0.8894,
        1.042,
        1.220,
        1.425,
        1.662,
        1.936,
        2.252,
        2.615,
        3.032,
        3.511,
        4.060,
        4.688,
        5.406,
        6.225,
        7.159,
        8.223,
        9.432,
        10.80,
        12.36,
        14.13,
        16.12,
        18.38,
        20.92,
        23.80,
        27.03,
        30.67,
        34.76,
        39.35,
        44.49,
        50.26,
        56.71,
        63.93,
        71.98,
        80.97,
        90.98,
        102.1,
        114.5,
        128.3,
        143.6,
        160.6,
        179.4,
        200.2,
        223.3,
        248.8,
        276.9,
        307.9,
        342.1,
        379.8,
        421.3,
        466.9,
        517.0,
        572.0,
        632.3,
        698.5,
        770.9,
        850.2,
        937.0,
        1032.0,
        1146.6,
        1272.0,
        1408.1,
        1556.7,
        1716.9,
        1890.3,
        2077.6,
        2279.6,
        2496.7,
        2729.8,
        2980.0,
        3247.8,
        3534.1,
        3839.8,
        4164.8,
        4510.5,
        4876.9,
        5265.1,
        5675.2,
        6107.8,
        6566.2,
        7054.7,
        7575.3,
        8129.4,
        8719.2,
        9346.5,
        10013.0,
        10722.0,
        11474.0,
        12272.0,
        13119.0,
        14017.0,
        14969.0,
        15977.0,
        17044.0,
        18173.0,
        19367.0,
        20630.0,
        21964.0,
        23373.0,
        24861.0,
        26430.0,
        28086.0,
        29831.0,
        31671.0,
        33608.0,
        35649.0,
        37796.0,
        40055.0,
        42430.0,
        44927.0,
        47551.0,
        50307.0,
        53200.0,
        56236.0,
        59422.0,
        62762.0,
        66264.0,
        69934.0,
        73777.0,
        77802.0,
        82015.0,
        86423.0,
        91034.0,
        95855.0,
        100890.0,
        106160.0,
        111660.0,
        117400.0,
        123400.0,
        129650.0,
        136170.0,
        142980.0,
        150070.0,
        157460.0,
        165160.0,
        173180.0,
        181530.0,
        190220.0,
        199260.0,
        208670.0,
        218450.0,
        228610.0,
        239180.0,
        250160.0,
        261560.0,
        273400.0,
        285700.0,
        298450.0,
        311690.0,
        325420.0,
        339650.0,
        354410.0,
        369710.0,
        385560.0,
        401980.0,
        418980.0,
        436590.0,
        454810.0,
        473670.0,
        493170.0,
        513350.0,
        534220.0,
        555800.0,
        578090.0,
        601130.0,
        624940.0,
        649530.0,
        674920.0,
        701130.0,
        728190.0,
        756110.0,
        784920.0,
        814630.0,
        845280.0,
        876880.0,
        909450.0,
        943020.0,
        977610.0,
        1013250.0,
        1049940.0,
        1087740.0,
        1087740.0,
    ],
    dtype=float,
)

HumidityKind = Literal["mixing_ratio", "specific_humidity"]
VisibilityMethod = Literal["RUC", "FSL"]


def _saturation_vapor_pressure_hpa(temperature_c: ArrayLike) -> np.ndarray:
    """水面饱和水汽压（hPa），李社宏（1994）公式。"""

    td = as_float_array(temperature_c)
    temperature_k = _T0 + td
    exponent = ((_L0 + _CL * _T0) * (temperature_k - _T0)) / (_RW_CAL * _T0 * temperature_k)
    return _E0 * np.exp(exponent) * (_T0 / temperature_k) ** (_CL / _RW_CAL)


def saturation_vapor_pressure(
    temperature: ArrayLike,
    *,
    temperature_unit: str = "C",
    output_pressure_unit: str = "hPa",
) -> ArrayOrScalar:
    """计算水面饱和水汽压。

    参数
    ----
    temperature:
        温度（对饱和水汽压而言即空气温度；旧接口常传入露点以得实际水汽压）。
    temperature_unit:
        温度单位，默认 ``C``。支持 ``C`` / ``degC`` / ``celsius``、``K`` / ``kelvin``、
        ``F`` / ``fahrenheit`` 等别名。
    output_pressure_unit:
        输出气压单位，默认 ``hPa``。支持 ``hPa`` / ``mb``、``Pa``、``kPa``、``atm``。

    返回
    ----
    float 或 ndarray
        饱和水汽压，单位由 ``output_pressure_unit`` 指定。

    备注
    ----
    公式针对水面，未区分冰面。0 °C 时结果为 6.1078 hPa。
    """

    temperature_c = from_kelvin(to_kelvin(temperature, temperature_unit), "C")
    pressure_hpa = _saturation_vapor_pressure_hpa(temperature_c)
    pressure = from_pascal(to_pascal(pressure_hpa, "hPa"), output_pressure_unit)
    return restore_shape(pressure, temperature)


def condensation_temperature(
    pressure: ArrayLike,
    temperature: ArrayLike,
    dewpoint: ArrayLike,
    *,
    pressure_unit: str = "hPa",
    temperature_unit: str = "C",
    output_temperature_unit: str = "C",
) -> ArrayOrScalar:
    """根据温度和露点计算抬升凝结高度上的温度。

    参数
    ----
    pressure:
        气块起始气压。
    temperature:
        气块温度。
    dewpoint:
        气块露点温度。
    pressure_unit:
        气压单位，默认 ``hPa``。
    temperature_unit:
        温度与露点单位，默认 ``C``。
    output_temperature_unit:
        输出温度单位，默认 ``C``。

    返回
    ----
    float 或 ndarray
        抬升凝结温度。

    备注
    ----
    采用李社宏（1994）迭代格式。网格插值导致露点高于温度或出现负混合比时
    可能不收敛，此时返回最后一次迭代值。
    """

    pressure_hpa = from_pascal(to_pascal(pressure, pressure_unit), "hPa")
    temperature_c = from_kelvin(to_kelvin(temperature, temperature_unit), "C")
    dewpoint_c = from_kelvin(to_kelvin(dewpoint, temperature_unit), "C")
    result_c = _condensation_temperature_c(pressure_hpa, temperature_c, dewpoint_c)
    result = from_kelvin(to_kelvin(result_c, "C"), output_temperature_unit)
    return restore_shape(result, pressure, temperature, dewpoint)


def _condensation_temperature_c(
    pressure_hpa: ArrayLike,
    temperature_c: ArrayLike,
    dewpoint_c: ArrayLike,
) -> np.ndarray:
    """抬升凝结温度的向量化迭代（内部，摄氏度 / hPa）。"""

    p8, t_c, td_c = np.broadcast_arrays(
        as_float_array(pressure_hpa),
        as_float_array(temperature_c),
        as_float_array(dewpoint_c),
    )
    cpd = 0.2403
    cpv = 0.445
    rd = 6.85578 * 0.01
    rw = 11.017874 * 0.01
    t0 = _T0
    temperature_k = t0 + t_c
    dewpoint_k = t0 + td_c
    etd = _saturation_vapor_pressure_hpa(td_c)
    mixing_ratio = (rd / rw) * etd / (p8 - etd)
    ml = (cpd * (1.0 + cpv * mixing_ratio / cpd)) / (rd * (1.0 + mixing_ratio / (rd / rw)))
    z0 = np.power(temperature_k, ml) / etd
    out = dewpoint_k.copy()
    step = np.full(out.shape, 10.0, dtype=float)
    for _ in range(10_000):
        z_val = np.power(out, ml) / _saturation_vapor_pressure_hpa(out - t0)
        residual = np.abs(z_val - z0)
        done = (residual <= 10.0) | (step < 1.0e-8)
        if bool(np.all(done)):
            break
        too_small = (~done) & (z_val < z0)
        too_large = (~done) & (z_val >= z0)
        out = np.where(too_small, out - step, out)
        out = np.where(too_large, out + step - step / 5.0, out)
        step = np.where(too_large, step / 5.0, step)
    return out - t0


def relative_humidity_from_dewpoint(
    temperature: ArrayLike,
    dewpoint: ArrayLike,
    *,
    temperature_unit: str = "C",
    output_humidity_unit: str = "%",
) -> ArrayOrScalar:
    """由温度和露点计算相对湿度。

    参数
    ----
    temperature:
        空气温度。
    dewpoint:
        露点温度。
    temperature_unit:
        二者的温度单位，默认 ``C``。
    output_humidity_unit:
        输出相对湿度单位，默认 ``%``。亦可为 ``fraction``。

    返回
    ----
    float 或 ndarray
        相对湿度。公式取自 Dutton，潜热随温度线性变化。
    """

    temperature_k = to_kelvin(temperature, temperature_unit)
    dewpoint_k = to_kelvin(dewpoint, temperature_unit)
    latent = 597.3 - 0.57 * (temperature_k - 273.15)
    rh_fraction = np.exp((latent / _GC_CAL) * (1.0 / temperature_k - 1.0 / dewpoint_k))
    rh = from_rh_fraction(rh_fraction, output_humidity_unit)
    return restore_shape(rh, temperature, dewpoint)


def dewpoint_from_relative_humidity(
    temperature: ArrayLike,
    relative_humidity: ArrayLike,
    *,
    temperature_unit: str = "C",
    humidity_unit: str = "%",
    output_temperature_unit: str = "C",
) -> ArrayOrScalar:
    """由温度和相对湿度计算露点。

    参数
    ----
    temperature:
        空气温度。
    relative_humidity:
        相对湿度。
    temperature_unit:
        温度单位，默认 ``C``。
    humidity_unit:
        相对湿度单位，默认 ``%``。
    output_temperature_unit:
        输出露点单位，默认 ``C``。

    返回
    ----
    float 或 ndarray
        露点温度。相对湿度不大于 0 时对应元素为 ``nan``。

    备注
    ----
    旧版 ``__main__`` 示例：18 °C、相对湿度 46.5% 时露点约为 6.30 °C。
    """

    temperature_k = to_kelvin(temperature, temperature_unit)
    rh_fraction = to_rh_fraction(relative_humidity, humidity_unit)
    temperature_b, rh_b = np.broadcast_arrays(temperature_k, rh_fraction)
    latent = (597.3 - 0.57 * (temperature_b - 273.15)) / _GC_CAL
    safe = rh_b > 0.0
    dewpoint_k = np.full(temperature_b.shape, np.nan, dtype=float)
    dewpoint_k = np.where(
        safe,
        temperature_b * latent / (latent - temperature_b * np.log(np.clip(rh_b, 1.0e-15, None))),
        dewpoint_k,
    )
    dewpoint = from_kelvin(dewpoint_k, output_temperature_unit)
    return restore_shape(dewpoint, temperature, relative_humidity)


def relative_humidity_from_mixing_ratio(
    temperature: ArrayLike,
    mixing_ratio: ArrayLike,
    pressure: ArrayLike,
    *,
    temperature_unit: str = "C",
    mixing_ratio_unit: str = "kg/kg",
    pressure_unit: str = "hPa",
    output_humidity_unit: str = "%",
) -> ArrayOrScalar:
    """由温度、混合比和气压计算相对湿度（NCL ``relhum`` 查表）。

    参数
    ----
    temperature:
        空气温度。
    mixing_ratio:
        水汽混合比。
    pressure:
        气压。
    temperature_unit:
        温度单位，默认 ``C``。
    mixing_ratio_unit:
        混合比单位，默认 ``kg/kg``。支持 ``g/kg``。
    pressure_unit:
        气压单位，默认 ``hPa``。
    output_humidity_unit:
        输出相对湿度单位，默认 ``%``。

    返回
    ----
    float 或 ndarray
        相对湿度。允许大于 100%；小于 0 时截断为 0.0001%。

    备注
    ----
    饱和水汽压由 173.16–375.16 K 查表线性内插得到，与 NCL 一致。
    """

    temperature_k = np.clip(to_kelvin(temperature, temperature_unit), 173.16, 375.16)
    mixing = to_kgkg(mixing_ratio, mixing_ratio_unit)
    pressure_pa = to_pascal(pressure, pressure_unit)
    index = np.floor(temperature_k - 173.16).astype(int)
    index = np.clip(index, 0, _RELHUM_ES_TABLE.size - 2)
    t_low = 173.16 + index
    es_pa = ((t_low + 1.0 - temperature_k) * _RELHUM_ES_TABLE[index]
             + (temperature_k - t_low) * _RELHUM_ES_TABLE[index + 1]) * 0.1
    rh_percent = (mixing * (pressure_pa - 0.378 * es_pa) / (0.622 * es_pa)) * 100.0
    rh_percent = np.where(rh_percent < 0.0, 0.0001, rh_percent)
    rh = from_rh_fraction(rh_percent / 100.0, output_humidity_unit)
    return restore_shape(rh, temperature, mixing_ratio, pressure)


def mixing_ratio_from_relative_humidity(
    pressure: ArrayLike,
    temperature: ArrayLike,
    relative_humidity: ArrayLike,
    *,
    pressure_unit: str = "hPa",
    temperature_unit: str = "C",
    humidity_unit: str = "%",
    output_humidity_unit: str = "kg/kg",
) -> ArrayOrScalar:
    """由气压、温度和相对湿度计算水汽混合比。

    参数
    ----
    pressure:
        气压。
    temperature:
        空气温度。
    relative_humidity:
        相对湿度。
    pressure_unit:
        气压单位，默认 ``hPa``。
    temperature_unit:
        温度单位，默认 ``C``。
    humidity_unit:
        相对湿度单位，默认 ``%``。
    output_humidity_unit:
        输出混合比单位，默认 ``kg/kg``。

    返回
    ----
    float 或 ndarray
        混合比。饱和水汽压采用 Tetens 公式（NCL ``mixhum_ptrh``）。
    """

    mixing = _tetens_mixing_ratio(pressure, temperature, relative_humidity,
                                  pressure_unit, temperature_unit, humidity_unit)
    return restore_shape(from_kgkg(mixing, output_humidity_unit),
                         pressure, temperature, relative_humidity)


def specific_humidity_from_relative_humidity(
    pressure: ArrayLike,
    temperature: ArrayLike,
    relative_humidity: ArrayLike,
    *,
    pressure_unit: str = "hPa",
    temperature_unit: str = "C",
    humidity_unit: str = "%",
    output_humidity_unit: str = "kg/kg",
) -> ArrayOrScalar:
    """由气压、温度和相对湿度计算比湿。

    参数含义同 :func:`mixing_ratio_from_relative_humidity`。比湿 ``q = w / (1 + w)``。
    """

    mixing = _tetens_mixing_ratio(pressure, temperature, relative_humidity,
                                  pressure_unit, temperature_unit, humidity_unit)
    specific = mixing / (1.0 + mixing)
    return restore_shape(from_kgkg(specific, output_humidity_unit),
                         pressure, temperature, relative_humidity)


def _tetens_mixing_ratio(
    pressure: ArrayLike,
    temperature: ArrayLike,
    relative_humidity: ArrayLike,
    pressure_unit: str,
    temperature_unit: str,
    humidity_unit: str,
) -> np.ndarray:
    """Tetens 饱和水汽压下的混合比（kg/kg）。"""

    pressure_hpa = from_pascal(to_pascal(pressure, pressure_unit), "hPa")
    temperature_k = to_kelvin(temperature, temperature_unit)
    rh_fraction = to_rh_fraction(relative_humidity, humidity_unit)
    est = _ES0 * np.exp((_TETENS_A * (temperature_k - _T0)) / (temperature_k - _TETENS_B))
    qst = (_EP * est) / (pressure_hpa - _ONEMEP * est)
    return qst * rh_fraction


def convert_humidity(
    value: ArrayLike,
    *,
    from_quantity: HumidityKind = "mixing_ratio",
    to_quantity: HumidityKind = "specific_humidity",
    humidity_unit: str = "kg/kg",
    output_humidity_unit: str = "kg/kg",
) -> ArrayOrScalar:
    """在混合比与比湿之间转换，并同时换算质量单位。

    参数
    ----
    value:
        混合比或比湿。
    from_quantity:
        输入物理量：``mixing_ratio``（混合比）或 ``specific_humidity``（比湿）。
    to_quantity:
        输出物理量，取值同上。
    humidity_unit:
        输入质量单位，默认 ``kg/kg``。支持 ``g/kg``、``mg/kg``。
    output_humidity_unit:
        输出质量单位，默认 ``kg/kg``。

    返回
    ----
    float 或 ndarray
        转换后的水汽含量。关系为 ``q = w / (1 + w)``、``w = q / (1 - q)``。
    """

    if from_quantity not in ("mixing_ratio", "specific_humidity"):
        raise ValueError("from_quantity 必须是 'mixing_ratio' 或 'specific_humidity'")
    if to_quantity not in ("mixing_ratio", "specific_humidity"):
        raise ValueError("to_quantity 必须是 'mixing_ratio' 或 'specific_humidity'")

    amount = to_kgkg(value, humidity_unit)
    if from_quantity != to_quantity:
        if from_quantity == "mixing_ratio":
            amount = amount / (1.0 + amount)
        else:
            amount = amount / (1.0 - amount)
    return restore_shape(from_kgkg(amount, output_humidity_unit), value)


def visibility(
    relative_humidity: ArrayLike,
    temperature: ArrayLike,
    method: VisibilityMethod = "RUC",
    *,
    humidity_unit: str = "%",
    temperature_unit: str = "C",
    output_distance_unit: str = "km",
) -> ArrayOrScalar:
    """用相对湿度和温度估算能见度。

    参数
    ----
    relative_humidity:
        相对湿度。
    temperature:
        空气温度。
    method:
        ``RUC``（Rapid Update Cycle 指数衰减）或 ``FSL``
        （Forecast Systems Laboratory 露点差公式）。
    humidity_unit:
        相对湿度单位，默认 ``%``。
    temperature_unit:
        温度单位，默认 ``C``。
    output_distance_unit:
        输出距离单位，默认 ``km``。

    返回
    ----
    float 或 ndarray
        能见度。

    备注
    ----
    旧实现把 RUC 结果又乘以 1000，与“单位为 km”的说明矛盾，会得到上万公里的
    荒谬值。此处按 60 km 乘指数衰减给出千米级能见度，再换算到输出单位。
    """

    rh_fraction = to_rh_fraction(relative_humidity, humidity_unit)
    temperature_k = to_kelvin(temperature, temperature_unit)
    rh_percent = from_rh_fraction(rh_fraction, "%")
    if method == "RUC":
        decrement = np.minimum(80.0, rh_fraction - 0.15)
        vis_m = 60.0 * np.exp(-2.5 * decrement) * 1000.0
    elif method == "FSL":
        dewpoint_k = to_kelvin(
            dewpoint_from_relative_humidity(
                temperature_k,
                rh_percent,
                temperature_unit="K",
                humidity_unit="%",
                output_temperature_unit="K",
            ),
            "K",
        )
        vis_km = 6000.0 * (temperature_k - dewpoint_k) / np.power(rh_percent, 1.75)
        vis_m = vis_km * 1000.0
    else:
        raise ValueError("method 必须是 'RUC' 或 'FSL'")
    return restore_shape(from_meters(vis_m, output_distance_unit), relative_humidity, temperature)
