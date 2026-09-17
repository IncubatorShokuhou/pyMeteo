"""NCL 函数名兼容层：固定 NCL 单位，转调现代 ``pymeteo`` API。

本模块**不**实现独立物理公式。各封装只翻译参数与单位约定，再调用
``pymeteo.thermo`` / ``pymeteo.wind`` 中的同名能力。需要灵活单位时请使用
顶层现代函数，不要从 ``pymeteo`` 顶层导入这些 NCL 名字。
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

import pymeteo.thermo as _thermo
import pymeteo.wind as _wind
from pymeteo.thermo import HumidityKind
from pymeteo.units import ArrayOrScalar, restore_shape

__all__ = [
    "dewtemp_trh",
    "mixhum_convert",
    "mixhum_ptrh",
    "relhum",
    "relhum_ttd",
    "wind_component",
    "wind_direction",
    "wind_speed",
]


def _scalar_number(value: ArrayLike, name: str) -> float:
    array = np.asarray(value)
    if array.size != 1:
        raise ValueError(f"{name} 必须是标量")
    return float(array.item())


def _iounit_to_mass_unit(flag: int, name: str) -> str:
    if flag == 0:
        return "kg/kg"
    if flag == 1:
        return "g/kg"
    raise ValueError(f"{name} 必须是 0（kg/kg）或 1（g/kg）")


def dewtemp_trh(tk: ArrayLike, rh: ArrayLike) -> ArrayOrScalar:
    """由温度和相对湿度计算露点（NCL ``dewtemp_trh`` 兼容封装）。

    这是 NCL 兼容薄封装，本身不含独立物理公式。NCL 约定：``tk`` 为开尔文、
    ``rh`` 为百分数，返回露点开尔文。灵活单位请用
    :func:`pymeteo.dewpoint_from_relative_humidity`。
    """

    return _thermo.dewpoint_from_relative_humidity(
        tk,
        rh,
        temperature_unit="K",
        humidity_unit="%",
        output_temperature_unit="K",
    )


def relhum_ttd(t: ArrayLike, td: ArrayLike, opt: ArrayLike) -> ArrayOrScalar:
    """由温度和露点计算相对湿度（NCL ``relhum_ttd`` 兼容封装）。

    这是 NCL 兼容薄封装。NCL 约定：``t`` / ``td`` 为开尔文；``opt=0`` 返回
    百分数，``opt=1`` 返回 0–1 小数。灵活单位请用
    :func:`pymeteo.relative_humidity_from_dewpoint`。
    """

    option = int(_scalar_number(opt, "opt"))
    if option == 0:
        output_humidity_unit = "%"
    elif option == 1:
        output_humidity_unit = "fraction"
    else:
        raise ValueError("opt 必须是 0（百分数）或 1（小数）")
    return _thermo.relative_humidity_from_dewpoint(
        t,
        td,
        temperature_unit="K",
        output_humidity_unit=output_humidity_unit,
    )


def relhum(t: ArrayLike, w: ArrayLike, p: ArrayLike) -> ArrayOrScalar:
    """由温度、混合比和气压计算相对湿度（NCL ``relhum`` 兼容封装）。

    这是 NCL 兼容薄封装。NCL 约定：``t`` 为开尔文、``w`` 为 kg/kg、``p`` 为
    帕斯卡，返回百分数。灵活单位请用
    :func:`pymeteo.relative_humidity_from_mixing_ratio`。
    """

    return _thermo.relative_humidity_from_mixing_ratio(
        t,
        w,
        p,
        temperature_unit="K",
        mixing_ratio_unit="kg/kg",
        pressure_unit="Pa",
        output_humidity_unit="%",
    )


def mixhum_ptrh(
    p: ArrayLike,
    tk: ArrayLike,
    rh: ArrayLike,
    iswit: ArrayLike,
) -> ArrayOrScalar:
    """由气压、温度和相对湿度计算混合比或比湿（NCL ``mixhum_ptrh`` 兼容封装）。

    这是 NCL 兼容薄封装。NCL 约定：``p`` 为 hPa、``tk`` 为开尔文、``rh`` 为
    百分数。``iswit`` 绝对值为 1 返回混合比、为 2 返回比湿；正数单位 kg/kg，
    负数单位 g/kg。灵活单位请用
    :func:`pymeteo.mixing_ratio_from_relative_humidity` 或
    :func:`pymeteo.specific_humidity_from_relative_humidity`。
    """

    flag = int(_scalar_number(iswit, "iswit"))
    kind = abs(flag)
    output_humidity_unit = "g/kg" if flag < 0 else "kg/kg"
    kwargs = {
        "pressure_unit": "hPa",
        "temperature_unit": "K",
        "humidity_unit": "%",
        "output_humidity_unit": output_humidity_unit,
    }
    if kind == 1:
        return _thermo.mixing_ratio_from_relative_humidity(p, tk, rh, **kwargs)
    if kind == 2:
        return _thermo.specific_humidity_from_relative_humidity(p, tk, rh, **kwargs)
    raise ValueError("iswit 的绝对值必须是 1（混合比）或 2（比湿）")


def mixhum_convert(wq: ArrayLike, wqType: ArrayLike, iounit: ArrayLike) -> ArrayOrScalar:
    """混合比与比湿互换（NCL ``mixhum_convert`` 兼容封装）。

    这是 NCL 兼容薄封装。``wqType`` 为 ``"w"`` / ``"W"`` 表示输入混合比、输出
    比湿；为 ``"q"`` / ``"Q"`` 则相反。``iounit`` 为长度 2 的整数序列：0 表示
    kg/kg，1 表示 g/kg，分别对应输入与输出。灵活单位请用
    :func:`pymeteo.convert_humidity`。
    """

    kind = str(np.asarray(wqType).item()).strip().lower()
    from_quantity: HumidityKind
    to_quantity: HumidityKind
    if kind == "w":
        from_quantity, to_quantity = "mixing_ratio", "specific_humidity"
    elif kind == "q":
        from_quantity, to_quantity = "specific_humidity", "mixing_ratio"
    else:
        raise ValueError("wqType 必须是 'w'/'W'（混合比→比湿）或 'q'/'Q'（比湿→混合比）")

    flags = np.asarray(iounit, dtype=int).reshape(-1)
    if flags.size != 2:
        raise ValueError("iounit 必须是长度为 2 的整数序列，如 (0, 0) 或 (1, 1)")
    humidity_unit = _iounit_to_mass_unit(int(flags[0]), "iounit[0]")
    output_humidity_unit = _iounit_to_mass_unit(int(flags[1]), "iounit[1]")
    return _thermo.convert_humidity(
        wq,
        from_quantity=from_quantity,
        to_quantity=to_quantity,
        humidity_unit=humidity_unit,
        output_humidity_unit=output_humidity_unit,
    )


def wind_speed(u: ArrayLike, v: ArrayLike) -> ArrayOrScalar:
    """由 u、v 分量计算风速（NCL ``wind_speed`` 兼容封装）。

    这是 NCL 兼容薄封装。NCL 约定：分量与风速均为 m/s。灵活单位请用
    :func:`pymeteo.wind_speed`。
    """

    return _wind.wind_speed(u, v, speed_unit="m/s", output_speed_unit="m/s")


def wind_direction(u: ArrayLike, v: ArrayLike, opt: ArrayLike = 0) -> ArrayOrScalar:
    """由 u、v 分量计算气象风向（NCL ``wind_direction`` 兼容封装）。

    这是 NCL 兼容薄封装。风向为来向、单位度。静风时 ``opt=0`` 返回 0，
    ``opt=1`` 返回 ``nan``，其它标量则作为静风填充值。现代接口见
    :func:`pymeteo.wind_direction`。
    """

    direction = _wind.wind_direction(u, v, speed_unit="m/s")
    option = _scalar_number(opt, "opt")
    speed = _wind.wind_speed(u, v, speed_unit="m/s")
    calm = np.asarray(speed) == 0.0
    if option == 0.0 or not np.any(calm):
        return direction
    fill = np.nan if option == 1.0 else option
    filled = np.where(calm, fill, np.asarray(direction, dtype=float))
    return restore_shape(filled, u, v)


def wind_component(
    wspd: ArrayLike,
    wdir: ArrayLike,
    opt: ArrayLike = 0,
) -> tuple[ArrayOrScalar, ArrayOrScalar]:
    """由风速和气象风向计算 u、v（NCL ``wind_component`` 兼容封装）。

    这是 NCL 兼容薄封装。NCL 中 ``opt`` 未使用，保留以兼容位置参数。返回
    ``(u, v)`` 元组（Python 惯用形式，而非 NCL 左侧长度为 2 的堆叠数组）。
    灵活单位请用 :func:`pymeteo.uv_from_speed_direction`。
    """

    _ = opt
    return _wind.uv_from_speed_direction(
        wspd, wdir, speed_unit="m/s", output_speed_unit="m/s"
    )
