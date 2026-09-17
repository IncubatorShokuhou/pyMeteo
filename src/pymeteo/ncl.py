"""NCL 函数名兼容层：固定 NCL 单位，转调现代 ``pymeteo`` API。

本模块**不**实现独立物理公式。各封装只翻译参数与单位约定，再调用
``pymeteo.thermo`` / ``pymeteo.wind`` / ``pymeteo.dynamics`` 中的同名能力。
需要灵活单位时请使用顶层现代函数，不要从 ``pymeteo`` 顶层导入这些 NCL 名字。
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

import pymeteo.dynamics as _dynamics
import pymeteo.thermo as _thermo
import pymeteo.wind as _wind
from pymeteo.thermo import HumidityKind
from pymeteo.units import ArrayOrScalar, as_float_array, restore_shape, to_rh_fraction

__all__ = [
    "coriolis_param",
    "dewtemp_trh",
    "lclvl",
    "mixhum_convert",
    "mixhum_ptd",
    "mixhum_ptrh",
    "omega_to_w",
    "pot_temp",
    "pot_temp_equiv",
    "relhum",
    "relhum_ttd",
    "temp_virtual",
    "vapor_pres_rh",
    "w_to_omega",
    "wetbulb_stull",
    "wind_component",
    "wind_direction",
    "wind_speed",
]

_NCL_TEMP_UNIT = {0: "C", 1: "K", 2: "F"}
_NCL_MASS_UNIT = {0: "kg/kg", 1: "g/kg"}


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


def mixhum_ptd(p: ArrayLike, tdk: ArrayLike, iswit: ArrayLike) -> ArrayOrScalar:
    """由气压和露点计算混合比或比湿（NCL ``mixhum_ptd`` 兼容封装）。

    这是 NCL 兼容薄封装。NCL 约定：``p`` 为 Pa、``tdk`` 为开尔文露点。
    ``iswit`` 与 ``mixhum_ptrh`` 相同：绝对值为 1 混合比、为 2 比湿；正号
    kg/kg，负号 g/kg。灵活单位请用 :func:`pymeteo.mixing_ratio_from_dewpoint`。
    """

    flag = int(_scalar_number(iswit, "iswit"))
    kind = abs(flag)
    output_humidity_unit = "g/kg" if flag < 0 else "kg/kg"
    mixing = _thermo.mixing_ratio_from_dewpoint(
        p,
        tdk,
        pressure_unit="Pa",
        temperature_unit="K",
        output_humidity_unit="kg/kg",
    )
    if kind == 1:
        return _thermo.convert_humidity(
            mixing,
            from_quantity="mixing_ratio",
            to_quantity="mixing_ratio",
            humidity_unit="kg/kg",
            output_humidity_unit=output_humidity_unit,
        )
    if kind == 2:
        return _thermo.convert_humidity(
            mixing,
            from_quantity="mixing_ratio",
            to_quantity="specific_humidity",
            humidity_unit="kg/kg",
            output_humidity_unit=output_humidity_unit,
        )
    raise ValueError("iswit 的绝对值必须是 1（混合比）或 2（比湿）")


def vapor_pres_rh(rh: ArrayLike, es: ArrayLike) -> ArrayOrScalar:
    """由相对湿度和饱和水汽压计算水汽压（NCL ``vapor_pres_rh`` 兼容封装）。

    这是 NCL 兼容薄封装：``e = (RH%/100) · e_s``，``es`` 与返回值单位相同。
    由温度求水汽压请用 :func:`pymeteo.vapor_pressure_from_relative_humidity`。
    """

    rh_fraction = to_rh_fraction(rh, "%")
    return restore_shape(as_float_array(es) * as_float_array(rh_fraction), rh, es)


def pot_temp(p: ArrayLike, t: ArrayLike, dim: ArrayLike = -1, opt: bool = False) -> ArrayOrScalar:
    """计算位温（NCL ``pot_temp`` 兼容封装）。

    这是 NCL 兼容薄封装。``p`` 为 Pa、``t`` 为 K，返回 K。``dim`` / ``opt``
    仅为兼容 NCL 签名，NumPy 广播下忽略。灵活单位请用
    :func:`pymeteo.potential_temperature`。
    """

    _ = dim, opt
    return _thermo.potential_temperature(
        p,
        t,
        pressure_unit="Pa",
        temperature_unit="K",
        output_temperature_unit="K",
    )


def pot_temp_equiv(
    p: ArrayLike,
    t: ArrayLike,
    w: ArrayLike,
    dim: ArrayLike = -1,
    humVarType: str = "r",
) -> ArrayOrScalar:
    """计算相当位温（NCL ``pot_temp_equiv`` 兼容封装，内部仍用 Bolton 1980）。

    这是薄封装。NCL 约定：``p`` 为 Pa、``t`` 为 K；``w`` 的含义由
    ``humVarType`` 决定：``"r"`` / ``"w"`` 为混合比 kg/kg，``"q"`` 为比湿
    kg/kg，``"rh"`` 为相对湿度百分数。``dim`` 仅为兼容签名，忽略。
    物理上转调 Bolton 式 (43)（含 LCL），比 NCL 6.4 无 LCL 近似更接近
    ``pot_temp_equiv_tlcl``。灵活单位请用
    :func:`pymeteo.equivalent_potential_temperature`。
    """

    _ = dim
    kind = str(humVarType).strip().lower()
    if kind in {"r", "w"}:
        dewpoint = _dewpoint_k_from_mixing_ratio(p, t, w)
    elif kind == "q":
        mixing = _thermo.convert_humidity(
            w,
            from_quantity="specific_humidity",
            to_quantity="mixing_ratio",
            humidity_unit="kg/kg",
            output_humidity_unit="kg/kg",
        )
        dewpoint = _dewpoint_k_from_mixing_ratio(p, t, mixing)
    elif kind == "rh":
        dewpoint = _thermo.dewpoint_from_relative_humidity(
            t,
            w,
            temperature_unit="K",
            humidity_unit="%",
            output_temperature_unit="K",
        )
    else:
        raise ValueError("humVarType 必须是 'r'/'w'（混合比）、'q'（比湿）或 'rh'（相对湿度 %）")
    return _thermo.equivalent_potential_temperature(
        p,
        t,
        dewpoint,
        pressure_unit="Pa",
        temperature_unit="K",
        output_temperature_unit="K",
    )


def _dewpoint_k_from_mixing_ratio(
    pressure_pa: ArrayLike, temperature_k: ArrayLike, mixing: ArrayLike
) -> ArrayOrScalar:
    """由混合比（kg/kg）反演开尔文露点，供 NCL ``pot_temp_equiv`` 转调。"""

    vapor_hpa = _thermo.vapor_pressure_from_mixing_ratio(
        pressure_pa,
        mixing,
        pressure_unit="Pa",
        mixing_ratio_unit="kg/kg",
        output_pressure_unit="hPa",
    )
    dewpoint_c = _thermo._dewpoint_c_from_vapor_pressure_hpa(vapor_hpa)
    return restore_shape(dewpoint_c + 273.15, pressure_pa, temperature_k, mixing)


def temp_virtual(t: ArrayLike, w: ArrayLike, iounit: ArrayLike) -> ArrayOrScalar:
    """计算虚温（NCL ``temp_virtual`` 兼容封装）。

    这是 NCL 兼容薄封装。``iounit`` 长度为 3：下标 0 为输入温度（0=°C、
    1=K、2=°F），下标 1 为混合比（0=kg/kg、1=g/kg），下标 2 为输出温度。
    现代实现用 ``T_v = T (1 + r/ε)/(1 + r)``，与 NCL 文档中 ``T(1+0.61 r)``
    近似略有差别。灵活单位请用 :func:`pymeteo.virtual_temperature`。
    """

    flags = np.asarray(iounit, dtype=int).reshape(-1)
    if flags.size != 3:
        raise ValueError("iounit 必须是长度为 3 的整数序列")
    try:
        temperature_unit = _NCL_TEMP_UNIT[int(flags[0])]
        mixing_ratio_unit = _NCL_MASS_UNIT[int(flags[1])]
        output_temperature_unit = _NCL_TEMP_UNIT[int(flags[2])]
    except KeyError as exc:
        raise ValueError("iounit 取值：温度 0/1/2（C/K/F），混合比 0/1（kg/kg、g/kg）") from exc
    return _thermo.virtual_temperature(
        t,
        w,
        temperature_unit=temperature_unit,
        mixing_ratio_unit=mixing_ratio_unit,
        output_temperature_unit=output_temperature_unit,
    )


def wetbulb_stull(
    t: ArrayLike, rh: ArrayLike, iounit: ArrayLike, opt: bool = False
) -> ArrayOrScalar:
    """海平面湿球温度（NCL ``wetbulb_stull`` 兼容封装）。

    这是 NCL 兼容薄封装。``rh`` 为百分数。``iounit`` 长度为 2：下标 0 为
    输入温度、下标 1 为输出温度（0=°C、1=K、2=°F）。``opt`` 未使用。
    灵活单位请用 :func:`pymeteo.wet_bulb_temperature`。
    """

    _ = opt
    flags = np.asarray(iounit, dtype=int).reshape(-1)
    if flags.size != 2:
        raise ValueError("iounit 必须是长度为 2 的整数序列")
    try:
        temperature_unit = _NCL_TEMP_UNIT[int(flags[0])]
        output_temperature_unit = _NCL_TEMP_UNIT[int(flags[1])]
    except KeyError as exc:
        raise ValueError("iounit 取值：0=°C、1=K、2=°F") from exc
    return _thermo.wet_bulb_temperature(
        t,
        rh,
        temperature_unit=temperature_unit,
        humidity_unit="%",
        output_temperature_unit=output_temperature_unit,
    )


def lclvl(p: ArrayLike, tk: ArrayLike, tdk: ArrayLike) -> ArrayOrScalar:
    """计算抬升凝结高度气压（NCL ``lclvl`` 兼容封装）。

    这是 NCL 兼容薄封装。``p`` 为 hPa，``tk`` / ``tdk`` 为开尔文，返回 LCL
    气压（hPa）。现代接口同时给出 LCL 温度，见
    :func:`pymeteo.lifting_condensation_level`。
    """

    pressure, _temperature = _thermo.lifting_condensation_level(
        p,
        tk,
        tdk,
        pressure_unit="hPa",
        temperature_unit="K",
        output_pressure_unit="hPa",
        output_temperature_unit="K",
    )
    return pressure


def coriolis_param(lat: ArrayLike) -> ArrayOrScalar:
    """计算科里奥利参数（NCL ``coriolis_param`` 兼容封装）。

    这是 NCL 兼容薄封装。纬度单位为度，返回 s⁻¹。灵活单位请用
    :func:`pymeteo.coriolis_parameter`。
    """

    return _dynamics.coriolis_parameter(lat, latitude_unit="deg")


def omega_to_w(omega: ArrayLike, p: ArrayLike, t: ArrayLike) -> ArrayOrScalar:
    """ω（Pa/s）转为 w（m/s）（NCL ``omega_to_w`` 兼容封装）。

    这是 NCL 兼容薄封装。参数顺序与 NCL 相同：``(omega, p, t)``，``p`` 为
    Pa、``t`` 为 K。现代 API 为 ``(omega, temperature, pressure)``。
    """

    return _dynamics.omega_to_w(
        omega,
        t,
        p,
        omega_unit="Pa/s",
        temperature_unit="K",
        pressure_unit="Pa",
        output_speed_unit="m/s",
    )


def w_to_omega(w: ArrayLike, p: ArrayLike, t: ArrayLike) -> ArrayOrScalar:
    """w（m/s）转为 ω（Pa/s）（NCL ``w_to_omega`` 兼容封装）。

    这是 NCL 兼容薄封装。参数顺序与 NCL 相同：``(w, p, t)``。
    """

    return _dynamics.w_to_omega(
        w,
        t,
        p,
        speed_unit="m/s",
        temperature_unit="K",
        pressure_unit="Pa",
        output_omega_unit="Pa/s",
    )
