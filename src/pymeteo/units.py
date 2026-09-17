"""单位字符串解析与换算。

本模块把公开 API 接受的单位别名规范为内部计算量纲，再按需换算回调用方
指定的输出单位。内部统一采用：

* 温度：开尔文（K）
* 气压：帕斯卡（Pa）
* 风速：米每秒（m/s）
* 相对湿度：小数（0–1）
* 混合比 / 比湿：千克每千克（kg/kg）
* 距离：米（m）
* 角度：度（deg），三角函数计算前再转为弧度
"""

from typing import Any, Dict, Union

import numpy as np

try:
    from numpy.typing import ArrayLike
except ImportError:  # NumPy 1.19 (last line that fully supports CPython 3.6)
    ArrayLike = Any  # type: ignore[misc,assignment]

ArrayOrScalar = Union[float, np.ndarray]


class UnitError(ValueError):
    """无法识别或不适用于当前物理量的单位字符串。"""


def as_float_array(value: ArrayLike) -> np.ndarray:
    """将输入转为 ``float`` 类型的 ``ndarray``。"""

    return np.asarray(value, dtype=float)


def restore_shape(result: ArrayLike, *originals: ArrayLike) -> ArrayOrScalar:
    """若全部原始输入都是标量，则把结果还原为 Python ``float``。"""

    arr = np.asarray(result, dtype=float)
    if all(np.ndim(np.asarray(item)) == 0 for item in originals):
        return float(arr)
    return arr


def _strip_unit(unit: str) -> str:
    """去掉空白、下划线、百分号以外的装饰字符，便于匹配别名。"""

    text = unit.strip().lower()
    text = text.replace(" ", "").replace("_", "").replace("°", "")
    text = text.replace("℃", "c").replace("℉", "f")
    return text


def _lookup(unit: str, table: Dict[str, str], kind: str) -> str:
    key = _strip_unit(unit)
    if key not in table:
        supported = ", ".join(sorted(set(table)))
        raise UnitError(f"无法识别的{kind}单位 {unit!r}。支持：{supported}")
    return table[key]


_TEMPERATURE_ALIASES = {
    "c": "C",
    "degc": "C",
    "celsius": "C",
    "centigrade": "C",
    "k": "K",
    "kelvin": "K",
    "f": "F",
    "degf": "F",
    "fahrenheit": "F",
}

_PRESSURE_ALIASES = {
    "pa": "Pa",
    "pascal": "Pa",
    "pascals": "Pa",
    "hpa": "hPa",
    "mb": "hPa",
    "mbar": "hPa",
    "millibar": "hPa",
    "millibars": "hPa",
    "kpa": "kPa",
    "atm": "atm",
}

_SPEED_ALIASES = {
    "m/s": "m/s",
    "ms-1": "m/s",
    "m/sec": "m/s",
    "mps": "m/s",
    "kt": "kt",
    "knot": "kt",
    "knots": "kt",
    "kn": "kt",
    "km/h": "km/h",
    "kmh": "km/h",
    "kph": "km/h",
    "mph": "mph",
}

_RH_ALIASES = {
    "%": "%",
    "percent": "%",
    "percentage": "%",
    "fraction": "fraction",
    "1": "fraction",
    "ratio": "fraction",
}

_MASS_HUMIDITY_ALIASES = {
    "kg/kg": "kg/kg",
    "kgkg-1": "kg/kg",
    "g/g": "kg/kg",
    "g/kg": "g/kg",
    "gkg-1": "g/kg",
    "mg/kg": "mg/kg",
}

_DISTANCE_ALIASES = {
    "m": "m",
    "meter": "m",
    "meters": "m",
    "metre": "m",
    "metres": "m",
    "km": "km",
    "kilometer": "km",
    "kilometers": "km",
    "kilometre": "km",
    "kilometres": "km",
    "cm": "cm",
    "mi": "mi",
    "mile": "mi",
    "miles": "mi",
    "ft": "ft",
    "foot": "ft",
    "feet": "ft",
    "nmi": "nmi",
    "nm": "nmi",
    "nauticalmile": "nmi",
    "nauticalmiles": "nmi",
}

_ANGLE_ALIASES = {
    "deg": "deg",
    "degree": "deg",
    "degrees": "deg",
    "rad": "rad",
    "radian": "rad",
    "radians": "rad",
}

_KNOT_TO_MPS = 1852.0 / 3600.0
_MILE_TO_M = 1609.344
_NMI_TO_M = 1852.0
_ATM_TO_PA = 101325.0


def canonical_temperature_unit(unit: str) -> str:
    """返回温度单位的规范名：``C`` / ``K`` / ``F``。"""

    return _lookup(unit, _TEMPERATURE_ALIASES, "温度")


def canonical_pressure_unit(unit: str) -> str:
    """返回气压单位的规范名。"""

    return _lookup(unit, _PRESSURE_ALIASES, "气压")


def canonical_speed_unit(unit: str) -> str:
    """返回风速单位的规范名。"""

    return _lookup(unit, _SPEED_ALIASES, "风速")


def canonical_rh_unit(unit: str) -> str:
    """返回相对湿度单位的规范名：``%`` 或 ``fraction``。"""

    return _lookup(unit, _RH_ALIASES, "相对湿度")


def canonical_mass_humidity_unit(unit: str) -> str:
    """返回混合比 / 比湿单位的规范名。"""

    return _lookup(unit, _MASS_HUMIDITY_ALIASES, "水汽质量混合比")


def canonical_distance_unit(unit: str) -> str:
    """返回距离单位的规范名。"""

    return _lookup(unit, _DISTANCE_ALIASES, "距离")


def canonical_angle_unit(unit: str) -> str:
    """返回角度单位的规范名：``deg`` 或 ``rad``。"""

    return _lookup(unit, _ANGLE_ALIASES, "角度")


def to_kelvin(value: ArrayLike, unit: str) -> np.ndarray:
    """把温度换算为开尔文。"""

    arr = as_float_array(value)
    kind = canonical_temperature_unit(unit)
    if kind == "K":
        return arr
    if kind == "C":
        return arr + 273.15
    return (arr + 459.67) * (5.0 / 9.0)


def from_kelvin(value: ArrayLike, unit: str) -> np.ndarray:
    """把开尔文温度换算为指定单位。"""

    arr = as_float_array(value)
    kind = canonical_temperature_unit(unit)
    if kind == "K":
        return arr
    if kind == "C":
        return arr - 273.15
    return arr * (9.0 / 5.0) - 459.67


def to_pascal(value: ArrayLike, unit: str) -> np.ndarray:
    """把气压换算为帕斯卡。"""

    arr = as_float_array(value)
    kind = canonical_pressure_unit(unit)
    if kind == "Pa":
        return arr
    if kind == "hPa":
        return arr * 100.0
    if kind == "kPa":
        return arr * 1000.0
    return arr * _ATM_TO_PA


def from_pascal(value: ArrayLike, unit: str) -> np.ndarray:
    """把帕斯卡换算为指定气压单位。"""

    arr = as_float_array(value)
    kind = canonical_pressure_unit(unit)
    if kind == "Pa":
        return arr
    if kind == "hPa":
        return arr / 100.0
    if kind == "kPa":
        return arr / 1000.0
    return arr / _ATM_TO_PA


def to_mps(value: ArrayLike, unit: str) -> np.ndarray:
    """把风速换算为米每秒。"""

    arr = as_float_array(value)
    kind = canonical_speed_unit(unit)
    if kind == "m/s":
        return arr
    if kind == "kt":
        return arr * _KNOT_TO_MPS
    if kind == "km/h":
        return arr / 3.6
    return arr * _MILE_TO_M / 3600.0


def from_mps(value: ArrayLike, unit: str) -> np.ndarray:
    """把米每秒换算为指定风速单位。"""

    arr = as_float_array(value)
    kind = canonical_speed_unit(unit)
    if kind == "m/s":
        return arr
    if kind == "kt":
        return arr / _KNOT_TO_MPS
    if kind == "km/h":
        return arr * 3.6
    return arr * 3600.0 / _MILE_TO_M


def to_rh_fraction(value: ArrayLike, unit: str) -> np.ndarray:
    """把相对湿度换算为 0–1 小数。"""

    arr = as_float_array(value)
    kind = canonical_rh_unit(unit)
    if kind == "fraction":
        return arr
    return arr / 100.0


def from_rh_fraction(value: ArrayLike, unit: str) -> np.ndarray:
    """把相对湿度小数换算为指定单位。"""

    arr = as_float_array(value)
    kind = canonical_rh_unit(unit)
    if kind == "fraction":
        return arr
    return arr * 100.0


def to_kgkg(value: ArrayLike, unit: str) -> np.ndarray:
    """把混合比或比湿换算为 kg/kg。"""

    arr = as_float_array(value)
    kind = canonical_mass_humidity_unit(unit)
    if kind == "kg/kg":
        return arr
    if kind == "g/kg":
        return arr / 1000.0
    return arr / 1.0e6


def from_kgkg(value: ArrayLike, unit: str) -> np.ndarray:
    """把 kg/kg 换算为指定水汽质量单位。"""

    arr = as_float_array(value)
    kind = canonical_mass_humidity_unit(unit)
    if kind == "kg/kg":
        return arr
    if kind == "g/kg":
        return arr * 1000.0
    return arr * 1.0e6


def to_meters(value: ArrayLike, unit: str) -> np.ndarray:
    """把距离换算为米。"""

    arr = as_float_array(value)
    kind = canonical_distance_unit(unit)
    factors = {"m": 1.0, "km": 1000.0, "cm": 0.01, "mi": _MILE_TO_M, "nmi": _NMI_TO_M, "ft": 0.3048}
    return arr * factors[kind]


def from_meters(value: ArrayLike, unit: str) -> np.ndarray:
    """把米换算为指定距离单位。"""

    arr = as_float_array(value)
    kind = canonical_distance_unit(unit)
    factors = {"m": 1.0, "km": 1000.0, "cm": 0.01, "mi": _MILE_TO_M, "nmi": _NMI_TO_M, "ft": 0.3048}
    return arr / factors[kind]


def to_radians(value: ArrayLike, unit: str = "deg") -> np.ndarray:
    """把角度换算为弧度。"""

    arr = as_float_array(value)
    kind = canonical_angle_unit(unit)
    if kind == "rad":
        return arr
    return np.deg2rad(arr)
