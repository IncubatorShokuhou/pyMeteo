"""对照 NCL 官网文档例题打印的数值做回归。

黄金值硬编码自公开页面（及页面上引用的 Wallace & Hobbs / Stull 印刷值），
运行时不依赖 MetPy、NCL 库或 Pint。``tests/test_ncl.py`` 只测封装接线；
本文件专门钉住官网例题，公式若被“改进”偏离文档数字应失败。
"""

import numpy as np
import pytest

import pymeteo as pm
from pymeteo import (
    dewpoint_from_relative_humidity,
    lifting_condensation_level,
    mixing_ratio_from_relative_humidity,
    potential_temperature,
    relative_humidity_from_dewpoint,
    specific_humidity_from_relative_humidity,
    wet_bulb_temperature,
)

# NCL dewtemp_trh / mixhum_ptrh Example 1 共用 Wallace & Hobbs p.74 输入
_TK_18C = 18.0 + 273.15
_RH_WH = 46.5
_P_HPA_WH = 1000.0


def test_dewtemp_trh_ncl_example_1() -> None:
    # NCL dewtemp_trh Example 1
    # https://www.ncl.ucar.edu/Document/Functions/Built-in/dewtemp_trh.shtml
    # 输入：tk = 18+273.15 K，rh = 46.5%
    # NCL 注释：td = 6.3 C；Wallace & Hobbs 书中约 6.4 °C
    ncl_k = pm.ncl.dewtemp_trh(_TK_18C, _RH_WH)
    ncl_c = ncl_k - 273.15
    modern_c = dewpoint_from_relative_humidity(
        18.0, _RH_WH, temperature_unit="C", humidity_unit="%", output_temperature_unit="C"
    )
    modern_k = dewpoint_from_relative_humidity(
        _TK_18C,
        _RH_WH,
        temperature_unit="K",
        humidity_unit="%",
        output_temperature_unit="K",
    )
    assert ncl_c == pytest.approx(modern_c)
    assert ncl_k == pytest.approx(modern_k)
    # 本库已给出约 6.299 °C，须钉在 NCL 印刷的 6.3 °C 附近（0.05 °C）
    assert ncl_c == pytest.approx(6.3, abs=0.05)
    # 与书中约 6.4 °C 允许约 0.15 °C（NCL 自身也写成 6.3 而非 6.4）
    assert ncl_c == pytest.approx(6.4, abs=0.15)


def test_mixhum_ptrh_ncl_example_1() -> None:
    # NCL mixhum_ptrh Example 1
    # https://www.ncl.ucar.edu/Document/Functions/Built-in/mixhum_ptrh.shtml
    # 输入：p = 1000 hPa，tk = 18+273.15 K，rh = 46.5%
    # iswit=1  → 0.006018462 kg/kg（混合比）
    # iswit=-1 → 6.018462 g/kg
    # iswit=2  → 0.005982457 kg/kg（比湿）
    # iswit=-2 → 5.982456 g/kg
    # NCL 最后一位受打印舍入影响（iswit=-2 印刷 5.982456 g/kg，
    # 而 iswit=2 的 0.005982457 kg/kg ×1000 = 5.982457 g/kg，差 1e-6 g/kg），
    # 相对 1e-6 仍覆盖；本库 Tetens 实现已达约 3e-7。
    mix_kg = pm.ncl.mixhum_ptrh(_P_HPA_WH, _TK_18C, _RH_WH, 1)
    mix_g = pm.ncl.mixhum_ptrh(_P_HPA_WH, _TK_18C, _RH_WH, -1)
    q_kg = pm.ncl.mixhum_ptrh(_P_HPA_WH, _TK_18C, _RH_WH, 2)
    q_g = pm.ncl.mixhum_ptrh(_P_HPA_WH, _TK_18C, _RH_WH, -2)

    assert mix_kg == pytest.approx(0.006018462, rel=1e-6)
    assert mix_g == pytest.approx(6.018462, rel=1e-6)
    assert q_kg == pytest.approx(0.005982457, rel=1e-6)
    assert q_g == pytest.approx(5.982456, rel=1e-6)

    kwargs = dict(pressure_unit="hPa", temperature_unit="K", humidity_unit="%")
    assert mix_kg == pytest.approx(
        mixing_ratio_from_relative_humidity(
            _P_HPA_WH, _TK_18C, _RH_WH, output_humidity_unit="kg/kg", **kwargs
        )
    )
    assert q_kg == pytest.approx(
        specific_humidity_from_relative_humidity(
            _P_HPA_WH, _TK_18C, _RH_WH, output_humidity_unit="kg/kg", **kwargs
        )
    )


# NCL pot_temp Example 1 表：p(Pa), t(K), pot(K)
# https://www.ncl.ucar.edu/Document/Functions/Contributed/pot_temp.shtml
# NCL 印刷三位小数；Poisson κ=0.286，abs=1e-3 K 对齐打印精度。
_POT_TEMP_NCL_ROWS = [
    (100800.0, 302.45, 301.762),
    (100000.0, 301.25, 301.25),
    (95000.0, 296.65, 301.034),
    (90000.0, 294.05, 303.045),
    (85000.0, 291.55, 305.421),
    (80000.0, 289.05, 308.098),
    (50000.0, 268.65, 327.553),
    (20000.0, 220.75, 349.789),
]


@pytest.mark.parametrize("pressure_pa, temperature_k, theta_ncl", _POT_TEMP_NCL_ROWS)
def test_pot_temp_ncl_example_1_table(
    pressure_pa: float, temperature_k: float, theta_ncl: float
) -> None:
    ncl = pm.ncl.pot_temp(pressure_pa, temperature_k)
    modern = potential_temperature(
        pressure_pa,
        temperature_k,
        pressure_unit="Pa",
        temperature_unit="K",
        output_temperature_unit="K",
    )
    assert ncl == pytest.approx(modern)
    assert ncl == pytest.approx(theta_ncl, abs=1e-3)


def test_wetbulb_stull_ncl_example_1() -> None:
    # NCL wetbulb_stull Example 1 / Stull (2011)
    # https://www.ncl.ucar.edu/Document/Functions/Contributed/wetbulb_stull.shtml
    # 输入：T = 20 °C，RH = 50%；NCL 打印 TW00 = 13.69934 °C
    ncl_c = pm.ncl.wetbulb_stull(20.0, 50.0, (0, 0), False)
    modern_c = wet_bulb_temperature(
        20.0, 50.0, temperature_unit="C", humidity_unit="%", output_temperature_unit="C"
    )
    ncl_from_k = pm.ncl.wetbulb_stull(20.0 + 273.15, 50.0, (1, 0))
    assert ncl_c == pytest.approx(modern_c)
    assert ncl_c == pytest.approx(13.69934, rel=1e-5)
    assert ncl_c == pytest.approx(13.69934, abs=1e-3)
    assert ncl_from_k == pytest.approx(13.69934, rel=1e-5)


def test_lclvl_ncl_example_1_wallace_hobbs() -> None:
    # NCL lclvl Example 1
    # https://www.ncl.ucar.edu/Document/Functions/Built-in/lclvl.shtml
    # 输入：p = 1000 hPa，t = 15 °C，td = 4 °C
    # NCL 注释：plcl = 848.6 hPa；Wallace & Hobbs 习题约 848 hPa
    # NCL 用 Stipanuk (1973) 迭代；本库 lclvl 转调 Bolton (1980) 式 (22)
    # 再沿干绝热求气压，当前约 847.12 hPa，与 848.6 差约 1.5 hPa，
    # 属经验公式差异而非实现错误，故对官网 848.6 用 abs=2.0 hPa。
    tk = 15.0 + 273.15
    tdk = 4.0 + 273.15
    ncl_p = pm.ncl.lclvl(1000.0, tk, tdk)
    modern_p, _t_lcl = lifting_condensation_level(
        1000.0, tk, tdk, pressure_unit="hPa", temperature_unit="K", output_pressure_unit="hPa"
    )
    assert ncl_p == pytest.approx(modern_p)
    assert ncl_p == pytest.approx(848.6, abs=2.0)
    assert ncl_p == pytest.approx(848.0, abs=1.0)


def test_ncl_example_2_sounding_spot_check() -> None:
    # NCL dewtemp_trh Example 2 / mixhum_ptrh Example 2 探空前几层
    # https://www.ncl.ucar.edu/Document/Functions/Built-in/dewtemp_trh.shtml
    # https://www.ncl.ucar.edu/Document/Functions/Built-in/mixhum_ptrh.shtml
    # 官网未给出整表黄金值；此处抽 3 层做有限性、RH 往返与混合比量级检查。
    pressure_hpa = np.array([1008.0, 1000.0, 850.0])
    temperature_c = np.array([29.3, 28.1, 18.4])
    relative_humidity = np.array([75.0, 60.0, 90.5])
    temperature_k = temperature_c + 273.15

    dewpoint_k = pm.ncl.dewtemp_trh(temperature_k, relative_humidity)
    mixing_gkg = pm.ncl.mixhum_ptrh(pressure_hpa, temperature_k, relative_humidity, -1)
    rh_round_trip = relative_humidity_from_dewpoint(
        temperature_c, dewpoint_k - 273.15, temperature_unit="C", output_humidity_unit="%"
    )

    assert np.all(np.isfinite(dewpoint_k))
    assert np.all(np.isfinite(mixing_gkg))
    assert np.all(dewpoint_k < temperature_k)
    # 对流层常见水汽混合比量级（g/kg）；饱和暖湿近地面可到十几 g/kg
    assert np.all(mixing_gkg > 1.0)
    assert np.all(mixing_gkg < 40.0)
    np.testing.assert_allclose(rh_round_trip, relative_humidity, rtol=1e-10)
