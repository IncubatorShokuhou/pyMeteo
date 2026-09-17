"""水汽热力学公式，对照旧实现与 README 示例。"""

import numpy as np
import pytest

from pymeteo import (
    condensation_temperature,
    convert_humidity,
    dewpoint_from_relative_humidity,
    equivalent_potential_temperature,
    lifting_condensation_level,
    mixing_ratio_from_dewpoint,
    mixing_ratio_from_relative_humidity,
    parcel_temperature_at_pressure,
    potential_temperature,
    relative_humidity_from_dewpoint,
    relative_humidity_from_mixing_ratio,
    saturation_mixing_ratio,
    saturation_vapor_pressure,
    specific_humidity_from_relative_humidity,
    vapor_pressure_from_mixing_ratio,
    vapor_pressure_from_relative_humidity,
    virtual_temperature,
    visibility,
    wet_bulb_temperature,
)


def test_saturation_vapor_pressure_known_values() -> None:
    assert saturation_vapor_pressure(0.0, temperature_unit="C") == pytest.approx(6.1078)
    assert saturation_vapor_pressure(20.0, temperature_unit="C") == pytest.approx(
        23.364614454202506
    )
    assert saturation_vapor_pressure(-10.0, temperature_unit="C") == pytest.approx(
        2.8622003456927323
    )
    pa = saturation_vapor_pressure(0.0, temperature_unit="C", output_pressure_unit="Pa")
    assert pa == pytest.approx(610.78)


def test_dewpoint_from_relative_humidity_readme_sample() -> None:
    dewpoint = dewpoint_from_relative_humidity(
        18.0, 46.5, temperature_unit="C", humidity_unit="%", output_temperature_unit="C"
    )
    assert dewpoint == pytest.approx(6.299777153836203, rel=1e-10)
    kelvin = dewpoint_from_relative_humidity(
        18.0 + 273.15,
        46.5,
        temperature_unit="K",
        humidity_unit="%",
        output_temperature_unit="K",
    )
    assert kelvin == pytest.approx(279.4497771538362, rel=1e-10)


def test_relative_humidity_dewpoint_round_trip() -> None:
    temperature = 18.0
    dewpoint = dewpoint_from_relative_humidity(temperature, 46.5)
    rh = relative_humidity_from_dewpoint(temperature, dewpoint)
    assert rh == pytest.approx(46.5, rel=1e-10)


def test_condensation_temperature_matches_legacy() -> None:
    result = condensation_temperature(850.0, 16.6, 0.6, pressure_unit="hPa", temperature_unit="C")
    assert result == pytest.approx(-2.699374080000041, rel=1e-6, abs=1e-5)


def test_relative_humidity_from_mixing_ratio_legacy() -> None:
    rh = relative_humidity_from_mixing_ratio(
        18.0,
        0.006,
        1000.0,
        temperature_unit="C",
        mixing_ratio_unit="kg/kg",
        pressure_unit="hPa",
        output_humidity_unit="%",
    )
    assert rh == pytest.approx(46.42262853199607, rel=1e-8)


def test_mixing_ratio_from_relative_humidity_legacy() -> None:
    mixing = mixing_ratio_from_relative_humidity(
        1000.0,
        18.0,
        46.5,
        pressure_unit="hPa",
        temperature_unit="C",
        humidity_unit="%",
        output_humidity_unit="kg/kg",
    )
    assert mixing == pytest.approx(0.006018464086616119, rel=1e-12)
    grams = mixing_ratio_from_relative_humidity(
        1000.0,
        18.0,
        46.5,
        output_humidity_unit="g/kg",
    )
    assert grams == pytest.approx(6.018464086616119, rel=1e-12)


def test_specific_humidity_from_relative_humidity() -> None:
    q_gkg = specific_humidity_from_relative_humidity(
        1000.0,
        18.0,
        46.5,
        output_humidity_unit="g/kg",
    )
    assert q_gkg == pytest.approx(5.98245887274087, rel=1e-12)


def test_convert_humidity_mixing_ratio_specific_humidity() -> None:
    specific = convert_humidity(
        0.01,
        from_quantity="mixing_ratio",
        to_quantity="specific_humidity",
        humidity_unit="kg/kg",
    )
    assert specific == pytest.approx(0.01 / 1.01)
    mixing = convert_humidity(
        specific,
        from_quantity="specific_humidity",
        to_quantity="mixing_ratio",
    )
    assert mixing == pytest.approx(0.01)
    grams = convert_humidity(
        10.0,
        from_quantity="mixing_ratio",
        to_quantity="specific_humidity",
        humidity_unit="g/kg",
        output_humidity_unit="g/kg",
    )
    assert grams == pytest.approx(9.900990099009901)


def test_convert_humidity_invalid_specific_humidity_is_nan() -> None:
    result = convert_humidity(1.0, from_quantity="specific_humidity", to_quantity="mixing_ratio")
    assert np.isnan(result)


def test_visibility_ruc_is_kilometre_scale() -> None:
    vis_km = visibility(80.0, 18.0, method="RUC", humidity_unit="%", temperature_unit="C")
    assert vis_km == pytest.approx(11.814700512251644, rel=1e-6)
    vis_m = visibility(80.0, 18.0, method="RUC", output_distance_unit="m")
    assert vis_m == pytest.approx(11814.700512251644, rel=1e-6)


def test_visibility_fsl_legacy_km() -> None:
    vis = visibility(80.0, 18.0, method="FSL")
    assert vis == pytest.approx(9.840102017412084, rel=1e-6)


def test_zero_relative_humidity_dewpoint_is_nan() -> None:
    result = dewpoint_from_relative_humidity(20.0, 0.0)
    assert np.isnan(result)


def test_thermo_arrays() -> None:
    e = saturation_vapor_pressure(np.array([0.0, 20.0]))
    np.testing.assert_allclose(e, [6.1078, 23.364614454202506])


def test_saturation_mixing_ratio_equals_ws_from_es() -> None:
    pressure = 1000.0
    temperature = 20.0
    es = saturation_vapor_pressure(temperature, temperature_unit="C")
    expected = 0.622 * es / (pressure - es)
    ws = saturation_mixing_ratio(pressure, temperature, output_humidity_unit="kg/kg")
    assert ws == pytest.approx(expected)
    grams = saturation_mixing_ratio(pressure, temperature, output_humidity_unit="g/kg")
    assert grams == pytest.approx(expected * 1000.0)


def test_mixing_ratio_from_dewpoint_wallace_hobbs_order() -> None:
    # Wallace & Hobbs：1000 hPa、Td=6.4 °C 时混合比约 6 g/kg
    mixing = mixing_ratio_from_dewpoint(
        1000.0, 6.4, pressure_unit="hPa", temperature_unit="C", output_humidity_unit="g/kg"
    )
    assert mixing == pytest.approx(6.0, abs=0.15)


def test_vapor_pressure_mixing_ratio_round_trip() -> None:
    pressure = 1000.0
    mixing = mixing_ratio_from_dewpoint(pressure, 10.0)
    vapor = vapor_pressure_from_mixing_ratio(pressure, mixing)
    es = saturation_vapor_pressure(10.0)
    assert vapor == pytest.approx(es, rel=1e-10)


def test_vapor_pressure_from_relative_humidity_half_saturation() -> None:
    es = saturation_vapor_pressure(15.0)
    vapor = vapor_pressure_from_relative_humidity(15.0, 50.0, humidity_unit="%")
    assert vapor == pytest.approx(0.5 * es)


def test_potential_temperature_at_1000_hpa_is_temperature() -> None:
    theta = potential_temperature(1000.0, 20.0, temperature_unit="C", output_temperature_unit="C")
    assert theta == pytest.approx(20.0)
    theta_k = potential_temperature(
        100000.0, 293.15, pressure_unit="Pa", temperature_unit="K", output_temperature_unit="K"
    )
    assert theta_k == pytest.approx(293.15)


def test_potential_temperature_poisson_850_hpa() -> None:
    theta = potential_temperature(850.0, 10.0, temperature_unit="C", output_temperature_unit="K")
    expected = (10.0 + 273.15) * (1000.0 / 850.0) ** 0.286
    assert theta == pytest.approx(expected)


def test_potential_temperature_ncl_documented_point() -> None:
    # NCL pot_temp 文档：p=100800 Pa、T=302.45 K → θ≈301.762 K
    theta = potential_temperature(
        100800.0, 302.45, pressure_unit="Pa", temperature_unit="K", output_temperature_unit="K"
    )
    assert theta == pytest.approx(301.762, abs=1e-3)


def test_equivalent_potential_temperature_dry_near_potential() -> None:
    theta = potential_temperature(1000.0, 20.0, output_temperature_unit="K")
    theta_e = equivalent_potential_temperature(1000.0, 20.0, -40.0, output_temperature_unit="K")
    assert theta_e == pytest.approx(theta, abs=0.5)


def test_equivalent_potential_temperature_saturated_exceeds_theta() -> None:
    theta = potential_temperature(1000.0, 20.0, output_temperature_unit="K")
    theta_e = equivalent_potential_temperature(1000.0, 20.0, 20.0, output_temperature_unit="K")
    assert theta_e > theta + 10.0


def test_equivalent_potential_temperature_bolton_formula() -> None:
    pressure = 1000.0
    temperature_k = 27.0 + 273.15
    dewpoint_k = 22.0 + 273.15
    es = saturation_vapor_pressure(22.0)
    mixing = 0.622 * es / (pressure - es)
    t_lcl = 1.0 / (1.0 / (dewpoint_k - 56.0) + np.log(temperature_k / dewpoint_k) / 800.0) + 56.0
    expected = temperature_k * np.exp((3376.0 / t_lcl - 2.54) * mixing * (1.0 + 0.81 * mixing))
    theta_e = equivalent_potential_temperature(
        pressure, 27.0, 22.0, temperature_unit="C", output_temperature_unit="K"
    )
    assert theta_e == pytest.approx(expected, rel=1e-10)


def test_virtual_temperature_exact_formula() -> None:
    mixing = 0.0135
    temperature_k = 20.0 + 273.15
    expected_c = temperature_k * (1.0 + mixing / 0.622) / (1.0 + mixing) - 273.15
    tv = virtual_temperature(20.0, mixing, mixing_ratio_unit="kg/kg")
    assert tv == pytest.approx(expected_c)
    tv_g = virtual_temperature(20.0, 13.5, mixing_ratio_unit="g/kg")
    assert tv_g == pytest.approx(expected_c)


def test_wet_bulb_stull_documented_example() -> None:
    # Stull (2011) 与 NCL wetbulb_stull：20 °C、50% → 约 13.7 °C
    tw = wet_bulb_temperature(20.0, 50.0)
    assert tw == pytest.approx(13.699341968988136, rel=1e-8)
    tw_k = wet_bulb_temperature(
        20.0 + 273.15, 0.5, temperature_unit="K", humidity_unit="fraction", output_temperature_unit="K"
    )
    assert tw_k == pytest.approx(13.699341968988136 + 273.15, rel=1e-8)


def test_lifting_condensation_level_saturated_is_starting_point() -> None:
    p_lcl, t_lcl = lifting_condensation_level(950.0, 12.0, 12.0)
    assert p_lcl == pytest.approx(950.0, rel=1e-6)
    assert t_lcl == pytest.approx(12.0, abs=1e-4)


def test_lifting_condensation_level_wallace_hobbs() -> None:
    # Wallace & Hobbs / NCL lclvl：1000 hPa、15 °C、Td=4 °C → 约 848 hPa
    p_lcl, t_lcl = lifting_condensation_level(1000.0, 15.0, 4.0)
    assert p_lcl == pytest.approx(848.0, abs=3.0)
    # 干绝热抬升中露点随气压下降，T_LCL 低于起始 T 与 Td
    assert t_lcl < 15.0
    assert t_lcl < 4.0


def test_parcel_temperature_dry_adiabatic_when_lcl_above_target() -> None:
    # 极干：LCL 远高于 500 hPa，500 hPa 气块温度应为干绝热
    parcel = parcel_temperature_at_pressure(
        1000.0, 20.0, -40.0, 500.0, output_temperature_unit="K"
    )
    expected = (20.0 + 273.15) * (500.0 / 1000.0) ** 0.286
    assert parcel == pytest.approx(expected, rel=1e-8)
