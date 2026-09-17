"""水汽热力学公式，对照旧实现与 README 示例。"""

import numpy as np
import pytest

from pymeteo import (
    condensation_temperature,
    convert_humidity,
    dewpoint_from_relative_humidity,
    mixing_ratio_from_relative_humidity,
    relative_humidity_from_dewpoint,
    relative_humidity_from_mixing_ratio,
    saturation_vapor_pressure,
    specific_humidity_from_relative_humidity,
    visibility,
)


def test_saturation_vapor_pressure_known_values() -> None:
    assert saturation_vapor_pressure(0.0, temperature_unit="C") == pytest.approx(6.1078)
    assert saturation_vapor_pressure(20.0, temperature_unit="C") == pytest.approx(23.364614454202506)
    assert saturation_vapor_pressure(-10.0, temperature_unit="C") == pytest.approx(2.8622003456927323)
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


def test_visibility_ruc_is_kilometre_scale() -> None:
    vis_km = visibility(80.0, 18.0, method="RUC", humidity_unit="%", temperature_unit="C")
    assert vis_km == pytest.approx(11.814700512251644, rel=1e-6)
    vis_m = visibility(
        80.0, 18.0, method="RUC", output_distance_unit="m"
    )
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
