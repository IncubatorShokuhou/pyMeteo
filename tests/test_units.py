"""单位换算往返与别名识别。"""

import numpy as np
import pytest

from pymeteo.units import (
    UnitError,
    from_kelvin,
    from_kgkg,
    from_meters,
    from_mps,
    from_pascal,
    from_rh_fraction,
    to_kelvin,
    to_kgkg,
    to_meters,
    to_mps,
    to_pascal,
    to_radians,
    to_rh_fraction,
)


def test_temperature_round_trip_c_k_f() -> None:
    celsius = 18.0
    kelvin = to_kelvin(celsius, "C")
    assert kelvin == pytest.approx(291.15)
    assert from_kelvin(kelvin, "degC") == pytest.approx(18.0)
    fahrenheit = from_kelvin(kelvin, "F")
    assert fahrenheit == pytest.approx(64.4)
    assert to_kelvin(fahrenheit, "fahrenheit") == pytest.approx(291.15)


def test_pressure_aliases_hpa_mb_pa() -> None:
    assert to_pascal(1013.25, "hPa") == pytest.approx(101325.0)
    assert to_pascal(1013.25, "mb") == pytest.approx(101325.0)
    assert from_pascal(101325.0, "kPa") == pytest.approx(101.325)


def test_speed_knots_and_kmh() -> None:
    mps = to_mps(1.0, "kt")
    assert mps == pytest.approx(1852.0 / 3600.0)
    assert from_mps(mps, "knots") == pytest.approx(1.0)
    assert from_mps(10.0, "km/h") == pytest.approx(36.0)


def test_humidity_and_distance_round_trips() -> None:
    assert to_rh_fraction(46.5, "%") == pytest.approx(0.465)
    assert from_rh_fraction(0.465, "percent") == pytest.approx(46.5)
    assert to_kgkg(6.0, "g/kg") == pytest.approx(0.006)
    assert from_kgkg(0.006, "g/kg") == pytest.approx(6.0)
    assert to_meters(1.0, "km") == pytest.approx(1000.0)
    assert from_meters(1852.0, "nmi") == pytest.approx(1.0)
    assert from_meters(1852.0, "nm") == pytest.approx(1.0)


def test_angle_degrees_to_radians() -> None:
    assert to_radians(180.0, "deg") == pytest.approx(np.pi)


def test_unknown_unit_raises() -> None:
    with pytest.raises(UnitError):
        to_kelvin(1.0, "rankine")


def test_array_conversion() -> None:
    result = to_kelvin(np.array([0.0, 100.0]), "celsius")
    np.testing.assert_allclose(result, [273.15, 373.15])
