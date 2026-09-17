"""几何距离、重力与海平面气压。"""

import math

import numpy as np
import pytest

from pymeteo import earth_distance, gravity, height_thickness, sea_level_pressure


def test_earth_distance_equator_one_degree() -> None:
    distance_m = earth_distance(0.0, 0.0, 0.0, 1.0, output_distance_unit="m")
    assert distance_m == pytest.approx(111319.491, rel=1e-5)
    distance_km = earth_distance(0.0, 0.0, 0.0, 1.0, output_distance_unit="km")
    assert distance_km == pytest.approx(distance_m / 1000.0)


def test_earth_distance_meridian_one_degree() -> None:
    distance_m = earth_distance(0.0, 0.0, 1.0, 0.0, output_distance_unit="m")
    assert distance_m == pytest.approx(110574.389, rel=1e-5)


def test_earth_distance_identical_points_zero() -> None:
    assert earth_distance(40.0, 116.0, 40.0, 116.0) == pytest.approx(0.0)


def test_earth_distance_array() -> None:
    result = earth_distance(
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 1.0]),
        np.array([1.0, 0.0]),
        output_distance_unit="km",
    )
    assert result.shape == (2,)
    assert result[0] == pytest.approx(111.319, rel=1e-4)


def test_gravity_latitude_in_degrees() -> None:
    assert gravity(0.0) == pytest.approx(9.7803)
    expected_45 = 9.7803 * (
        1.0
        + 0.0053024 * math.sin(math.radians(45.0)) ** 2
        - 0.000005 * math.sin(math.radians(90.0)) ** 2
    )
    assert gravity(45.0, latitude_unit="deg") == pytest.approx(expected_45)
    assert gravity(math.radians(45.0), latitude_unit="rad") == pytest.approx(expected_45)


def test_sea_level_pressure_legacy() -> None:
    slp = sea_level_pressure(1000.0, 100.0, 20.0, 18.0, pressure_unit="hPa", temperature_unit="C")
    assert slp == pytest.approx(1011.7583630997631)


def test_height_thickness_hypsometric_1000_to_500() -> None:
    rd = 287.058
    g = 9.80665
    t_k = 273.15
    expected_m = rd * t_k / g * math.log(1000.0 / 500.0)
    thickness = height_thickness(1000.0, 500.0, 0.0, temperature_unit="C", output_distance_unit="m")
    assert thickness == pytest.approx(expected_m)
    thickness_km = height_thickness(
        1000.0, 500.0, 273.15, temperature_unit="K", output_distance_unit="km"
    )
    assert thickness_km == pytest.approx(expected_m / 1000.0)
