"""风场恒等式与气象风向惯例。"""

import numpy as np
import pytest

from pymeteo import (
    bulk_wind_shear,
    uv_from_speed_direction,
    wind_components,
    wind_direction,
    wind_speed,
)


def test_wind_speed_from_components() -> None:
    assert wind_speed(3.0, 4.0) == pytest.approx(5.0)


def test_wind_direction_legacy() -> None:
    assert wind_direction(3.0, 4.0) == pytest.approx(216.86989764584402)
    assert wind_direction(0.0, 0.0) == pytest.approx(0.0)


def test_uv_from_meteorological_direction() -> None:
    u, v = uv_from_speed_direction(5.0, 90.0)
    assert u == pytest.approx(-5.0)
    assert v == pytest.approx(0.0, abs=1e-12)
    u0, v0 = uv_from_speed_direction(5.0, 0.0)
    assert u0 == pytest.approx(0.0, abs=1e-12)
    assert v0 == pytest.approx(-5.0)
    u180, v180 = uv_from_speed_direction(5.0, 180.0)
    assert u180 == pytest.approx(0.0, abs=1e-12)
    assert v180 == pytest.approx(5.0)


def test_round_trip_speed_direction() -> None:
    speed = 12.0
    direction = 247.5
    u, v = wind_components(speed, direction)
    assert wind_speed(u, v) == pytest.approx(speed)
    assert wind_direction(u, v) == pytest.approx(direction)


def test_knot_unit_conversion_on_speed() -> None:
    speed_kt = wind_speed(10.0, 0.0, speed_unit="m/s", output_speed_unit="kt")
    assert speed_kt == pytest.approx(10.0 * 3600.0 / 1852.0)


def test_wind_arrays() -> None:
    speed = wind_speed(np.array([3.0, 0.0]), np.array([4.0, 5.0]))
    np.testing.assert_allclose(speed, [5.0, 5.0])


def test_bulk_wind_shear_vector_difference() -> None:
    shear = bulk_wind_shear(0.0, 0.0, 3.0, 4.0)
    assert shear == pytest.approx(5.0)
    shear_kt = bulk_wind_shear(10.0, 0.0, 0.0, 0.0, speed_unit="m/s", output_speed_unit="kt")
    assert shear_kt == pytest.approx(10.0 * 3600.0 / 1852.0)
    # 同向同速 → 0
    assert bulk_wind_shear(5.0, 5.0, 5.0, 5.0) == pytest.approx(0.0)
