"""科里奥利参数与 ω ↔ w 往返。"""

import math

import numpy as np
import pytest

from pymeteo import coriolis_parameter, omega_to_w, w_to_omega


def test_coriolis_parameter_equator_and_poles() -> None:
    omega = 7.292e-5
    assert coriolis_parameter(0.0) == pytest.approx(0.0)
    assert coriolis_parameter(90.0) == pytest.approx(2.0 * omega)
    assert coriolis_parameter(-90.0) == pytest.approx(-2.0 * omega)


def test_coriolis_parameter_45_and_ncl_35() -> None:
    omega = 7.292e-5
    expected_45 = 2.0 * omega * math.sin(math.radians(45.0))
    assert coriolis_parameter(45.0, latitude_unit="deg") == pytest.approx(expected_45)
    assert coriolis_parameter(math.radians(45.0), latitude_unit="rad") == pytest.approx(expected_45)
    # NCL 文档：35° → 8.365038e-5 s⁻¹
    assert coriolis_parameter(35.0) == pytest.approx(8.365038e-5, rel=1e-6)


def test_omega_w_round_trip() -> None:
    omega = 0.1
    temperature = 0.0
    pressure = 850.0
    w = omega_to_w(omega, temperature, pressure)
    back = w_to_omega(w, temperature, pressure)
    assert back == pytest.approx(omega)
    # 下沉 ω>0 → w<0
    assert w < 0.0


def test_omega_to_w_matches_hydrostatic_ideal_gas() -> None:
    rd = 287.058
    g = 9.80665
    omega = 0.1
    t_k = 273.15
    p_pa = 85000.0
    expected = -omega / ((p_pa / (rd * t_k)) * g)
    w = omega_to_w(
        omega,
        t_k,
        p_pa,
        temperature_unit="K",
        pressure_unit="Pa",
        output_speed_unit="m/s",
    )
    assert w == pytest.approx(expected)


def test_omega_hpa_unit_and_arrays() -> None:
    w_pa = omega_to_w(0.1, 0.0, 1000.0, omega_unit="Pa/s")
    w_hpa = omega_to_w(0.001, 0.0, 1000.0, omega_unit="hPa/s")
    assert w_hpa == pytest.approx(w_pa)
    result = omega_to_w(np.array([0.1, -0.1]), np.array([0.0, 10.0]), np.array([1000.0, 850.0]))
    assert result.shape == (2,)
    assert result[0] < 0.0
    assert result[1] > 0.0
