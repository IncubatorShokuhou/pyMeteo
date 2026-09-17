"""稳定度指数：对照旧公式与文献订正后的 SWEAT。"""

import numpy as np
import pytest

from pymeteo import (
    a_index,
    k_index,
    layer_temperature_difference,
    lifted_index,
    lifted_index_from_surface,
    parcel_temperature_at_pressure,
    showalter_index,
    sweat_index,
    temperature_dewpoint_depression,
    total_totals_index,
)


def test_k_index_legacy_sample() -> None:
    result = k_index(16.6, 0.6, 7.0, -2.0, -15.9, temperature_unit="C")
    assert result == pytest.approx(24.1)


def test_a_index_legacy_sample() -> None:
    result = a_index(16.6, 0.6, 7.0, -2.0, -15.9, -20.0, temperature_unit="C")
    assert result == pytest.approx(3.4)


def test_layer_diffs() -> None:
    assert temperature_dewpoint_depression(7.0, -2.0) == pytest.approx(9.0)
    assert layer_temperature_difference(16.6, -15.9) == pytest.approx(32.5)


def test_total_totals_from_relative_humidity() -> None:
    tt = total_totals_index(
        18.0,
        -15.9,
        relative_humidity_850=46.5,
        temperature_unit="C",
        humidity_unit="%",
    )
    assert tt == pytest.approx(56.09977715383616, rel=1e-10)


def test_total_totals_from_dewpoint_matches_kelvin_inputs() -> None:
    dewpoint = 6.299777153836203
    from_c = total_totals_index(18.0, -15.9, dewpoint_850=dewpoint, temperature_unit="C")
    from_k = total_totals_index(
        18.0 + 273.15,
        -15.9 + 273.15,
        dewpoint_850=dewpoint + 273.15,
        temperature_unit="K",
    )
    assert from_c == pytest.approx(from_k)


def test_showalter_legacy_sample() -> None:
    result = showalter_index(16.6, 0.6, -15.9, temperature_unit="C")
    assert result == pytest.approx(1.1, abs=1e-3)


def test_showalter_array() -> None:
    result = showalter_index(
        np.array([16.6, 16.6]),
        np.array([0.6, 0.6]),
        np.array([-15.9, -15.9]),
    )
    np.testing.assert_allclose(result, [1.1, 1.1], atol=1e-3)


def test_sweat_uses_celsius_dewpoint_and_knots() -> None:
    result = sweat_index(
        18.0,
        -15.9,
        10.0,
        5.0,
        20.0,
        10.0,
        relative_humidity_850=46.5,
        temperature_unit="C",
        humidity_unit="%",
        speed_unit="m/s",
    )
    # 12*6.299777 + 20*(56.099777-49) + 2*21.573 + 43.146，切变项因风向差为 0 而置零
    assert result == pytest.approx(304.524, abs=0.05)
    # 旧实现用开尔文露点会得到约 3609，此处必须远小于该错误量级
    assert result < 500.0


def test_sweat_dewpoint_argument_matches_relative_humidity() -> None:
    dewpoint = 6.299777153836203
    from_td = sweat_index(
        18.0,
        -15.9,
        10.0,
        5.0,
        20.0,
        10.0,
        dewpoint_850=dewpoint,
        temperature_unit="C",
        speed_unit="m/s",
    )
    from_rh = sweat_index(
        18.0,
        -15.9,
        10.0,
        5.0,
        20.0,
        10.0,
        relative_humidity_850=46.5,
        temperature_unit="C",
        humidity_unit="%",
        speed_unit="m/s",
    )
    assert from_td == pytest.approx(from_rh, rel=1e-10)


def test_sweat_shear_term_when_veering_and_strong_winds() -> None:
    # 850 hPa 风向 180°、500 hPa 风向 250°，风速均 ≥ 15 kt，切变项应保留
    result = sweat_index(
        18.0,
        -15.9,
        0.0,
        20.0,
        30.0 * np.sin(np.deg2rad(70.0)),
        30.0 * np.cos(np.deg2rad(70.0)),
        dewpoint_850=6.299777153836203,
        temperature_unit="C",
        speed_unit="kt",
    )
    tt = 18.0 + 6.299777153836203 - 2.0 * (-15.9)
    expected = (
        12.0 * 6.299777153836203
        + 20.0 * (tt - 49.0)
        + 2.0 * 20.0
        + 30.0
        + 125.0 * (np.sin(np.deg2rad(70.0)) + 0.2)
    )
    assert result == pytest.approx(expected, rel=1e-6)


def test_k_index_independent_of_c_or_k() -> None:
    in_c = k_index(16.6, 0.6, 7.0, -2.0, -15.9, temperature_unit="C")
    in_k = k_index(
        16.6 + 273.15,
        0.6 + 273.15,
        7.0 + 273.15,
        -2.0 + 273.15,
        -15.9 + 273.15,
        temperature_unit="K",
    )
    assert in_c == pytest.approx(in_k)


def test_a_index_independent_of_c_or_k() -> None:
    in_c = a_index(16.6, 0.6, 7.0, -2.0, -15.9, -20.0, temperature_unit="C")
    in_k = a_index(
        16.6 + 273.15,
        0.6 + 273.15,
        7.0 + 273.15,
        -2.0 + 273.15,
        -15.9 + 273.15,
        -20.0 + 273.15,
        temperature_unit="K",
    )
    assert in_c == pytest.approx(in_k)


def test_lifted_index_is_environment_minus_parcel() -> None:
    assert lifted_index(-15.0, -10.0) == pytest.approx(-5.0)
    assert lifted_index(-15.0, -20.0) == pytest.approx(5.0)
    in_k = lifted_index(-15.0 + 273.15, -10.0 + 273.15, temperature_unit="K")
    assert in_k == pytest.approx(-5.0)


def test_lifted_index_from_surface_dry_matches_poisson() -> None:
    t500 = -20.0
    li = lifted_index_from_surface(1000.0, 20.0, -40.0, t500)
    parcel_500 = (20.0 + 273.15) * (500.0 / 1000.0) ** 0.286 - 273.15
    assert li == pytest.approx(t500 - parcel_500, rel=1e-8)


def test_lifted_index_from_surface_moist_warmer_than_dry() -> None:
    dry = lifted_index_from_surface(1000.0, 20.0, -40.0, -20.0)
    moist = lifted_index_from_surface(1000.0, 20.0, 18.0, -20.0)
    # 湿气块到达 500 hPa 更暖，LI 更小（更不稳定）
    assert moist < dry


def test_lifted_index_from_surface_matches_helper() -> None:
    parcel = parcel_temperature_at_pressure(1000.0, 25.0, 15.0, 500.0)
    li = lifted_index_from_surface(1000.0, 25.0, 15.0, -10.0)
    assert li == pytest.approx(lifted_index(-10.0, parcel))


def test_lifted_index_from_surface_pressure_unit_does_not_shift_500hpa() -> None:
    li_hpa = lifted_index_from_surface(1000.0, 20.0, -40.0, -20.0, pressure_unit="hPa")
    li_pa = lifted_index_from_surface(100000.0, 20.0, -40.0, -20.0, pressure_unit="Pa")
    assert li_hpa == pytest.approx(li_pa)
