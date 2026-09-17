"""稳定度指数：对照旧公式与文献订正后的 SWEAT。"""

import numpy as np
import pytest

from pymeteo import (
    a_index,
    k_index,
    layer_temperature_difference,
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
