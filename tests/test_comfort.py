"""热指数与风寒：对照 NWS 公布公式。"""

import numpy as np
import pytest

from pymeteo import heat_index, wind_chill


def test_heat_index_nws_90f_60_percent() -> None:
    # Rothfusz 全式在 90 °F、60% 时约 99.7 °F（常被四舍五入为 100 °F）
    hi = heat_index(90.0, 60.0, temperature_unit="F", output_temperature_unit="F")
    expected = (
        -42.379
        + 2.04901523 * 90.0
        + 10.14333127 * 60.0
        - 0.22475541 * 90.0 * 60.0
        - 0.00683783 * 90.0**2
        - 0.05481717 * 60.0**2
        + 0.00122874 * 90.0**2 * 60.0
        + 0.00085282 * 90.0 * 60.0**2
        - 0.00000199 * 90.0**2 * 60.0**2
    )
    assert hi == pytest.approx(expected)
    assert hi == pytest.approx(100.0, abs=0.4)


def test_heat_index_below_80f_uses_simple_average() -> None:
    t_f = 70.0
    rh = 40.0
    simple = 0.5 * (t_f + 61.0 + ((t_f - 68.0) * 1.2) + (rh * 0.094))
    expected = 0.5 * (simple + t_f)
    hi = heat_index(t_f, rh, temperature_unit="F", output_temperature_unit="F")
    assert hi == pytest.approx(expected)


def test_heat_index_celsius_default() -> None:
    hi_f = heat_index(90.0, 60.0, temperature_unit="F", output_temperature_unit="F")
    hi_c = heat_index(
        (90.0 - 32.0) * 5.0 / 9.0, 60.0, temperature_unit="C", output_temperature_unit="C"
    )
    assert hi_c == pytest.approx((hi_f - 32.0) * 5.0 / 9.0, rel=1e-6)


def test_heat_index_low_humidity_adjustment() -> None:
    t_f = 90.0
    rh = 10.0
    roth = (
        -42.379
        + 2.04901523 * t_f
        + 10.14333127 * rh
        - 0.22475541 * t_f * rh
        - 0.00683783 * t_f**2
        - 0.05481717 * rh**2
        + 0.00122874 * t_f**2 * rh
        + 0.00085282 * t_f * rh**2
        - 0.00000199 * t_f**2 * rh**2
    )
    adj = ((13.0 - rh) / 4.0) * np.sqrt((17.0 - abs(t_f - 95.0)) / 17.0)
    expected = roth - adj
    hi = heat_index(t_f, rh, temperature_unit="F", output_temperature_unit="F")
    assert hi == pytest.approx(expected)


def test_wind_chill_0f_10mph() -> None:
    wc = wind_chill(0.0, 10.0, temperature_unit="F", speed_unit="mph", output_temperature_unit="F")
    expected = 35.74 + 0.6215 * 0.0 - 35.75 * (10.0**0.16) + 0.4275 * 0.0 * (10.0**0.16)
    assert wc == pytest.approx(expected)
    assert wc == pytest.approx(-16.0, abs=0.1)


def test_wind_chill_celsius_and_array() -> None:
    wc = wind_chill(
        np.array([0.0, 10.0]),
        np.array([5.0, 5.0]),
        temperature_unit="C",
        speed_unit="m/s",
    )
    assert wc.shape == (2,)
    assert wc[0] < 0.0
