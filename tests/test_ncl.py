"""NCL 兼容薄封装：单位约定应与调用现代 API 等价。"""

import numpy as np
import pytest

import pymeteo as pm
from pymeteo import (
    convert_humidity,
    dewpoint_from_relative_humidity,
    mixing_ratio_from_relative_humidity,
    relative_humidity_from_dewpoint,
    relative_humidity_from_mixing_ratio,
    specific_humidity_from_relative_humidity,
    uv_from_speed_direction,
    wind_speed,
)


def test_ncl_import_paths() -> None:
    from pymeteo.ncl import dewtemp_trh

    tk = 18.0 + 273.15
    rh = 46.5
    via_from = dewtemp_trh(tk, rh)
    via_attr = pm.ncl.dewtemp_trh(tk, rh)
    assert via_from == pytest.approx(via_attr)
    assert not hasattr(pm, "dewtemp_trh")
    assert "dewtemp_trh" not in pm.__all__
    assert "ncl" in pm.__all__


def test_dewtemp_trh_matches_modern_kelvin_percent() -> None:
    tk = 18.0 + 273.15
    rh = 46.5
    ncl = pm.ncl.dewtemp_trh(tk, rh)
    modern = dewpoint_from_relative_humidity(
        tk,
        rh,
        temperature_unit="K",
        humidity_unit="%",
        output_temperature_unit="K",
    )
    assert ncl == pytest.approx(modern)
    assert ncl == pytest.approx(279.4497771538362, rel=1e-10)


def test_relhum_ttd_opt_percent_and_fraction() -> None:
    t = 18.0 + 273.15
    td = 6.299777153836203 + 273.15
    percent = pm.ncl.relhum_ttd(t, td, 0)
    fraction = pm.ncl.relhum_ttd(t, td, 1)
    modern_percent = relative_humidity_from_dewpoint(
        t, td, temperature_unit="K", output_humidity_unit="%"
    )
    modern_fraction = relative_humidity_from_dewpoint(
        t, td, temperature_unit="K", output_humidity_unit="fraction"
    )
    assert percent == pytest.approx(modern_percent)
    assert fraction == pytest.approx(modern_fraction)
    assert percent == pytest.approx(46.5, rel=1e-10)
    assert fraction == pytest.approx(0.465, rel=1e-10)


def test_relhum_ncl_units_kelvin_kgkg_pascal() -> None:
    t = 18.0 + 273.15
    w = 0.006
    p = 1000.0 * 100.0
    ncl = pm.ncl.relhum(t, w, p)
    modern = relative_humidity_from_mixing_ratio(
        t,
        w,
        p,
        temperature_unit="K",
        mixing_ratio_unit="kg/kg",
        pressure_unit="Pa",
        output_humidity_unit="%",
    )
    assert ncl == pytest.approx(modern)
    assert ncl == pytest.approx(46.42262853199607, rel=1e-8)


def test_mixhum_ptrh_iswit_signs() -> None:
    p = 1000.0
    tk = 18.0 + 273.15
    rh = 46.5
    kwargs = dict(pressure_unit="hPa", temperature_unit="K", humidity_unit="%")

    mix_kg = pm.ncl.mixhum_ptrh(p, tk, rh, 1)
    mix_g = pm.ncl.mixhum_ptrh(p, tk, rh, -1)
    q_kg = pm.ncl.mixhum_ptrh(p, tk, rh, 2)
    q_g = pm.ncl.mixhum_ptrh(p, tk, rh, -2)

    assert mix_kg == pytest.approx(
        mixing_ratio_from_relative_humidity(p, tk, rh, output_humidity_unit="kg/kg", **kwargs)
    )
    assert mix_g == pytest.approx(
        mixing_ratio_from_relative_humidity(p, tk, rh, output_humidity_unit="g/kg", **kwargs)
    )
    assert q_kg == pytest.approx(
        specific_humidity_from_relative_humidity(p, tk, rh, output_humidity_unit="kg/kg", **kwargs)
    )
    assert q_g == pytest.approx(
        specific_humidity_from_relative_humidity(p, tk, rh, output_humidity_unit="g/kg", **kwargs)
    )
    assert mix_g == pytest.approx(mix_kg * 1000.0)
    assert q_g == pytest.approx(q_kg * 1000.0)


def test_mixhum_convert_ncl_wqtype_and_iounit() -> None:
    mixing_gkg = 15.2
    specific_gkg = pm.ncl.mixhum_convert(mixing_gkg, "w", (1, 1))
    modern = convert_humidity(
        mixing_gkg,
        from_quantity="mixing_ratio",
        to_quantity="specific_humidity",
        humidity_unit="g/kg",
        output_humidity_unit="g/kg",
    )
    assert specific_gkg == pytest.approx(modern)
    back = pm.ncl.mixhum_convert(specific_gkg, "q", (1, 1))
    assert back == pytest.approx(mixing_gkg)

    mixing_kg = 0.01
    specific_kg = pm.ncl.mixhum_convert(mixing_kg, "W", [0, 0])
    assert specific_kg == pytest.approx(
        convert_humidity(
            mixing_kg,
            from_quantity="mixing_ratio",
            to_quantity="specific_humidity",
            humidity_unit="kg/kg",
            output_humidity_unit="kg/kg",
        )
    )
    grams_out = pm.ncl.mixhum_convert(mixing_kg, "w", (0, 1))
    assert grams_out == pytest.approx(specific_kg * 1000.0)


def test_wind_speed_matches_modern_mps() -> None:
    ncl = pm.ncl.wind_speed(3.0, 4.0)
    modern = wind_speed(3.0, 4.0, speed_unit="m/s", output_speed_unit="m/s")
    assert ncl == pytest.approx(modern)
    assert ncl == pytest.approx(5.0)


def test_wind_direction_and_component_match_modern() -> None:
    direction = pm.ncl.wind_direction(3.0, 4.0)
    assert direction == pytest.approx(pm.wind_direction(3.0, 4.0))
    u, v = pm.ncl.wind_component(12.0, 247.5)
    mu, mv = uv_from_speed_direction(12.0, 247.5, speed_unit="m/s", output_speed_unit="m/s")
    assert u == pytest.approx(mu)
    assert v == pytest.approx(mv)
    calm_zero = pm.ncl.wind_direction(0.0, 0.0, 0)
    calm_missing = pm.ncl.wind_direction(0.0, 0.0, 1)
    calm_custom = pm.ncl.wind_direction(0.0, 0.0, -999)
    assert calm_zero == pytest.approx(0.0)
    assert np.isnan(calm_missing)
    assert calm_custom == pytest.approx(-999.0)


def test_ncl_wrappers_broadcast_arrays() -> None:
    tk = np.array([18.0, 20.0]) + 273.15
    rh = np.array([46.5, 80.0])
    dew = pm.ncl.dewtemp_trh(tk, rh)
    modern = dewpoint_from_relative_humidity(
        tk,
        rh,
        temperature_unit="K",
        humidity_unit="%",
        output_temperature_unit="K",
    )
    np.testing.assert_allclose(dew, modern)
    speed = pm.ncl.wind_speed(np.array([3.0, 0.0]), np.array([4.0, 5.0]))
    np.testing.assert_allclose(speed, [5.0, 5.0])


def test_ncl_all_exports_expected_names() -> None:
    expected = {
        "dewtemp_trh",
        "mixhum_convert",
        "mixhum_ptrh",
        "relhum",
        "relhum_ttd",
        "wind_component",
        "wind_direction",
        "wind_speed",
    }
    assert set(pm.ncl.__all__) == expected
    for name in expected:
        assert hasattr(pm.ncl, name)
    # 与现代 API 撞名的 wind_* 仍在顶层；NCL 专有名字不得再导出
    ncl_only = expected - {"wind_speed", "wind_direction"}
    for name in ncl_only:
        assert name not in pm.__all__
        assert not hasattr(pm, name)
