"""NCL 兼容薄封装：单位约定应与调用现代 API 等价。"""

import numpy as np
import pytest

import pymeteo as pm
from pymeteo import (
    convert_humidity,
    dewpoint_from_relative_humidity,
    mixing_ratio_from_dewpoint,
    mixing_ratio_from_relative_humidity,
    potential_temperature,
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
        "coriolis_param",
        "dewtemp_trh",
        "lclvl",
        "mixhum_convert",
        "mixhum_ptd",
        "mixhum_ptrh",
        "omega_to_w",
        "pot_temp",
        "pot_temp_equiv",
        "relhum",
        "relhum_ttd",
        "temp_virtual",
        "vapor_pres_rh",
        "w_to_omega",
        "wetbulb_stull",
        "wind_component",
        "wind_direction",
        "wind_speed",
    }
    assert set(pm.ncl.__all__) == expected
    for name in expected:
        assert hasattr(pm.ncl, name)
    # 与现代 API 撞名的 wind_* / omega_* 仍在顶层；NCL 专有名字不得再导出
    ncl_only = expected - {"wind_speed", "wind_direction", "omega_to_w", "w_to_omega"}
    for name in ncl_only:
        assert name not in pm.__all__
        assert not hasattr(pm, name)


def test_mixhum_ptd_matches_modern() -> None:
    p = 100000.0
    tdk = 6.4 + 273.15
    mix_kg = pm.ncl.mixhum_ptd(p, tdk, 1)
    mix_g = pm.ncl.mixhum_ptd(p, tdk, -1)
    modern = mixing_ratio_from_dewpoint(
        p, tdk, pressure_unit="Pa", temperature_unit="K", output_humidity_unit="kg/kg"
    )
    assert mix_kg == pytest.approx(modern)
    assert mix_g == pytest.approx(mix_kg * 1000.0)
    q_kg = pm.ncl.mixhum_ptd(p, tdk, 2)
    assert q_kg == pytest.approx(mix_kg / (1.0 + mix_kg))


def test_pot_temp_and_equiv_ncl_units() -> None:
    theta = pm.ncl.pot_temp(85000.0, 10.0 + 273.15)
    modern = potential_temperature(
        85000.0, 10.0 + 273.15, pressure_unit="Pa", temperature_unit="K", output_temperature_unit="K"
    )
    assert theta == pytest.approx(modern)
    from pymeteo import equivalent_potential_temperature, relative_humidity_from_dewpoint

    td = 0.0 + 273.15
    t = 15.0 + 273.15
    p = 100000.0
    mixing = mixing_ratio_from_dewpoint(
        p, td, pressure_unit="Pa", temperature_unit="K", output_humidity_unit="kg/kg"
    )
    theta_e = pm.ncl.pot_temp_equiv(p, t, mixing, -1, "r")
    modern_e = equivalent_potential_temperature(
        p, t, td, pressure_unit="Pa", temperature_unit="K", output_temperature_unit="K"
    )
    assert theta_e == pytest.approx(modern_e, rel=1e-4)
    rh = relative_humidity_from_dewpoint(t, td, temperature_unit="K", output_humidity_unit="%")
    theta_e_rh = pm.ncl.pot_temp_equiv(p, t, rh, -1, "rh")
    assert theta_e_rh == pytest.approx(modern_e, rel=1e-4)


def test_temp_virtual_iounit_and_wetbulb_lcl_coriolis_omega() -> None:
    from pymeteo import (
        coriolis_parameter,
        lifting_condensation_level,
        virtual_temperature,
        wet_bulb_temperature,
    )

    tv = pm.ncl.temp_virtual(20.0, 13.5, (0, 1, 0))
    modern_tv = virtual_temperature(
        20.0, 13.5, temperature_unit="C", mixing_ratio_unit="g/kg", output_temperature_unit="C"
    )
    assert tv == pytest.approx(modern_tv)
    tw = pm.ncl.wetbulb_stull(20.0, 50.0, (0, 0), False)
    assert tw == pytest.approx(
        wet_bulb_temperature(20.0, 50.0, temperature_unit="C", humidity_unit="%")
    )
    tw_k = pm.ncl.wetbulb_stull(20.0 + 273.15, 50.0, (1, 1))
    assert tw_k == pytest.approx(tw + 273.15, rel=1e-8)
    plcl = pm.ncl.lclvl(1000.0, 15.0 + 273.15, 4.0 + 273.15)
    modern_p, _ = lifting_condensation_level(
        1000.0, 15.0 + 273.15, 4.0 + 273.15, temperature_unit="K"
    )
    assert plcl == pytest.approx(modern_p)
    assert pm.ncl.coriolis_param(35.0) == pytest.approx(
        coriolis_parameter(35.0, latitude_unit="deg")
    )
    omega = 0.05
    p = 85000.0
    t = 273.15
    w = pm.ncl.omega_to_w(omega, p, t)
    back = pm.ncl.w_to_omega(w, p, t)
    assert back == pytest.approx(omega)
    from pymeteo import omega_to_w as modern_omega_to_w

    modern_w = modern_omega_to_w(
        omega, t, p, omega_unit="Pa/s", temperature_unit="K", pressure_unit="Pa"
    )
    assert w == pytest.approx(modern_w)
    vapor = pm.ncl.vapor_pres_rh(50.0, 23.37)
    assert vapor == pytest.approx(11.685)
