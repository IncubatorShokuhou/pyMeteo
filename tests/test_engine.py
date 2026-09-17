"""MeteoEngine 知识面：对照真实公开函数名。"""

import math

import pytest

import pymeteo
from pymeteo.engine import MeteoEngine


def test_engine_is_not_on_classic_root() -> None:
    assert "MeteoEngine" not in pymeteo.__all__
    assert not hasattr(pymeteo, "MeteoEngine")


def test_list_functions_covers_public_callables() -> None:
    engine = MeteoEngine()
    listed = engine.list_functions()
    names = {item["name"] for item in listed}
    public = {
        name
        for name in pymeteo.__all__
        if name not in {"UnitError", "__version__", "ncl"}
    }
    assert names == public
    thermo = engine.list_functions(module="thermo")
    assert {item["name"] for item in thermo} <= names
    assert all(item["module"] == "thermo" for item in thermo)
    assert "dewpoint_from_relative_humidity" in {item["name"] for item in thermo}


def test_list_functions_rejects_unknown_module() -> None:
    engine = MeteoEngine()
    with pytest.raises(ValueError, match="module"):
        engine.list_functions(module="adengine")


def test_explain_real_function_signature() -> None:
    engine = MeteoEngine()
    info = engine.explain("dewpoint_from_relative_humidity")
    assert info["name"] == "dewpoint_from_relative_humidity"
    assert info["module"] == "thermo"
    assert "temperature" in info["signature"]
    assert "relative_humidity" in info["signature"]
    assert info["defaults"]["temperature_unit"] == "C"
    assert info["defaults"]["humidity_unit"] == "%"
    assert "callable" in info
    assert "dewpoint" in info["doc"].lower() or "露点" in info["doc"]


def test_explain_unknown_name() -> None:
    engine = MeteoEngine()
    with pytest.raises(ValueError, match="Unknown"):
        engine.explain("IForest")


def test_recommend_dewpoint_and_showalter() -> None:
    engine = MeteoEngine()
    dew = engine.recommend("dewpoint from relative humidity in percent")
    names = [item["name"] for item in dew]
    assert "dewpoint_from_relative_humidity" in names
    assert dew[0]["why"]
    si = engine.recommend("沙氏指数")
    assert "showalter_index" in [item["name"] for item in si]
    structured = engine.recommend(need={"quantity": "great-circle distance"})
    assert "earth_distance" in [item["name"] for item in structured]


def test_ncl_lookup_dewtemp() -> None:
    engine = MeteoEngine()
    hit = engine.ncl_lookup("dewtemp_trh")
    assert hit["ncl"] == "dewtemp_trh"
    assert hit["pymeteo"] == "dewpoint_from_relative_humidity"
    assert "K" in hit["units"]
    assert hit["callable"].endswith("dewtemp_trh")
    missing = engine.ncl_lookup("showalter_index")
    assert "error" in missing


def test_unit_help_aliases() -> None:
    engine = MeteoEngine()
    all_kinds = engine.unit_help()
    assert "temperature" in all_kinds
    assert "C" in all_kinds["temperature"]["canonical"]
    assert "mb" in all_kinds["pressure"]["aliases"]
    temp = engine.unit_help("temperature")
    assert temp["kind"] == "temperature"
    omega = engine.unit_help("omega")
    assert "Pa/s" in omega["canonical"]


def test_run_calc_scalar_dewpoint() -> None:
    engine = MeteoEngine()
    result = engine.run_calc(
        "dewpoint_from_relative_humidity",
        args=[18.0, 46.5],
        kwargs={"temperature_unit": "C", "humidity_unit": "%"},
    )
    assert result["ok"] is True
    assert result["value"] == pytest.approx(6.3, abs=0.05)


def test_run_calc_rejects_unknown_and_huge_payloads() -> None:
    engine = MeteoEngine()
    bad = engine.run_calc("eval", args=["1+1"])
    assert bad["ok"] is False
    huge = engine.run_calc("wind_speed", args=[list(range(20000)), list(range(20000))])
    assert huge["ok"] is False


def test_run_calc_serializes_tuple() -> None:
    engine = MeteoEngine()
    result = engine.run_calc("uv_from_speed_direction", args=[10.0, 90.0])
    assert result["ok"] is True
    assert isinstance(result["value"], list)
    assert len(result["value"]) == 2
    assert result["value"][0] == pytest.approx(-10.0)
    assert math.isfinite(result["value"][1])
