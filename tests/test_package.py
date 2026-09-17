"""包导入面与破坏性变更。"""

import pymeteo


def test_public_api_exports() -> None:
    expected = {
        "UnitError",
        "__version__",
        "a_index",
        "condensation_temperature",
        "convert_humidity",
        "dewpoint_from_relative_humidity",
        "earth_distance",
        "gravity",
        "k_index",
        "layer_temperature_difference",
        "mixing_ratio_from_relative_humidity",
        "relative_humidity_from_dewpoint",
        "relative_humidity_from_mixing_ratio",
        "saturation_vapor_pressure",
        "sea_level_pressure",
        "showalter_index",
        "specific_humidity_from_relative_humidity",
        "sweat_index",
        "temperature_dewpoint_depression",
        "total_totals_index",
        "uv_from_speed_direction",
        "visibility",
        "wind_components",
        "wind_direction",
        "wind_speed",
    }
    assert set(pymeteo.__all__) == expected
    for name in expected:
        if name == "__version__":
            assert isinstance(pymeteo.__version__, str)
        else:
            assert hasattr(pymeteo, name)


def test_version_present() -> None:
    assert pymeteo.__version__ == "2.0.0"


def test_old_single_file_names_are_gone() -> None:
    for legacy in (
        "showalter",
        "E_WATER",
        "Tc",
        "K",
        "A",
        "TT",
        "ws",
        "wd",
        "SWEAT_calculate",
        "relhum",
        "dewtemp_trh",
        "mixhum_ptrh",
        "mixhum_convert",
    ):
        assert not hasattr(pymeteo, legacy)
        assert legacy not in pymeteo.__all__
