"""包导入面与破坏性变更。"""

import pymeteo


def test_public_api_exports() -> None:
    for name in (
        "showalter_index",
        "k_index",
        "a_index",
        "total_totals_index",
        "sweat_index",
        "saturation_vapor_pressure",
        "condensation_temperature",
        "relative_humidity_from_dewpoint",
        "dewpoint_from_relative_humidity",
        "relative_humidity_from_mixing_ratio",
        "mixing_ratio_from_relative_humidity",
        "convert_humidity",
        "wind_speed",
        "wind_direction",
        "wind_components",
        "uv_from_speed_direction",
        "earth_distance",
        "gravity",
        "sea_level_pressure",
    ):
        assert name in pymeteo.__all__
        assert callable(getattr(pymeteo, name))


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
