"""包导入面与破坏性变更。"""

from importlib.metadata import metadata, version

import pymeteo


def test_public_api_exports() -> None:
    expected = {
        "UnitError",
        "__version__",
        "a_index",
        "bulk_wind_shear",
        "condensation_temperature",
        "convert_humidity",
        "coriolis_parameter",
        "dewpoint_from_relative_humidity",
        "earth_distance",
        "equivalent_potential_temperature",
        "gravity",
        "heat_index",
        "height_thickness",
        "k_index",
        "layer_temperature_difference",
        "lifting_condensation_level",
        "lifted_index",
        "lifted_index_from_surface",
        "mixing_ratio_from_dewpoint",
        "mixing_ratio_from_relative_humidity",
        "ncl",
        "omega_to_w",
        "parcel_temperature_at_pressure",
        "potential_temperature",
        "relative_humidity_from_dewpoint",
        "relative_humidity_from_mixing_ratio",
        "saturation_mixing_ratio",
        "saturation_vapor_pressure",
        "sea_level_pressure",
        "showalter_index",
        "specific_humidity_from_relative_humidity",
        "sweat_index",
        "temperature_dewpoint_depression",
        "total_totals_index",
        "uv_from_speed_direction",
        "vapor_pressure_from_mixing_ratio",
        "vapor_pressure_from_relative_humidity",
        "virtual_temperature",
        "visibility",
        "w_to_omega",
        "wet_bulb_temperature",
        "wind_chill",
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
    assert pymeteo.__version__ == "2.2.2"


def test_pypi_distribution_name() -> None:
    meta = metadata("py-meteo")
    assert meta["Name"].replace("_", "-").lower() == "py-meteo"
    assert version("py-meteo") == "2.2.2"


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
