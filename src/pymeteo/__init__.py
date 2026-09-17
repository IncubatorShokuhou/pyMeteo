"""pymeteo：带显式单位参数的轻量气象诊断函数库。"""

from pymeteo import ncl
from pymeteo.geo import earth_distance, gravity, sea_level_pressure
from pymeteo.indices import (
    a_index,
    k_index,
    layer_temperature_difference,
    showalter_index,
    sweat_index,
    temperature_dewpoint_depression,
    total_totals_index,
)
from pymeteo.thermo import (
    condensation_temperature,
    convert_humidity,
    dewpoint_from_relative_humidity,
    mixing_ratio_from_relative_humidity,
    relative_humidity_from_dewpoint,
    relative_humidity_from_mixing_ratio,
    saturation_vapor_pressure,
    specific_humidity_from_relative_humidity,
    visibility,
)
from pymeteo.units import UnitError
from pymeteo.wind import (
    uv_from_speed_direction,
    wind_components,
    wind_direction,
    wind_speed,
)

__version__ = "2.1.0"

__all__ = [
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
    "ncl",
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
]
