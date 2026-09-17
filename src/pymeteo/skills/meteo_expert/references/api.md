# Public API (root)

Re-exported at `import pymeteo as pm`. Source of truth is `pymeteo.__all__`. NCL shims are listed in `ncl.md` and are **not** in this table.

## thermo

| Function | Typical use |
|----------|-------------|
| `saturation_vapor_pressure` | Liquid es(T) |
| `condensation_temperature` | LCL temperature, 李社宏 iterator |
| `relative_humidity_from_dewpoint` | T, Td → RH |
| `dewpoint_from_relative_humidity` | T, RH → Td |
| `relative_humidity_from_mixing_ratio` | T, r, p → RH (NCL table) |
| `mixing_ratio_from_relative_humidity` | p, T, RH → r (Tetens) |
| `specific_humidity_from_relative_humidity` | p, T, RH → q |
| `convert_humidity` | mixing ratio ↔ specific humidity |
| `visibility` | RH, T → visibility (`RUC` or `FSL`) |
| `saturation_mixing_ratio` | p, T → rs |
| `mixing_ratio_from_dewpoint` | p, Td → r |
| `vapor_pressure_from_mixing_ratio` | p, r → e |
| `vapor_pressure_from_relative_humidity` | T, RH → e |
| `potential_temperature` | p, T → θ (output default K) |
| `equivalent_potential_temperature` | p, T, Td → θe (Bolton) |
| `virtual_temperature` | T, r → Tv |
| `wet_bulb_temperature` | T, RH → Tw (Stull, sea level) |
| `lifting_condensation_level` | p, T, Td → (p_LCL, T_LCL) |
| `parcel_temperature_at_pressure` | lift parcel to a target p |

## indices

| Function | Typical use |
|----------|-------------|
| `temperature_dewpoint_depression` | T − Td |
| `layer_temperature_difference` | T_lower − T_upper |
| `k_index` | K |
| `a_index` | A |
| `total_totals_index` | TT |
| `showalter_index` | SI |
| `sweat_index` | SWEAT |
| `lifted_index` | Tenv(500) − Tparcel(500) |
| `lifted_index_from_surface` | surface parcel LI |

## wind

| Function | Typical use |
|----------|-------------|
| `wind_speed` | (u, v) → speed |
| `wind_direction` | (u, v) → from-direction deg |
| `uv_from_speed_direction` | speed, dir → (u, v) |
| `wind_components` | alias of `uv_from_speed_direction` |
| `bulk_wind_shear` | two-level vector shear |

## geo

| Function | Typical use |
|----------|-------------|
| `earth_distance` | WGS84 geodesic (default km) |
| `gravity` | g(latitude) |
| `sea_level_pressure` | station p → MSLP |
| `height_thickness` | hypsometric ΔZ (default m) |

## dynamics

| Function | Typical use |
|----------|-------------|
| `coriolis_parameter` | f = 2Ω sinφ |
| `omega_to_w` | ω → w (modern order T then p) |
| `w_to_omega` | w → ω |

## comfort

| Function | Typical use |
|----------|-------------|
| `heat_index` | NWS heat index |
| `wind_chill` | NWS 2001 wind chill |

Also at the root: `UnitError`, `__version__`, and the `ncl` submodule.
