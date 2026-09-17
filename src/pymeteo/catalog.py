"""Public-function catalog for agents: list, explain, recommend, NCL, units.

Introspects the live ``pymeteo`` package. This is knowledge/planning only;
it does not orchestrate multi-step sounding workflows.
"""

import inspect
import math

import numpy as np

import pymeteo
from pymeteo import ncl as ncl_mod
from pymeteo.units import unit_alias_tables

_PUBLIC_SKIP = frozenset({"UnitError", "__version__", "ncl"})
_MAX_RUN_CALC_VALUES = 10000

_MODULE_TAGS = {
    "thermo": (
        "thermo",
        "thermodynamic",
        "vapour",
        "vapor",
        "humidity",
        "dewpoint",
        "mixing",
        "theta",
        "lcl",
        "wet-bulb",
        "virtual",
        "visibility",
        "热力学",
        "露点",
        "湿度",
        "混合比",
        "位温",
        "湿球",
    ),
    "indices": (
        "index",
        "indices",
        "stability",
        "showalter",
        "sweat",
        "lifted",
        "thunderstorm",
        "指数",
        "沙氏",
        "稳定度",
        "对流",
    ),
    "wind": (
        "wind",
        "u-v",
        "uv",
        "shear",
        "direction",
        "风",
        "风向",
        "风速",
        "切变",
    ),
    "geo": (
        "geo",
        "distance",
        "great-circle",
        "gravity",
        "sea-level",
        "thickness",
        "距离",
        "重力",
        "海平面",
        "厚度",
    ),
    "dynamics": (
        "dynamics",
        "coriolis",
        "omega",
        "vertical",
        "科氏",
        "垂直速度",
    ),
    "comfort": (
        "comfort",
        "heat-index",
        "wind-chill",
        "体感",
        "热指数",
        "风寒",
    ),
}

# Phrase → function boosts for natural-language recommend().
_PHRASES = (
    ("dewpoint from relative humidity", "dewpoint_from_relative_humidity"),
    ("dew point", "dewpoint_from_relative_humidity"),
    ("dewpoint", "dewpoint_from_relative_humidity"),
    ("露点", "dewpoint_from_relative_humidity"),
    ("relative humidity from dewpoint", "relative_humidity_from_dewpoint"),
    ("relative humidity from mixing", "relative_humidity_from_mixing_ratio"),
    ("mixing ratio from dewpoint", "mixing_ratio_from_dewpoint"),
    ("mixing ratio from rh", "mixing_ratio_from_relative_humidity"),
    ("specific humidity", "specific_humidity_from_relative_humidity"),
    ("saturation vapor", "saturation_vapor_pressure"),
    ("saturation mixing", "saturation_mixing_ratio"),
    ("potential temperature", "potential_temperature"),
    ("equivalent potential", "equivalent_potential_temperature"),
    ("wet bulb", "wet_bulb_temperature"),
    ("wet-bulb", "wet_bulb_temperature"),
    ("virtual temperature", "virtual_temperature"),
    ("lifting condensation", "lifting_condensation_level"),
    ("condensation temperature", "condensation_temperature"),
    ("parcel temperature", "parcel_temperature_at_pressure"),
    ("visibility", "visibility"),
    ("showalter", "showalter_index"),
    ("沙氏", "showalter_index"),
    ("k index", "k_index"),
    ("k-index", "k_index"),
    ("sweat", "sweat_index"),
    ("lifted index from surface", "lifted_index_from_surface"),
    ("lifted index", "lifted_index"),
    ("抬升指数", "lifted_index_from_surface"),
    ("total totals", "total_totals_index"),
    ("a index", "a_index"),
    ("depression", "temperature_dewpoint_depression"),
    ("layer temperature", "layer_temperature_difference"),
    ("wind speed", "wind_speed"),
    ("wind direction", "wind_direction"),
    ("wind shear", "bulk_wind_shear"),
    ("bulk shear", "bulk_wind_shear"),
    ("u v", "uv_from_speed_direction"),
    ("u/v", "uv_from_speed_direction"),
    ("great-circle", "earth_distance"),
    ("great circle", "earth_distance"),
    ("distance", "earth_distance"),
    ("sea-level", "sea_level_pressure"),
    ("sea level", "sea_level_pressure"),
    ("thickness", "height_thickness"),
    ("gravity", "gravity"),
    ("coriolis", "coriolis_parameter"),
    ("科氏", "coriolis_parameter"),
    ("omega", "omega_to_w"),
    ("heat index", "heat_index"),
    ("wind chill", "wind_chill"),
    ("风寒", "wind_chill"),
    ("热指数", "heat_index"),
)

_NCL_MAP = {
    "dewtemp_trh": {
        "pymeteo": "dewpoint_from_relative_humidity",
        "units": "tk in K, rh in %, dewpoint returned in K",
        "notes": "Modern API defaults to C and %; NCL is K in / K out.",
    },
    "relhum_ttd": {
        "pymeteo": "relative_humidity_from_dewpoint",
        "units": "t and td in K; opt=0 → %, opt=1 → fraction",
        "notes": "Modern API takes temperature_unit and output_humidity_unit.",
    },
    "relhum": {
        "pymeteo": "relative_humidity_from_mixing_ratio",
        "units": "t in K, w in kg/kg, p in Pa, result in %",
        "notes": "RH > 100 is allowed; RH < 0 is clipped to 0.0001%. Lookup table 173.16–375.16 K.",
    },
    "mixhum_ptrh": {
        "pymeteo": "mixing_ratio_from_relative_humidity",
        "units": "p in hPa (not Pa), tk in K, rh in %; |iswit|=1 mixing ratio, 2 specific humidity; negative → g/kg",
        "notes": "iswit=2 maps to specific_humidity_from_relative_humidity.",
    },
    "mixhum_ptd": {
        "pymeteo": "mixing_ratio_from_dewpoint",
        "units": "p in Pa (not hPa), tdk in K; iswit as mixhum_ptrh",
        "notes": "Pressure unit differs from mixhum_ptrh.",
    },
    "mixhum_convert": {
        "pymeteo": "convert_humidity",
        "units": "wqType w/W mixing→specific, q/Q reverse; iounit=(in,out) 0=kg/kg 1=g/kg",
        "notes": None,
    },
    "vapor_pres_rh": {
        "pymeteo": "vapor_pressure_from_relative_humidity",
        "units": "rh in %; es and the result share a unit. NCL form is RH/100 · e_s (no temperature).",
        "notes": "Modern vapor_pressure_from_relative_humidity needs temperature; NCL wrapper does not.",
    },
    "pot_temp": {
        "pymeteo": "potential_temperature",
        "units": "p in Pa, t in K, result in K. dim/opt ignored",
        "notes": "Modern defaults: pressure_unit=hPa, temperature_unit=C, output K.",
    },
    "pot_temp_equiv": {
        "pymeteo": "equivalent_potential_temperature",
        "units": "p in Pa, t in K; humVarType r/w mixing kg/kg, q specific humidity, rh in %",
        "notes": "Bolton with LCL; closer to NCL pot_temp_equiv_tlcl than the no-LCL approx.",
    },
    "temp_virtual": {
        "pymeteo": "virtual_temperature",
        "units": "iounit length 3: T 0=C/1=K/2=F, mixing 0=kg/kg 1=g/kg, output T",
        "notes": "Uses T(1+r/ε)/(1+r), not T(1+0.61 r).",
    },
    "wetbulb_stull": {
        "pymeteo": "wet_bulb_temperature",
        "units": "rh in %; iounit length 2 (0=C, 1=K, 2=F). opt unused",
        "notes": "Sea-level Stull 2011 only.",
    },
    "lclvl": {
        "pymeteo": "lifting_condensation_level",
        "units": "p in hPa, tk/tdk in K; returns LCL pressure only (hPa)",
        "notes": "Modern function also returns LCL temperature.",
    },
    "wind_speed": {
        "pymeteo": "wind_speed",
        "units": "m/s in and out",
        "notes": None,
    },
    "wind_direction": {
        "pymeteo": "wind_direction",
        "units": "from-direction degrees; calm: opt=0 → 0, opt=1 → nan, other scalar → fill",
        "notes": "Modern wind_direction always returns 0 for calm.",
    },
    "wind_component": {
        "pymeteo": "uv_from_speed_direction",
        "units": "m/s, from-direction → (u, v) tuple. NCL opt unused",
        "notes": "Python returns a tuple, not an NCL stacked array.",
    },
    "coriolis_param": {
        "pymeteo": "coriolis_parameter",
        "units": "latitude in degrees → s^-1",
        "notes": None,
    },
    "omega_to_w": {
        "pymeteo": "omega_to_w",
        "units": "NCL order (omega, p, t) with p in Pa, t in K → m/s",
        "notes": "Modern order is (omega, temperature, pressure); modern pressure default is hPa.",
    },
    "w_to_omega": {
        "pymeteo": "w_to_omega",
        "units": "NCL order (w, p, t)",
        "notes": "Modern order is (w, temperature, pressure).",
    },
}

_KIND_DEFAULTS = {
    "temperature": "C (potential temperature output defaults to K)",
    "pressure": "hPa (NCL pot_temp / mixhum_ptd / omega use Pa)",
    "speed": "m/s",
    "humidity": "%",
    "mass_humidity": "kg/kg",
    "distance": "km for earth_distance; m for height_thickness",
    "angle": "deg",
    "omega": "Pa/s",
}


def public_function_names():
    """Return sorted public callable names re-exported at the package root."""

    names = []
    for name in pymeteo.__all__:
        if name in _PUBLIC_SKIP:
            continue
        names.append(name)
    return sorted(names)


def _module_of(name):
    obj = getattr(pymeteo, name)
    module = getattr(obj, "__module__", "")
    if module.startswith("pymeteo."):
        return module.split(".")[1]
    return module


def _first_paragraph(doc):
    if not doc:
        return ""
    lines = []
    for line in inspect.cleandoc(doc).splitlines():
        if not line.strip():
            break
        lines.append(line.strip())
    return " ".join(lines)


def _jsonable_default(value):
    if value is None or isinstance(value, (bool, str)):
        return value
    if isinstance(value, (int, float)):
        if isinstance(value, float) and not math.isfinite(value):
            return None
        return value
    return repr(value)


def _record(name):
    func = getattr(pymeteo, name)
    module = _module_of(name)
    doc = inspect.getdoc(func) or ""
    try:
        signature = str(inspect.signature(func))
    except (TypeError, ValueError):
        signature = "(...)"
    defaults = {}
    try:
        sig = inspect.signature(func)
        for key, param in sig.parameters.items():
            if param.default is not inspect.Parameter.empty:
                defaults[key] = _jsonable_default(param.default)
    except (TypeError, ValueError):
        pass
    tags = list(_MODULE_TAGS.get(module, ()))
    tags.extend(part for part in name.split("_") if part)
    return {
        "name": name,
        "module": module,
        "signature": name + signature,
        "defaults": defaults,
        "doc": doc,
        "summary": _first_paragraph(doc),
        "callable": "pymeteo." + name,
        "tags": tags,
    }


def list_functions(module=None):
    """List public functions, optionally filtered by module name."""

    wanted = None
    if module:
        wanted = module.split(".")[-1]
        known = set(_MODULE_TAGS)
        if wanted not in known:
            raise ValueError("unknown module %r; expected one of %s" % (module, ", ".join(sorted(known))))
    records = []
    for name in public_function_names():
        rec = _record(name)
        if wanted is not None and rec["module"] != wanted:
            continue
        records.append(rec)
    return records


def explain_function(name):
    """Return signature, docstring, and defaults for a public or NCL name."""

    key = name.strip()
    if key.startswith("pymeteo."):
        key = key.split(".", 1)[1]
    if key.startswith("ncl."):
        ncl_name = key.split(".", 1)[1]
        return _explain_ncl(ncl_name)
    if key in pymeteo.__all__ and key not in _PUBLIC_SKIP:
        return _record(key)
    if key in ncl_mod.__all__:
        return _explain_ncl(key)
    raise ValueError("Unknown function: %s" % name)


def _explain_ncl(ncl_name):
    if ncl_name not in ncl_mod.__all__:
        raise ValueError("Unknown function: ncl.%s" % ncl_name)
    func = getattr(ncl_mod, ncl_name)
    doc = inspect.getdoc(func) or ""
    try:
        signature = str(inspect.signature(func))
    except (TypeError, ValueError):
        signature = "(...)"
    mapped = dict(_NCL_MAP.get(ncl_name, {}))
    mapped.update(
        {
            "name": ncl_name,
            "module": "ncl",
            "signature": ncl_name + signature,
            "doc": doc,
            "summary": _first_paragraph(doc),
            "callable": "pymeteo.ncl." + ncl_name,
            "ncl": ncl_name,
        }
    )
    if "pymeteo" not in mapped:
        mapped["pymeteo"] = ncl_name
    return mapped


def recommend(query="", need=None, top_k=5):
    """Score public functions against a natural-language or structured need."""

    parts = []
    if query:
        parts.append(str(query))
    if isinstance(need, dict):
        for key in ("quantity", "inputs", "module", "ncl"):
            value = need.get(key)
            if value:
                if isinstance(value, (list, tuple)):
                    parts.append(" ".join(str(item) for item in value))
                else:
                    parts.append(str(value))
        extra = need.get("text") or need.get("query")
        if extra:
            parts.append(str(extra))
    text = " ".join(parts).strip()
    if not text:
        raise ValueError("recommend() needs a query string or a need dict")
    text_l = text.lower()
    haystack = text_l + " " + text

    scores = {}
    reasons = {}
    for rec in list_functions():
        name = rec["name"]
        score = 0
        why = []
        lowered = name.replace("_", " ")
        if name in haystack or lowered in text_l:
            score += 12
            why.append("name match")
        for token in name.split("_"):
            if len(token) > 2 and token in text_l:
                score += 2
        for tag in rec["tags"]:
            needle = tag.lower()
            if needle and needle in haystack:
                score += 3
                why.append("keyword %s" % tag)
        scores[name] = score
        reasons[name] = why

    for phrase, name in _PHRASES:
        if phrase.lower() in haystack or phrase in haystack:
            scores[name] = scores.get(name, 0) + 20
            reasons.setdefault(name, []).append("phrase %s" % phrase)

    for ncl_name, meta in _NCL_MAP.items():
        if ncl_name.lower() in text_l:
            modern = meta["pymeteo"]
            scores[modern] = scores.get(modern, 0) + 18
            reasons.setdefault(modern, []).append("NCL %s" % ncl_name)

    if isinstance(need, dict) and need.get("module"):
        module = str(need["module"]).split(".")[-1]
        for rec in list_functions():
            if rec["module"] == module:
                scores[rec["name"]] = scores.get(rec["name"], 0) + 4

    ranked = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    out = []
    for name, score in ranked:
        if score <= 0:
            continue
        rec = _record(name)
        why = []
        seen = set()
        for item in reasons.get(name, []):
            if item not in seen:
                seen.add(item)
                why.append(item)
        out.append(
            {
                "name": name,
                "module": rec["module"],
                "score": score,
                "why": why or ["token overlap"],
                "summary": rec["summary"],
                "callable": rec["callable"],
            }
        )
        if len(out) >= top_k:
            break
    return out


def ncl_lookup(name):
    """Map an NCL builtin name to the modern pymeteo callable."""

    key = name.strip()
    if key.startswith("pymeteo.ncl."):
        key = key.split(".")[-1]
    if key.startswith("ncl."):
        key = key.split(".", 1)[1]
    meta = _NCL_MAP.get(key)
    if meta is None or key not in ncl_mod.__all__:
        return {
            "error": "No NCL shim named %r. Showalter, K, SWEAT, heat index, and wind chill are modern-only."
            % name
        }
    result = {
        "ncl": key,
        "pymeteo": meta["pymeteo"],
        "units": meta["units"],
        "callable": "pymeteo.ncl." + key,
        "modern_callable": "pymeteo." + meta["pymeteo"],
    }
    if meta.get("notes"):
        result["notes"] = meta["notes"]
    return result


def _kind_payload(kind, alias_map, extra_canonical=None):
    canonical = sorted(set(alias_map.values()))
    if extra_canonical:
        for item in extra_canonical:
            if item not in canonical:
                canonical.append(item)
    return {
        "kind": kind,
        "canonical": canonical,
        "aliases": sorted(alias_map.keys()),
        "default": _KIND_DEFAULTS.get(kind, ""),
    }


def unit_help(kind=None):
    """Describe string unit aliases and defaults."""

    tables = unit_alias_tables()
    omega = {
        "pa/s": "Pa/s",
        "pas-1": "Pa/s",
        "pa s-1": "Pa/s",
        "hpa/s": "hPa/s",
        "mb/s": "hPa/s",
    }
    payloads = {}
    for key, alias_map in tables.items():
        payloads[key] = _kind_payload(key, alias_map)
    payloads["omega"] = _kind_payload("omega", omega, extra_canonical=["Pa/s", "hPa/s"])
    if kind is None or kind == "":
        return payloads
    key = kind.strip().lower().replace(" ", "_")
    aliases = {
        "rh": "humidity",
        "relative_humidity": "humidity",
        "mixing_ratio": "mass_humidity",
        "specific_humidity": "mass_humidity",
        "temp": "temperature",
        "wind": "speed",
        "w": "omega",
    }
    key = aliases.get(key, key)
    if key not in payloads:
        raise ValueError("unknown unit kind %r" % kind)
    return payloads[key]


def _count_values(value):
    if isinstance(value, list):
        total = 0
        for item in value:
            total += _count_values(item)
        return total
    return 1


def _ensure_run_calc_value(value):
    if isinstance(value, (bool, str)) or value is None:
        return
    if isinstance(value, (int, float)):
        return
    if isinstance(value, list):
        for item in value:
            _ensure_run_calc_value(item)
        return
    raise TypeError("run_calc only accepts JSON numbers, strings, bools, null, and lists")


def _to_jsonable(value):
    if isinstance(value, tuple):
        return [_to_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return _to_jsonable(value.item())
        return _to_jsonable(value.tolist())
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return _to_jsonable(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, list):
        return [_to_jsonable(item) for item in value]
    return value


def resolve_callable(name):
    """Resolve a public or NCL name to a Python callable."""

    key = name.strip()
    if key.startswith("pymeteo.ncl."):
        key = "ncl." + key.split(".")[-1]
    elif key.startswith("pymeteo."):
        key = key.split(".", 1)[1]
    if key.startswith("ncl."):
        ncl_name = key.split(".", 1)[1]
        if ncl_name not in ncl_mod.__all__:
            raise ValueError("Unknown function: %s" % name)
        return getattr(ncl_mod, ncl_name)
    if key in pymeteo.__all__ and key not in _PUBLIC_SKIP:
        return getattr(pymeteo, key)
    if key in ncl_mod.__all__:
        return getattr(ncl_mod, key)
    raise ValueError("Unknown function: %s" % name)


def run_calc(name, args=None, kwargs=None):
    """Call a whitelisted function with JSON-serializable scalars/lists.

    Does not execute arbitrary strings. Arrays larger than
    ``_MAX_RUN_CALC_VALUES`` values are rejected.
    """

    args = [] if args is None else args
    kwargs = {} if kwargs is None else kwargs
    if not isinstance(args, list):
        return {"ok": False, "error": "args must be a JSON list"}
    if not isinstance(kwargs, dict):
        return {"ok": False, "error": "kwargs must be a JSON object"}
    try:
        _ensure_run_calc_value(args)
        for key, value in kwargs.items():
            if not isinstance(key, str) or not key.isidentifier():
                raise TypeError("invalid kwarg name")
            _ensure_run_calc_value(value)
        n_values = _count_values(args) + sum(_count_values(v) for v in kwargs.values())
        if n_values > _MAX_RUN_CALC_VALUES:
            raise ValueError("payload has %d values; max is %d" % (n_values, _MAX_RUN_CALC_VALUES))
        func = resolve_callable(name)
    except (TypeError, ValueError) as exc:
        return {"ok": False, "error": str(exc)}
    try:
        result = func(*args, **kwargs)
    except Exception as exc:
        return {"ok": False, "error": "%s: %s" % (type(exc).__name__, exc)}
    return {"ok": True, "name": name, "value": _to_jsonable(result)}
