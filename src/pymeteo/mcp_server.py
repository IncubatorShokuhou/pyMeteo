"""pymeteo MCP server: stateless knowledge/planning tools.

Usage::

    python -m pymeteo.mcp_server

Importing this module is safe when the optional ``mcp`` extra is missing.
``main()`` prints an install hint and returns 1 in that case. Core pymeteo
stays on Python 3.7; the ``mcp`` package itself needs a newer Python.
"""

import importlib.util
import json
import sys

from pymeteo.engine import MeteoEngine

_engine = None


def _check_mcp():
    """Return FastMCP if the optional extra is installed, else None.

    Probe by importing. ``find_spec("mcp.server.fastmcp")`` can succeed on
    mcp 2.x even though the module raises and FastMCP was renamed.
    """

    if importlib.util.find_spec("mcp") is None:
        return None
    try:
        from mcp.server.fastmcp import FastMCP
    except (ModuleNotFoundError, ImportError):
        return None
    return FastMCP


def _get_engine():
    global _engine
    if _engine is None:
        _engine = MeteoEngine()
    return _engine


def _to_json(obj):
    return json.dumps(obj, indent=2, default=str, ensure_ascii=False)


def _parse_json_object_or_list(text, empty):
    if text is None or text == "":
        return empty
    try:
        parsed = json.loads(text)
    except (json.JSONDecodeError, TypeError) as exc:
        raise ValueError(f"Invalid JSON: {exc}") from exc
    return parsed


def list_functions(module=""):
    """List pymeteo public functions.

    Args:
        module: Optional module filter: thermo, indices, wind, geo, dynamics, comfort.
    """

    try:
        return _to_json(_get_engine().list_functions(module=module or None))
    except ValueError as exc:
        return _to_json({"error": str(exc)})


def explain_function(name):
    """Explain a pymeteo function or NCL shim: signature, defaults, docstring."""

    try:
        return _to_json(_get_engine().explain(name))
    except ValueError as exc:
        return _to_json({"error": str(exc)})


def recommend(query="", need_json=""):
    """Recommend pymeteo functions for a meteorological question.

    Args:
        query: Natural-language need, e.g. "dewpoint from RH".
        need_json: Optional JSON object with quantity/inputs/module/ncl.
    """

    need = None
    if need_json:
        try:
            need = _parse_json_object_or_list(need_json, None)
        except ValueError as exc:
            return _to_json({"error": str(exc)})
        if need is not None and not isinstance(need, dict):
            return _to_json({"error": "need_json must be a JSON object"})
    try:
        return _to_json(_get_engine().recommend(query=query, need=need))
    except ValueError as exc:
        return _to_json({"error": str(exc)})


def ncl_lookup(name):
    """Map an NCL builtin name to the modern pymeteo callable and unit notes."""

    return _to_json(_get_engine().ncl_lookup(name))


def unit_help(kind=""):
    """Describe string unit aliases and defaults (temperature, pressure, ...)."""

    try:
        return _to_json(_get_engine().unit_help(kind or None))
    except ValueError as exc:
        return _to_json({"error": str(exc)})


def run_calc(name, args_json="[]", kwargs_json="{}"):
    """Call a whitelisted pymeteo function with JSON scalars/lists.

    Safe: no eval, no imports, no file I/O. Rejects large arrays.
    Use this for point checks, not for arbitrary code execution.

    Args:
        name: Public function name, or ncl.<name>.
        args_json: JSON list of positional arguments.
        kwargs_json: JSON object of keyword arguments (unit strings, etc.).
    """

    try:
        args = _parse_json_object_or_list(args_json, [])
        kwargs = _parse_json_object_or_list(kwargs_json, {})
    except ValueError as exc:
        return _to_json({"ok": False, "error": str(exc)})
    return _to_json(_get_engine().run_calc(name, args=args, kwargs=kwargs))


_TOOL_FUNCTIONS = (
    list_functions,
    explain_function,
    recommend,
    ncl_lookup,
    unit_help,
    run_calc,
)


def main():
    """Entry point for ``python -m pymeteo.mcp_server`` and ``pymeteo mcp serve``."""

    FastMCP = _check_mcp()
    if FastMCP is None:
        print(
            "pymeteo MCP server requires the 'mcp' package. "
            "Install with: pip install 'pymeteo-kit[mcp]'",
            file=sys.stderr,
        )
        return 1

    mcp = FastMCP("pymeteo")
    for fn in _TOOL_FUNCTIONS:
        mcp.tool()(fn)
    mcp.run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
