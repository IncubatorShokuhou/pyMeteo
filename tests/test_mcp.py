"""MCP 模块：无 mcp 包时导入安全；有 mcp 时注册工具。"""

import json

import pytest

from pymeteo import mcp_server


def test_mcp_module_imports_without_running_server() -> None:
    assert mcp_server._check_mcp is not None
    payload = json.loads(mcp_server.list_functions())
    names = {item["name"] for item in payload}
    assert "showalter_index" in names
    explained = json.loads(mcp_server.explain_function("wind_speed"))
    assert explained["module"] == "wind"
    rec = json.loads(mcp_server.recommend("Coriolis parameter"))
    assert rec[0]["name"] == "coriolis_parameter"
    ncl = json.loads(mcp_server.ncl_lookup("pot_temp"))
    assert ncl["pymeteo"] == "potential_temperature"
    units = json.loads(mcp_server.unit_help("speed"))
    assert "kt" in units["aliases"]
    calc = json.loads(mcp_server.run_calc("wind_speed", args_json="[3, 4]"))
    assert calc["ok"] is True
    assert calc["value"] == 5.0


def test_mcp_import_does_not_exit_when_mcp_missing(monkeypatch) -> None:
    monkeypatch.setattr(mcp_server, "_check_mcp", lambda: None)
    assert mcp_server.main.__name__ == "main"
    assert mcp_server._check_mcp() is None


@pytest.mark.skipif(mcp_server._check_mcp() is None, reason="mcp extra not installed")
def test_fastmcp_tool_registration_smoke() -> None:
    FastMCP = mcp_server._check_mcp()
    server = FastMCP("pymeteo-test")
    for fn in mcp_server._TOOL_FUNCTIONS:
        server.tool()(fn)
    tools = getattr(server, "_tool_manager", None)
    if tools is not None:
        names = set(getattr(tools, "_tools", {}))
        if names:
            assert "list_functions" in names
            assert "recommend" in names
            assert "run_calc" in names
