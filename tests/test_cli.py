"""pymeteo CLI：skill 安装路径与 info 输出。"""

from pathlib import Path

import pymeteo
from pymeteo.cli import main


def test_info_prints_version_and_paths(capsys) -> None:
    assert main(["info"]) == 0
    out = capsys.readouterr().out
    assert pymeteo.__version__ in out
    assert "Classic API: OK" in out
    assert "meteo-expert" in out
    assert "MCP extra:" in out


def test_install_skill_project(tmp_path: Path, monkeypatch, capsys) -> None:
    monkeypatch.chdir(tmp_path)
    assert main(["install-skill", "--project"]) == 0
    skill = tmp_path / "skills" / "meteo-expert" / "SKILL.md"
    assert skill.is_file()
    out = capsys.readouterr().out
    assert "meteo-expert" in out
    assert str(skill) in out


def test_install_skill_target(tmp_path: Path, capsys) -> None:
    target = tmp_path / "agent-skills"
    assert main(["install-skill", "--target", str(target)]) == 0
    assert (target / "meteo-expert" / "SKILL.md").is_file()
    assert (target / "meteo-expert" / "references" / "units.md").is_file()


def test_install_skill_list(capsys) -> None:
    assert main(["install-skill", "--list"]) == 0
    assert "meteo-expert" in capsys.readouterr().out


def test_help_without_args(capsys) -> None:
    assert main([]) == 0
    out = capsys.readouterr().out
    assert "install-skill" in out
    assert "mcp" in out


def test_mcp_serve_without_extra_exits_1(monkeypatch) -> None:
    import pymeteo.mcp_server as mcp_server

    monkeypatch.setattr(mcp_server, "_check_mcp", lambda: None)
    assert main(["mcp", "serve"]) == 1
