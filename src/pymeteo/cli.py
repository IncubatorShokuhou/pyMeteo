"""Command-line interface: skill install, MCP serve, and diagnostics."""

import argparse
import importlib.util
import shutil
import sys
from pathlib import Path

from pymeteo.skills import packaged_skill_names, run_install


def _mcp_available():
    if importlib.util.find_spec("mcp") is None:
        return False
    try:
        return importlib.util.find_spec("mcp.server.fastmcp") is not None
    except ModuleNotFoundError:
        return False


def _cmd_install_skill(args):
    return run_install(
        target=args.target,
        project=args.project,
        skill=args.skill,
        list_skills=args.list_skills,
    )


def _cmd_info(args):
    import pymeteo
    from pymeteo.engine import MeteoEngine

    engine = MeteoEngine()
    listed = engine.list_functions()
    counts = {}
    for item in listed:
        counts[item["module"]] = counts.get(item["module"], 0) + 1

    classic_ok = True
    try:
        from pymeteo import dewpoint_from_relative_humidity  # noqa: F401
    except Exception:
        classic_ok = False

    claude_dir = Path.home() / ".claude"
    codex_dir = Path.home() / ".codex"
    user_skill = claude_dir / "skills" / "meteo-expert" / "SKILL.md"
    project_skill = Path.cwd() / "skills" / "meteo-expert" / "SKILL.md"
    user_installed = user_skill.is_file()
    project_installed = project_skill.is_file()

    agents = []
    if shutil.which("claude") is not None or (Path.home() / ".claude.json").exists() or claude_dir.is_dir():
        agents.append("Claude Code")
    if codex_dir.is_dir():
        agents.append("Codex")

    print("pymeteo version: %s" % pymeteo.__version__)
    parts = []
    for key in ("thermo", "indices", "wind", "geo", "dynamics", "comfort"):
        if counts.get(key):
            parts.append("%d %s" % (counts[key], key))
    print("Functions: %d total (%s)" % (len(listed), ", ".join(parts)))
    print("Classic API: %s" % ("OK" if classic_ok else "ERROR"))
    print("MeteoEngine: OK")
    if _mcp_available():
        print("MCP extra: OK (run: pymeteo mcp serve)")
    else:
        print("MCP extra: NOT INSTALLED (install: pip install 'pymeteo-kit[mcp]')")

    if user_installed and project_installed:
        print("meteo-expert skill: INSTALLED (user-global) at %s" % user_skill)
        print("                    INSTALLED (project) at %s" % project_skill)
    elif user_installed:
        print("meteo-expert skill: INSTALLED (user-global) at %s" % user_skill)
        if "Codex" in agents:
            print("  Codex does not read ~/.claude/skills/. Use: pymeteo install-skill --project")
    elif project_installed:
        print("meteo-expert skill: INSTALLED (project) at %s" % project_skill)
        if "Claude Code" in agents:
            print("  For a user-global Claude Code install: pymeteo install-skill")
    else:
        print("meteo-expert skill: NOT INSTALLED")
        print("  Claude Code: pymeteo install-skill")
        print("  Codex:       pymeteo install-skill --project")
    if agents:
        print("Detected agents: %s" % ", ".join(agents))
    print("Packaged skills: %s" % ", ".join(packaged_skill_names()))
    return 0


def _cmd_mcp_serve(args):
    from pymeteo import mcp_server

    return mcp_server.main()


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="pymeteo",
        description="pymeteo CLI for the meteo-expert skill, MCP server, and diagnostics.",
    )
    sub = parser.add_subparsers(dest="command")

    skill_p = sub.add_parser(
        "install-skill",
        help="Copy the meteo-expert skill into a Claude Code or Codex skill directory.",
    )
    skill_p.add_argument(
        "--target",
        type=Path,
        default=None,
        help="Custom target directory (parent of meteo-expert/). Overrides --project.",
    )
    skill_p.add_argument(
        "--project",
        action="store_true",
        help="Install into ./skills/meteo-expert in the current working directory.",
    )
    skill_p.add_argument(
        "--skill",
        default="meteo-expert",
        help="Packaged skill name (default: meteo-expert).",
    )
    skill_p.add_argument(
        "--list",
        action="store_true",
        dest="list_skills",
        help="List packaged skills and exit.",
    )
    skill_p.set_defaults(func=_cmd_install_skill)

    info_p = sub.add_parser("info", help="Print version, function counts, and agent-ready paths.")
    info_p.set_defaults(func=_cmd_info)

    mcp_p = sub.add_parser("mcp", help="MCP server commands.")
    mcp_sub = mcp_p.add_subparsers(dest="mcp_command")
    serve_p = mcp_sub.add_parser("serve", help="Run the pymeteo MCP server.")
    serve_p.set_defaults(func=_cmd_mcp_serve)

    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        return 0
    if args.command == "mcp" and not getattr(args, "mcp_command", None):
        mcp_p.print_help()
        return 0
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
