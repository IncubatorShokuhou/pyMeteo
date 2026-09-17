"""Packaged agent skills (Markdown) and the installer that copies them.

The skill lives at ``pymeteo/skills/meteo_expert/`` inside the wheel
(underscore: Python import path). Claude Code and Codex look for the
hyphenated directory ``meteo-expert`` from the SKILL.md frontmatter.
"""

import os
import shutil
import sys
from pathlib import Path

_INSTALL_DIRNAME_MAP = {
    "meteo_expert": "meteo-expert",
}

__all__ = ["get_skill_path", "install", "packaged_skill_names"]


def _normalize_to_package_name(skill):
    reverse = {value: key for key, value in _INSTALL_DIRNAME_MAP.items()}
    return reverse.get(skill, skill.replace("-", "_"))


def _install_dirname(skill_pkg):
    return _INSTALL_DIRNAME_MAP.get(skill_pkg, skill_pkg.replace("_", "-"))


def packaged_skill_names():
    """Return hyphenated skill identifiers shipped in this package."""

    return [_install_dirname(name) for name in sorted(_INSTALL_DIRNAME_MAP)]


def get_skill_path(skill_name="meteo-expert"):
    """Return the directory that contains SKILL.md inside the installed package."""

    pkg_name = _normalize_to_package_name(skill_name)
    path = Path(__file__).resolve().parent / pkg_name
    return path


def _ignored_name(name):
    return name in {"__pycache__", "__init__.py"} or name.endswith(".pyc")


def _copytree_overwrite(src, dst):
    """Copy a directory tree, overwriting files. Python 3.7 compatible."""

    if not os.path.isdir(dst):
        os.makedirs(dst)
    for name in os.listdir(src):
        if _ignored_name(name):
            continue
        source = os.path.join(src, name)
        target = os.path.join(dst, name)
        if os.path.isdir(source):
            _copytree_overwrite(source, target)
        else:
            parent = os.path.dirname(target)
            if parent and not os.path.isdir(parent):
                os.makedirs(parent)
            shutil.copy2(source, target)


def install(target_dir, skill_name="meteo-expert"):
    """Copy a packaged skill (SKILL.md + references/) into ``target_dir``."""

    pkg_name = _normalize_to_package_name(skill_name)
    source_dir = get_skill_path(pkg_name)
    source_skill = source_dir / "SKILL.md"
    if not source_skill.is_file():
        raise FileNotFoundError(
            f"Packaged skill not found: {source_skill}. Reinstalling pymeteo-kit may fix this."
        )
    dest_dir = Path(target_dir).expanduser().resolve() / _install_dirname(pkg_name)
    _copytree_overwrite(str(source_dir), str(dest_dir))
    return dest_dir / "SKILL.md"


def run_install(target=None, project=False, skill="meteo-expert", list_skills=False):
    """Shared install path for ``pymeteo install-skill``. Returns an exit code."""

    if list_skills:
        print("Available skills:")
        for pkg_name, install_name in _INSTALL_DIRNAME_MAP.items():
            source = get_skill_path(pkg_name) / "SKILL.md"
            marker = "ok" if source.is_file() else "MISSING"
            print(f"  {install_name} ({marker})")
        return 0

    if target is not None:
        resolved_target = Path(target)
        mode = "custom"
    elif project:
        resolved_target = Path.cwd() / "skills"
        mode = "project"
    else:
        resolved_target = Path.home() / ".claude" / "skills"
        mode = "user-global"

    try:
        dest = install(resolved_target, skill_name=skill)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    canonical = _install_dirname(_normalize_to_package_name(skill))
    print(f"Installed {canonical} skill to: {dest}")
    if mode == "user-global":
        print("Claude Code picks this up from ~/.claude/skills/ on the next session.")
    elif mode == "project":
        print(
            f"Project-local skill is at ./skills/{canonical}/ "
            "(Codex / Claude Code in this directory)."
        )
    else:
        print("Restart the agent session if it caches skills at startup.")
    return 0
