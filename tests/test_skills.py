"""meteo-expert skill 打包与公开 API 一致性。"""

import re
import zipfile
from pathlib import Path

import pytest

import pymeteo
from pymeteo.ncl import __all__ as NCL_NAMES
from pymeteo.skills import get_skill_path, install, packaged_skill_names

_REPO_ROOT = Path(__file__).resolve().parents[1]

_SKILL_VOCAB = {
    "meteo-expert",
    "meteo_expert",
    "pymeteo",
    "pymeteo-kit",
    "temperature_unit",
    "output_temperature_unit",
    "humidity_unit",
    "output_humidity_unit",
    "pressure_unit",
    "output_pressure_unit",
    "mixing_ratio_unit",
    "speed_unit",
    "output_speed_unit",
    "output_distance_unit",
    "output_omega_unit",
    "angle_unit",
    "latitude_unit",
    "height_unit",
    "from_quantity",
    "to_quantity",
    "humVarType",
    "wqType",
    "iounit",
    "iswit",
    "dewpoint_850",
    "relative_humidity_850",
    "temperature_850",
    "temperature_700",
    "temperature_500",
    "dewpoint_700",
    "dewpoint_500",
    "lapse_rate",
    "station_pressure",
    "station_height",
    "mean_temperature",
    "pressure_bottom",
    "pressure_top",
    "mcp_server",
    "list_functions",
    "explain_function",
    "ncl_lookup",
    "unit_help",
    "run_calc",
    "MeteoEngine",
    "install_skill",
    "__all__",
    "method",
}


def test_packaged_skill_tree_exists() -> None:
    names = packaged_skill_names()
    assert "meteo-expert" in names
    root = get_skill_path("meteo-expert")
    assert (root / "SKILL.md").is_file()
    text = (root / "SKILL.md").read_text(encoding="utf-8")
    assert text.startswith("---")
    assert "name: meteo-expert" in text
    refs = root / "references"
    for name in ("workflow.md", "units.md", "pitfalls.md", "ncl.md"):
        assert (refs / name).is_file(), name


def test_install_skill_copies_tree(tmp_path: Path) -> None:
    dest = install(tmp_path, skill_name="meteo-expert")
    assert dest.is_file()
    assert dest.parent.name == "meteo-expert"
    assert (dest.parent / "references" / "workflow.md").is_file()
    assert not (dest.parent / "__init__.py").exists()


def test_skill_backticks_match_public_api() -> None:
    root = get_skill_path("meteo-expert")
    texts = [path.read_text(encoding="utf-8") for path in root.rglob("*.md")]
    names = set()
    for text in texts:
        names.update(re.findall(r"`([A-Za-z_][A-Za-z0-9_]*)`", text))
    public = set(pymeteo.__all__) | set(NCL_NAMES)
    public.update(
        {
            "thermo",
            "indices",
            "wind",
            "geo",
            "dynamics",
            "comfort",
            "ncl",
            "units",
            "UnitError",
        }
    )
    unknown = []
    for name in sorted(names):
        if name in public or name in _SKILL_VOCAB:
            continue
        if "_" not in name:
            continue
        unknown.append(name)
    assert unknown == []


@pytest.mark.skipif(
    __import__("sys").version_info < (3, 8),
    reason="hatchling requires Python 3.8+; 3.7 installs a prebuilt wheel",
)
def test_wheel_includes_skill_markdown(tmp_path: Path) -> None:
    import subprocess
    import sys

    subprocess.check_call(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            str(_REPO_ROOT),
            "-w",
            str(tmp_path),
            "--no-deps",
        ],
    )
    wheels = list(tmp_path.glob("pymeteo_kit-*.whl"))
    assert len(wheels) == 1
    with zipfile.ZipFile(wheels[0]) as zf:
        names = zf.namelist()
    assert "pymeteo/skills/meteo_expert/SKILL.md" in names
    assert "pymeteo/skills/meteo_expert/references/ncl.md" in names
