"""包导入面与破坏性变更。"""

import ast
import compileall
import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

import pymeteo

try:
    from importlib.metadata import metadata, version
except ImportError:  # Python < 3.8
    from importlib_metadata import metadata, version

_REPO_ROOT = Path(__file__).resolve().parents[1]


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
    assert pymeteo.__version__ == "2.3.0"


def test_pypi_distribution_name() -> None:
    meta = metadata("pymeteo-kit")
    assert meta["Name"] == "pymeteo-kit"
    assert version("pymeteo-kit") == "2.3.0"
    assert meta["Requires-Python"] == ">=3.6"


def test_packaging_declares_python_36() -> None:
    text = (_REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'requires-python = ">=3.6"' in text
    assert "Programming Language :: Python :: 3.6" in text
    assert "numpy>=1.19,<1.20" in text
    assert "numpy>=2.1" in text
    assert 'core-metadata-version = "2.1"' in text
    assert 'requires = ["hatchling>=1.18"]' in text


def test_hatchling_wheel_maps_src_pymeteo() -> None:
    text = (_REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'name = "pymeteo-kit"' in text
    assert "[tool.hatch.build.targets.wheel]" in text
    assert 'packages = ["src/pymeteo"]' in text


def test_package_compiles() -> None:
    src = str(_REPO_ROOT / "src" / "pymeteo")
    assert compileall.compile_dir(src, quiet=1, force=True)


_PEP585_BASES = frozenset({"list", "dict", "tuple", "set", "frozenset", "type"})


def _assert_annotation_is_py36(node: ast.AST, filename: str, lineno: int) -> None:
    for child in ast.walk(node):
        if isinstance(child, ast.BinOp) and isinstance(child.op, ast.BitOr):
            raise AssertionError(
                f"{filename}:{lineno}: PEP 604 `X | Y` in a type annotation is invalid at runtime on 3.6"
            )
        if isinstance(child, ast.Subscript) and isinstance(child.value, ast.Name):
            if child.value.id in _PEP585_BASES:
                raise AssertionError(
                    f"{filename}:{lineno}: PEP 585 `{child.value.id}[...]` is not subscriptable on 3.6"
                )


_NamedExpr = getattr(ast, "NamedExpr", None)


def test_sources_avoid_syntax_newer_than_36() -> None:
    src_root = _REPO_ROOT / "src" / "pymeteo"
    for path in sorted(src_root.glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        filename = str(path.relative_to(_REPO_ROOT))
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module == "__future__":
                names = {alias.name for alias in node.names}
                if "annotations" in names:
                    raise AssertionError(
                        f"{filename}:{node.lineno}: `from __future__ import annotations` is 3.7+"
                    )
            if _NamedExpr is not None and isinstance(node, _NamedExpr):
                raise AssertionError(f"{filename}:{node.lineno}: walrus `:=` is 3.8+")
            if type(node).__name__ == "Match":
                raise AssertionError(f"{filename}:{node.lineno}: match/case is 3.10+")
            if isinstance(node, ast.FunctionDef):
                args = list(node.args.args) + list(node.args.kwonlyargs)
                args.extend(getattr(node.args, "posonlyargs", []))
                if getattr(node.args, "vararg", None) is not None:
                    args.append(node.args.vararg)
                if getattr(node.args, "kwarg", None) is not None:
                    args.append(node.args.kwarg)
                for arg in args:
                    if arg.annotation is not None:
                        _assert_annotation_is_py36(arg.annotation, filename, arg.lineno)
                if node.returns is not None:
                    _assert_annotation_is_py36(node.returns, filename, node.lineno)
                if getattr(node.args, "posonlyargs", None):
                    raise AssertionError(
                        f"{filename}:{node.lineno}: positional-only `/` is 3.8+"
                    )
            if isinstance(node, ast.AnnAssign) and node.annotation is not None:
                _assert_annotation_is_py36(node.annotation, filename, node.lineno)



@pytest.mark.skipif(
    sys.version_info < (3, 8),
    reason="hatchling requires Python 3.8+; 3.6/3.7 install a prebuilt wheel",
)
def test_built_wheel_ships_pymeteo_import_package(tmp_path: Path) -> None:
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
    assert len(wheels) == 1, wheels
    with zipfile.ZipFile(wheels[0]) as zf:
        names = zf.namelist()
        meta = zf.read(next(n for n in names if n.endswith(".dist-info/METADATA"))).decode()
    assert "pymeteo/__init__.py" in names
    assert not any(n.startswith("pymeteo_kit/") for n in names)
    assert "Name: pymeteo-kit" in meta
    assert "Requires-Python: >=3.6" in meta
    assert "Metadata-Version: 2.1" in meta
    assert f"Version: {pymeteo.__version__}" in meta


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
