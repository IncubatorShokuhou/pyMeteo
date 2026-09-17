# Install

Python 3.7+. Runtime dependency is NumPy.

```bash
pip install pymeteo-kit
```

```python
import pymeteo as pm
```

The PyPI name is `pymeteo-kit`; the import is `pymeteo`.

From a git checkout:

```bash
pip install -e ".[dev]"
pytest
```

Optional agent extras: `pymeteo install-skill`, `pip install "pymeteo-kit[mcp]"`, `pymeteo mcp serve`. See [Agent](Agent.md).

[中文](Install_zh.md)
