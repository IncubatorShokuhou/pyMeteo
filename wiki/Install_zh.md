# 安装

Python 3.7+。运行时只依赖 NumPy。

```bash
pip install pymeteo-kit
```

```python
import pymeteo as pm
```

PyPI 上的名字是 `pymeteo-kit`，导入名是 `pymeteo`。

从仓库安装：

```bash
pip install -e ".[dev]"
pytest
```

给 agent 用的可选命令：`pymeteo install-skill`，`pip install "pymeteo-kit[mcp]"`，`pymeteo mcp serve`。见 [Agent](Agent_zh.md)。

[English](Install.md)
