# Agent

[English](Agent.md)

pymeteo 带了一个 `meteo-expert` skill，以及可选的 MCP 服务，方便编程 agent 选函数、核对单位字符串、对照 NCL 名字。原来的 `import pymeteo as pm` 写法不变。

```bash
pip install pymeteo-kit
pymeteo install-skill              # Claude Code → ~/.claude/skills/meteo-expert
pymeteo install-skill --project    # Codex → ./skills/meteo-expert
pip install "pymeteo-kit[mcp]"
pymeteo mcp serve                  # python -m pymeteo.mcp_server
pymeteo info
```

不用 MCP 时，在 Python 里：

```python
from pymeteo.engine import MeteoEngine

engine = MeteoEngine()
engine.recommend("dewpoint from RH")
engine.explain("showalter_index")
engine.ncl_lookup("dewtemp_trh")
```

核心安装仍是 Python 3.7+、只依赖 NumPy。`mcp` extra 跟着 `mcp` 包自己的 Python 要求走。
