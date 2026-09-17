# Agent

[中文](Agent_zh.md)

pymeteo ships a `meteo-expert` skill and an optional MCP server so coding agents can pick a function, check unit strings, and map NCL names. The classic `import pymeteo as pm` calls do not change.

```bash
pip install pymeteo-kit
pymeteo install-skill              # Claude Code → ~/.claude/skills/meteo-expert
pymeteo install-skill --project    # Codex → ./skills/meteo-expert
pip install "pymeteo-kit[mcp]"
pymeteo mcp serve                  # python -m pymeteo.mcp_server
pymeteo info
```

From Python, without MCP:

```python
from pymeteo.engine import MeteoEngine

engine = MeteoEngine()
engine.recommend("dewpoint from RH")
engine.explain("showalter_index")
engine.ncl_lookup("dewtemp_trh")
```

Core install stays Python 3.7+ and NumPy-only. The `mcp` extra follows whatever Python the `mcp` package itself requires.
