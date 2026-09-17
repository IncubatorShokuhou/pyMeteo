"""Lightweight planning/knowledge API for agents.

``MeteoEngine`` mirrors the MCP tools without requiring the ``mcp`` extra.
It does not keep a session or run a multi-step sounding workflow.
"""

from pymeteo import catalog


class MeteoEngine:
    """Pure-Python catalog over the classic pymeteo functions."""

    def list_functions(self, module=None):
        """List public functions. ``module`` is thermo/indices/wind/geo/dynamics/comfort."""

        return catalog.list_functions(module=module)

    def explain(self, name):
        """Explain a public function or NCL shim by name."""

        return catalog.explain_function(name)

    def recommend(self, query="", need=None, top_k=5):
        """Rank functions for a meteorological question."""

        return catalog.recommend(query=query, need=need, top_k=top_k)

    def ncl_lookup(self, name):
        """Map an NCL builtin name to a pymeteo callable plus unit notes."""

        return catalog.ncl_lookup(name)

    def unit_help(self, kind=None):
        """Return unit aliases and defaults. ``kind`` is optional."""

        return catalog.unit_help(kind)

    def run_calc(self, name, args=None, kwargs=None):
        """Call a whitelisted function with JSON-serializable scalars/lists."""

        return catalog.run_calc(name, args=args, kwargs=kwargs)
