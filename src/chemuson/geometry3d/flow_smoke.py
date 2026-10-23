"""Verificación empaquetada del flujo normal de generación de conformeros 3D.

Ejecuta la misma ruta de producción que la GUI (``conformer_3d_for_graph`` ->
``rdkit_backend.generate_conformer`` -> ``rdkit_safe.conformer_3d_stream_isolated``)
sin inyectar escenas ni fixtures: construye un MolGraph 2D y exige geometrías 3D
reales publicadas por el servicio.
"""

from __future__ import annotations

import sys
import time
from typing import Any

from chemuson.core.model import MolGraph
from chemuson.geometry3d.service import conformer_3d_for_graph


def _ethanol_graph() -> MolGraph:
    graph = MolGraph()
    carbon_one = graph.add_atom("C", 0.0, 0.0).id
    carbon_two = graph.add_atom("C", 1.5, 0.0).id
    oxygen = graph.add_atom("O", 3.0, 0.0).id
    graph.add_bond(carbon_one, carbon_two, order=1)
    graph.add_bond(carbon_two, oxygen, order=1)
    return graph


def run_smoke(*, timeout_s: float = 30.0) -> dict[str, Any]:
    """Genera un conformero por el flujo normal y reporta el diagnóstico."""
    started = time.monotonic()
    events: list[dict[str, Any]] = []

    def on_progress(event: dict[str, Any]) -> None:
        events.append(dict(event))

    result = conformer_3d_for_graph(
        _ethanol_graph(), timeout_s=timeout_s, force=True, on_progress=on_progress
    )
    positions = {
        atom_id: [round(float(value), 6) for value in coordinates]
        for atom_id, coordinates in sorted(result.atom_positions.items())
    }
    ok = bool(positions) and len(positions) == 3 and not result.message
    return {
        "ok": ok,
        "entry": "geometry3d.service.conformer_3d_for_graph",
        "frozen": bool(getattr(sys, "frozen", False)),
        "atom_position_count": len(positions),
        "positions": {str(key): value for key, value in positions.items()},
        "source": result.source,
        "method": result.method,
        "energy": result.energy,
        "message": result.message,
        "progress_stages": [str(event.get("stage", "")) for event in events],
        "elapsed_s": round(time.monotonic() - started, 3),
        "diagnostics": dict(result.metadata or {}),
    }
