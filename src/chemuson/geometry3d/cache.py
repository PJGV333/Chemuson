from __future__ import annotations

"""Claves de caché para conformaciones/optimizaciones 3D."""

import hashlib

from chemuson.core.model import MolGraph
from chemuson.geometry3d.model import OptimizationSettings
from chemuson.geometry3d.service import chemical_graph_hash


def cache_key_for_3d(graph: MolGraph, settings: OptimizationSettings, backend: str) -> str:
    """Clave estable que separa backend, force field, iteraciones y semilla."""
    payload = "|".join(
        (
            chemical_graph_hash(graph),
            str(backend),
            str(settings.forcefield.value),
            str(int(settings.max_iters)),
            str(int(settings.steps_per_update)),
            str(int(settings.seed)),
        )
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()
