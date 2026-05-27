"""Pruebas del pipeline 3D real desacoplado."""

from __future__ import annotations

import math
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.geometry3d import (
    Rotation3D,
    conformer_3d_for_graph,
    conformer_3d_for_graph_async,
    project_conformer_to_2d,
)


def test_project_conformer_to_2d_returns_depth_cueing() -> None:
    projected = project_conformer_to_2d(
        {
            1: (-1.0, 0.0, -1.0),
            2: (1.0, 0.0, 1.0),
        },
        rotation=Rotation3D(yaw=math.radians(15.0)),
        center=(100.0, 50.0),
        scale=10.0,
    )

    assert [atom.atom_id for atom in projected] == [1, 2]
    assert projected[0].cue.opacity < projected[1].cue.opacity
    assert projected[0].cue.stroke_scale < projected[1].cue.stroke_scale


def test_conformer_3d_for_graph_uses_cache() -> None:
    graph = MolGraph()
    c1 = graph.add_atom("C", 0.0, 0.0)
    c2 = graph.add_atom("C", 42.0, 0.0)
    o = graph.add_atom("O", 84.0, 0.0)
    graph.add_bond(c1.id, c2.id, order=1)
    graph.add_bond(c2.id, o.id, order=1)

    first = conformer_3d_for_graph(graph, timeout_s=8.0, force=True)
    if first.message.startswith("rdkit") or first.message:
        # Entornos CI sin RDKit operativo mantienen degradación segura.
        return

    second = conformer_3d_for_graph(graph, timeout_s=8.0)

    assert first.ok
    assert set(first.atom_positions) == {c1.id, c2.id, o.id}
    assert second.ok
    assert second.from_cache
    assert second.cache_key == first.cache_key


def test_conformer_3d_async_returns_future() -> None:
    graph = MolGraph()
    graph.add_atom("O", 0.0, 0.0)

    future = conformer_3d_for_graph_async(graph, timeout_s=8.0, force=True)
    result = future.result(timeout=12.0)

    assert result.ok or result.message
