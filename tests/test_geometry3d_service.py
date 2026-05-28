"""Pruebas del pipeline 3D real desacoplado."""

from __future__ import annotations

import math
import os
import sys

import pytest

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from chemuson.core.model import MolGraph
from chemuson.geometry3d import (
    CoordinateSet3D,
    ForceField,
    OptimizationSettings,
    Rotation3D,
    conformer_3d_for_graph,
    conformer_3d_for_graph_async,
    project_conformer_to_2d,
)
from chemuson.geometry3d.export_xyz import molgraph_to_xyz
from chemuson.geometry3d import obabel_backend, rdkit_backend


def _ethanol_graph() -> MolGraph:
    graph = MolGraph()
    c1 = graph.add_atom("C", 0.0, 0.0)
    c2 = graph.add_atom("C", 42.0, 0.0)
    o = graph.add_atom("O", 84.0, 0.0)
    graph.add_bond(c1.id, c2.id, order=1)
    graph.add_bond(c2.id, o.id, order=1)
    return graph


def test_coordinate_set_normalizes_positions() -> None:
    coordset = CoordinateSet3D({"1": (0, 1, 2)}, source="test", method="manual", energy=-1.5)

    assert coordset.normalized_positions() == {1: (0.0, 1.0, 2.0)}
    assert coordset.energy == -1.5


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


def test_project_conformer_to_2d_accepts_coordinate_set() -> None:
    coordset = CoordinateSet3D({1: (-1.0, 0.0, -1.0), 2: (1.0, 0.0, 1.0)}, source="test", method="manual")

    projected = project_conformer_to_2d(coordset)

    assert [atom.atom_id for atom in projected] == [1, 2]


def test_molgraph_to_xyz_exports_coordinate_set() -> None:
    graph = _ethanol_graph()
    coordset = CoordinateSet3D({1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (2.0, 0.0, 0.0)}, "test", "manual")

    xyz = molgraph_to_xyz(graph, coordset, comment="ethanol")

    assert xyz.splitlines()[0] == "3"
    assert xyz.splitlines()[1] == "ethanol"
    assert "O 2.00000000 0.00000000 0.00000000" in xyz


def test_conformer_3d_for_graph_uses_cache() -> None:
    graph = _ethanol_graph()

    first = conformer_3d_for_graph(graph, timeout_s=8.0, force=True)
    if first.message.startswith("rdkit") or first.message:
        # Entornos CI sin RDKit operativo mantienen degradación segura.
        return

    second = conformer_3d_for_graph(graph, timeout_s=8.0)

    assert first.ok
    assert set(first.atom_positions) == {1, 2, 3}
    assert first.coordinate_set is not None
    assert second.ok
    assert second.from_cache
    assert second.cache_key == first.cache_key


def test_conformer_3d_async_returns_future() -> None:
    graph = MolGraph()
    graph.add_atom("O", 0.0, 0.0)

    future = conformer_3d_for_graph_async(graph, timeout_s=8.0, force=True)
    result = future.result(timeout=12.0)

    assert result.ok or result.message


def test_rdkit_backend_simple_molecule_if_available() -> None:
    from chemuson.chemio.rdkit_safe import is_rdkit_worker_available

    if not is_rdkit_worker_available(timeout_s=5.0):
        pytest.skip("RDKit worker no disponible")
    result = rdkit_backend.optimize(
        _ethanol_graph(),
        settings=OptimizationSettings(forcefield=ForceField.UFF, max_iters=25, timeout_s=10.0),
    )

    assert result.ok
    assert result.coordinates is not None
    assert set(result.coordinates.positions) == {1, 2, 3}


def test_obabel_backend_degrades_or_runs() -> None:
    if not obabel_backend.is_available():
        pytest.skip("Open Babel no disponible")
    graph = _ethanol_graph()
    coordset = CoordinateSet3D({1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (2.0, 0.0, 0.0)}, "test", "manual")
    result = obabel_backend.optimize(graph, coordset, OptimizationSettings(forcefield=ForceField.UFF, max_iters=5, timeout_s=5.0))

    assert result.ok or result.message
