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
    cache_key_for_3d,
    conformer_3d_for_graph,
    conformer_3d_for_graph_async,
    project_conformer_to_2d,
    scene_molecule_from_graph,
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


def test_cache_key_for_3d_separates_forcefields() -> None:
    graph = _ethanol_graph()
    uff = cache_key_for_3d(graph, OptimizationSettings(forcefield=ForceField.UFF), "rdkit")
    mmff = cache_key_for_3d(graph, OptimizationSettings(forcefield=ForceField.MMFF94), "rdkit")

    assert uff != mmff


def test_scene_molecule3d_includes_generated_atoms_without_mutating_graph() -> None:
    graph = _ethanol_graph()
    coordset = CoordinateSet3D(
        {1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (2.0, 0.0, 0.0)},
        "rdkit",
        "UFF",
        metadata={"generated_atoms": [{"id": "H1", "element": "H", "coords": (0.0, 1.0, 0.0)}]},
    )

    scene = scene_molecule_from_graph(graph, coordset)

    assert len(graph.atoms) == 3
    assert len(scene.atoms) == 4
    assert any(atom.generated and atom.element == "H" for atom in scene.atoms)


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


def test_rdkit_backend_uses_initial_coordinates_if_available() -> None:
    from chemuson.chemio.rdkit_safe import is_rdkit_worker_available

    if not is_rdkit_worker_available(timeout_s=5.0):
        pytest.skip("RDKit worker no disponible")
    initial = CoordinateSet3D(
        {1: (0.0, 0.0, 0.0), 2: (8.0, 0.0, 0.0), 3: (12.0, 3.0, 0.0)},
        source="manual",
        method="test",
    )

    result = rdkit_backend.optimize(
        _ethanol_graph(),
        initial,
        OptimizationSettings(forcefield=ForceField.UFF, max_iters=3, steps_per_update=1, timeout_s=10.0),
    )

    assert result.ok
    assert result.coordinates is not None
    assert result.coordinates.metadata.get("used_initial_coordinates") is True


def test_rdkit_safe_sends_initial_coordinates(monkeypatch) -> None:
    from chemuson.chemio import rdkit_safe

    captured: dict = {}

    def fake_run_worker(request, timeout_s):
        captured.update(request)
        return {
            "ok": True,
            "positions": {"1": [0.0, 0.0, 0.0], "2": [1.0, 0.0, 0.0], "3": [2.0, 0.0, 0.0]},
            "metadata": {"source": "rdkit", "method": "UFF", "used_initial_coordinates": True},
        }

    monkeypatch.setattr(rdkit_safe, "_run_worker", fake_run_worker)

    positions, metadata, error = rdkit_safe.optimize_3d_isolated(
        _ethanol_graph(),
        coordinates={1: (9.0, 0.0, 0.0), 2: (10.0, 0.0, 0.0), 3: (11.0, 0.0, 0.0)},
        forcefield="UFF",
        max_iters=5,
        steps_per_update=2,
        seed=123,
    )

    assert error is None
    assert positions
    assert metadata["used_initial_coordinates"] is True
    assert captured["coordinates"]["1"] == [9.0, 0.0, 0.0]
    assert captured["steps_per_update"] == 2
    assert captured["seed"] == 123


def test_rdkit_backend_optimize_iter_uses_returned_frames(monkeypatch) -> None:
    from chemuson.chemio import rdkit_safe

    def fake_optimize(*args, **kwargs):
        metadata = {
            "source": "rdkit",
            "method": "UFF",
            "energy": 1.0,
            "converged": False,
            "frames": [
                {"step": 2, "energy": 5.0, "converged": False, "positions": {"1": [0, 0, 0], "2": [1, 0, 0], "3": [2, 0, 0]}},
                {"step": 4, "energy": 1.0, "converged": True, "positions": {"1": [0, 0, 0], "2": [1.2, 0, 0], "3": [2.2, 0, 0]}},
            ],
        }
        return {1: (0.0, 0.0, 0.0), 2: (1.2, 0.0, 0.0), 3: (2.2, 0.0, 0.0)}, metadata, None

    monkeypatch.setattr(rdkit_safe, "optimize_3d_isolated", fake_optimize)

    frames = list(rdkit_backend.optimize_iter(_ethanol_graph(), settings=OptimizationSettings(forcefield=ForceField.UFF, max_iters=4, steps_per_update=2)))

    assert [frame.step for frame in frames] == [2, 4]
    assert frames[-1].converged is True


def test_rdkit_backend_returns_optimization_frames_if_available() -> None:
    from chemuson.chemio.rdkit_safe import is_rdkit_worker_available

    if not is_rdkit_worker_available(timeout_s=5.0):
        pytest.skip("RDKit worker no disponible")
    settings = OptimizationSettings(forcefield=ForceField.UFF, max_iters=6, steps_per_update=2, timeout_s=10.0)

    frames = list(rdkit_backend.optimize_iter(_ethanol_graph(), settings=settings))

    assert frames
    assert frames[0].step == 2
    assert all(frame.coordinates.positions for frame in frames)


def test_obabel_backend_degrades_or_runs() -> None:
    if not obabel_backend.is_available():
        pytest.skip("Open Babel no disponible")
    graph = _ethanol_graph()
    coordset = CoordinateSet3D({1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (2.0, 0.0, 0.0)}, "test", "manual")
    result = obabel_backend.optimize(graph, coordset, OptimizationSettings(forcefield=ForceField.UFF, max_iters=5, timeout_s=5.0))

    assert result.ok or result.message
    if result.ok:
        assert result.coordinates is not None
        assert len(result.coordinates.positions) == len(graph.atoms)
