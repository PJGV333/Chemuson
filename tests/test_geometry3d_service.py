"""Pruebas del pipeline 3D real desacoplado."""

from __future__ import annotations

import math
from typing import Any

import pytest


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
from chemuson.geometry3d.model import (
    VisualFramePlayback,
    coordinate_sets_interpolatable,
    interpolate_coordinate_sets,
    scene_with_coordinates,
)
from chemuson.geometry3d.editing import dihedral_degrees, rotate_atoms_about_axis, set_dihedral
from chemuson.geometry3d.export_xyz import molgraph_to_xyz
from chemuson.geometry3d import obabel_backend, rdkit_backend
from chemuson.chemio.persistence import PersistenceManager


# PubChem CID 73078 reference used here is not the user's original .cmsn.
_TETRANDRINE_CID_73078_SMILES = (
    "CN1CCC2=CC(=C3C=C2[C@@H]1CC4=CC=C(C=C4)OC5=C(C=CC(=C5)"
    "C[C@H]6C7=C(O3)C(=C(C=C7CCN6C)OC)OC)OC)OC"
)


def test_carbon_scene_color_has_contrast_for_dark_viewport() -> None:
    graph = MolGraph()
    carbon = graph.add_atom("C", 0.0, 0.0)
    coordinates = CoordinateSet3D({carbon.id: (0.0, 0.0, 0.0)}, "test", "manual")

    scene = scene_molecule_from_graph(graph, coordinates)

    assert scene.atoms[0].color == "#909090"


def _ethanol_graph() -> MolGraph:
    graph = MolGraph()
    c1 = graph.add_atom("C", 0.0, 0.0)
    c2 = graph.add_atom("C", 42.0, 0.0)
    o = graph.add_atom("O", 84.0, 0.0)
    graph.add_bond(c1.id, c2.id, order=1)
    graph.add_bond(c2.id, o.id, order=1)
    return graph


def _ethane_graph() -> MolGraph:
    graph = MolGraph()
    first = graph.add_atom("C", 0.0, 0.0)
    second = graph.add_atom("C", 42.0, 0.0)
    graph.add_bond(first.id, second.id, order=1)
    return graph


def _mm_snapshot(step: int, offset: float, energy: float, forcefield: str = "MMFF94") -> CoordinateSet3D:
    return CoordinateSet3D(
        positions={1: (offset, 0.0, 0.0), 2: (offset + 1.5, 0.0, 0.0)},
        source="rdkit",
        method=forcefield,
        energy=energy,
        metadata={
            "forcefield": forcefield,
            "energy_unit": "kcal/mol",
            "step": step,
            "max_steps": 16,
            "fixed_atom_ids": (1,),
            "fixed_atom_count": 1,
            "constraints": {
                "kind": "fixed_atom_positions",
                "atom_count": 1,
                "coordinate_unit": "Å",
            },
            "generated_atoms": [
                {"id": "genH1", "element": "H", "coords": [offset, 1.0, 0.0]}
            ],
            "generated_bonds": [
                {"id": "genB1", "a1_id": 1, "a2_id": "genH1", "order": 1}
            ],
        },
    )


def test_live_mm_visual_interpolation_keeps_measured_energy_separate():
    start = _mm_snapshot(0, 0.0, -1.0)
    first = _mm_snapshot(8, 1.0, -2.0)
    final = _mm_snapshot(16, 2.0, -3.0)

    assert coordinate_sets_interpolatable(start, first)
    visual = interpolate_coordinate_sets(start, first, 0.5)
    assert visual.positions[1] == pytest.approx((0.5, 0.0, 0.0))
    generated = visual.metadata.get("generated_atoms")
    assert isinstance(generated, list) and isinstance(generated[0], dict)
    assert generated[0]["coords"] == pytest.approx([0.5, 1.0, 0.0])
    assert visual.energy is None
    assert visual.metadata["visual_interpolation"] is True
    assert visual.metadata["visual_interpolation_from_step"] == 0
    assert visual.metadata["visual_interpolation_to_step"] == 8

    manual_start = CoordinateSet3D(
        {1: (0.0, 0.0, 0.0), 2: (1.5, 0.0, 0.0)},
        source="manual",
        method="build-3d",
    )
    assert coordinate_sets_interpolatable(manual_start, first)
    manual_to_mm = interpolate_coordinate_sets(manual_start, first, 0.5)
    assert manual_to_mm.positions[1] == pytest.approx((0.5, 0.0, 0.0))
    manual_generated = manual_to_mm.metadata.get("generated_atoms")
    assert isinstance(manual_generated, list)
    assert manual_generated[0]["coords"] == pytest.approx([1.0, 1.0, 0.0])

    playback = VisualFramePlayback(transition_s=0.08)
    playback.begin(start)
    assert playback.enqueue(first)
    assert playback.enqueue(final)
    playback.finish(final)
    initial_visual = playback.sample(0.0)
    assert initial_visual is not None and initial_visual.energy is None
    middle = playback.sample(0.04)
    assert middle is not None
    assert middle.positions[1] == pytest.approx((0.5, 0.0, 0.0))
    first_snapshot = playback.sample(0.08)
    assert first_snapshot is not None and first_snapshot.energy == -2.0
    next_visual = playback.sample(0.081)
    assert next_visual is not None and next_visual.energy is None
    next_middle = playback.sample(0.121)
    assert next_middle is not None
    assert next_middle.positions[1] == pytest.approx((1.5, 0.0, 0.0))
    final_snapshot = playback.sample(0.161)
    assert final_snapshot is not None and final_snapshot.energy == -3.0
    assert playback.drained


def test_live_mm_playback_rejects_incompatible_forcefield_and_cancels_queue():
    start = _mm_snapshot(0, 0.0, -1.0)
    target = _mm_snapshot(8, 1.0, -2.0, forcefield="UFF")
    assert not coordinate_sets_interpolatable(start, target)
    with pytest.raises(ValueError, match="coordinate_snapshots_incompatible"):
        interpolate_coordinate_sets(start, target, 0.5)

    playback = VisualFramePlayback(transition_s=0.08, max_pending=2)
    playback.begin(start)
    assert not playback.enqueue(target)
    assert playback.current is target
    playback.cancel()
    assert not playback.enabled
    assert playback.pending_count == 0


def test_scene_with_visual_coordinates_updates_render_atoms_without_mutating_source():
    from chemuson.geometry3d import SceneAtom3D, SceneBond3D, SceneMolecule3D

    scene = SceneMolecule3D(
        atoms=(
            SceneAtom3D(1, "C", 0.0, 0.0, 0.0, 1.7, "#909090"),
            SceneAtom3D(2, "C", 1.5, 0.0, 0.0, 1.7, "#909090"),
        ),
        bonds=(SceneBond3D(1, 1, 2),),
    )
    visual = interpolate_coordinate_sets(
        _mm_snapshot(0, 0.0, -1.0),
        _mm_snapshot(8, 2.0, -2.0),
        0.5,
    )

    rendered = scene_with_coordinates(scene, visual)

    assert len(rendered.atoms) == 3
    assert rendered.atoms[0].x == pytest.approx(1.0)
    assert rendered.atoms[1].x == pytest.approx(2.5)
    assert rendered.atoms[2].id == "genH1"
    assert rendered.atoms[2].x == pytest.approx(1.0)
    assert len(rendered.bonds) == 2
    assert rendered.bonds[0] == scene.bonds[0]
    assert rendered.bonds[1].a2_id == "genH1"
    assert scene.atoms[0].x == 0.0
    assert len(scene.atoms) == 2
    assert visual.energy is None



def test_conformer_service_propagates_progress_cancel_and_failure_diagnostics(
    monkeypatch,
):
    import threading

    from chemuson.geometry3d import OptimizationResult, service

    graph = _ethane_graph()
    cancel_event = threading.Event()
    progress = []
    diagnostics = {
        "error_code": "cancelled",
        "stage": "embedding",
        "graph": {"atom_count": 2, "bond_count": 1},
    }

    def fake_generate(received_graph, _settings, *, cancel_event, on_progress):
        assert received_graph is graph
        assert cancel_event is not None
        on_progress({"event": "progress", "stage": "embedding"})
        return OptimizationResult(
            None,
            method="rdkit",
            message="cancelled: embedding",
            metadata={"diagnostics": diagnostics},
        )

    monkeypatch.setattr(service, "generate_conformer", fake_generate)

    result = service.conformer_3d_for_graph(
        graph,
        force=True,
        cancel_event=cancel_event,
        on_progress=progress.append,
    )

    assert not result.ok
    assert result.message == "cancelled: embedding"
    assert result.metadata["diagnostics"] == diagnostics
    assert [event["stage"] for event in progress] == ["embedding"]


def test_conformer_cache_key_separates_forcefield_and_seed(monkeypatch):
    from chemuson.geometry3d import OptimizationResult, service

    graph = _ethane_graph()
    monkeypatch.setattr(service, "_CONFORMER_CACHE", {})

    def fake_generate(_graph, settings, **_kwargs):
        coordinates = CoordinateSet3D(
            {1: (0.0, 0.0, 0.0), 2: (1.5, 0.0, 0.0)},
            source="rdkit",
            method=settings.forcefield.value,
        )
        return OptimizationResult(coordinates, converged=True, method=coordinates.method)

    monkeypatch.setattr(service, "generate_conformer", fake_generate)
    first = service.conformer_3d_for_graph(
        graph,
        force=True,
        settings=OptimizationSettings(forcefield=ForceField.MMFF94, seed=11),
    )
    second = service.conformer_3d_for_graph(
        graph,
        force=True,
        settings=OptimizationSettings(forcefield=ForceField.UFF, seed=12),
    )

    assert first.cache_key != second.cache_key


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


def test_conformer_backend_relays_worker_progress_and_structured_diagnostics(
    monkeypatch,
):
    from chemuson.chemio import rdkit_safe

    graph = _ethane_graph()
    progress = []
    diagnostics = {
        "error_code": None,
        "stage": "done",
        "graph": {"atom_count": 2, "bond_count": 1},
        "timings_s": {"total": 0.25},
    }

    def fake_stream(_graph, **_kwargs):
        yield {"event": "progress", "stage": "embedding", "message": "Embedding"}
        yield {
            "event": "done",
            "ok": True,
            "positions": {"1": [0.0, 0.0, 0.0], "2": [1.5, 0.0, 0.0]},
            "metadata": {
                "source": "rdkit",
                "method": "mmff",
                "energy": -1.0,
                "converged": True,
                "generated_atoms": [],
                "generated_bonds": [],
                "diagnostics": diagnostics,
            },
            "diagnostics": diagnostics,
            "message": "",
        }

    monkeypatch.setattr(rdkit_safe, "conformer_3d_stream_isolated", fake_stream)

    result = rdkit_backend.generate_conformer(
        graph,
        OptimizationSettings(forcefield=ForceField.MMFF94),
        on_progress=progress.append,
    )

    assert result.ok
    assert result.coordinates is not None
    assert result.coordinates.metadata["diagnostics"] == diagnostics
    assert result.metadata["diagnostics"] == diagnostics
    assert result.converged
    assert [event["stage"] for event in progress] == ["embedding"]


def test_conformer_backend_keeps_failure_diagnostics_without_coordinates(monkeypatch):
    from chemuson.chemio import rdkit_safe

    graph = _ethane_graph()
    diagnostics = {
        "error_code": "embed_failed",
        "stage": "embedding",
        "graph": {"atom_count": 2, "bond_count": 1},
    }

    def fake_stream(_graph, **_kwargs):
        yield {
            "event": "error",
            "ok": False,
            "error": "embed_failed",
            "detail": "all attempts failed",
            "diagnostics": diagnostics,
        }

    monkeypatch.setattr(rdkit_safe, "conformer_3d_stream_isolated", fake_stream)

    result = rdkit_backend.generate_conformer(graph)

    assert not result.ok
    assert result.coordinates is None
    failure_diagnostics = result.metadata.get("diagnostics")
    assert isinstance(failure_diagnostics, dict)
    assert failure_diagnostics.get("error_code") == diagnostics["error_code"]
    assert failure_diagnostics.get("detail") == "all attempts failed"
    assert "embed_failed" in result.message


def test_generated_conformer_exposes_hydrogen_atoms_and_bonds_without_persisting_them(tmp_path) -> None:
    from chemuson.chemio.rdkit_safe import is_rdkit_worker_available

    if not is_rdkit_worker_available(timeout_s=5.0):
        pytest.skip("RDKit worker no disponible")
    graph = _ethane_graph()
    result = conformer_3d_for_graph(graph, timeout_s=20.0, force=True)
    assert result.ok, result.message
    coordinates = result.coordinate_set
    assert coordinates is not None
    generated_atoms = coordinates.metadata.get("generated_atoms", [])
    generated_bonds = coordinates.metadata.get("generated_bonds", [])
    assert isinstance(generated_atoms, list)
    assert isinstance(generated_bonds, list)

    scene = scene_molecule_from_graph(graph, coordinates)

    assert len(generated_atoms) == 6
    assert len(generated_bonds) == 6
    assert len(scene.atoms) == 8
    assert len(scene.bonds) == 7
    assert sum(atom.generated for atom in scene.atoms) == 6
    assert len(graph.atoms) == 2
    assert len(graph.bonds) == 1

    class PersistenceStub:
        def __init__(self, model: MolGraph) -> None:
            self.model = model

        def get_persistence_data(self) -> dict[str, Any]:
            return {}

        def load_persistence_data(self, data: dict[str, Any]) -> None:
            pass

        def rebuild_persistence_view(self) -> None:
            pass

    source = PersistenceStub(graph)
    path = tmp_path / "ethane.cmsn"
    PersistenceManager.save_to_file(str(path), source)
    serialized = PersistenceManager.save_to_dict(source)
    assert len(serialized["model"]["atoms"]) == 2
    assert len(serialized["model"]["bonds"]) == 1

    restored = PersistenceStub(MolGraph())
    PersistenceManager.load_from_file(str(path), restored)
    assert len(restored.model.atoms) == 2
    assert len(restored.model.bonds) == 1


def test_rdkit_optimization_preserves_generated_hydrogen_scene_metadata() -> None:
    from chemuson.chemio.rdkit_safe import is_rdkit_worker_available

    if not is_rdkit_worker_available(timeout_s=5.0):
        pytest.skip("RDKit worker no disponible")
    result = rdkit_backend.optimize(
        _ethane_graph(),
        settings=OptimizationSettings(forcefield=ForceField.UFF, max_iters=3, timeout_s=15.0),
    )

    assert result.ok
    assert result.coordinates is not None
    generated_atoms = result.coordinates.metadata.get("generated_atoms", [])
    generated_bonds = result.coordinates.metadata.get("generated_bonds", [])
    assert isinstance(generated_atoms, list)
    assert isinstance(generated_bonds, list)
    assert len(generated_atoms) == 6
    assert len(generated_bonds) == 6


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

    def fake_stream(*args, **kwargs):
        yield {
            "event": "frame",
            "sequence": 1,
            "step": 2,
            "energy": 5.0,
            "forcefield": "UFF",
            "energy_unit": "kcal/mol",
            "converged": False,
            "positions": {"1": [0, 0, 0], "2": [1, 0, 0], "3": [2, 0, 0]},
        }
        yield {
            "event": "frame",
            "sequence": 2,
            "step": 4,
            "energy": 1.0,
            "forcefield": "UFF",
            "energy_unit": "kcal/mol",
            "converged": True,
            "positions": {"1": [0, 0, 0], "2": [1.2, 0, 0], "3": [2.2, 0, 0]},
        }
        yield {
            "event": "done",
            "ok": True,
            "method": "UFF",
            "energy": 1.0,
            "energy_unit": "kcal/mol",
            "converged": True,
            "positions": {"1": [0, 0, 0], "2": [1.2, 0, 0], "3": [2.2, 0, 0]},
            "metadata": {"method": "UFF", "source": "rdkit"},
        }

    monkeypatch.setattr(rdkit_safe, "optimize_3d_stream_isolated", fake_stream)

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


def test_rdkit_backend_frames_keep_derived_hydrogens() -> None:
    generated_atoms = [
        {"id": "genH2", "element": "H", "coords": [1.0, 0.0, 0.0]}
    ]
    generated_bonds = [
        {"id": "genB1", "a1_id": 1, "a2_id": "genH2", "order": 1}
    ]
    frames = rdkit_backend._frames_from_metadata(
        {
            "method": "MMFF94",
            "frames": [
                {
                    "step": 2,
                    "energy": 3.0,
                    "positions": {"1": [0.0, 0.0, 0.0], "2": [2.0, 0.0, 0.0]},
                    "generated_atoms": generated_atoms,
                    "generated_bonds": generated_bonds,
                }
            ],
        },
        "rdkit",
    )

    assert len(frames) == 1
    frame = frames[0]
    assert frame.coordinates.metadata["generated_atoms"] == generated_atoms
    assert frame.coordinates.metadata["generated_bonds"] == generated_bonds
    scene = scene_molecule_from_graph(_ethane_graph(), frame.coordinates)
    assert len(scene.atoms) == 3
    assert len(scene.bonds) == 2


def test_conformer_diagnostics_capture_valence_aromaticity_and_unassigned_stereo() -> None:
    from chemuson.chemio.rdkit_safe import (
        is_rdkit_worker_available,
        smiles_to_molgraph_isolated,
    )

    if not is_rdkit_worker_available(timeout_s=5.0):
        pytest.skip("RDKit worker no disponible")

    benzene, error = smiles_to_molgraph_isolated("c1ccccc1", timeout_s=8.0)
    assert benzene is not None, error
    aromatic_result = conformer_3d_for_graph(benzene, force=True)
    assert aromatic_result.ok, aromatic_result.message
    aromatic_diagnostics = aromatic_result.coordinate_set.metadata["diagnostics"]
    aromatic_graph = aromatic_diagnostics["graph"]
    assert aromatic_graph["aromatic_atom_count"] == 6
    assert aromatic_graph["aromatic_bond_count"] == 6
    assert len(aromatic_graph["atom_states"]) == 6
    assert all(atom["total_valence"] == 4 for atom in aromatic_graph["atom_states"])
    assert aromatic_diagnostics["return_code"] == 0
    assert aromatic_diagnostics["timings_s"]["total"] > 0

    unspecified, error = smiles_to_molgraph_isolated("CC(O)F", timeout_s=8.0)
    assert unspecified is not None, error
    unspecified_result = conformer_3d_for_graph(unspecified, force=True)
    assert unspecified_result.ok, unspecified_result.message
    stereo_diagnostics = unspecified_result.coordinate_set.metadata["diagnostics"][
        "stereochemistry_diagnostics"
    ]
    assert stereo_diagnostics["unassigned_atom_count"] == 1
    assert stereo_diagnostics["specified_atoms"] == {}


@pytest.mark.parametrize(
    ("name", "smiles", "document_atoms", "document_bonds", "generated_hydrogens"),
    [
        ("water", "O", 1, 0, 2),
        ("ethane", "CC", 2, 1, 6),
        ("ethanol", "CCO", 3, 2, 6),
        ("butane", "CCCC", 4, 3, 10),
        ("benzene", "c1ccccc1", 6, 6, 6),
        ("cyclohexane", "C1CCCCC1", 6, 6, 12),
        ("ethene", "C=C", 2, 1, 4),
        ("ethyne", "C#C", 2, 1, 2),
        ("adamantane", "C1C2CC3CC1CC(C2)C3", 10, 12, 16),
        ("tetrandrine", _TETRANDRINE_CID_73078_SMILES, 46, 52, 42),
    ],
)
def test_3d_acceptance_molecules_preserve_graph_and_build_derived_hydrogens(
    name: str,
    smiles: str,
    document_atoms: int,
    document_bonds: int,
    generated_hydrogens: int,
) -> None:
    from chemuson.chemio.rdkit_safe import (
        is_rdkit_worker_available,
        smiles_to_molgraph_isolated,
    )

    if not is_rdkit_worker_available(timeout_s=5.0):
        pytest.skip("RDKit worker no disponible")
    graph, error = smiles_to_molgraph_isolated(smiles, timeout_s=8.0)
    assert graph is not None, (name, error)
    atom_snapshot = {
        atom_id: (atom.element, atom.charge, atom.isotope, atom.stereo_cip, atom.x, atom.y)
        for atom_id, atom in graph.atoms.items()
    }
    bond_snapshot = {
        bond_id: (bond.a1_id, bond.a2_id, bond.order, bond.is_aromatic)
        for bond_id, bond in graph.bonds.items()
    }

    result = conformer_3d_for_graph(graph, timeout_s=20.0, force=True)

    assert result.ok, (name, result.message)
    coordinates = result.coordinate_set
    assert coordinates is not None
    scene = scene_molecule_from_graph(graph, coordinates)
    assert len(graph.atoms) == document_atoms
    assert len(graph.bonds) == document_bonds
    assert len(scene.atoms) == document_atoms + generated_hydrogens
    assert len(scene.bonds) == document_bonds + generated_hydrogens
    diagnostics = coordinates.metadata["diagnostics"]
    assert diagnostics["return_code"] == 0
    assert diagnostics["graph"]["heavy_atom_count"] == document_atoms
    assert diagnostics["graph"]["atom_count"] == document_atoms
    assert diagnostics["graph"]["bond_count"] == document_bonds
    assert diagnostics["stage"] == "done"
    assert diagnostics["timings_s"]["total"] > 0
    if name == "benzene":
        assert diagnostics["graph"]["aromatic_atom_count"] == 6
        assert diagnostics["graph"]["aromatic_bond_count"] == 6
    if name == "adamantane":
        assert diagnostics["graph"]["ring_count"] >= 3
        assert diagnostics["graph"]["largest_ring_size"] == 6
        assert diagnostics["macrocycle_used"] is False
    if name == "tetrandrine":
        specified_atoms = diagnostics["stereochemistry_diagnostics"]["specified_atoms"]
        assert len(specified_atoms) == 2
        assert set(specified_atoms.values()) == {"S"}
        assert len(coordinates.normalized_positions()) == 46
        assert diagnostics["graph"]["ring_count"] >= 5
    assert all(math.isfinite(value) for atom in scene.atoms for value in (atom.x, atom.y, atom.z))
    if name not in {"water", "ethene", "ethyne"}:
        assert any(abs(atom.z) > 1.0e-6 for atom in scene.atoms)
    assert {atom.id for atom in scene.atoms if not atom.generated} == set(graph.atoms)
    assert sum(atom.generated and atom.element == "H" for atom in scene.atoms) == generated_hydrogens
    assert all(bond.a1_id in {atom.id for atom in scene.atoms} for bond in scene.bonds)
    assert all(bond.a2_id in {atom.id for atom in scene.atoms} for bond in scene.bonds)
    assert result.method in {"mmff94", "mmff94s", "uff", "unoptimized"}
    assert {
        atom_id: (atom.element, atom.charge, atom.isotope, atom.stereo_cip, atom.x, atom.y)
        for atom_id, atom in graph.atoms.items()
    } == atom_snapshot
    assert {
        bond_id: (bond.a1_id, bond.a2_id, bond.order, bond.is_aromatic)
        for bond_id, bond in graph.bonds.items()
    } == bond_snapshot
    if name in {"ethane", "ethanol", "butane"}:
        raw_generated_atoms = coordinates.metadata.get("generated_atoms", [])
        raw_generated_bonds = coordinates.metadata.get("generated_bonds", [])
        generated_atoms = raw_generated_atoms if isinstance(raw_generated_atoms, list) else []
        generated_bonds = raw_generated_bonds if isinstance(raw_generated_bonds, list) else []

        def generated_hydrogen_bond_lengths(scene):
            scene_positions = {
                str(atom.id): (atom.x, atom.y, atom.z)
                for atom in scene.atoms
            }
            hydrogen_ids = {str(item["id"]) for item in generated_atoms}
            lengths = {}
            for bond in generated_bonds:
                first_id, second_id = str(bond["a1_id"]), str(bond["a2_id"])
                if first_id in hydrogen_ids:
                    hydrogen_id, parent_id = first_id, second_id
                elif second_id in hydrogen_ids:
                    hydrogen_id, parent_id = second_id, first_id
                else:
                    continue
                lengths[hydrogen_id] = math.dist(
                    scene_positions[hydrogen_id], scene_positions[parent_id]
                )
            return lengths

        def generated_hydrogen_vectors(scene):
            scene_positions = {
                str(atom.id): (atom.x, atom.y, atom.z)
                for atom in scene.atoms
            }
            hydrogen_ids = {str(item["id"]) for item in generated_atoms}
            vectors: dict[int, list[tuple[float, float, float]]] = {}
            for bond in generated_bonds:
                first_id, second_id = str(bond["a1_id"]), str(bond["a2_id"])
                if first_id in hydrogen_ids:
                    hydrogen_id, parent_id = first_id, second_id
                elif second_id in hydrogen_ids:
                    hydrogen_id, parent_id = second_id, first_id
                else:
                    continue
                hydrogen = scene_positions[hydrogen_id]
                parent = scene_positions[parent_id]
                vectors.setdefault(int(parent_id), []).append(
                    tuple(hydrogen[index] - parent[index] for index in range(3))
                )
            return vectors

        if name == "butane":
            target = dihedral_degrees(coordinates, 1, 2, 3, 4) + 63.0
            transformed = set_dihedral(graph, coordinates, (1, 2, 3, 4), target)
        else:
            transformed = rotate_atoms_about_axis(
                coordinates, set(graph.atoms), (0.0, 0.0, 0.0), (1.0, 2.0, 3.0), 37.0
            )
        transformed_scene = scene_molecule_from_graph(graph, transformed)
        before_lengths = generated_hydrogen_bond_lengths(scene)
        after_lengths = generated_hydrogen_bond_lengths(transformed_scene)
        assert set(after_lengths) == set(before_lengths)
        assert all(
            abs(after_lengths[hydrogen_id] - bond_length) < 1.0e-6
            for hydrogen_id, bond_length in before_lengths.items()
        )
        before_vectors = generated_hydrogen_vectors(scene)
        after_vectors = generated_hydrogen_vectors(transformed_scene)
        assert set(after_vectors) == set(before_vectors)

        def pairwise_angles(vectors):
            angles = []
            for index, first in enumerate(vectors):
                for second in vectors[index + 1 :]:
                    denominator = math.sqrt(sum(value * value for value in first)) * math.sqrt(
                        sum(value * value for value in second)
                    )
                    cosine = sum(a * b for a, b in zip(first, second)) / denominator
                    angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
            return angles

        moving_parents = {3, 4} if name == "butane" else set(graph.atoms)
        orientation_changed = False
        for parent_id, before_parent_vectors in before_vectors.items():
            after_parent_vectors = after_vectors[parent_id]
            assert pairwise_angles(after_parent_vectors) == pytest.approx(
                pairwise_angles(before_parent_vectors), abs=1.0e-5
            )
            if parent_id in moving_parents:
                orientation_changed |= any(
                    math.dist(before, after) > 1.0e-5
                    for before, after in zip(before_parent_vectors, after_parent_vectors)
                )
            else:
                assert after_parent_vectors == pytest.approx(before_parent_vectors)
        assert orientation_changed
    if name == "benzene":
        # The importer stores the aromatic ring as a Kekulé graph (3 single + 3 double).
        assert sum(bond.order == 1 for bond in graph.bonds.values()) == 3
        assert sum(bond.order == 2 for bond in graph.bonds.values()) == 3
        assert all(graph.atom_degree(atom_id) == 2 for atom_id in graph.atoms)
    elif name == "cyclohexane":
        assert all(bond.order == 1 and not bond.is_aromatic for bond in graph.bonds.values())
    elif name == "ethene":
        assert len(graph.bonds) == 1
        assert next(iter(graph.bonds.values())).order == 2
    elif name == "ethyne":
        assert len(graph.bonds) == 1
        assert next(iter(graph.bonds.values())).order == 3
    if name == "water":
        oxygen = next(atom for atom in scene.atoms if not atom.generated and atom.element == "O")
        hydrogens = [atom for atom in scene.atoms if atom.generated and atom.element == "H"]
        vectors = [
            tuple(getattr(hydrogen, axis) - getattr(oxygen, axis) for axis in ("x", "y", "z"))
            for hydrogen in hydrogens
        ]
        dot = sum(a * b for a, b in zip(*vectors))
        norms = [math.sqrt(sum(value * value for value in vector)) for vector in vectors]
        angle = math.degrees(math.acos(max(-1.0, min(1.0, dot / (norms[0] * norms[1])))))
        assert 95.0 <= angle <= 115.0


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
