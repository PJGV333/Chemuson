from __future__ import annotations

import math
import threading
import time

import pytest

from chemuson.chemio import rdkit_safe
from chemuson.core.model import MolGraph


def _butane() -> tuple[MolGraph, dict[int, tuple[float, float, float]]]:
    graph = MolGraph()
    atoms = [graph.add_atom("C", index * 40.0, 0.0) for index in range(4)]
    for first, second in zip(atoms, atoms[1:]):
        graph.add_bond(first.id, second.id, order=1)
    coordinates = {
        atoms[0].id: (0.0, 0.0, 0.0),
        atoms[1].id: (2.5, 0.0, 0.0),
        atoms[2].id: (4.8, 1.4, 0.0),
        atoms[3].id: (7.1, 0.3, 2.0),
    }
    return graph, coordinates


def test_worker_stream_yields_frame_before_process_finishes(tmp_path, monkeypatch):
    worker = tmp_path / "stream_worker.py"
    worker.write_text(
        "import json, time\n"
        "print(json.dumps({'event':'frame','sequence':1}), flush=True)\n"
        "time.sleep(0.5)\n"
        "print(json.dumps({'event':'done','ok':True}), flush=True)\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(rdkit_safe, "_worker_path", lambda: worker)
    started = time.monotonic()

    events = rdkit_safe._run_worker_stream({}, timeout_s=3.0)
    first = next(events)
    first_elapsed = time.monotonic() - started

    assert first == {"event": "frame", "sequence": 1}
    assert first_elapsed < 0.4
    assert next(events) == {"event": "done", "ok": True}
    assert time.monotonic() - started >= 0.4


def test_worker_stream_cancels_child_process(tmp_path, monkeypatch):
    worker = tmp_path / "cancel_worker.py"
    worker.write_text(
        "import json, time\n"
        "print(json.dumps({'event':'progress','stage':'embedding'}), flush=True)\n"
        "time.sleep(10)\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(rdkit_safe, "_worker_path", lambda: worker)
    cancel_event = threading.Event()
    events = rdkit_safe._run_worker_stream(
        {}, timeout_s=15.0, cancel_event=cancel_event
    )

    assert next(events)["stage"] == "embedding"
    started = time.monotonic()
    cancel_event.set()
    cancelled = next(events)

    assert cancelled["event"] == "cancelled"
    assert cancelled["stage"] == "embedding"
    assert time.monotonic() - started < 1.0


def test_conformer_stream_reports_real_stage_progress_and_diagnostics():
    graph, _coordinates = _butane()
    stream = getattr(rdkit_safe, "conformer_3d_stream_isolated")
    events = list(stream(graph, timeout_s=20.0, seed=1234))

    progress = [event for event in events if event.get("event") == "progress"]
    terminal = events[-1]
    diagnostics = terminal["diagnostics"]

    assert terminal["event"] == "done", terminal
    assert terminal["ok"] is True
    assert [event["stage"] for event in progress] == [
        "payload_received",
        "sanitization",
        "hydrogen_stereo",
        "embedding",
        "minimization",
        "geometry_validation",
        "publication",
    ]
    assert all("percent" not in event for event in progress)
    assert diagnostics["error_code"] is None
    assert diagnostics["stage"] == "done"
    assert diagnostics["graph"]["atom_count"] == 4
    assert diagnostics["graph"]["bond_count"] == 3
    assert diagnostics["macrocycle_used"] is False
    assert diagnostics["timings_s"]
    assert all(value >= 0.0 for value in diagnostics["timings_s"].values())
    attempts = diagnostics["attempts"]
    assert attempts
    assert len({attempt["seed"] for attempt in attempts}) == len(attempts)
    assert set(terminal["positions"]) == {"1", "2", "3", "4"}


def test_conformer_preserves_embedded_geometry_when_forcefields_are_unavailable(
    monkeypatch,
):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from chemuson.chemio import _rdkit_worker

    graph, _coordinates = _butane()
    request = rdkit_safe._graph_request_payload(graph, [], "graph_conformer3d")
    monkeypatch.setattr(AllChem, "MMFFHasAllMoleculeParams", lambda _mol: False)
    monkeypatch.setattr(AllChem, "UFFGetMoleculeForceField", lambda _mol: None)

    result = _rdkit_worker._handle_graph_conformer3d_mode(Chem, request)

    assert result["ok"] is True
    assert set(result["positions"]) == {"1", "2", "3", "4"}
    assert result["metadata"]["method"] == "unoptimized"
    assert result["metadata"]["energy"] is None
    assert result["diagnostics"]["degradation_code"] == "missing_forcefield_parameters"


def test_conformer_retries_with_distinct_seeds_and_enables_macrocycle_params(
    monkeypatch,
):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from chemuson.chemio import _rdkit_worker

    graph = MolGraph()
    atoms = [graph.add_atom("C", index * 40.0, 0.0) for index in range(8)]
    for index, atom in enumerate(atoms):
        graph.add_bond(atom.id, atoms[(index + 1) % len(atoms)].id)
    request = rdkit_safe._graph_request_payload(graph, [], "graph_conformer3d")
    original_embed = AllChem.EmbedMolecule
    parameters = []
    calls = 0

    def fail_once(mol, params):
        nonlocal calls
        calls += 1
        parameters.append(params)
        if calls == 1:
            return -1
        return original_embed(mol, params)

    monkeypatch.setattr(AllChem, "EmbedMolecule", fail_once)

    result = _rdkit_worker._handle_graph_conformer3d_mode(Chem, request)

    attempts = result["diagnostics"]["attempts"]
    assert result["ok"] is True
    assert result["diagnostics"]["macrocycle_used"] is True
    assert len(attempts) == 2
    assert attempts[0]["seed"] != attempts[1]["seed"]
    assert all(params.useMacrocycleTorsions for params in parameters)
    assert all(params.useMacrocycle14config for params in parameters)


@pytest.mark.parametrize("forcefield", ["MMFF94", "MMFF94s", "UFF"])
def test_rdkit_stream_emits_progressive_butane_frames_with_fixed_atom(forcefield):
    from chemuson.chemio.rdkit_safe import is_rdkit_worker_available

    if not is_rdkit_worker_available(timeout_s=5.0):
        pytest.skip("RDKit worker no disponible")
    graph, initial = _butane()
    events = list(
        rdkit_safe.optimize_3d_stream_isolated(
            graph,
            coordinates=initial,
            forcefield=forcefield,
            max_iters=80,
            steps_per_update=2,
            fixed_atom_ids={1},
            timeout_s=20.0,
        )
    )

    terminal = events[-1]
    if terminal.get("error") == "missing_mmff_parameters":
        pytest.skip(f"{forcefield} sin parámetros para el butano de prueba")
    frames = [event for event in events if event.get("event") == "frame"]

    assert terminal["event"] == "done", terminal
    assert terminal["ok"] is True
    assert terminal["method"] == forcefield
    assert len(frames) >= 2
    assert [frame["sequence"] for frame in frames] == sorted(
        frame["sequence"] for frame in frames
    )
    assert all(math.isfinite(frame["energy"]) for frame in frames)
    assert all(frame["positions"]["1"] == list(initial[1]) for frame in frames)
    assert any(
        frame["positions"]["2"] != list(initial[2]) for frame in frames
    )
    assert any(
        not math.isclose(first["energy"], second["energy"], abs_tol=1e-10)
        for first, second in zip(frames, frames[1:])
    )
    assert terminal["energy_unit"] == "kcal/mol"

    from chemuson.geometry3d.model import scene_molecule_from_graph
    from chemuson.geometry3d.rdkit_backend import _frame_from_event

    for event in frames:
        live_frame = _frame_from_event(event, source="rdkit")
        assert live_frame is not None
        live_scene = scene_molecule_from_graph(graph, live_frame.coordinates)
        assert len(live_scene.atoms) == 14
        assert len(live_scene.bonds) == 13


def test_live_mm_after_butane_dihedral_emits_bound_derived_hydrogens():
    from chemuson.geometry3d import CoordinateSet3D, conformer_3d_for_graph
    from chemuson.geometry3d.editing import dihedral_degrees, set_dihedral
    from chemuson.geometry3d.model import scene_molecule_from_graph

    if not rdkit_safe.is_rdkit_worker_available(timeout_s=5.0):
        pytest.skip("RDKit worker no disponible")
    graph, _initial_positions = _butane()
    conformer = conformer_3d_for_graph(graph, timeout_s=20.0, force=True)
    assert conformer.ok
    original = conformer.coordinate_set
    assert original is not None
    transformed = set_dihedral(graph, original, (1, 2, 3, 4), 60.0)
    assert dihedral_degrees(transformed, 1, 2, 3, 4) == pytest.approx(60.0, abs=1.0e-5)

    events = list(
        rdkit_safe.optimize_3d_stream_isolated(
            graph,
            coordinates=transformed.normalized_positions(),
            forcefield="MMFF94",
            max_iters=20,
            steps_per_update=4,
            timeout_s=20.0,
        )
    )

    frames = [event for event in events if event.get("event") == "frame"]
    assert events[-1]["event"] == "done", events[-1]
    assert frames
    for event in frames:
        coordinate_set = CoordinateSet3D(
            {int(atom_id): tuple(position) for atom_id, position in event["positions"].items()},
            source="rdkit",
            method=event["forcefield"],
            metadata={
                "generated_atoms": event["generated_atoms"],
                "generated_bonds": event["generated_bonds"],
            },
        )
        scene = scene_molecule_from_graph(graph, coordinate_set)
        assert len(scene.atoms) == 14
        assert len(scene.bonds) == 13
        positions = {str(atom.id): (atom.x, atom.y, atom.z) for atom in scene.atoms}
        hydrogen_ids = {str(atom.id) for atom in scene.atoms if atom.generated}
        assert len(hydrogen_ids) == 10
        for bond in scene.bonds:
            first_id, second_id = str(bond.a1_id), str(bond.a2_id)
            if first_id in hydrogen_ids:
                hydrogen_id, parent_id = first_id, second_id
            elif second_id in hydrogen_ids:
                hydrogen_id, parent_id = second_id, first_id
            else:
                continue
            assert 0.7 <= math.dist(positions[hydrogen_id], positions[parent_id]) <= 1.3


def test_rdkit_backend_stream_returns_last_frame_and_emits_callbacks():
    from chemuson.geometry3d import CoordinateSet3D, ForceField, OptimizationSettings
    from chemuson.geometry3d.rdkit_backend import optimize_stream

    graph, initial = _butane()
    initial_coordinates = CoordinateSet3D(initial, source="test", method="manual")
    callbacks = []
    result = optimize_stream(
        graph,
        initial_coordinates,
        OptimizationSettings(
            forcefield=ForceField.MMFF94,
            max_iters=80,
            steps_per_update=2,
            timeout_s=20.0,
        ),
        fixed_atom_ids=(1,),
        on_frame=callbacks.append,
    )

    assert result.ok, result.message
    assert result.energy is not None and math.isfinite(result.energy)
    assert len(callbacks) >= 2
    assert tuple(callbacks) == result.frames
    assert result.coordinates is not None
    assert result.coordinates.positions == callbacks[-1].coordinates.positions
    assert all(frame.coordinates.positions[1] == initial[1] for frame in callbacks)


def test_worker_stream_times_out_and_terminates_child(tmp_path, monkeypatch):
    worker = tmp_path / "timeout_worker.py"
    worker.write_text(
        "import json, time\n"
        "print(json.dumps({'event':'progress','stage':'minimization'}), flush=True)\n"
        "time.sleep(10)\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(rdkit_safe, "_worker_path", lambda: worker)
    started = time.monotonic()

    events = list(rdkit_safe._run_worker_stream({}, timeout_s=0.15))

    assert events[-1]["error"] == "timeout"
    assert events[-1]["stage"] == "minimization"
    assert time.monotonic() - started < 1.0


def test_worker_stream_reports_nonzero_exit_after_progress(tmp_path, monkeypatch):
    worker = tmp_path / "failed_worker.py"
    worker.write_text(
        "import json, sys\n"
        "print(json.dumps({'event':'frame','sequence':1}), flush=True)\n"
        "sys.exit(7)\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(rdkit_safe, "_worker_path", lambda: worker)

    events = list(rdkit_safe._run_worker_stream({}, timeout_s=3.0))

    assert events[0]["event"] == "frame"
    assert events[-1]["event"] == "error"
    assert events[-1]["error"] == "worker_exit_code"
    assert events[-1]["returncode"] == 7
