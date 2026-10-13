from __future__ import annotations

import math
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

from chemuson.chemio import rdkit_io
from chemuson.chemio.rdkit_io import (
    molfile_to_molgraph,
    molgraph_to_molfile,
    molgraph_to_smiles,
    smiles_to_molgraph,
)
from chemuson.core.model import MolGraph
from chemuson.gui.controllers import file_controller as file_controller_module
from chemuson.gui.controllers.file_controller import FileController


S_SMILES = "C[C@H](F)CO"
R_SMILES = "C[C@@H](F)CO"
MOL_FILTER = "Archivo MOL (*.mol)"


def _canonical(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def _cip(smiles: str) -> tuple[str, ...]:
    molecule = Chem.MolFromSmiles(smiles)
    assert molecule is not None
    return tuple(
        sorted(
            label
            for _index, label in Chem.FindMolChiralCenters(
                molecule,
                includeUnassigned=False,
                useLegacyImplementation=False,
            )
        )
    )


def _molblock_geometry(molblock: str) -> tuple[list[tuple[float, float]], list[float]]:
    lines = molblock.splitlines()
    counts_index = next(index for index, line in enumerate(lines) if "V2000" in line)
    atom_count = int(lines[counts_index][:3])
    bond_count = int(lines[counts_index][3:6])
    atom_lines = lines[counts_index + 1 : counts_index + 1 + atom_count]
    coordinates = [(float(line.split()[0]), float(line.split()[1])) for line in atom_lines]
    bond_lines = lines[counts_index + 1 + atom_count :][:bond_count]
    lengths = []
    for line in bond_lines:
        tokens = line.split()
        first, second = int(tokens[0]) - 1, int(tokens[1]) - 1
        lengths.append(math.dist(coordinates[first], coordinates[second]))
    return coordinates, lengths


def _geometry_graph() -> MolGraph:
    graph = MolGraph()
    atoms = [
        graph.add_atom("C", 100.0, 200.0),
        graph.add_atom("C", 140.0, 200.0),
        graph.add_atom("O", 160.0, 230.0, is_explicit=True),
        graph.add_atom("F", 190.0, 230.0, is_explicit=True),
    ]
    for first, second in zip(atoms, atoms[1:], strict=False):
        graph.add_bond(first.id, second.id)
    return graph


def _save_context(tmp_path: Path):
    class UndoStack:
        def __init__(self) -> None:
            self.clean = False

        def setClean(self) -> None:
            self.clean = True

    class Autosave:
        def __init__(self) -> None:
            self.cleaned = False

        def cleanup_after_save(self) -> None:
            self.cleaned = True

    class TabManager:
        def __init__(self) -> None:
            self.path = None
            self.autosave = Autosave()

        def file_path_for(self, _canvas):
            return self.path

        def autosave_manager_for(self, _canvas):
            return self.autosave

        def set_canvas_file_path(self, _canvas, path):
            self.path = path

    canvas = SimpleNamespace(graph=_geometry_graph(), undo_stack=UndoStack())
    tab_manager = TabManager()
    recent: list[str] = []
    statuses: list[str] = []
    context = SimpleNamespace(
        parent=None,
        canvas=canvas,
        tab_manager=tab_manager,
        add_recent_file=recent.append,
        show_status=statuses.append,
    )
    return context, tab_manager, recent, statuses


def _run_open_probe(path: Path, expected_smiles: str, tmp_path: Path) -> dict[str, str]:
    probe = Path(__file__).parent / "probes" / "p0_open_mol.py"
    env = os.environ.copy()
    project_paths = [str(Path(__file__).resolve().parents[1] / "src")]
    chem_packages = Path(__file__).resolve().parents[1] / "chem/lib/python3.14/site-packages"
    if chem_packages.is_dir():
        project_paths.insert(0, str(chem_packages))
    if env.get("PYTHONPATH"):
        project_paths.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(project_paths)
    env["QT_QPA_PLATFORM"] = "offscreen"
    env["XDG_CONFIG_HOME"] = str(tmp_path / "config")
    result = subprocess.run(
        [sys.executable, str(probe), str(path), expected_smiles],
        capture_output=True,
        text=True,
        timeout=30.0,
        env=env,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return {"stdout": result.stdout, "stderr": result.stderr}


def test_mol_export_uses_same_uniform_scale_in_rdkit_and_fallback(monkeypatch) -> None:
    for backend in ("fallback", "rdkit"):
        graph = _geometry_graph()
        original_coordinates = {atom.id: (atom.x, atom.y) for atom in graph.atoms.values()}
        original_lengths = [
            math.dist(
                original_coordinates[bond.a1_id],
                original_coordinates[bond.a2_id],
            )
            for bond in graph.bonds.values()
        ]
        original_area = (
            (graph.atoms[2].x - graph.atoms[1].x) * (graph.atoms[3].y - graph.atoms[1].y)
            - (graph.atoms[2].y - graph.atoms[1].y) * (graph.atoms[3].x - graph.atoms[1].x)
        )
        if backend == "fallback":
            monkeypatch.setattr(rdkit_io, "_rdkit_available", lambda: False)
        else:
            monkeypatch.setattr(rdkit_io, "Chem", Chem)
            monkeypatch.setattr(rdkit_io, "AllChem", AllChem)
            monkeypatch.setattr(rdkit_io, "_rdkit_available", lambda: True)

        molblock = molgraph_to_molfile(graph)
        coordinates, lengths = _molblock_geometry(molblock)
        assert sorted(lengths)[len(lengths) // 2] == pytest.approx(1.5, abs=2e-4)
        assert graph.atoms[1].x == original_coordinates[1][0]
        assert graph.atoms[1].y == original_coordinates[1][1]
        assert {atom.id: (atom.x, atom.y) for atom in graph.atoms.values()} == original_coordinates
        assert original_area > 0
        exported_area = (
            (coordinates[1][0] - coordinates[0][0]) * (coordinates[2][1] - coordinates[0][1])
            - (coordinates[1][1] - coordinates[0][1]) * (coordinates[2][0] - coordinates[0][0])
        )
        assert exported_area > 0, "uniform positive scaling must not reflect the drawing"
        ratios = [new / old for new, old in zip(lengths, original_lengths, strict=True)]
        assert max(ratios) - min(ratios) < 2e-5


def test_s_r_mol_roundtrip_preserves_identity_cip_and_graph_data() -> None:
    expected_formula = "C3H7FO"
    observed_canonical = []
    for smiles, expected_cip in ((S_SMILES, ("S",)), (R_SMILES, ("R",))):
        graph = smiles_to_molgraph(smiles)
        source_coordinates = {atom.id: (atom.x, atom.y) for atom in graph.atoms.values()}
        molblock = molgraph_to_molfile(graph)
        _coordinates, bond_lengths = _molblock_geometry(molblock)
        assert len(graph.atoms) == 5
        assert len(graph.bonds) == 4
        assert sorted(bond_lengths)[len(bond_lengths) // 2] == pytest.approx(1.5, abs=2e-4)
        assert {atom.id: (atom.x, atom.y) for atom in graph.atoms.values()} == source_coordinates

        restored = molfile_to_molgraph(molblock)
        assert len(restored.atoms) == 5
        assert len(restored.bonds) == 4
        assert sorted(atom.element for atom in restored.atoms.values()) == ["C", "C", "C", "F", "O"]
        assert all(atom.charge == 0 for atom in restored.atoms.values())
        assert all(atom.isotope in (None, 0) for atom in restored.atoms.values())
        restored_smiles = molgraph_to_smiles(restored)
        restored_molecule = Chem.MolFromSmiles(restored_smiles)
        assert restored_molecule is not None
        assert rdMolDescriptors.CalcMolFormula(restored_molecule) == expected_formula
        assert Chem.MolToSmiles(restored_molecule, canonical=True, isomericSmiles=True) == _canonical(smiles)
        graph_cip = tuple(sorted(atom.stereo_cip for atom in restored.atoms.values() if atom.stereo_cip))
        assert graph_cip == expected_cip
        assert _cip(restored_smiles) == expected_cip
        observed_canonical.append(_canonical(restored_smiles))
    assert observed_canonical[0] != observed_canonical[1]


def test_fallback_wedge_roundtrip_recovers_cip_without_mutating_source(monkeypatch) -> None:
    monkeypatch.setattr(rdkit_io, "_rdkit_available", lambda: False)
    for smiles, expected_cip in ((S_SMILES, ("S",)), (R_SMILES, ("R",))):
        graph = smiles_to_molgraph(smiles)
        for atom in graph.atoms.values():
            atom.stereo_cip = None
        source_coordinates = {atom.id: (atom.x, atom.y) for atom in graph.atoms.values()}

        molblock = molgraph_to_molfile(graph)
        restored = molfile_to_molgraph(molblock)

        assert {atom.id: (atom.x, atom.y) for atom in graph.atoms.values()} == source_coordinates
        assert tuple(sorted(atom.stereo_cip for atom in restored.atoms.values() if atom.stereo_cip)) == expected_cip
        assert _canonical(molgraph_to_smiles(restored)) == _canonical(smiles)


def test_rdkit_generated_mol_opens_through_chemuson_file_controller(tmp_path) -> None:
    molecule = Chem.MolFromSmiles(S_SMILES)
    assert molecule is not None
    AllChem.Compute2DCoords(molecule)
    rdkit_path = tmp_path / "rdkit-generated.mol"
    rdkit_path.write_text(Chem.MolToMolBlock(molecule), encoding="utf-8")

    _run_open_probe(rdkit_path, S_SMILES, tmp_path)


def test_chemuson_exported_mol_opens_through_chemuson_file_controller(tmp_path) -> None:
    molblock = molgraph_to_molfile(smiles_to_molgraph(S_SMILES))
    chemuson_path = tmp_path / "chemuson-export.mol"
    chemuson_path.write_text(molblock, encoding="utf-8")

    _run_open_probe(chemuson_path, S_SMILES, tmp_path)


@pytest.mark.parametrize(
    ("entered_name", "expected_name"),
    (("molecule", "molecule.mol"), ("molecule.mol", "molecule.mol")),
)
def test_save_as_mol_adds_extension_once_and_updates_document_after_write(
    tmp_path,
    monkeypatch,
    entered_name: str,
    expected_name: str,
) -> None:
    context, tab_manager, recent, statuses = _save_context(tmp_path)
    selected_path = tmp_path / entered_name
    monkeypatch.setattr(
        file_controller_module.QFileDialog,
        "getSaveFileName",
        lambda *_args, **_kwargs: (str(selected_path), MOL_FILTER),
    )

    FileController().save_file(context)

    expected_path = tmp_path / expected_name
    assert expected_path.is_file()
    assert "V2000" in expected_path.read_text(encoding="utf-8")
    assert tab_manager.path == str(expected_path)
    assert context.canvas.undo_stack.clean
    assert tab_manager.autosave.cleaned
    assert recent == [str(expected_path)]
    assert statuses == [f"Guardado: {expected_path}"]


def test_save_as_cmsn_keeps_native_extension_and_writer(monkeypatch, tmp_path) -> None:
    context, tab_manager, recent, _statuses = _save_context(tmp_path)
    selected_path = tmp_path / "native-document"
    writes: list[str] = []

    def fake_save(path, _canvas):
        writes.append(path)
        Path(path).write_text("native cmsn payload", encoding="utf-8")

    monkeypatch.setattr(file_controller_module.PersistenceManager, "save_to_file", fake_save)
    monkeypatch.setattr(
        file_controller_module.QFileDialog,
        "getSaveFileName",
        lambda *_args, **_kwargs: (str(selected_path), "Archivo de ChemUSON (*.cmsn)"),
    )

    FileController().save_file(context)

    expected_path = tmp_path / "native-document.cmsn"
    assert writes == [str(expected_path)]
    assert expected_path.read_text(encoding="utf-8") == "native cmsn payload"
    assert tab_manager.path == str(expected_path)
    assert context.canvas.undo_stack.clean
    assert recent == [str(expected_path)]


def test_mol_write_failure_is_reported_without_marking_document_saved(monkeypatch, tmp_path) -> None:
    context, tab_manager, recent, statuses = _save_context(tmp_path)
    selected_path = tmp_path / "not-written"
    errors: list[str] = []
    monkeypatch.setattr(
        file_controller_module.QFileDialog,
        "getSaveFileName",
        lambda *_args, **_kwargs: (str(selected_path), MOL_FILTER),
    )
    monkeypatch.setattr(
        file_controller_module,
        "molgraph_to_molfile",
        lambda _graph: (_ for _ in ()).throw(OSError("disk full")),
    )
    monkeypatch.setattr(
        file_controller_module.QMessageBox,
        "critical",
        lambda _parent, _title, message: errors.append(message),
    )

    FileController().save_file(context)

    assert errors == ["No se pudo guardar:\ndisk full"]
    assert not selected_path.with_suffix(".mol").exists()
    assert tab_manager.path is None
    assert not context.canvas.undo_stack.clean
    assert not tab_manager.autosave.cleaned
    assert recent == []
    assert statuses == []


def test_sdf_is_rejected_instead_of_saving_a_mol_block(monkeypatch, tmp_path) -> None:
    context, tab_manager, recent, statuses = _save_context(tmp_path)
    selected_path = tmp_path / "misleading.sdf"
    errors: list[str] = []
    monkeypatch.setattr(
        file_controller_module.QFileDialog,
        "getSaveFileName",
        lambda *_args, **_kwargs: (str(selected_path), "Todos los archivos (*.*)"),
    )
    monkeypatch.setattr(
        file_controller_module.QMessageBox,
        "critical",
        lambda _parent, _title, message: errors.append(message),
    )

    FileController().save_file(context)

    assert errors == [
        "SDF multi-registro no está soportado. Guarda la estructura como MOL (.mol)."
    ]
    assert not selected_path.exists()
    assert tab_manager.path is None
    assert not context.canvas.undo_stack.clean
    assert recent == []
    assert statuses == []


def test_mol_open_filter_does_not_advertise_unimplemented_sdf() -> None:
    source = (Path(__file__).resolve().parents[1] / "src/chemuson/gui/main_window.py").read_text(
        encoding="utf-8"
    )
    open_handler = source.split("def _on_file_open(self)", maxsplit=1)[1].split(
        "def _load_file_into_canvas", maxsplit=1
    )[0]
    assert "Archivos MOL (*.mol)" in open_handler
    assert "*.sdf" not in open_handler


def test_sdf_import_is_rejected_before_clearing_canvas(tmp_path) -> None:
    path = tmp_path / "multi-record.sdf"
    path.write_text("not an SDF", encoding="utf-8")

    class Canvas:
        cleared = False

        def clear_canvas(self):
            self.cleared = True

    canvas = Canvas()
    with pytest.raises(ValueError, match="SDF multi-registro no está soportado"):
        FileController().load_file_into_canvas(str(path), canvas)
    assert not canvas.cleared


def test_mol_read_failure_is_shown_by_open_file_workflow(monkeypatch, tmp_path) -> None:
    path = tmp_path / "invalid.mol"
    path.write_text("not a CTAB", encoding="utf-8")
    errors: list[str] = []

    class UndoStack:
        def setClean(self) -> None:
            pass

    class Canvas:
        undo_stack = UndoStack()

        def clear_canvas(self) -> None:
            pass

        def _insert_molgraph(self, _graph) -> None:
            pass

    class Tabs:
        def count(self) -> int:
            return 1

    class TabManager:
        def set_canvas_file_path(self, _canvas, _path) -> None:
            pass

        def discard_canvas(self, _canvas) -> bool:
            return True

    canvas = Canvas()
    context = SimpleNamespace(
        parent=None,
        canvas=canvas,
        tabs=Tabs(),
        tab_manager=TabManager(),
        create_document_tab=lambda _make_current: canvas,
        apply_toolbar_defaults=lambda _canvas: None,
        set_active_canvas=lambda _canvas: None,
        before_canvas_discard=lambda _canvas: None,
        add_recent_file=lambda _path: None,
        refresh_recent_menu=lambda: None,
        update_total_charge_indicator=lambda: None,
        show_status=lambda _message: None,
    )
    monkeypatch.setattr(
        file_controller_module.QMessageBox,
        "critical",
        lambda _parent, _title, message: errors.append(message),
    )

    FileController().open_file_path(context, str(path))

    assert errors and errors[0].startswith("No se pudo abrir el archivo:")
    assert "Mol inválido" in errors[0]
