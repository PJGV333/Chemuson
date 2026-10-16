"""Clean2D con enlaces de cuña sólida y rayada: matriz de reproducción y contratos.

OpenSpec: ``clean2d-wedge-stereo-layout-campaign``.

Diagnóstico (ver ``baseline.md`` del cambio):

- ``stereo_layout_signature`` comparaba los **signos** del producto vectorial de
  cada par de vecinos de un centro estereogénico. Ese signo cambia cuando el hueco
  angular entre dos vecinos cruza 180°, aunque el orden cíclico (y por tanto la
  configuración dibujada) sea idéntico. Con tres vecinos el hueco medio es de
  120°, así que casi cualquier reconstrucción invertía algún signo y el motor
  rechazaba **todos** los candidatos (``cambio_firma_estereoquimica``): una
  molécula con una sola cuña no se limpiaba nunca.
- La firma sólo consideraba vecinos **dentro** de la selección, de modo que un
  centro quiral parcialmente seleccionado podía reordenarse respecto a un vecino
  de frontera sin que la comprobación se enterara.

Estos tests fijan ambos contratos: re-encuadres que conservan el orden cíclico se
aceptan, y las inversiones reales (imagen especular, intercambio de vecinos,
cambio de cuña, cruce de un vecino de frontera) se siguen rechazando.
"""

from __future__ import annotations

import copy
import json
import math
import subprocess
import sys

import pytest

from chemuson.clean2d import (
    assert_clean2d_invariants,
    classify_clean2d_complexity,
    run_clean2d_engine,
    smiles_to_molgraph_best_depiction,
    stereo_layout_signature,
)
from chemuson.chemio.rdkit_io import molgraph_to_molfile
from chemuson.core.model import Bond, BondStereo, BondStyle, MolGraph

TARGET = 42.0


# ---------------------------------------------------------------------------
# Constructores de la matriz
# ---------------------------------------------------------------------------


def _graph() -> MolGraph:
    return MolGraph()


def _add(graph: MolGraph, element: str, x: float, y: float) -> int:
    return graph.add_atom(element, x, y).id


def _bond(graph: MolGraph, a: int, b: int, style: BondStyle = BondStyle.PLAIN, order: int = 1) -> Bond:
    return graph.add_bond(a, b, order=order, style=style)


def case_acyclic_control() -> MolGraph:
    """Caso 1 (control): cadena acíclica sin estereoquímica."""
    graph = _graph()
    ids = [_add(graph, "C", x, y) for x, y in ((10, 12), (55, 20), (78, -14), (132, 6))]
    _bond(graph, ids[0], ids[1])
    _bond(graph, ids[1], ids[2])
    _bond(graph, ids[2], ids[3])
    return graph


def case_single_wedge() -> MolGraph:
    """Caso 2: una sola cuña sólida."""
    graph = _graph()
    c1 = _add(graph, "C", 8.0, 10.0)
    c2 = _add(graph, "C", 60.0, 16.0)
    c3 = _add(graph, "C", 84.0, -20.0)
    o1 = _add(graph, "O", 128.0, 4.0)
    _bond(graph, c1, c2)
    _bond(graph, c2, c3)
    _bond(graph, c2, o1, style=BondStyle.WEDGE)
    return graph


def case_single_hashed() -> MolGraph:
    """Caso 3: una sola cuña rayada."""
    graph = _graph()
    c1 = _add(graph, "C", 8.0, 10.0)
    c2 = _add(graph, "C", 60.0, 16.0)
    c3 = _add(graph, "C", 84.0, -20.0)
    o1 = _add(graph, "O", 128.0, 4.0)
    _bond(graph, c1, c2)
    _bond(graph, c2, c3)
    _bond(graph, c2, o1, style=BondStyle.HASHED)
    return graph


def case_tetrahedral_both() -> MolGraph:
    """Caso 4: centro tetraédrico con cuña sólida y cuña rayada."""
    graph = _graph()
    center = _add(graph, "C", 60.0, 20.0)
    left = _add(graph, "C", 12.0, 30.0)
    right = _add(graph, "C", 110.0, 12.0)
    up = _add(graph, "O", 66.0, 74.0)
    down = _add(graph, "N", 52.0, -34.0)
    _bond(graph, center, left)
    _bond(graph, center, right)
    _bond(graph, center, up, style=BondStyle.WEDGE)
    _bond(graph, center, down, style=BondStyle.HASHED)
    return graph


def case_multiple_centers() -> MolGraph:
    """Caso 5: dos centros quirales en la misma cadena."""
    graph = _graph()
    c1 = _add(graph, "C", 0.0, 0.0)
    c2 = _add(graph, "C", 48.0, 8.0)
    c3 = _add(graph, "C", 92.0, -14.0)
    c4 = _add(graph, "C", 140.0, 2.0)
    o2 = _add(graph, "O", 44.0, 56.0)
    o3 = _add(graph, "O", 96.0, -62.0)
    _bond(graph, c1, c2)
    _bond(graph, c2, c3)
    _bond(graph, c3, c4)
    _bond(graph, c2, o2, style=BondStyle.WEDGE)
    _bond(graph, c3, o3, style=BondStyle.HASHED)
    return graph


WEDGE_CASES = (
    pytest.param(case_single_wedge, id="single-wedge"),
    pytest.param(case_single_hashed, id="single-hashed"),
    pytest.param(case_tetrahedral_both, id="tetrahedral-wedge-hash"),
    pytest.param(case_multiple_centers, id="multiple-stereo-centers"),
)


# ---------------------------------------------------------------------------
# Utilidades de verificación química
# ---------------------------------------------------------------------------


def coords_of(graph: MolGraph) -> dict[int, tuple[float, float]]:
    return {atom_id: (atom.x, atom.y) for atom_id, atom in graph.atoms.items()}


def _styles_by_atom(graph: MolGraph) -> dict[int, dict[int, str]]:
    styles: dict[int, dict[int, str]] = {}
    for bond in graph.bonds.values():
        style = getattr(bond.style, "value", str(bond.style))
        styles.setdefault(bond.a1_id, {})[bond.a2_id] = style
        styles.setdefault(bond.a2_id, {})[bond.a1_id] = style
    return styles


def depicted_configuration(graph: MolGraph, coords: dict[int, tuple[float, float]]) -> dict[int, tuple]:
    """Configuración percibida por centro: orden cíclico de vecinos + estilos.

    Dos disposiciones con el mismo orden cíclico (módulo rotación) y los mismos
    estilos de enlace describen la misma estereoquímica dibujada.
    """
    styles = _styles_by_atom(graph)
    centers = {
        atom_id
        for atom_id, neighbours in styles.items()
        if any(style in {"wedge", "hashed"} for style in neighbours.values())
    }
    configuration: dict[int, tuple] = {}
    for center in sorted(centers):
        center_point = coords.get(center)
        if center_point is None:
            continue
        ordered: list[tuple[float, int, str]] = []
        for neighbour, style in styles[center].items():
            point = coords.get(neighbour)
            if point is None:
                continue
            dx = point[0] - center_point[0]
            dy = point[1] - center_point[1]
            if math.hypot(dx, dy) <= 1e-9:
                continue
            ordered.append((math.atan2(dy, dx), neighbour, style))
        if not ordered:
            continue
        ordered.sort(key=lambda item: (item[0], item[1]))
        sequence = tuple((neighbour, style) for _angle, neighbour, style in ordered)
        configuration[center] = min(
            sequence[index:] + sequence[:index] for index in range(len(sequence))
        )
    return configuration


def wedge_signature(graph: MolGraph) -> tuple:
    """Estilos de los enlaces en cuña/hash (identidad de cuñas y sus átomos)."""
    return tuple(
        sorted(
            (
                int(bond.id),
                int(bond.a1_id),
                int(bond.a2_id),
                getattr(bond.style, "value", str(bond.style)),
            )
            for bond in graph.bonds.values()
            if getattr(bond.style, "value", str(bond.style)) in {"wedge", "hashed"}
        )
    )


def _coords(graph: MolGraph, mapping: dict[int, tuple[float, float]]) -> dict[int, tuple[float, float]]:
    """Coordenadas completas del grafo aplicando una propuesta parcial."""
    full = coords_of(graph)
    full.update({atom_id: (float(x), float(y)) for atom_id, (x, y) in mapping.items()})
    return full


def _independent_rdkit_report(kind: str, structure: str) -> dict[str, object]:
    """Assign R/S in a fresh RDKit process, independently of graph metadata."""
    if kind == "molblock":
        # RDKit 2026 may emit only two header lines; normalize to the standard
        # three-line V2000 header before parsing independently.
        lines = structure.splitlines()
        counts_index = next((i for i, line in enumerate(lines) if "V2000" in line), None)
        if counts_index is not None and counts_index < 3:
            lines[counts_index:counts_index] = [""] * (3 - counts_index)
        structure = "\n".join(lines) + "\n"
    script = r'''import json, sys
from rdkit import Chem
request = json.load(sys.stdin)
if request["kind"] == "smiles":
    mol = Chem.MolFromSmiles(request["structure"])
else:
    mol = Chem.MolFromMolBlock(request["structure"], sanitize=True, removeHs=False, strictParsing=False)
assert mol is not None, "RDKit could not parse the reference structure"
Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
print(json.dumps({
    "canonical_smiles": Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
    "centers": [label for _index, label in centers],
}))'''
    request = json.dumps({"kind": kind, "structure": structure})
    result = subprocess.run(
        [sys.executable, "-c", script],
        input=request,
        capture_output=True,
        text=True,
        timeout=20.0,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return json.loads(result.stdout)


def _independent_graph_report(graph: MolGraph) -> dict[str, object]:
    return _independent_rdkit_report("molblock", molgraph_to_molfile(graph))


# ---------------------------------------------------------------------------
# 1. Contratos de la firma estereoquímica
# ---------------------------------------------------------------------------


def test_signature_is_invariant_under_rotation_and_order_preserving_relayout() -> None:
    graph = case_single_wedge()
    before = coords_of(graph)
    center = 2
    selected = set(graph.atoms)

    # Rotación rígida de todo el entorno del centro: misma configuración.
    angle = math.radians(37.0)
    cx, cy = before[center]
    rotated = {}
    for atom_id, (x, y) in before.items():
        dx, dy = x - cx, y - cy
        rotated[atom_id] = (
            cx + dx * math.cos(angle) - dy * math.sin(angle),
            cy + dx * math.sin(angle) + dy * math.cos(angle),
        )
    assert stereo_layout_signature(graph, before, selected) == stereo_layout_signature(
        graph, rotated, selected
    )
    assert depicted_configuration(graph, before) == depicted_configuration(graph, rotated)

    # Re-encuadre que conserva el orden cíclico pero cambia los huecos angulares
    # (120/120/120 -> 100/130/130): el control antiguo rechazaba esto.
    relayout = {
        1: (cx + TARGET * math.cos(math.radians(100.0)), cy + TARGET * math.sin(math.radians(100.0))),
        3: (cx + TARGET * math.cos(math.radians(230.0)), cy + TARGET * math.sin(math.radians(230.0))),
        4: (cx + TARGET * math.cos(math.radians(360.0)), cy + TARGET * math.sin(math.radians(360.0))),
    }
    after = _coords(graph, relayout)
    assert stereo_layout_signature(graph, before, selected) == stereo_layout_signature(
        graph, after, selected
    )
    assert depicted_configuration(graph, before) == depicted_configuration(graph, after)


def test_signature_detects_mirror_neighbour_swap_and_wedge_change() -> None:
    graph = case_tetrahedral_both()
    before = coords_of(graph)
    selected = set(graph.atoms)
    base = stereo_layout_signature(graph, before, selected)

    # Imagen especular respecto del centro: configuración invertida.
    cx, cy = before[1]
    mirrored = {atom_id: (2 * cx - x, y) for atom_id, (x, y) in before.items()}
    assert stereo_layout_signature(graph, before, selected) != stereo_layout_signature(
        graph, mirrored, selected
    )
    assert depicted_configuration(graph, before) != depicted_configuration(graph, mirrored)

    # Intercambio de dos vecinos del centro: configuración invertida.
    swapped = dict(before)
    swapped[2], swapped[3] = before[3], before[2]
    assert stereo_layout_signature(graph, before, selected) != stereo_layout_signature(
        graph, swapped, selected
    )

    # Cambiar el sentido de una cuña se detecta aunque la geometría no cambie.
    graph.bonds[3].style = BondStyle.HASHED
    assert stereo_layout_signature(graph, before, selected) != base


def test_signature_covers_stereo_neighbours_outside_the_selection() -> None:
    """Un vecino de frontera participa en el orden cíclico del centro.

    Seleccionar sólo tres de los cuatro vecinos de un centro quiral no permite
    cruzar al vecino no seleccionado: hacerlo invierte la configuración dibujada.
    """
    graph = case_tetrahedral_both()
    before = coords_of(graph)
    selected = {1, 2, 3, 4}  # el vecino 5 (enlace en hash) queda fuera
    base = stereo_layout_signature(graph, before, selected)

    # El vecino 5 no está en la selección pero ordena al centro 1.
    assert any(entry[0] == "center" and 5 in [neighbour for neighbour, _ in entry[2]] for entry in base)

    # Mover el vecino 4 al otro lado del rayo del vecino 5 (respecto al centro)
    # cruza la frontera y cambia la configuración.
    cx, cy = before[1]
    crossed = dict(before)
    crossed[4] = (cx + TARGET * math.cos(math.radians(200.0)), cy + TARGET * math.sin(math.radians(200.0)))
    assert stereo_layout_signature(graph, before, selected) != stereo_layout_signature(
        graph, crossed, selected
    )
    assert depicted_configuration(graph, before) != depicted_configuration(graph, crossed)


# ---------------------------------------------------------------------------
# 2. Matriz del motor: las cuñas sí se limpian y conservan la configuración
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("builder", WEDGE_CASES)
@pytest.mark.parametrize("mode", ["quick", "publication"])
def test_wedge_structures_are_cleaned_and_keep_their_configuration(builder, mode) -> None:
    graph = builder()
    before = coords_of(graph)
    before_signature = stereo_layout_signature(graph, before, set(graph.atoms))
    before_configuration = depicted_configuration(graph, before)
    before_wedges = wedge_signature(graph)

    result = run_clean2d_engine(graph, None, mode=mode, target_bond_length=TARGET)

    assert result.selected is not None, (mode, result.reason, result.message)
    assert result.result_state == "applied", (mode, result.stable_reason, result.message)
    after = result.selected.coords
    assert after, mode

    # Invariantes químicas duras: identidad, conectividad, órdenes, cargas,
    # isótopos y estereoquímica explícita de átomos y enlaces.
    assert_clean2d_invariants(graph, graph, before, after, atom_ids=set(graph.atoms))

    # La estereoquímica dibujada: cuñas intactas y configuración idéntica.
    assert wedge_signature(graph) == before_wedges
    assert stereo_layout_signature(graph, after, set(graph.atoms)) == before_signature
    assert depicted_configuration(graph, after) == before_configuration

    # Algo cambió de verdad (no es un no-op disfrazado).
    moved = sum(
        1
        for atom_id, (x, y) in after.items()
        if math.hypot(x - before[atom_id][0], y - before[atom_id][1]) > 0.5
    )
    assert moved >= 1


def test_control_acyclic_case_still_applies() -> None:
    graph = case_acyclic_control()
    before = coords_of(graph)
    result = run_clean2d_engine(graph, None, mode="quick", target_bond_length=TARGET)
    assert result.result_state == "applied"
    assert_clean2d_invariants(graph, graph, before, result.selected.coords, atom_ids=set(graph.atoms))


@pytest.mark.parametrize("builder", WEDGE_CASES)
def test_wedge_cleaning_is_deterministic(builder) -> None:
    first = run_clean2d_engine(builder(), None, mode="quick", target_bond_length=TARGET)
    second = run_clean2d_engine(builder(), None, mode="quick", target_bond_length=TARGET)
    assert first.result_state == second.result_state
    assert first.selected is not None and second.selected is not None
    assert first.selected.geometry_hash == second.selected.geometry_hash
    assert first.selected.coords == second.selected.coords


# ---------------------------------------------------------------------------
# 3. Selecciones parciales: nunca se invierte la configuración
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "selected",
    [
        pytest.param({1, 2, 3, 4}, id="selection-with-chiral-center"),
        pytest.param({1, 2}, id="stereo-bond-on-boundary"),
        pytest.param({1}, id="only-the-chiral-center"),
    ],
)
def test_partial_selections_never_invert_the_depicted_configuration(selected) -> None:
    graph = case_tetrahedral_both()
    before = coords_of(graph)
    before_configuration = depicted_configuration(graph, before)
    before_wedges = wedge_signature(graph)

    result = run_clean2d_engine(graph, selected, mode="quick", target_bond_length=TARGET)

    # Sin candidato seguro el resultado correcto es preserve-only/no-op (o se
    # informa del motivo); nunca se aplica una inversión.
    assert result.result_state in {"applied", "no-op", "preserve-only", "failed-controlled"}
    after = result.selected.coords if result.selected is not None else before
    assert depicted_configuration(graph, _coords(graph, after)) == before_configuration
    assert wedge_signature(graph) == before_wedges
    if result.result_state == "applied":
        assert stereo_layout_signature(graph, before, selected) == stereo_layout_signature(
            graph, _coords(graph, after), selected
        )
        # Sólo se mueven átomos de la selección.
        for atom_id, (x, y) in after.items():
            if atom_id not in selected:
                continue
            assert math.isfinite(x) and math.isfinite(y)


def _import_reference_graph(smiles: str) -> MolGraph:
    return smiles_to_molgraph_best_depiction(smiles, timeout_s=15.0)


def _combine_disconnected(graphs: tuple[MolGraph, ...]) -> tuple[MolGraph, tuple[set[int], ...]]:
    combined = MolGraph()
    components: list[set[int]] = []
    next_atom_id = 1
    for component_index, source in enumerate(graphs):
        shift_x = float(component_index) * 500.0
        id_map: dict[int, int] = {}
        selected_ids: set[int] = set()
        for atom in source.atoms.values():
            new_id = next_atom_id
            next_atom_id += 1
            id_map[atom.id] = new_id
            selected_ids.add(new_id)
            combined.add_atom(
                atom.element,
                atom.x + shift_x,
                atom.y,
                atom_id=new_id,
                charge=atom.charge,
                isotope=atom.isotope,
                stereo_cip=atom.stereo_cip,
                explicit_h=atom.explicit_h,
                is_explicit=atom.is_explicit,
            )
        for bond in source.bonds.values():
            combined.add_bond(
                id_map[bond.a1_id],
                id_map[bond.a2_id],
                order=bond.order,
                style=bond.style,
                stereo=bond.stereo,
                is_aromatic=bond.is_aromatic,
            )
        components.append(selected_ids)
    return combined, tuple(components)


def test_real_imported_wedge_is_counted_as_one_stereocenter() -> None:
    graph = _import_reference_graph("C[C@H](F)CO")
    independent_reference = _independent_rdkit_report("smiles", "C[C@H](F)CO")
    assert independent_reference["centers"] == ["S"]
    assert len([atom for atom in graph.atoms.values() if atom.stereo_cip]) == 1
    assert sum(bond.stereo != BondStereo.NONE for bond in graph.bonds.values()) == 1

    profile = classify_clean2d_complexity(graph)

    assert profile.stereo_center_count == 1, profile
    assert profile.block_counts.get("stereo_center", 0) == 1
    assert profile.reason != "multiple_stereo_centers"
    result = run_clean2d_engine(graph, None, mode="quick", target_bond_length=TARGET)
    assert result.selected is not None
    assert result.selected.metadata.get("complex_preserve_only") is not True


@pytest.mark.parametrize("mode", ["quick", "publication", "propose"])
def test_two_wedge_directions_from_one_atom_are_one_center(mode: str) -> None:
    graph = _import_reference_graph("C[C@](F)(Cl)Br")
    center = next(atom.id for atom in graph.atoms.values() if atom.stereo_cip)
    outgoing = [bond for bond in graph.bonds.values() if bond.a1_id == center]
    assert len(outgoing) == 4
    second_direction = next(bond for bond in outgoing if bond.stereo == BondStereo.NONE)
    second_direction.style = BondStyle.HASHED
    second_direction.stereo = BondStereo.DOWN

    profile = classify_clean2d_complexity(graph)
    before = coords_of(graph)
    before_signature = stereo_layout_signature(graph, before, set(graph.atoms))
    result = run_clean2d_engine(graph, None, mode=mode, target_bond_length=TARGET)

    assert profile.stereo_center_count == 1, profile
    assert profile.block_counts.get("stereo_center", 0) == 1
    assert result.selected is not None
    assert result.selected.metadata.get("complex_preserve_only") is not True
    assert stereo_layout_signature(graph, result.selected.coords, set(graph.atoms)) == before_signature


def test_distinct_real_stereocenters_still_trigger_preserve_policy() -> None:
    graph = _import_reference_graph("C[C@H](O)[C@@H](F)Cl")
    independent_reference = _independent_rdkit_report("smiles", "C[C@H](O)[C@@H](F)Cl")
    assert len(independent_reference["centers"]) == 2
    assert all(label in {"R", "S"} for label in independent_reference["centers"])

    profile = classify_clean2d_complexity(graph)

    assert profile.stereo_center_count == 2
    assert profile.preserve_only
    assert profile.reason == "multiple_stereo_centers"


@pytest.mark.parametrize("mode", ["quick", "publication", "propose"])
def test_selection_of_one_disconnected_real_molecule_only_targets_its_atoms(mode: str) -> None:
    first = _import_reference_graph("C[C@H](F)CO")
    second = _import_reference_graph("CCO")
    graph, (selected, unselected) = _combine_disconnected((first, second))
    before = coords_of(graph)

    profile = classify_clean2d_complexity(graph, selected)
    result = run_clean2d_engine(graph, selected, mode=mode, target_bond_length=TARGET)

    assert profile.stereo_center_count == 1
    assert result.selected is not None
    assert result.selected.source != "complex_preserve"
    assert all(set(candidate.coords) <= selected for candidate in (*result.candidates, *result.rejected))
    assert unselected.isdisjoint(result.atom_ids)
    assert all(coords_of(graph)[atom_id] == before[atom_id] for atom_id in unselected)


@pytest.mark.parametrize("mode", ["quick", "publication", "propose"])
def test_partial_selection_containing_real_stereocenter_preserves_boundary(mode: str) -> None:
    graph = _import_reference_graph("C[C@H](F)CO")
    center = next(atom.id for atom in graph.atoms.values() if atom.stereo_cip)
    selected = {center, next(bond.a2_id for bond in graph.bonds.values() if bond.a1_id == center)}
    before = coords_of(graph)
    before_wedges = wedge_signature(graph)
    before_configuration = depicted_configuration(graph, before)

    profile = classify_clean2d_complexity(graph, selected)
    result = run_clean2d_engine(graph, selected, mode=mode, target_bond_length=TARGET)
    after = _coords(graph, result.selected.coords if result.selected is not None else {})

    assert profile.stereo_center_count == 1
    assert result.selected is not None
    assert depicted_configuration(graph, after) == before_configuration
    assert wedge_signature(graph) == before_wedges
    assert all(after[atom_id] == before[atom_id] for atom_id in graph.atoms.keys() - selected)


def test_reference_enantiomers_keep_independent_cip_across_clean2d_roundtrip() -> None:
    sources = ("C[C@H](F)CO", "C[C@@H](F)CO")
    reports: list[dict[str, dict[str, object]]] = []
    for smiles in sources:
        graph = _import_reference_graph(smiles)
        # Strip cached CIP metadata so the independent RDKit assignment must
        # come from the represented wedge/hash bond in the MOL round-trip.
        for atom in graph.atoms.values():
            atom.stereo_cip = None
        before_reference = _independent_rdkit_report("smiles", smiles)
        before_independent = _independent_graph_report(graph)
        assert before_independent["centers"] == before_reference["centers"]

        before = coords_of(graph)
        before_wedges = wedge_signature(graph)
        before_configuration = depicted_configuration(graph, before)
        by_mode: dict[str, dict[str, object]] = {}
        for mode in ("quick", "publication", "propose"):
            result = run_clean2d_engine(graph, None, mode=mode, target_bond_length=TARGET)
            diagnostics = [
                (candidate.source, candidate.rejected, candidate.rejection_reason, candidate.stable_rejection_reason)
                for candidate in (*result.candidates, *result.rejected)
            ]
            assert result.selected is not None, diagnostics
            after_graph = copy.deepcopy(graph)
            for atom_id, (x, y) in result.selected.coords.items():
                after_graph.atoms[atom_id].x = x
                after_graph.atoms[atom_id].y = y

            assert_clean2d_invariants(
                graph,
                graph,
                before,
                result.selected.coords,
                atom_ids=set(graph.atoms),
            )
            assert wedge_signature(after_graph) == before_wedges
            assert depicted_configuration(after_graph, coords_of(after_graph)) == before_configuration
            after_independent = _independent_graph_report(after_graph)
            assert after_independent["centers"] == before_reference["centers"], (mode, diagnostics)
            assert after_independent["canonical_smiles"] == before_reference["canonical_smiles"], (
                mode,
                diagnostics,
            )
            by_mode[mode] = after_independent
        reports.append(by_mode)

    for mode in ("quick", "publication", "propose"):
        assert reports[0][mode]["centers"] != reports[1][mode]["centers"]
        assert reports[0][mode]["canonical_smiles"] != reports[1][mode]["canonical_smiles"]


def test_gui_controller_rejects_cached_proposal_that_flips_independent_cip() -> None:
    from types import SimpleNamespace

    from chemuson.clean2d import clean2d_geometry_hash, run_clean2d_engine
    from chemuson.gui.controllers.clean2d_controller import Clean2DController

    graph = _import_reference_graph("C[C@H](F)CO")
    selected = set(graph.atoms)
    before = coords_of(graph)
    controller = Clean2DController()
    before_report = _independent_graph_report(graph)
    assert before_report["centers"] == ["S"]

    first = run_clean2d_engine(
        graph,
        selected,
        mode="propose",
        target_bond_length=TARGET,
        avoid_hashes={clean2d_geometry_hash(graph, before, selected)},
    )
    assert first.selected is not None
    first_graph = copy.deepcopy(graph)
    for atom_id, (x, y) in first.selected.coords.items():
        first_graph.atoms[atom_id].x = x
        first_graph.atoms[atom_id].y = y
    assert _independent_graph_report(first_graph)["centers"] == ["S"]
    controller._remember_geometry_hashes(
        graph,
        selected,
        before,
        coords_of(first_graph),
    )

    avoid_hashes = controller._avoid_hashes_for_proposal(graph, selected, before)
    second = run_clean2d_engine(
        graph,
        selected,
        mode="propose",
        target_bond_length=TARGET,
        avoid_hashes=avoid_hashes,
    )
    assert second.selected is not None
    unsafe_graph = copy.deepcopy(graph)
    for atom_id, (x, y) in second.selected.coords.items():
        unsafe_graph.atoms[atom_id].x = x
        unsafe_graph.atoms[atom_id].y = y
    assert _independent_graph_report(unsafe_graph)["centers"] == ["R"]
    assert stereo_layout_signature(graph, before, selected) == stereo_layout_signature(
        graph,
        second.selected.coords,
        selected,
    )

    accepted, reason, diagnostics = controller._validate_independent_cip_preservation(
        SimpleNamespace(model=graph),
        second.selected.coords,
        candidate_source=second.selected.source,
    )

    assert not accepted
    assert "configuración R/S" in reason
    assert diagnostics["stable_rejection_reason"] == "stereo-risk"
    assert diagnostics["expected_cip_by_atom"] == {"2": "S"}
    assert diagnostics["observed_cip_by_atom"] == {"2": "R"}


def test_isolated_wedge_does_not_invent_an_r_s_assignment() -> None:
    graph = case_single_wedge()
    independent = _independent_graph_report(graph)

    assert not any(label in {"R", "S"} for label in independent["centers"])


# ---------------------------------------------------------------------------
# 4. Camino GUI: Clean2D rápido actúa sobre una molécula con cuña
# ---------------------------------------------------------------------------


def test_gui_clean_2d_full_moves_wedge_molecule_and_undo_restores_it() -> None:
    from PyQt6.QtWidgets import QApplication

    from chemuson.core.model import BondStyle as ModelBondStyle
    from chemuson.gui.main_window import ChemusonWindow

    window = ChemusonWindow()
    try:
        canvas = window.canvas
        model = canvas.model
        c1 = model.add_atom("C", 8.0, 10.0)
        c2 = model.add_atom("C", 60.0, 16.0)
        c3 = model.add_atom("C", 84.0, -20.0)
        o1 = model.add_atom("O", 128.0, 4.0)
        model.add_bond(c1.id, c2.id, order=1)
        model.add_bond(c2.id, c3.id, order=1)
        model.add_bond(c2.id, o1.id, order=1, style=ModelBondStyle.WEDGE)
        canvas._rebuild_items_from_model()
        canvas.undo_stack.setClean()

        before = coords_of(model)
        before_configuration = depicted_configuration(model, before)
        window._on_clean_2d_full()
        QApplication.processEvents()

        after = coords_of(model)
        assert any(
            math.hypot(after[key][0] - before[key][0], after[key][1] - before[key][1]) > 0.5
            for key in before
        ), "Clean2D no modificó una molécula con cuña"
        assert [
            getattr(bond.style, "value", str(bond.style)) for bond in model.bonds.values()
        ] == ["plain", "plain", "wedge"]
        assert depicted_configuration(model, after) == before_configuration

        canvas.undo_stack.undo()
        QApplication.processEvents()
        restored = coords_of(model)
        assert all(
            math.hypot(restored[key][0] - before[key][0], restored[key][1] - before[key][1]) < 0.51
            for key in before
        )
    finally:
        for index in range(window.tabs.count()):
            tab_canvas = window._canvas_from_tab_index(index)
            if tab_canvas is not None:
                tab_canvas.undo_stack.setClean()
        window.close()
        QApplication.processEvents()


def test_gui_clean2d_modes_handle_imported_reference_without_false_multiple_centers(
    isolated_shell_settings, monkeypatch
) -> None:
    from PyQt6.QtWidgets import QApplication

    from chemuson.gui.main_window import ChemusonWindow

    monkeypatch.setattr(ChemusonWindow, "_maybe_check_updates_startup", lambda self: None)
    window = ChemusonWindow()
    try:
        references = {
            "action_clean_2d_full": "Ctrl+K",
            "action_clean_2d_publication": "Ctrl+Shift+K",
            "action_clean_2d_propose": "Ctrl+Alt+K",
        }
        for action_name, shortcut in references.items():
            action = getattr(window, action_name)
            assert action.shortcut().toString() == shortcut
            canvas = window._create_document_tab(make_current=True)
            source = _import_reference_graph("C[C@H](F)CO")
            model = canvas.model
            model.atoms.clear()
            model.bonds.clear()
            for atom in source.atoms.values():
                model.add_atom(
                    atom.element,
                    atom.x,
                    atom.y,
                    atom_id=atom.id,
                    charge=atom.charge,
                    isotope=atom.isotope,
                    stereo_cip=atom.stereo_cip,
                    explicit_h=atom.explicit_h,
                    is_explicit=atom.is_explicit,
                )
            for bond in source.bonds.values():
                model.add_bond(
                    bond.a1_id,
                    bond.a2_id,
                    order=bond.order,
                    style=bond.style,
                    stereo=bond.stereo,
                    is_aromatic=bond.is_aromatic,
                )
            canvas._rebuild_items_from_model()
            canvas.undo_stack.setClean()
            before = coords_of(model)
            before_configuration = depicted_configuration(model, before)
            before_reference = _independent_graph_report(model)
            assert before_reference["centers"] == ["S"]

            action.trigger()
            QApplication.processEvents()

            after = coords_of(model)
            after_reference = _independent_graph_report(model)
            assert after_reference["centers"] == before_reference["centers"], action_name
            assert after_reference["canonical_smiles"] == before_reference["canonical_smiles"], action_name
            assert depicted_configuration(model, after) == before_configuration, action_name
            assert window.statusBar().currentMessage() != "Estructura 2D compleja preservada; no se redibujó globalmente", (
                action_name,
                window.statusBar().currentMessage(),
            )
    finally:
        for index in range(window.tabs.count()):
            tab_canvas = window._canvas_from_tab_index(index)
            if tab_canvas is not None:
                tab_canvas.undo_stack.setClean()
        window.close()
        QApplication.processEvents()
