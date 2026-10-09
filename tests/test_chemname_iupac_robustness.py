"""Reference-backed regression tests for the bounded IUPAC robustness campaign."""

from __future__ import annotations

import pytest

from chemuson.chemio.rdkit_io import smiles_to_molgraph
from chemuson.chemname import NameOptions, iupac_name
from chemuson.chemname.errors import ChemNameNotSupported
from chemuson.core.model import MolGraph


def _primary_amide() -> MolGraph:
    graph = MolGraph()
    methyl = graph.add_atom("C", 0.0, 0.0)
    carbonyl = graph.add_atom("C", 1.0, 0.0)
    oxygen = graph.add_atom("O", 1.5, 0.8)
    nitrogen = graph.add_atom("N", 1.5, -0.8)
    graph.add_bond(methyl.id, carbonyl.id, order=1)
    graph.add_bond(carbonyl.id, oxygen.id, order=2)
    graph.add_bond(carbonyl.id, nitrogen.id, order=1)
    return graph


def _aryl_ketone(
    *,
    stereogenic_aryl_branch: bool = False,
    multiple_parent_attachments: bool = False,
) -> MolGraph:
    graph = MolGraph()
    ring = [graph.add_atom("C", float(index), 0.0) for index in range(6)]
    for index in range(6):
        graph.add_bond(
            ring[index].id,
            ring[(index + 1) % 6].id,
            order=1,
            is_aromatic=True,
        )

    methylene = graph.add_atom("C", -1.0, 0.0)
    carbonyl = graph.add_atom("C", -2.0, 0.0)
    oxygen = graph.add_atom("O", -2.5, 0.8)
    terminal_methyl = graph.add_atom("C", -3.0, 0.0)
    graph.add_bond(ring[0].id, methylene.id, order=1)
    if multiple_parent_attachments:
        graph.add_bond(ring[3].id, methylene.id, order=1)
    graph.add_bond(methylene.id, carbonyl.id, order=1)
    graph.add_bond(carbonyl.id, oxygen.id, order=2)
    graph.add_bond(carbonyl.id, terminal_methyl.id, order=1)

    amino = graph.add_atom("N", 4.0, 1.0)
    graph.add_bond(ring[4].id, amino.id, order=1)
    if stereogenic_aryl_branch:
        branch = graph.add_atom("C", 1.0, 1.0)
        branch.stereo_cip = "R"
        branch_methyl = graph.add_atom("C", 0.5, 2.0)
        branch_oxygen = graph.add_atom("O", 1.5, 2.0)
        graph.add_bond(ring[1].id, branch.id, order=1)
        graph.add_bond(branch.id, branch_methyl.id, order=1)
        graph.add_bond(branch.id, branch_oxygen.id, order=1)
    else:
        ring_methyl = graph.add_atom("C", 1.0, 1.0)
        graph.add_bond(ring[1].id, ring_methyl.id, order=1)
    return graph


def test_primary_amide_uses_requested_systematic_name_without_amino_prefix() -> None:
    assert iupac_name(_primary_amide()) == "ethanamide"


def test_supported_multisubstituted_phenyl_retains_groups_and_locants() -> None:
    assert iupac_name(_aryl_ketone()) == "1-(5-amino-2-methylphenyl)propan-2-one"


def test_stereogenic_aryl_branch_fails_closed_instead_of_losing_stereo() -> None:
    graph = _aryl_ketone(stereogenic_aryl_branch=True)
    assert iupac_name(graph) == "N/D"
    with pytest.raises(ChemNameNotSupported):
        iupac_name(graph, NameOptions(return_nd_on_fail=False))


def test_isotopically_modified_aryl_branch_fails_closed() -> None:
    graph = _aryl_ketone()
    methyl = next(
        atom for atom in graph.atoms.values() if atom.element == "C" and atom.x == 1.0 and atom.y == 1.0
    )
    methyl.isotope = 13
    assert iupac_name(graph) == "N/D"


def test_charged_aryl_branch_fails_closed() -> None:
    graph = _aryl_ketone()
    methyl = next(
        atom for atom in graph.atoms.values() if atom.element == "C" and atom.x == 1.0 and atom.y == 1.0
    )
    methyl.charge = 1
    assert iupac_name(graph) == "N/D"


def test_multiple_ring_to_parent_connections_fail_closed() -> None:
    assert iupac_name(_aryl_ketone(multiple_parent_attachments=True)) == "N/D"


@pytest.mark.parametrize(
    "smiles",
    [
        "CCO.CC",
        "CC([13CH3])CC",
        "CC[C@H](O)C",
        "CC[C@@H](O)C",
        "CC([CH2+])CC",
        "CC[NH3+]",
        "Cc1ccc([13CH3])cc1",
        "Nc1ccc([NH3+])cc1",
        "[13c]1ccccc1",
        "Cc1ccc(C[C@H](O)C)cc1",
        "O=C(O)c1ccc([13CH3])cc1",
    ],
    ids=[
        "disconnected-fragment",
        "isotopic-linear-branch",
        "stereo-up",
        "stereo-down",
        "charged-carbon-branch",
        "unsupported-ethylammonium",
        "isotopic-ring-substituent",
        "charged-ring-substituent",
        "isotopic-ring-parent",
        "stereogenic-ring-substituent",
        "isotopic-benzoic-acid-decoration",
    ],
)
def test_unrepresented_components_and_annotations_fail_closed(smiles: str) -> None:
    graph = smiles_to_molgraph(smiles)
    assert iupac_name(graph) == "N/D"
    with pytest.raises(ChemNameNotSupported):
        iupac_name(graph, NameOptions(return_nd_on_fail=False))
