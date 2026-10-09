from __future__ import annotations

import os
import subprocess
import sys

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from chemuson.chemio.rdkit_io import molfile_to_molgraph, molgraph_to_molfile, molgraph_to_smiles
from chemuson.chemio.rdkit_safe import (
    molgraph_to_smiles_isolated,
    smiles_to_molgraph_isolated,
    text_to_molblock,
)
from chemuson.clean2d import smiles_to_molgraph_best_depiction
from chemuson.core.model import BondStyle


TETRANDRINE_SMILES = "CN1CCC2=CC(=C3C=C2C1CC4=CC=C(C=C4)OC5=C(C=CC(=C5)CC6C7=C(O3)C(=C(C=C7CCN6C)OC)OC)OC)OC"
VANCOMYCIN_SMILES = "C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)[C@H]([C@H](C(=O)N[C@H](C(=O)N[C@H]5C(=O)N[C@@H]7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)[C@H](NC(=O)[C@H]([C@@H](C1=CC(=C(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)[C@@H](CC(C)C)NC)O)Cl)CO)O)O)(C)N)O"


def _import(smiles: str):
    graph, error = smiles_to_molgraph_isolated(smiles, timeout_s=12.0)
    assert error is None, error
    assert graph is not None
    return graph


def _export(graph) -> str:
    smiles, error = molgraph_to_smiles_isolated(graph, timeout_s=12.0)
    assert error is None, error
    assert smiles
    return smiles


def _rdkit_identity(smiles: str) -> tuple[str, str, tuple[str, ...], tuple[str, ...]]:
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, f"RDKit rejected {smiles!r}"
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    centers = tuple(
        sorted(
            label
            for _idx, label in Chem.FindMolChiralCenters(
                mol,
                includeUnassigned=True,
                useLegacyImplementation=False,
            )
        )
    )
    ez = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        stereo = bond.GetStereo()
        if stereo in {
            Chem.BondStereo.STEREOE,
            getattr(Chem.BondStereo, "STEREOTRANS", Chem.BondStereo.STEREOE),
        }:
            ez.append("E")
        elif stereo in {
            Chem.BondStereo.STEREOZ,
            getattr(Chem.BondStereo, "STEREOCIS", Chem.BondStereo.STEREOZ),
        }:
            ez.append("Z")
        else:
            ez.append("UNSPECIFIED")
    return canonical, formula, centers, tuple(sorted(ez))


def _assert_graph_atom_bond_identity(smiles: str, graph, *, assert_cip_metadata: bool = True) -> None:
    reference = Chem.MolFromSmiles(smiles)
    assert reference is not None
    graph_atoms = [graph.atoms[atom_id] for atom_id in sorted(graph.atoms)]
    assert len(graph_atoms) == reference.GetNumAtoms()
    for rd_atom, graph_atom in zip(reference.GetAtoms(), graph_atoms, strict=True):
        assert graph_atom.element == rd_atom.GetSymbol()
        assert graph_atom.charge == rd_atom.GetFormalCharge()
        assert (graph_atom.isotope or 0) == rd_atom.GetIsotope()

    expected_bonds = {
        (
            min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            int(round(bond.GetBondTypeAsDouble())),
            bool(bond.GetIsAromatic()),
        )
        for bond in reference.GetBonds()
    }
    graph_index = {atom.id: index for index, atom in enumerate(graph_atoms)}
    actual_bonds = {
        (
            min(graph_index[bond.a1_id], graph_index[bond.a2_id]),
            max(graph_index[bond.a1_id], graph_index[bond.a2_id]),
            int(bond.order),
            bool(bond.is_aromatic),
        )
        for bond in graph.bonds.values()
    }
    assert actual_bonds == expected_bonds

    expected_cip = {
        index: label
        for index, label in Chem.FindMolChiralCenters(
            reference,
            includeUnassigned=False,
            useLegacyImplementation=False,
        )
    }
    actual_cip = {
        index: atom.stereo_cip
        for index, atom in enumerate(graph_atoms)
        if atom.stereo_cip
    }
    if assert_cip_metadata:
        assert actual_cip == expected_cip


def _assert_smiles_roundtrip(smiles: str) -> None:
    graph = _import(smiles)
    _assert_graph_atom_bond_identity(smiles, graph)
    expected = _rdkit_identity(smiles)
    exported = _export(graph)
    assert _rdkit_identity(exported) == expected
    assert _export(graph) == exported, "ChemIO SMILES serialization is not deterministic"


def _assert_mol_sdf_roundtrip(smiles: str) -> None:
    graph = _import(smiles)
    expected = _rdkit_identity(smiles)
    molblock = molgraph_to_molfile(graph)
    assert molgraph_to_molfile(graph) == molblock, "ChemIO MOL serialization is not deterministic"
    for record in (molblock, f"{molblock.rstrip()}\n> <chemuson-case>\nstereo-roundtrip\n\n$$$$\n"):
        reimported = molfile_to_molgraph(record)
        _assert_graph_atom_bond_identity(smiles, reimported, assert_cip_metadata=False)
        if expected[2]:
            assert _wedge_hash_count(reimported) >= 1
        assert _rdkit_identity(_export(reimported)) == expected


def test_chiral_smiles_import_creates_wedge_or_hash() -> None:
    source = "C[C@H](O)F"
    graph = smiles_to_molgraph_best_depiction(source, timeout_s=12.0)
    assert _wedge_hash_count(graph) >= 1
    assert _rdkit_identity(_export(graph)) == _rdkit_identity(source)


def test_amino_acid_chiral_smiles_import_creates_wedge_or_hash() -> None:
    source = "N[C@@H](C)C(=O)O"
    graph = smiles_to_molgraph_best_depiction(source, timeout_s=12.0)
    assert _wedge_hash_count(graph) >= 1
    assert _rdkit_identity(_export(graph)) == _rdkit_identity(source)


def test_mol_stereo_bond_keeps_its_ctab_endpoint_orientation() -> None:
    source = "C[C@H](O)F"
    response = text_to_molblock(source, fmt="smiles", timeout_s=8.0)
    assert response.get("ok")
    lines = str(response.get("molblock", "")).splitlines()
    counts_index = next(index for index, line in enumerate(lines) if "V2000" in line)
    atom_count = int(lines[counts_index][:3])
    bond_count = int(lines[counts_index][3:6])
    bond_start = counts_index + 1 + atom_count
    stereo_lines = [
        line
        for line in lines[bond_start : bond_start + bond_count]
        if int(line[9:12].strip() or "0") in {1, 6}
    ]
    assert len(stereo_lines) == 1
    ctab_a1, ctab_a2 = int(stereo_lines[0][:3]), int(stereo_lines[0][3:6])

    graph = molfile_to_molgraph("\n".join(lines))
    stereo_bond = next(bond for bond in graph.bonds.values() if bond.style in {BondStyle.WEDGE, BondStyle.HASHED})
    assert (stereo_bond.a1_id, stereo_bond.a2_id) == (ctab_a1, ctab_a2)
    assert _rdkit_identity(_export(graph)) == _rdkit_identity(source)


def test_opposite_enantiomers_remain_distinct_across_smiles_roundtrip() -> None:
    pairs = (
        ("C[C@H](O)F", "C[C@@H](O)F"),
        ("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"),
    )
    for first, second in pairs:
        first_identity = _rdkit_identity(first)
        second_identity = _rdkit_identity(second)
        assert first_identity[0] != second_identity[0]
        _assert_smiles_roundtrip(first)
        _assert_smiles_roundtrip(second)


def test_neighbor_order_variants_and_multicenter_smiles_roundtrip() -> None:
    equivalent_inputs = (
        ("C[C@H](O)F", "F[C@@H](O)C"),
        ("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O"),
    )
    for first, reordered in equivalent_inputs:
        assert _rdkit_identity(first) == _rdkit_identity(reordered)
        _assert_smiles_roundtrip(first)
        _assert_smiles_roundtrip(reordered)

    _assert_smiles_roundtrip("C[C@H](O)[C@@H](F)Cl")


def test_mol_sdf_roundtrips_preserve_stereo_identity_and_atom_bond_data() -> None:
    for smiles in (
        "C[C@H](O)F",
        "C[C@@H](O)F",
        "FC=CF",
        "N[C@@H](C)C(=O)O",
        "C[C@H](O)[C@@H](F)Cl",
        "[13CH3][C@H](O)F",
        "[NH3+][C@H](C)C(=O)[O-]",
    ):
        _assert_mol_sdf_roundtrip(smiles)


def test_rdkit_direct_and_isolated_conversion_paths_preserve_stereo() -> None:
    script = r'''
from rdkit import Chem
from chemuson.chemio import rdkit_io, rdkit_safe

for source in ("C[C@H](O)F", "C[C@@H](O)F", "F/C=C/F", "F/C=C\\F"):
    reference = Chem.MolFromSmiles(source)
    expected = Chem.MolToSmiles(reference, canonical=True, isomericSmiles=True)
    graph, error = rdkit_safe.smiles_to_molgraph_isolated(source, timeout_s=8.0)
    assert error is None, error
    direct, _mapping = rdkit_io.molgraph_to_rdkit_with_map(graph)
    assert Chem.MolToSmiles(direct, canonical=True, isomericSmiles=True) == expected
    molblock = rdkit_io.molgraph_to_molfile(graph)
    imported = rdkit_io.molfile_to_molgraph(molblock)
    exported, error = rdkit_safe.molgraph_to_smiles_isolated(imported, timeout_s=8.0)
    assert error is None, error
    actual = Chem.MolFromSmiles(exported)
    assert actual is not None
    assert Chem.MolToSmiles(actual, canonical=True, isomericSmiles=True) == expected

for source in ("F/C=C/F", "F/C=C\\F"):
    source_mol = Chem.MolFromSmiles(source)
    graph = rdkit_io.rdkit_to_molgraph(source_mol)
    assert any(bond.stereo_ez for bond in graph.bonds.values())
    exported = rdkit_io.molgraph_to_smiles(graph)
    assert Chem.MolToSmiles(Chem.MolFromSmiles(exported), canonical=True, isomericSmiles=True) == Chem.MolToSmiles(source_mol, canonical=True, isomericSmiles=True)
'''
    env = os.environ.copy()
    env["CHEMUSON_ENABLE_DIRECT_RDKIT"] = "1"
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        env=env,
        timeout=45.0,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_stereo_worker_failure_is_not_silently_exported(monkeypatch) -> None:
    graph = _import("C[C@H](O)F")
    monkeypatch.setattr(
        "chemuson.chemio.rdkit_safe.molgraph_to_smiles_isolated",
        lambda _graph, timeout_s=5.0: (None, "timeout"),
    )
    with pytest.raises(RuntimeError, match="Stereo-preserving SMILES export failed"):
        molgraph_to_smiles(graph)


def test_ez_alkene_roundtrips_and_opposites_remain_distinct() -> None:
    for pair in (("F/C=C/F", "F/C=C\\F"),):
        first, second = pair
        assert _rdkit_identity(first)[0] != _rdkit_identity(second)[0]
        _assert_smiles_roundtrip(first)
        _assert_smiles_roundtrip(second)
        _assert_mol_sdf_roundtrip(first)
        _assert_mol_sdf_roundtrip(second)


def test_achiral_and_unspecified_tetrahedral_centers_stay_unassigned() -> None:
    for smiles in ("CCO", "c1ccccc1", "CC(O)F", "FC=CF"):
        graph = _import(smiles)
        assert not any(atom.stereo_cip for atom in graph.atoms.values())
        assert _wedge_hash_count(graph) == 0
        expected = _rdkit_identity(smiles)
        assert expected[2].count("?") == (1 if smiles == "CC(O)F" else 0)
        if smiles == "FC=CF":
            assert expected[3] == ("UNSPECIFIED",)
        assert _rdkit_identity(_export(graph)) == expected


def test_tetrandrine_import_does_not_invent_visual_stereo() -> None:
    reference = Chem.MolFromSmiles(TETRANDRINE_SMILES)
    assert reference is not None
    centers = Chem.FindMolChiralCenters(
        reference,
        includeUnassigned=True,
        useLegacyImplementation=False,
    )
    assert centers and all(label == "?" for _idx, label in centers)
    assert not Chem.FindMolChiralCenters(
        reference,
        includeUnassigned=False,
        useLegacyImplementation=False,
    )

    graph = smiles_to_molgraph_best_depiction(TETRANDRINE_SMILES, timeout_s=30.0)
    assert _wedge_hash_count(graph) == 0
    assert not any(atom.stereo_cip for atom in graph.atoms.values())
    assert _rdkit_identity(_export(graph)) == _rdkit_identity(TETRANDRINE_SMILES)


def test_vancomycin_import_preserves_visual_stereo() -> None:
    graph = smiles_to_molgraph_best_depiction(VANCOMYCIN_SMILES, timeout_s=30.0)
    assert _wedge_hash_count(graph) >= 4
    assert _rdkit_identity(_export(graph)) == _rdkit_identity(VANCOMYCIN_SMILES)


def _wedge_hash_count(graph) -> int:
    return sum(
        1
        for bond in graph.bonds.values()
        if bond.style in {BondStyle.WEDGE, BondStyle.HASHED}
    )
