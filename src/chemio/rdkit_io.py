from __future__ import annotations

from typing import Dict, Tuple

from core.model import BondStyle, BondStereo, MolGraph

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
except Exception:  # pragma: no cover - optional dependency at runtime
    Chem = None
    AllChem = None
    rdMolDraw2D = None


def _require_rdkit():
    if Chem is None or AllChem is None or rdMolDraw2D is None:
        raise RuntimeError("RDKit no disponible")


def molgraph_to_rdkit(molgraph: MolGraph):
    _require_rdkit()
    rw = Chem.RWMol()
    id_map: Dict[int, int] = {}

    for atom in sorted(molgraph.atoms.values(), key=lambda a: a.id):
        rd_atom = Chem.Atom(atom.element)
        rd_atom.SetFormalCharge(atom.charge)
        if atom.isotope is not None:
            rd_atom.SetIsotope(atom.isotope)
        rd_idx = rw.AddAtom(rd_atom)
        id_map[atom.id] = rd_idx

    for bond in molgraph.bonds.values():
        if bond.is_aromatic:
            bond_type = Chem.BondType.AROMATIC
        elif bond.order == 2:
            bond_type = Chem.BondType.DOUBLE
        elif bond.order == 3:
            bond_type = Chem.BondType.TRIPLE
        else:
            bond_type = Chem.BondType.SINGLE
        rw.AddBond(id_map[bond.a1_id], id_map[bond.a2_id], bond_type)

    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for atom_id, idx in id_map.items():
        atom = molgraph.atoms[atom_id]
        conf.SetAtomPosition(idx, (atom.x, atom.y, 0.0))
    mol.AddConformer(conf, assignId=True)
    return mol


def molgraph_to_smiles(molgraph: MolGraph) -> str:
    mol = molgraph_to_rdkit(molgraph)
    return Chem.MolToSmiles(mol, canonical=True)


def molgraph_to_molfile(molgraph: MolGraph) -> str:
    mol = molgraph_to_rdkit(molgraph)
    return Chem.MolToMolBlock(mol)


def molgraph_to_svg(molgraph: MolGraph, size: Tuple[int, int] = (300, 200)) -> str:
    _require_rdkit()
    mol = molgraph_to_rdkit(molgraph)
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def molfile_to_molgraph(molfile: str) -> MolGraph:
    _require_rdkit()
    mol = Chem.MolFromMolBlock(molfile, sanitize=True)
    return rdkit_to_molgraph(mol)


def smiles_to_molgraph(smiles: str) -> MolGraph:
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    return rdkit_to_molgraph(mol)


def rdkit_to_molgraph(mol) -> MolGraph:
    _require_rdkit()
    if mol is None:
        raise ValueError("Mol inválido")
    if mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()

    graph = MolGraph()
    idx_map: Dict[int, int] = {}

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        new_atom = graph.add_atom(atom.GetSymbol(), pos.x, pos.y)
        new_atom.charge = atom.GetFormalCharge()
        if atom.GetIsotope():
            new_atom.isotope = atom.GetIsotope()
        idx_map[idx] = new_atom.id

    for bond in mol.GetBonds():
        order = 1
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            order = 2
        elif bond.GetBondType() == Chem.BondType.TRIPLE:
            order = 3
        graph.add_bond(
            idx_map[bond.GetBeginAtomIdx()],
            idx_map[bond.GetEndAtomIdx()],
            order,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=bond.GetIsAromatic(),
        )

    _scale_to_default(graph)
    return graph


def _scale_to_default(graph: MolGraph, target: float = 40.0) -> None:
    if not graph.bonds:
        return
    lengths = []
    for bond in graph.bonds.values():
        a1 = graph.atoms[bond.a1_id]
        a2 = graph.atoms[bond.a2_id]
        dx = a2.x - a1.x
        dy = a2.y - a1.y
        lengths.append((dx * dx + dy * dy) ** 0.5)
    avg = sum(lengths) / len(lengths)
    if avg <= 0:
        return
    scale = target / avg
    for atom in graph.atoms.values():
        atom.x *= scale
        atom.y *= scale
