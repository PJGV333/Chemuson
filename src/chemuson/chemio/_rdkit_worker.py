"""Worker aislado para tareas RDKit (invocado por `rdkit_safe.py`)."""

from __future__ import annotations

import json
import sys
from typing import Any


def _fail(message: str, **extra: Any) -> int:
    payload = {"ok": False, "error": message}
    payload.update(extra)
    sys.stdout.write(json.dumps(payload))
    return 0


def _bond_type_from_payload(Chem, bond_data: dict[str, Any]):
    """Convierte metadatos del enlace de Chemuson a tipo de enlace RDKit."""
    style = str(bond_data.get("style", "") or "").lower()
    if style == "coordination":
        return Chem.BondType.DATIVE
    if bool(bond_data.get("is_aromatic", False)):
        return Chem.BondType.AROMATIC
    order = int(bond_data.get("order", 1) or 1)
    if order == 2:
        return Chem.BondType.DOUBLE
    if order == 3:
        return Chem.BondType.TRIPLE
    return Chem.BondType.SINGLE


def _build_mol_from_graph_payload(Chem, request: dict[str, Any]):
    """Reconstruye un RDKit Mol desde payload serializado."""
    atoms = request.get("atoms", [])
    bonds = request.get("bonds", [])
    if not isinstance(atoms, list) or not isinstance(bonds, list):
        return None, {}, "invalid_graph_payload"

    rw = Chem.RWMol()
    id_map: dict[int, int] = {}
    periodic = Chem.GetPeriodicTable()

    for atom_data in atoms:
        atom_id = int(atom_data.get("id"))
        element = str(atom_data.get("element", "C"))
        atomic_number = 0
        try:
            atomic_number = int(periodic.GetAtomicNumber(element))
        except Exception:
            atomic_number = 0
        rd_atom = Chem.Atom(atomic_number if atomic_number > 0 else 0)
        if atomic_number <= 0:
            rd_atom.SetProp("atomLabel", element)
        rd_atom.SetFormalCharge(int(atom_data.get("formal_charge", 0) or 0))
        isotope = atom_data.get("isotope")
        if isotope not in {None, 0}:
            rd_atom.SetIsotope(int(isotope))
        radical = int(atom_data.get("radical_electrons", 0) or 0)
        if radical > 0:
            rd_atom.SetNumRadicalElectrons(radical)
        id_map[atom_id] = rw.AddAtom(rd_atom)

    for bond_data in bonds:
        a1_id = int(bond_data.get("a1_id"))
        a2_id = int(bond_data.get("a2_id"))
        if a1_id not in id_map or a2_id not in id_map:
            continue
        begin = a1_id
        end = a2_id
        if str(bond_data.get("style", "") or "").lower() == "coordination":
            donor = bond_data.get("donor_atom_id")
            try:
                donor_id = int(donor) if donor is not None else None
            except Exception:
                donor_id = None
            if donor_id == a2_id:
                begin, end = a2_id, a1_id
            elif donor_id == a1_id:
                begin, end = a1_id, a2_id
        bond_type = _bond_type_from_payload(Chem, bond_data)
        rw.AddBond(id_map[begin], id_map[end], bond_type)
        if bond_type == Chem.BondType.AROMATIC:
            rw.GetAtomWithIdx(id_map[a1_id]).SetIsAromatic(True)
            rw.GetAtomWithIdx(id_map[a2_id]).SetIsAromatic(True)
    return rw.GetMol(), id_map, ""


def _stereo_descriptors_for_chain(Chem, mol, id_map: dict[int, int], chain: list[int]) -> list[str]:
    """Extrae descriptores (R/S, E/Z) y los mapea por locante en la cadena."""
    locant_by_atom = {int(atom_id): idx + 1 for idx, atom_id in enumerate(chain)}
    idx_to_atom = {rd_idx: atom_id for atom_id, rd_idx in id_map.items()}
    descriptors: list[tuple[int, str]] = []

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    Chem.FindPotentialStereoBonds(mol)

    centers = Chem.FindMolChiralCenters(
        mol,
        includeUnassigned=False,
        useLegacyImplementation=False,
    )
    for rd_idx, label in centers:
        atom_id = idx_to_atom.get(int(rd_idx))
        if atom_id not in locant_by_atom:
            continue
        if label not in {"R", "S"}:
            continue
        loc = locant_by_atom[atom_id]
        descriptors.append((loc, f"{loc}{label}"))

    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        stereo = bond.GetStereo()
        if stereo not in {Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ}:
            continue
        a_atom = idx_to_atom.get(int(bond.GetBeginAtomIdx()))
        b_atom = idx_to_atom.get(int(bond.GetEndAtomIdx()))
        if a_atom not in locant_by_atom or b_atom not in locant_by_atom:
            continue
        loc = min(locant_by_atom[a_atom], locant_by_atom[b_atom])
        label = "E" if stereo == Chem.BondStereo.STEREOE else "Z"
        descriptors.append((loc, f"{loc}{label}"))

    descriptors.sort(key=lambda item: (item[0], item[1]))
    return [text for _loc, text in descriptors]


def _handle_text_mode(Chem, request: dict[str, Any]) -> dict[str, Any]:
    """Ruta simple para SMILES/molblock sin mapeo de locantes."""
    fmt = str(request.get("format", "smiles") or "smiles").lower()
    value = str(request.get("value", "") or "")
    if not value:
        return {"ok": False, "error": "empty_input"}
    if fmt == "molblock":
        mol = Chem.MolFromMolBlock(value, sanitize=False)
    else:
        mol = Chem.MolFromSmiles(value)
    if mol is None:
        return {"ok": False, "error": "invalid_input"}
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    Chem.FindPotentialStereoBonds(mol)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        stereo = bond.GetStereo()
        if stereo == Chem.BondStereo.STEREOE:
            bonds.append("E")
        elif stereo == Chem.BondStereo.STEREOZ:
            bonds.append("Z")
    return {
        "ok": True,
        "chiral_centers": [[int(idx), str(label)] for idx, label in centers],
        "double_bond_stereo": bonds,
    }


def main() -> int:
    try:
        request = json.loads(sys.stdin.read() or "{}")
    except Exception:
        return _fail("invalid_json")
    if not isinstance(request, dict):
        return _fail("invalid_request")

    try:
        from rdkit import Chem
    except Exception as exc:  # pragma: no cover - entorno sin rdkit
        return _fail("rdkit_unavailable", detail=str(exc))

    mode = str(request.get("mode", "graph") or "graph")
    if mode == "text":
        result = _handle_text_mode(Chem, request)
        sys.stdout.write(json.dumps(result))
        return 0

    mol, id_map, error = _build_mol_from_graph_payload(Chem, request)
    if mol is None:
        return _fail(error or "invalid_graph")
    chain = request.get("chain", [])
    if not isinstance(chain, list) or not chain:
        return _fail("missing_chain")
    try:
        descriptors = _stereo_descriptors_for_chain(Chem, mol, id_map, [int(x) for x in chain])
    except Exception as exc:
        return _fail("stereo_failed", detail=str(exc))
    sys.stdout.write(json.dumps({"ok": True, "descriptors": descriptors}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
