"""Worker aislado para tareas RDKit (invocado por `rdkit_safe.py`)."""

from __future__ import annotations

import importlib
import importlib.machinery
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, TextIO

try:
    _SRC_ROOT = Path(__file__).resolve().parents[2]
    if str(_SRC_ROOT) not in sys.path:
        sys.path.insert(0, str(_SRC_ROOT))
    from chemuson.core.model import ATOMIC_NUMBERS as _ATOMIC_NUMBERS
except Exception:  # pragma: no cover - fallback defensivo en worker aislado
    _ATOMIC_NUMBERS = {}


def _fail(message: str, *, output: TextIO | None = None, **extra: Any) -> int:
    payload = {"ok": False, "error": message}
    payload.update(extra)
    stream = output if output is not None else sys.stdout
    if stream is not None:
        stream.write(json.dumps(payload))
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


def _bond_priority(Chem, bond_type: Any) -> int:
    if bond_type == Chem.BondType.AROMATIC:
        return 5
    if bond_type == Chem.BondType.TRIPLE:
        return 4
    if bond_type == Chem.BondType.DOUBLE:
        return 3
    if bond_type == Chem.BondType.SINGLE:
        return 2
    if bond_type == Chem.BondType.DATIVE:
        return 1
    return 0


def _atomic_number_for_symbol(symbol: str) -> int:
    text = str(symbol or "").strip()
    if not text:
        return 0
    if text in _ATOMIC_NUMBERS:
        return int(_ATOMIC_NUMBERS[text])
    if len(text) == 1:
        normalized = text.upper()
    else:
        normalized = text[0].upper() + text[1:].lower()
    return int(_ATOMIC_NUMBERS.get(normalized, 0) or 0)


def _build_mol_from_graph_payload(Chem, request: dict[str, Any]):
    """Reconstruye un RDKit Mol desde payload serializado."""
    atoms = request.get("atoms", [])
    bonds = request.get("bonds", [])
    if not isinstance(atoms, list) or not isinstance(bonds, list):
        return None, {}, "invalid_graph_payload"

    rw = Chem.RWMol()
    id_map: dict[int, int] = {}
    for atom_data in atoms:
        atom_id = int(atom_data.get("id"))
        element = str(atom_data.get("element", "C"))
        atomic_number = _atomic_number_for_symbol(element)
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
        if atom_data.get("stereo_axial"):
            rd_atom.SetProp("_ChemusonStereoAxial", str(atom_data.get("stereo_axial")))
        if atom_data.get("stereo_helical"):
            rd_atom.SetProp("_ChemusonStereoHelical", str(atom_data.get("stereo_helical")))
        if atom_data.get("stereo_si_re"):
            rd_atom.SetProp("_ChemusonStereoSiRe", str(atom_data.get("stereo_si_re")))
        if atom_data.get("stereo_cip"):
            rd_atom.SetProp("_ChemusonStereoCIP", str(atom_data.get("stereo_cip")))
        try:
            rd_atom.SetDoubleProp("_Chemuson2DX", float(atom_data.get("x", 0.0) or 0.0))
            rd_atom.SetDoubleProp("_Chemuson2DY", float(atom_data.get("y", 0.0) or 0.0))
        except (TypeError, ValueError):
            pass
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
        begin_idx = id_map[begin]
        end_idx = id_map[end]
        if begin_idx == end_idx:
            continue

        rd_bond = rw.GetBondBetweenAtoms(begin_idx, end_idx)
        if rd_bond is None:
            try:
                rw.AddBond(begin_idx, end_idx, bond_type)
                rd_bond = rw.GetBondBetweenAtoms(begin_idx, end_idx)
            except Exception:
                # Defensive: avoid crashing worker on RDKit internal bond errors
                continue
        elif _bond_priority(Chem, bond_type) > _bond_priority(Chem, rd_bond.GetBondType()):
            rd_bond.SetBondType(bond_type)
        if rd_bond is not None:
            bond_id = bond_data.get("id")
            if bond_id is not None:
                rd_bond.SetProp("_ChemusonBondId", str(int(bond_id)))
            stereo_ez = str(bond_data.get("stereo_ez", "") or "").strip().upper()
            if stereo_ez:
                rd_bond.SetProp("_ChemusonStereoEZ", stereo_ez)
            style = str(bond_data.get("style", "") or "").lower()
            stereo = str(bond_data.get("stereo", "") or "").lower()
            if style == "wavy" or stereo == "either":
                rd_bond.SetProp("_ChemusonAmbiguousStereo", "ambiguous_bond_mark")
            elif stereo == "up" or (stereo in {"", "none"} and style == "wedge"):
                rd_bond.SetBondDir(Chem.BondDir.BEGINWEDGE)
            elif stereo == "down" or (stereo in {"", "none"} and style == "hashed"):
                rd_bond.SetBondDir(Chem.BondDir.BEGINDASH)
            if bond_data.get("stereo_axial"):
                rd_bond.SetProp("_ChemusonStereoAxial", str(bond_data.get("stereo_axial")))
            if bond_data.get("stereo_endo_exo"):
                rd_bond.SetProp("_ChemusonStereoEndoExo", str(bond_data.get("stereo_endo_exo")))
            if bond_data.get("stereo_helical"):
                rd_bond.SetProp("_ChemusonStereoHelical", str(bond_data.get("stereo_helical")))
        if bond_type == Chem.BondType.AROMATIC:
            rw.GetAtomWithIdx(id_map[a1_id]).SetIsAromatic(True)
            rw.GetAtomWithIdx(id_map[a2_id]).SetIsAromatic(True)
    mol = rw.GetMol()
    if any("x" in atom_data or "y" in atom_data for atom_data in atoms):
        conformer = Chem.Conformer(mol.GetNumAtoms())
        conformer.Set3D(False)
        for atom_data in atoms:
            atom_id = int(atom_data.get("id"))
            conformer.SetAtomPosition(
                id_map[atom_id],
                (float(atom_data.get("x", 0.0)), float(atom_data.get("y", 0.0)), 0.0),
            )
        mol.AddConformer(conformer, assignId=True)
    try:
        _apply_graph_stereo(Chem, mol, request, id_map)
    except Exception as exc:
        return None, {}, f"stereo_assignment_failed:{exc}"
    return mol, id_map, ""


def _apply_graph_stereo(Chem, mol, request: dict[str, Any], id_map: dict[int, int]) -> None:
    """Reconstruye sólo las asignaciones estéreo declaradas en el payload."""
    try:
        Chem.AssignChiralTypesFromBondDirs(mol, replaceExistingTags=True)
    except Exception:
        pass

    for atom_data in request.get("atoms", []):
        expected = str(atom_data.get("stereo_cip", "") or "").strip().upper()
        if not expected:
            continue
        if expected not in {"R", "S"}:
            raise ValueError(f"unsupported_tetrahedral_descriptor:{expected}")
        atom_id = int(atom_data.get("id"))
        rd_atom = mol.GetAtomWithIdx(id_map[atom_id])
        matched = False
        for tag in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            rd_atom.SetChiralTag(tag)
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            actual = rd_atom.GetProp("_CIPCode") if rd_atom.HasProp("_CIPCode") else ""
            if actual == expected:
                matched = True
                break
        if not matched:
            raise ValueError(f"cannot_preserve_tetrahedral_descriptor:{atom_id}:{expected}")

    stereo_bonds = [
        bond_data
        for bond_data in request.get("bonds", [])
        if str(bond_data.get("stereo_ez", "") or "").strip()
    ]
    if stereo_bonds:
        from rdkit.Chem import rdCIPLabeler

        rdCIPLabeler.AssignCIPLabels(mol)

    for bond_data in stereo_bonds:
        expected = str(bond_data.get("stereo_ez", "") or "").strip().upper()
        if expected not in {"E", "Z"}:
            raise ValueError(f"unsupported_double_bond_descriptor:{expected}")
        a1_id = int(bond_data.get("a1_id"))
        a2_id = int(bond_data.get("a2_id"))
        begin_idx = id_map[a1_id]
        end_idx = id_map[a2_id]
        rd_bond = mol.GetBondBetweenAtoms(begin_idx, end_idx)
        if rd_bond is None or rd_bond.GetBondType() != Chem.BondType.DOUBLE:
            raise ValueError(f"ez_annotation_not_double_bond:{a1_id}:{a2_id}")

        def _highest_priority_neighbor(atom_idx: int, other_idx: int) -> int:
            candidates = [
                neighbor
                for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
                if neighbor.GetIdx() != other_idx
            ]
            if not candidates:
                raise ValueError(f"ez_annotation_has_no_substituent:{atom_idx}")
            ranked = sorted(
                candidates,
                key=lambda neighbor: int(neighbor.GetProp("_CIPRank"))
                if neighbor.HasProp("_CIPRank")
                else -1,
                reverse=True,
            )
            if len(ranked) > 1 and (
                not ranked[0].HasProp("_CIPRank")
                or ranked[0].GetProp("_CIPRank") == ranked[1].GetProp("_CIPRank")
            ):
                raise ValueError(f"ez_annotation_has_tied_substituents:{atom_idx}")
            return int(ranked[0].GetIdx())

        first = _highest_priority_neighbor(begin_idx, end_idx)
        second = _highest_priority_neighbor(end_idx, begin_idx)
        rd_bond.SetStereoAtoms(first, second)
        rd_bond.SetStereo(Chem.BondStereo.STEREOE if expected == "E" else Chem.BondStereo.STEREOZ)

    if stereo_bonds:
        Chem.SetDoubleBondNeighborDirections(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)


def _infer_wedge_tags_from_2d(Chem, mol, wedge_participants: list[tuple[int, int]]) -> None:
    """Uses RDKit's molfile parser to interpret drawn wedge/hash directions."""
    if not wedge_participants or any(
        not atom.HasProp("_Chemuson2DX") or not atom.HasProp("_Chemuson2DY")
        for atom in mol.GetAtoms()
    ):
        return
    try:
        from rdkit.Geometry import Point3D

        candidate = Chem.Mol(mol)
        conformer = Chem.Conformer(candidate.GetNumAtoms())
        for atom in candidate.GetAtoms():
            conformer.SetAtomPosition(
                int(atom.GetIdx()),
                Point3D(
                    float(atom.GetDoubleProp("_Chemuson2DX")),
                    float(atom.GetDoubleProp("_Chemuson2DY")),
                    0.0,
                ),
            )
        candidate.AddConformer(conformer, assignId=True)
        molblock = Chem.MolToMolBlock(candidate)
        parsed = Chem.MolFromMolBlock(molblock, sanitize=False, removeHs=False)
        if parsed is None:
            return
        Chem.AssignStereochemistry(parsed, cleanIt=True, force=True)
        participants = {idx for pair in wedge_participants for idx in pair}
        for rd_idx in participants:
            parsed_atom = parsed.GetAtomWithIdx(int(rd_idx))
            if parsed_atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                continue
            original_atom = mol.GetAtomWithIdx(int(rd_idx))
            if original_atom.HasProp("_ChemusonStereoCIP"):
                continue
            original_atom.SetChiralTag(parsed_atom.GetChiralTag())
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception:
        return


def _apply_documented_tetrahedral_stereo(
    Chem,
    mol,
    id_map: dict[int, int],
) -> tuple[dict[str, str], str | None]:
    """Installs explicit CIP labels and wedge/hash stereo on an RDKit molecule."""
    wedge_participants: list[tuple[int, int]] = []
    for bond in mol.GetBonds():
        if bond.HasProp("_ChemusonAmbiguousStereo"):
            bond_id = bond.GetProp("_ChemusonBondId") if bond.HasProp("_ChemusonBondId") else str(bond.GetIdx())
            return {}, f"unsupported_stereochemistry:ambiguous_bond:{bond_id}"
        if bond.GetBondDir() in {Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH}:
            wedge_participants.append(
                (int(bond.GetBeginAtomIdx()), int(bond.GetEndAtomIdx()))
            )

    _infer_wedge_tags_from_2d(Chem, mol, wedge_participants)
    try:
        Chem.AssignChiralTypesFromBondDirs(mol, replaceExistingTags=True)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception as exc:
        return {}, f"unsupported_stereochemistry:wedge_assignment:{exc}"

    wedge_centers: set[int] = set()
    for first_idx, second_idx in wedge_participants:
        candidates = [
            idx
            for idx in (first_idx, second_idx)
            if mol.GetAtomWithIdx(idx).HasProp("_ChemusonStereoCIP")
            or mol.GetAtomWithIdx(idx).GetChiralTag()
            != Chem.ChiralType.CHI_UNSPECIFIED
        ]
        if len(candidates) != 1:
            return {}, "unsupported_stereochemistry:wedge_center_unresolved"
        center_idx = candidates[0]
        center = mol.GetAtomWithIdx(center_idx)
        label = (
            center.GetProp("_CIPCode")
            if center.HasProp("_CIPCode")
            else center.GetProp("_ChemusonStereoCIP")
            if center.HasProp("_ChemusonStereoCIP")
            else ""
        )
        if label not in {"R", "S"}:
            return {}, f"unsupported_stereochemistry:wedge_cip_unavailable:{center_idx}"
        wedge_centers.add(center_idx)

    expected: dict[str, str] = {}
    for atom_id, rd_idx in id_map.items():
        atom = mol.GetAtomWithIdx(int(rd_idx))
        requested = (
            atom.GetProp("_ChemusonStereoCIP").strip().upper()
            if atom.HasProp("_ChemusonStereoCIP")
            else ""
        )
        if requested and requested not in {"R", "S"}:
            return {}, f"unsupported_stereochemistry:cip:{atom_id}:{requested}"
        if requested:
            if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                selected_tag = None
                for tag in (
                    Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                    Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
                ):
                    candidate = Chem.Mol(mol)
                    candidate.GetAtomWithIdx(int(rd_idx)).SetChiralTag(tag)
                    Chem.AssignStereochemistry(candidate, cleanIt=True, force=True)
                    candidate_atom = candidate.GetAtomWithIdx(int(rd_idx))
                    if candidate_atom.HasProp("_CIPCode") and candidate_atom.GetProp("_CIPCode") == requested:
                        selected_tag = tag
                        break
                if selected_tag is None:
                    return {}, f"unsupported_stereochemistry:unresolvable_cip:{atom_id}:{requested}"
                atom.SetChiralTag(selected_tag)
                Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            observed = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else ""
            if observed != requested:
                return {}, f"stereochemistry_conflict:atom:{atom_id}:{requested}:{observed or 'unknown'}"
            expected[str(atom_id)] = requested
        elif int(rd_idx) in wedge_centers:
            if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                return {}, f"unsupported_stereochemistry:wedge_unresolved:{atom_id}"
            observed = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else ""
            if observed not in {"R", "S"}:
                return {}, f"unsupported_stereochemistry:wedge_cip_unavailable:{atom_id}"
            expected[str(atom_id)] = observed
    return expected, None


def _apply_documented_double_bond_stereo(
    Chem,
    mol,
) -> tuple[dict[str, str], str | None]:
    """Installs E/Z tags with RDKit's high-priority substituent references."""
    expected: dict[str, str] = {}
    try:
        potential_stereo = list(Chem.FindPotentialStereo(mol))
    except Exception as exc:
        return {}, f"unsupported_stereochemistry:double_bond_analysis:{exc}"
    for bond in mol.GetBonds():
        if bond.HasProp("_ChemusonAmbiguousStereo"):
            bond_id = bond.GetProp("_ChemusonBondId") if bond.HasProp("_ChemusonBondId") else str(bond.GetIdx())
            return {}, f"unsupported_stereochemistry:ambiguous_bond:{bond_id}"
        if not bond.HasProp("_ChemusonStereoEZ"):
            continue
        target = bond.GetProp("_ChemusonStereoEZ").strip().upper()
        bond_id = bond.GetProp("_ChemusonBondId") if bond.HasProp("_ChemusonBondId") else ""
        if target not in {"E", "Z"} or not bond_id:
            return {}, f"unsupported_stereochemistry:invalid_double_bond_label:{bond_id or bond.GetIdx()}"
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            return {}, f"unsupported_stereochemistry:ez_on_non_double_bond:{bond_id}"
        stereo_info = next(
            (
                info
                for info in potential_stereo
                if int(info.centeredOn) == int(bond.GetIdx())
                and info.type == Chem.StereoType.Bond_Double
            ),
            None,
        )
        controlling = list(stereo_info.controllingAtoms) if stereo_info is not None else []
        if len(controlling) < 4:
            return {}, f"unsupported_stereochemistry:ez_references_unavailable:{bond_id}"
        first_ref, second_ref = int(controlling[0]), int(controlling[2])
        atom_count = int(mol.GetNumAtoms())
        if not (0 <= first_ref < atom_count and 0 <= second_ref < atom_count):
            return {}, f"unsupported_stereochemistry:ez_references_unavailable:{bond_id}"
        try:
            bond.SetStereoAtoms(first_ref, second_ref)
            bond.SetStereo(
                Chem.BondStereo.STEREOE if target == "E" else Chem.BondStereo.STEREOZ
            )
        except Exception as exc:
            return {}, f"unsupported_stereochemistry:ez_assignment:{bond_id}:{exc}"
        expected[bond_id] = target
    return expected, None


def _stereochemistry_metadata(
    atom_targets: dict[str, str],
    bond_targets: dict[str, str],
) -> dict[str, dict[str, str]]:
    return {"atoms": dict(atom_targets), "bonds": dict(bond_targets)}


def _conformer_stereochemistry_diagnostics(
    Chem,
    mol,
    id_map: dict[int, int],
    atom_targets: dict[str, str],
    bond_targets: dict[str, str],
) -> dict[str, Any]:
    """Records specified and unspecified potential stereo without assigning new identity."""
    id_by_index = {int(index): int(atom_id) for atom_id, index in id_map.items()}
    labels = {
        int(index): str(label)
        for index, label in Chem.FindMolChiralCenters(
            mol,
            force=True,
            includeUnassigned=True,
            includeCIP=True,
            useLegacyImplementation=False,
        )
    }
    potential_atoms: list[dict[str, Any]] = []
    potential_bonds: list[dict[str, Any]] = []
    for stereo in Chem.FindPotentialStereo(mol):
        stereo_type = str(stereo.type)
        specified = str(stereo.specified).lower() == "specified"
        if stereo_type == "Atom_Tetrahedral":
            rd_index = int(stereo.centeredOn)
            atom_id = id_by_index.get(rd_index)
            if atom_id is not None:
                potential_atoms.append(
                    {
                        "atom_id": atom_id,
                        "label": labels.get(rd_index, "?"),
                        "descriptor": str(stereo.descriptor),
                        "specified": specified,
                    }
                )
        elif stereo_type == "Bond_Double":
            bond = mol.GetBondWithIdx(int(stereo.centeredOn))
            bond_id = bond.GetProp("_ChemusonBondId") if bond.HasProp("_ChemusonBondId") else str(bond.GetIdx())
            potential_bonds.append(
                {
                    "bond_id": bond_id,
                    "descriptor": str(stereo.descriptor),
                    "specified": specified,
                }
            )
    unassigned_atoms = [item for item in potential_atoms if not item["specified"]]
    unassigned_bonds = [item for item in potential_bonds if not item["specified"]]
    return {
        "atoms": dict(atom_targets),
        "bonds": dict(bond_targets),
        "specified_atoms": dict(atom_targets),
        "specified_bonds": dict(bond_targets),
        "potential_atom_centers": potential_atoms,
        "potential_bond_stereo": potential_bonds,
        "unassigned_atom_count": len(unassigned_atoms),
        "unassigned_bond_count": len(unassigned_bonds),
    }


def _rdkit_graph_state(Chem, mol, id_map: dict[int, int]) -> dict[str, Any]:
    """Captures RDKit-recognized connectivity, charge, valence, and aromaticity."""
    id_by_index = {int(index): int(atom_id) for atom_id, index in id_map.items()}
    atom_states = [
        {
            "atom_id": id_by_index.get(int(atom.GetIdx())),
            "element": str(atom.GetSymbol()),
            "formal_charge": int(atom.GetFormalCharge()),
            "total_valence": int(atom.GetTotalValence()),
            "aromatic": bool(atom.GetIsAromatic()),
        }
        for atom in mol.GetAtoms()
    ]
    bond_states = [
        {
            "bond_id": bond.GetProp("_ChemusonBondId") if bond.HasProp("_ChemusonBondId") else str(bond.GetIdx()),
            "atom_ids": [
                id_by_index.get(int(bond.GetBeginAtomIdx())),
                id_by_index.get(int(bond.GetEndAtomIdx())),
            ],
            "order": float(bond.GetBondTypeAsDouble()),
            "aromatic": bool(bond.GetIsAromatic()),
        }
        for bond in mol.GetBonds()
    ]
    return {
        "aromatic_atom_count": sum(bool(atom["aromatic"]) for atom in atom_states),
        "aromatic_bond_count": sum(bool(bond["aromatic"]) for bond in bond_states),
        "atom_states": atom_states,
        "bond_states": bond_states,
    }


def _validate_geometry_stereochemistry(
    Chem,
    mol,
    id_map: dict[int, int],
    targets: dict[str, dict[str, str]],
) -> str | None:
    """Rejects 3D frames that invert documented R/S or E/Z descriptors."""
    atom_targets = targets.get("atoms", {})
    bond_targets = targets.get("bonds", {})
    if not atom_targets and not bond_targets:
        return None
    try:
        observed = Chem.Mol(mol)
        Chem.AssignStereochemistryFrom3D(observed, replaceExistingTags=True)
        Chem.AssignStereochemistry(observed, cleanIt=True, force=True)
    except Exception as exc:
        return f"stereochemistry_unverifiable:{exc}"

    for raw_atom_id, target in atom_targets.items():
        atom_id = int(raw_atom_id)
        rd_idx = id_map.get(atom_id)
        if rd_idx is None:
            return f"stereochemistry_unverifiable:atom:{atom_id}"
        atom = observed.GetAtomWithIdx(int(rd_idx))
        actual = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else ""
        if actual != target:
            return f"stereochemistry_inversion:atom:{atom_id}:{target}:{actual or 'unknown'}"

    if bond_targets:
        conformer = mol.GetConformer()
        for bond in mol.GetBonds():
            if not bond.HasProp("_ChemusonBondId"):
                continue
            bond_id = bond.GetProp("_ChemusonBondId")
            target = bond_targets.get(bond_id)
            if target is None:
                continue
            references = list(bond.GetStereoAtoms())
            if len(references) != 2:
                return f"stereochemistry_unverifiable:bond:{bond_id}"
            begin = conformer.GetAtomPosition(int(bond.GetBeginAtomIdx()))
            end = conformer.GetAtomPosition(int(bond.GetEndAtomIdx()))
            first = conformer.GetAtomPosition(int(references[0]))
            second = conformer.GetAtomPosition(int(references[1]))
            axis = (float(end.x - begin.x), float(end.y - begin.y), float(end.z - begin.z))
            first_vec = (float(first.x - begin.x), float(first.y - begin.y), float(first.z - begin.z))
            second_vec = (float(second.x - end.x), float(second.y - end.y), float(second.z - end.z))
            axis_sq = sum(value * value for value in axis)
            if axis_sq <= 1e-12:
                return f"stereochemistry_unverifiable:bond:{bond_id}"
            first_scale = sum(a * b for a, b in zip(first_vec, axis)) / axis_sq
            second_scale = sum(a * b for a, b in zip(second_vec, axis)) / axis_sq
            first_projected = tuple(a - first_scale * b for a, b in zip(first_vec, axis))
            second_projected = tuple(a - second_scale * b for a, b in zip(second_vec, axis))
            dot = sum(a * b for a, b in zip(first_projected, second_projected))
            magnitude = math.sqrt(
                sum(value * value for value in first_projected)
                * sum(value * value for value in second_projected)
            )
            if magnitude <= 1e-12 or abs(dot) <= magnitude * 1e-6:
                return f"stereochemistry_unverifiable:bond:{bond_id}"
            actual = "Z" if dot > 0 else "E"
            if actual != target:
                return f"stereochemistry_inversion:bond:{bond_id}:{target}:{actual}"
    return None


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
        e_stereo = {Chem.BondStereo.STEREOE, getattr(Chem.BondStereo, "STEREOTRANS", Chem.BondStereo.STEREOE)}
        z_stereo = {Chem.BondStereo.STEREOZ, getattr(Chem.BondStereo, "STEREOCIS", Chem.BondStereo.STEREOZ)}
        if stereo not in e_stereo | z_stereo:
            continue
        a_atom = idx_to_atom.get(int(bond.GetBeginAtomIdx()))
        b_atom = idx_to_atom.get(int(bond.GetEndAtomIdx()))
        if a_atom not in locant_by_atom or b_atom not in locant_by_atom:
            continue
        loc = min(locant_by_atom[a_atom], locant_by_atom[b_atom])
        label = "E" if stereo in e_stereo else "Z"
        descriptors.append((loc, f"{loc}{label}"))

    descriptors.sort(key=lambda item: (item[0], item[1]))
    return [text for _loc, text in descriptors]


def _normalize_axial_label(label: Any) -> str | None:
    value = str(label or "").strip()
    lowered = value.lower().replace("-", "").replace("_", "")
    if lowered in {"ra", "r"}:
        return "R_a"
    if lowered in {"sa", "s"}:
        return "S_a"
    if value in {"R_a", "S_a"}:
        return value
    return None


def _normalize_helical_label(label: Any) -> str | None:
    value = str(label or "").strip().upper()
    if value in {"M", "P"}:
        return value
    return None


def _normalize_face_label(label: Any) -> str | None:
    value = str(label or "").strip().lower()
    if value in {"si", "re", "endo", "exo"}:
        return value
    return None


def _advanced_stereo_descriptors_for_chain(
    Chem,
    mol,
    id_map: dict[int, int],
    chain: list[int],
    request: dict[str, Any],
) -> list[str]:
    """Combina descriptores base con etiquetas avanzadas serializadas."""
    base = _stereo_descriptors_for_chain(Chem, mol, id_map, chain)
    locant_by_atom = {int(atom_id): idx + 1 for idx, atom_id in enumerate(chain)}
    advanced: list[str] = []

    for atom_data in request.get("atoms", []):
        try:
            atom_id = int(atom_data.get("id"))
        except Exception:
            continue
        loc = locant_by_atom.get(atom_id)
        if loc is None:
            continue
        helicity = _normalize_helical_label(atom_data.get("stereo_helical"))
        if helicity:
            advanced.append(helicity)
        axial = _normalize_axial_label(atom_data.get("stereo_axial"))
        if axial:
            advanced.append(f"{loc}{axial}")
        face = _normalize_face_label(atom_data.get("stereo_si_re"))
        if face in {"si", "re"}:
            advanced.append(f"{loc}{face}")

    for bond_data in request.get("bonds", []):
        try:
            a1_id = int(bond_data.get("a1_id"))
            a2_id = int(bond_data.get("a2_id"))
        except Exception:
            continue
        if a1_id not in locant_by_atom or a2_id not in locant_by_atom:
            continue
        loc = min(locant_by_atom[a1_id], locant_by_atom[a2_id])
        axial = _normalize_axial_label(bond_data.get("stereo_axial"))
        if axial:
            advanced.append(f"{loc}{axial}")
        face = _normalize_face_label(bond_data.get("stereo_endo_exo"))
        if face in {"endo", "exo"}:
            advanced.append(f"{loc}{face}")
        helicity = _normalize_helical_label(bond_data.get("stereo_helical"))
        if helicity:
            advanced.append(helicity)

    merged: list[str] = []
    seen: set[str] = set()
    for token in list(base) + advanced:
        text = str(token).strip()
        if not text or text in seen:
            continue
        seen.add(text)
        merged.append(text)
    return merged


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
    bonds: list[str] = []
    double_bond_descriptors: list[dict[str, Any]] = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        stereo = bond.GetStereo()
        if stereo in {
            Chem.BondStereo.STEREOE,
            getattr(Chem.BondStereo, "STEREOTRANS", Chem.BondStereo.STEREOE),
        }:
            label = "E"
        elif stereo in {
            Chem.BondStereo.STEREOZ,
            getattr(Chem.BondStereo, "STEREOCIS", Chem.BondStereo.STEREOZ),
        }:
            label = "Z"
        else:
            continue
        bonds.append(label)
        double_bond_descriptors.append(
            {
                "begin_atom_index": int(bond.GetBeginAtomIdx()),
                "end_atom_index": int(bond.GetEndAtomIdx()),
                "stereo": label,
            }
        )
    return {
        "ok": True,
        "chiral_centers": [[int(idx), str(label)] for idx, label in centers],
        "double_bond_stereo": bonds,
        "double_bond_descriptors": double_bond_descriptors,
        "stereo_metadata": _visual_stereo_metadata(Chem, mol),
    }


def _handle_to_molblock_mode(Chem, request: dict[str, Any]) -> dict[str, Any]:
    """Convierte entrada textual a molblock en proceso aislado."""
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
    try:
        if fmt != "molblock" and mol.GetNumConformers() == 0:
            from rdkit.Chem import AllChem

            AllChem.Compute2DCoords(mol)
        if fmt != "molblock":
            _assign_visual_stereo(Chem, mol)
        molblock = Chem.MolToMolBlock(mol)
    except Exception as exc:
        return {"ok": False, "error": "molblock_failed", "detail": str(exc)}
    return {"ok": True, "molblock": molblock, "metadata": _visual_stereo_metadata(Chem, mol)}


def _handle_smiles_depict_candidates_mode(Chem, request: dict[str, Any]) -> dict[str, Any]:
    smiles = str(request.get("smiles", "") or "").strip()
    target = float(request.get("target_bond_length", 40.0) or 40.0)
    if not smiles:
        return {"ok": False, "error": "empty_input"}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"ok": False, "error": "invalid_input"}

    candidates: list[dict[str, Any]] = []
    candidates.append(_depict_with_coordgen(Chem, mol))
    candidates.append(_depict_with_compute2d(Chem, mol))
    return {
        "ok": True,
        "candidates": candidates,
        "diagnostics": {
            "atom_count": int(mol.GetNumAtoms()),
            "bond_count": int(mol.GetNumBonds()),
            "target_bond_length": target,
        },
    }


def _depict_with_coordgen(Chem, mol) -> dict[str, Any]:
    source = "rdcoordgen"
    try:
        from rdkit.Chem import rdCoordGen

        candidate = Chem.Mol(mol)
        candidate.RemoveAllConformers()
        try:
            rdCoordGen.AddCoords(candidate)
        except AttributeError:
            rdCoordGen.AddCoords(candidate, clearConfs=True)
        if candidate.GetNumConformers() == 0:
            return {"source": source, "ok": False, "error": "no_conformer", "metadata": {"depictor": source}}
        _assign_visual_stereo(Chem, candidate)
        return {
            "source": source,
            "ok": True,
            "molblock": Chem.MolToMolBlock(candidate),
            "metadata": {"depictor": source, **_visual_stereo_metadata(Chem, candidate)},
        }
    except Exception as exc:
        return {"source": source, "ok": False, "error": "rdcoordgen_failed", "metadata": {"detail": str(exc)}}


def _depict_with_compute2d(Chem, mol) -> dict[str, Any]:
    source = "rddepictor_compute2d"
    try:
        from rdkit.Chem import AllChem

        candidate = Chem.Mol(mol)
        candidate.RemoveAllConformers()
        AllChem.Compute2DCoords(candidate)
        if candidate.GetNumConformers() == 0:
            return {"source": source, "ok": False, "error": "no_conformer", "metadata": {"depictor": source}}
        _assign_visual_stereo(Chem, candidate)
        return {
            "source": source,
            "ok": True,
            "molblock": Chem.MolToMolBlock(candidate),
            "metadata": {"depictor": source, **_visual_stereo_metadata(Chem, candidate)},
        }
    except Exception as exc:
        return {"source": source, "ok": False, "error": "compute2d_failed", "metadata": {"detail": str(exc)}}


def _assign_visual_stereo(Chem, mol) -> None:
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception:
        pass
    try:
        if mol.GetNumConformers() > 0:
            Chem.WedgeMolBonds(mol, mol.GetConformer())
    except Exception:
        pass


def _visual_stereo_metadata(Chem, mol) -> dict[str, Any]:
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    except Exception:
        centers = []
    chiral_centers = [
        {"atom_idx": int(atom_idx), "cip": str(label)}
        for atom_idx, label in centers
        if str(label) in {"R", "S"}
    ]
    tetrahedral_centers = [
        {"atom_index": int(atom_idx), "cip": str(label)}
        for atom_idx, label in centers
        if str(label).strip()
    ]
    stereo_bonds: list[dict[str, Any]] = []
    double_bond_stereo: list[dict[str, Any]] = []
    try:
        for bond in mol.GetBonds():
            stereo = bond.GetStereo()
            if stereo in {
                Chem.BondStereo.STEREOE,
                getattr(Chem.BondStereo, "STEREOTRANS", Chem.BondStereo.STEREOE),
            }:
                label = "E"
            elif stereo in {
                Chem.BondStereo.STEREOZ,
                getattr(Chem.BondStereo, "STEREOCIS", Chem.BondStereo.STEREOZ),
            }:
                label = "Z"
            else:
                continue
            double_bond_stereo.append(
                {
                    "begin_atom_index": int(bond.GetBeginAtomIdx()),
                    "end_atom_index": int(bond.GetEndAtomIdx()),
                    "stereo": label,
                }
            )
    except Exception:
        double_bond_stereo = []
    wedge_dirs = set()
    try:
        wedge_dirs = {Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH}
    except Exception:
        pass
    wedge_count = 0
    for bond in mol.GetBonds():
        try:
            if bond.GetBondDir() in wedge_dirs:
                wedge_count += 1
        except Exception:
            pass
        try:
            stereo = bond.GetStereo()
            if bond.GetBondType() == Chem.BondType.DOUBLE and stereo in {
                Chem.BondStereo.STEREOE,
                Chem.BondStereo.STEREOZ,
                getattr(Chem.BondStereo, "STEREOTRANS", Chem.BondStereo.STEREOE),
                getattr(Chem.BondStereo, "STEREOCIS", Chem.BondStereo.STEREOZ),
            }:
                stereo_bonds.append(
                    {
                        "begin_atom_idx": int(bond.GetBeginAtomIdx()),
                        "end_atom_idx": int(bond.GetEndAtomIdx()),
                        "stereo": "E"
                        if stereo in {
                            Chem.BondStereo.STEREOE,
                            getattr(Chem.BondStereo, "STEREOTRANS", Chem.BondStereo.STEREOE),
                        }
                        else "Z",
                    }
                )
        except Exception:
            pass
    return {
        "chiral_center_count": len(centers),
        "chiral_centers": chiral_centers,
        "stereo_bonds": stereo_bonds,
        "wedge_bond_count": wedge_count,
        "has_wedged_bonds": wedge_count > 0,
        "tetrahedral_centers": tetrahedral_centers,
        "double_bond_stereo": double_bond_stereo,
    }


def _rdkit_import_failure_code() -> str:
    try:
        spec = importlib.util.find_spec("rdkit")
    except Exception:
        spec = None
    return "rdkit_extension_import_failed" if spec is not None else "rdkit_unavailable"


def _handle_diagnostics_mode() -> dict[str, Any]:
    """Importa RDKit y verifica extensiones compiladas desde este worker."""
    frozen = bool(getattr(sys, "frozen", False))
    meipass = Path(str(getattr(sys, "_MEIPASS", ""))).resolve() if frozen else None
    base = {
        "python_executable": sys.executable,
        "worker_file": __file__,
        "worker_frozen": frozen,
        "worker_meipass": str(meipass) if meipass else "",
        "sys_path_head": [str(item) for item in sys.path[:6]],
    }
    try:
        import rdkit
    except Exception as exc:
        return {"ok": False, "error": _rdkit_import_failure_code(), "detail": str(exc), **base}

    module_names = (
        "rdkit.Chem.rdchem",
        "rdkit.Chem.rdMolDescriptors",
        "rdkit.Chem.rdDistGeom",
    )
    native_extensions: dict[str, dict[str, Any]] = {}
    try:
        for name in module_names:
            module = importlib.import_module(name)
            module_path = Path(str(getattr(module, "__file__", ""))).resolve()
            native = any(str(module_path).endswith(suffix) for suffix in importlib.machinery.EXTENSION_SUFFIXES)
            inside_bundle = bool(meipass and module_path.is_relative_to(meipass))
            native_extensions[name] = {
                "path": str(module_path),
                "exists": module_path.is_file(),
                "native": native,
                "inside_bundle": inside_bundle if frozen else None,
            }
    except Exception as exc:
        return {
            "ok": False,
            "error": "rdkit_extension_import_failed",
            "detail": str(exc),
            "native_extensions": native_extensions,
            **base,
        }

    invalid = [
        name
        for name, info in native_extensions.items()
        if not info["exists"] or not info["native"] or (frozen and not info["inside_bundle"])
    ]
    if invalid:
        return {
            "ok": False,
            "error": "rdkit_extension_not_packaged",
            "detail": ", ".join(invalid),
            "rdkit_version": str(getattr(rdkit, "__version__", "")),
            "native_extensions": native_extensions,
            **base,
        }
    return {
        "ok": True,
        "rdkit_version": str(getattr(rdkit, "__version__", "")),
        "native_extensions": native_extensions,
        **base,
    }


def _handle_graph_to_molblock_mode(Chem, request: dict[str, Any]) -> dict[str, Any]:
    """Convierte un grafo ChemIO a MOL preservando anotaciones stereo."""
    mol, _id_map, error = _build_mol_from_graph_payload(Chem, request)
    if mol is None:
        return {"ok": False, "error": error or "invalid_graph"}
    try:
        for stereo_info in Chem.FindPotentialStereo(mol, flagPossible=True):
            if (
                stereo_info.type == Chem.StereoType.Bond_Double
                and stereo_info.specified == Chem.StereoSpecified.Unspecified
            ):
                mol.GetBondWithIdx(int(stereo_info.centeredOn)).SetStereo(Chem.BondStereo.STEREOANY)
        return {"ok": True, "molblock": Chem.MolToMolBlock(mol)}
    except Exception as exc:
        return {"ok": False, "error": "molblock_failed", "detail": str(exc)}


def _handle_graph_to_smiles_mode(Chem, request: dict[str, Any]) -> dict[str, Any]:
    """Convierte un payload de grafo Chemuson a SMILES en proceso aislado."""
    mol, _id_map, error = _build_mol_from_graph_payload(Chem, request)
    if mol is None:
        return {"ok": False, "error": error or "invalid_graph"}
    try:
        smiles = Chem.MolToSmiles(mol, canonical=True)
    except Exception as exc:
        return {"ok": False, "error": "smiles_failed", "detail": str(exc)}
    return {"ok": True, "smiles": smiles}


def _handle_graph_to_inchi_mode(Chem, request: dict[str, Any]) -> dict[str, Any]:
    """Convierte un payload de grafo Chemuson a InChI en proceso aislado."""
    mol, _id_map, error = _build_mol_from_graph_payload(Chem, request)
    if mol is None:
        return {"ok": False, "error": error or "invalid_graph"}
    try:
        from rdkit.Chem.inchi import MolToInchi

        inchi = MolToInchi(mol)
    except Exception as exc:
        return {"ok": False, "error": "inchi_failed", "detail": str(exc)}
    return {"ok": True, "inchi": inchi}


def _handle_graph_descriptors_mode(Chem, request: dict[str, Any]) -> dict[str, Any]:
    """Calcula descriptores fisicoquímicos básicos desde payload de grafo."""
    mol, _id_map, error = _build_mol_from_graph_payload(Chem, request)
    if mol is None:
        return {"ok": False, "error": error or "invalid_graph"}
    try:
        from rdkit.Chem import rdMolDescriptors
    except Exception as exc:
        return {"ok": False, "error": "rdkit_extension_import_failed", "detail": str(exc)}
    try:
        from rdkit.Chem import Crippen, Descriptors, Lipinski
        Chem.SanitizeMol(mol)

        hbd = int(Lipinski.NumHDonors(mol))
        hba = int(Lipinski.NumHAcceptors(mol))
        molecular_weight = float(Descriptors.MolWt(mol))
        descriptors = {
            "logp": float(Crippen.MolLogP(mol)),
            "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
            "hbd": hbd,
            "hba": hba,
            "rotatable_bonds": int(Lipinski.NumRotatableBonds(mol)),
            "heavy_atoms": int(mol.GetNumHeavyAtoms()),
            "molecular_weight": molecular_weight,
            "lipinski_violations": int(
                (molecular_weight > 500.0)
                + (float(Crippen.MolLogP(mol)) > 5.0)
                + (hbd > 5)
                + (hba > 10)
            ),
        }
    except Exception as exc:
        return {"ok": False, "error": "descriptors_failed", "detail": str(exc)}
    return {"ok": True, "descriptors": descriptors}


def _generated_hydrogen_metadata(mol, conformer, original_atom_count: int, reverse_map: dict[int, int]) -> dict[str, list[dict[str, Any]]]:
    """Return non-document atoms/bonds appended by ``Chem.AddHs``."""
    generated_atoms: list[dict[str, Any]] = []
    generated_bonds: list[dict[str, Any]] = []
    generated_atom_ids: dict[int, str] = {}
    for atom in mol.GetAtoms():
        atom_index = int(atom.GetIdx())
        if atom_index < original_atom_count or atom.GetSymbol() != "H":
            continue
        position = conformer.GetAtomPosition(atom_index)
        generated_id = f"genH{atom_index}"
        generated_atom_ids[atom_index] = generated_id
        generated_atoms.append(
            {
                "id": generated_id,
                "element": "H",
                "coords": [float(position.x), float(position.y), float(position.z)],
            }
        )

    for bond in mol.GetBonds():
        start_index = int(bond.GetBeginAtomIdx())
        end_index = int(bond.GetEndAtomIdx())
        if start_index not in generated_atom_ids and end_index not in generated_atom_ids:
            continue
        start_id = reverse_map.get(start_index, generated_atom_ids.get(start_index))
        end_id = reverse_map.get(end_index, generated_atom_ids.get(end_index))
        if start_id is None or end_id is None:
            continue
        generated_bonds.append(
            {
                "id": f"genB{int(bond.GetIdx())}",
                "a1_id": start_id,
                "a2_id": end_id,
                "order": max(1, int(round(float(bond.GetBondTypeAsDouble())))),
            }
        )
    return {"generated_atoms": generated_atoms, "generated_bonds": generated_bonds}


def _handle_graph_conformer3d_mode(
    Chem,
    request: dict[str, Any],
    *,
    emit_progress: bool = False,
) -> dict[str, Any]:
    """Genera un conformero con intentos ETKDG acotados y diagnóstico por etapa."""
    from rdkit import rdBase
    from rdkit.Chem import AllChem

    started = time.perf_counter()
    stage_started = started
    stage = "payload_received"
    atoms = request.get("atoms", [])
    bonds = request.get("bonds", [])
    element_counts: dict[str, int] = {}
    if isinstance(atoms, list):
        for atom_data in atoms:
            if isinstance(atom_data, dict):
                element = str(atom_data.get("element", "?"))
                element_counts[element] = element_counts.get(element, 0) + 1
    diagnostics: dict[str, Any] = {
        "error_code": None,
        "degradation_code": None,
        "stage": stage,
        "graph": {
            "atom_count": len(atoms) if isinstance(atoms, list) else 0,
            "bond_count": len(bonds) if isinstance(bonds, list) else 0,
            "heavy_atom_count": sum(
                count for element, count in element_counts.items() if element != "H"
            ),
            "element_counts": dict(sorted(element_counts.items())),
            "formal_charge": sum(
                int(atom.get("formal_charge", 0) or 0)
                for atom in atoms
                if isinstance(atom, dict)
            ) if isinstance(atoms, list) else 0,
        },
        "rdkit": {"version": str(rdBase.rdkitVersion)},
        "return_code": None,
        "seed": int(request.get("seed", 0xC0FFEE) or 0xC0FFEE),
        "macrocycle_used": False,
        "attempts": [],
        "timings_s": {},
    }

    def enter_stage(next_stage: str, message: str, **extra: Any) -> None:
        nonlocal stage, stage_started
        now = time.perf_counter()
        diagnostics["timings_s"][stage] = max(0.0, now - stage_started)
        stage = next_stage
        stage_started = now
        diagnostics["stage"] = stage
        if emit_progress:
            _write_worker_event(
                {
                    "event": "progress",
                    "stage": stage,
                    "message": message,
                    "graph": diagnostics["graph"],
                    **extra,
                }
            )

    def fail(error_code: str, detail: str = "") -> dict[str, Any]:
        now = time.perf_counter()
        diagnostics["timings_s"][stage] = max(0.0, now - stage_started)
        diagnostics["timings_s"]["total"] = max(0.0, now - started)
        diagnostics["stage"] = stage
        diagnostics["error_code"] = error_code
        diagnostics["detail"] = str(detail or "")
        return {
            "ok": False,
            "error": error_code,
            "detail": str(detail or ""),
            "diagnostics": diagnostics,
        }

    if emit_progress:
        _write_worker_event(
            {
                "event": "progress",
                "stage": stage,
                "message": "Payload recibido; preparando grafo químico.",
                "graph": diagnostics["graph"],
            }
        )
    try:
        mol, id_map, error = _build_mol_from_graph_payload(Chem, request)
    except Exception as exc:
        return fail("invalid_graph", str(exc))
    if mol is None or not id_map:
        return fail("invalid_graph", error or "empty_graph")

    enter_stage("sanitization", "Validando valencias y sanitización química.")
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        return fail("invalid_graph", f"sanitize_failed:{exc}")
    if mol.GetNumAtoms() == 0:
        return fail("invalid_graph", "empty_graph")
    diagnostics["graph"].update(_rdkit_graph_state(Chem, mol, id_map))
    ring_sizes = [len(ring) for ring in mol.GetRingInfo().AtomRings()]
    macrocycle_used = any(size >= 8 for size in ring_sizes)
    diagnostics["graph"].update(
        {
            "ring_count": len(ring_sizes),
            "largest_ring_size": max(ring_sizes, default=0),
        }
    )
    diagnostics["macrocycle_used"] = macrocycle_used

    enter_stage("hydrogen_stereo", "Añadiendo H derivados y validando estereoquímica.")
    atom_targets, stereo_error = _apply_documented_tetrahedral_stereo(
        Chem, mol, id_map
    )
    if stereo_error:
        return fail("unsupported_stereochemistry", stereo_error)
    reverse_map = {rd_idx: atom_id for atom_id, rd_idx in id_map.items()}
    original_atom_count = int(mol.GetNumAtoms())
    try:
        mol_h = Chem.AddHs(mol, addCoords=True)
    except Exception as exc:
        return fail("invalid_graph", f"hydrogen_addition_failed:{exc}")
    bond_targets, stereo_error = _apply_documented_double_bond_stereo(Chem, mol_h)
    if stereo_error:
        return fail("unsupported_stereochemistry", stereo_error)
    stereo_targets = {"atoms": atom_targets, "bonds": bond_targets}
    stereo_diagnostics = _conformer_stereochemistry_diagnostics(
        Chem, mol_h, id_map, atom_targets, bond_targets
    )
    diagnostics["stereochemistry"] = _stereochemistry_metadata(
        atom_targets, bond_targets
    )
    diagnostics["stereochemistry_diagnostics"] = {
        key: value
        for key, value in stereo_diagnostics.items()
        if key not in {"atoms", "bonds"}
    }

    enter_stage("embedding", "Buscando una geometría ETKDG reproducible.")
    base_seed = abs(int(diagnostics["seed"])) % 2_147_483_646 + 1
    attempt_limit = 4
    embedded = False
    for attempt_index in range(attempt_limit):
        seed = (base_seed + attempt_index * 104_729 - 1) % 2_147_483_646 + 1
        random_coords = attempt_index == attempt_limit - 1
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        params.useRandomCoords = random_coords
        params.numThreads = 1
        params.maxIterations = 100
        if macrocycle_used:
            params.useMacrocycleTorsions = True
            params.useMacrocycle14config = True
        attempt_started = time.perf_counter()
        attempt_record: dict[str, Any] = {
            "index": attempt_index + 1,
            "seed": seed,
            "random_coords": random_coords,
            "macrocycle": macrocycle_used,
        }
        try:
            embed_code = int(AllChem.EmbedMolecule(mol_h, params))
        except Exception as exc:
            embed_code = -1
            attempt_record["detail"] = str(exc)
        attempt_record["return_code"] = embed_code
        attempt_record["elapsed_s"] = max(0.0, time.perf_counter() - attempt_started)
        diagnostics["attempts"].append(attempt_record)
        diagnostics["return_code"] = embed_code
        if embed_code == 0:
            embedded = True
            break
    if not embedded:
        return fail(
            "embed_failed",
            f"ETKDG agotó {attempt_limit} intentos; último retorno RDKit={diagnostics['return_code']}",
        )
    stereo_error = _validate_geometry_stereochemistry(
        Chem, mol_h, id_map, stereo_targets
    )
    if stereo_error:
        return fail("unsupported_stereochemistry", stereo_error)
    geometry_error = _validate_conformer_geometry(Chem, mol_h)
    if geometry_error:
        return fail("geometry_invalid", geometry_error)
    embedded_mol = Chem.Mol(mol_h)

    enter_stage("minimization", "Intentando minimización MMFF/UFF.")
    requested = str(request.get("forcefield", "MMFF94") or "MMFF94").upper()
    forcefield_attempts = (
        ("MMFF94S", "UFF") if requested == "MMFF94S"
        else ("UFF", "MMFF94") if requested == "UFF"
        else ("MMFF94", "UFF")
    )
    selected_forcefield = ""
    forcefield = None
    minimization_detail = ""
    for candidate in forcefield_attempts:
        try:
            if candidate.startswith("MMFF"):
                if not AllChem.MMFFHasAllMoleculeParams(mol_h):
                    minimization_detail = "missing_mmff_parameters"
                    continue
                variant = "MMFF94s" if candidate == "MMFF94S" else "MMFF94"
                properties = AllChem.MMFFGetMoleculeProperties(
                    mol_h, mmffVariant=variant
                )
                forcefield = AllChem.MMFFGetMoleculeForceField(mol_h, properties)
                selected_forcefield = variant.lower()
            else:
                forcefield = AllChem.UFFGetMoleculeForceField(mol_h)
                selected_forcefield = "uff"
            if forcefield is not None:
                break
            minimization_detail = "forcefield_unavailable"
        except Exception as exc:
            forcefield = None
            minimization_detail = f"forcefield_setup_failed:{exc}"
    energy = None
    method = "unoptimized"
    degradation_code = "missing_forcefield_parameters"
    converged = False
    if forcefield is not None:
        try:
            minimize_code = int(
                forcefield.Minimize(maxIts=max(1, int(request.get("max_iters", 200) or 200)))
            )
            candidate_energy = float(forcefield.CalcEnergy())
            if minimize_code == 0 and math.isfinite(candidate_energy):
                validation_error = _validate_geometry_stereochemistry(
                    Chem, mol_h, id_map, stereo_targets
                ) or _validate_conformer_geometry(Chem, mol_h)
                if validation_error:
                    minimization_detail = f"minimized_geometry_invalid:{validation_error}"
                    mol_h = embedded_mol
                    degradation_code = "minimization_geometry_invalid"
                else:
                    method = selected_forcefield
                    energy = candidate_energy
                    converged = True
                    degradation_code = None
            else:
                minimization_detail = (
                    "minimization_not_converged"
                    if minimize_code != 0
                    else "minimization_invalid_energy"
                )
                mol_h = embedded_mol
                degradation_code = minimization_detail
        except Exception as exc:
            minimization_detail = f"minimization_failed:{exc}"
            mol_h = embedded_mol
            degradation_code = "minimization_failed"
    if degradation_code:
        diagnostics["degradation_code"] = degradation_code
        diagnostics["minimization_detail"] = minimization_detail
    diagnostics["minimization"] = {
        "requested_forcefield": requested,
        "selected_forcefield": selected_forcefield or None,
        "converged": converged,
        "energy": energy,
        "detail": minimization_detail,
    }

    enter_stage("geometry_validation", "Verificando coordenadas y estereoquímica 3D.")
    stereo_error = _validate_geometry_stereochemistry(
        Chem, mol_h, id_map, stereo_targets
    )
    geometry_error = stereo_error or _validate_conformer_geometry(Chem, mol_h)
    if geometry_error:
        return fail("geometry_invalid", geometry_error)
    conformer = mol_h.GetConformer()
    positions: dict[str, list[float]] = {}
    for rd_idx, atom_id in reverse_map.items():
        position = conformer.GetAtomPosition(int(rd_idx))
        positions[str(atom_id)] = [
            float(position.x),
            float(position.y),
            float(position.z),
        ]
    generated = _generated_hydrogen_metadata(
        mol_h, conformer, original_atom_count, reverse_map
    )
    enter_stage("publication", "Conformero válido listo para publicar.")
    now = time.perf_counter()
    diagnostics["timings_s"][stage] = max(0.0, now - stage_started)
    diagnostics["timings_s"]["total"] = max(0.0, now - started)
    diagnostics["stage"] = "done"
    diagnostics["error_code"] = None
    metadata = {
        "source": "rdkit",
        "method": method,
        "energy": energy,
        "energy_unit": "kcal/mol" if energy is not None else None,
        "converged": converged,
        "stereochemistry": diagnostics["stereochemistry"],
        "diagnostics": diagnostics,
        **generated,
    }
    message = (
        f"{degradation_code}: {minimization_detail}"
        if degradation_code
        else ""
    )
    return {
        "ok": True,
        "positions": positions,
        "metadata": metadata,
        "diagnostics": diagnostics,
        "message": message,
    }


def _validate_conformer_geometry(Chem, mol) -> str | None:
    """Rejects non-finite, coincident, or implausibly stretched bonded coordinates."""
    if mol.GetNumConformers() == 0:
        return "missing_conformer"
    conformer = mol.GetConformer()
    for atom_index in range(mol.GetNumAtoms()):
        position = conformer.GetAtomPosition(atom_index)
        if not all(math.isfinite(float(value)) for value in (position.x, position.y, position.z)):
            return f"non_finite_coordinate:{atom_index}"
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DATIVE:
            continue
        first = conformer.GetAtomPosition(int(bond.GetBeginAtomIdx()))
        second = conformer.GetAtomPosition(int(bond.GetEndAtomIdx()))
        distance = math.dist(
            (float(first.x), float(first.y), float(first.z)),
            (float(second.x), float(second.y), float(second.z)),
        )
        if not math.isfinite(distance) or distance < 0.2 or distance > 4.0:
            return f"bond_length_out_of_range:{bond.GetIdx()}:{distance:.6g}"
    return None


def _handle_graph_optimize3d_mode(
    Chem,
    request: dict[str, Any],
    *,
    emit_frames: bool = False,
) -> dict[str, Any]:
    """Optimiza coordenadas 3D con UFF/MMFF94/MMFF94s."""
    mol, id_map, error = _build_mol_from_graph_payload(Chem, request)
    if mol is None:
        return {"ok": False, "error": error or "invalid_graph"}
    if mol.GetNumAtoms() == 0:
        return {"ok": False, "error": "empty_graph"}
    requested = str(request.get("forcefield", "MMFF94") or "MMFF94").upper()
    max_iters = max(1, int(request.get("max_iters", 200) or 200))
    steps_per_update = max(1, int(request.get("steps_per_update", max_iters) or max_iters))
    seed = int(request.get("seed", 0xC0FFEE) or 0xC0FFEE)
    try:
        Chem.SanitizeMol(mol)
        atom_targets, stereo_error = _apply_documented_tetrahedral_stereo(
            Chem, mol, id_map
        )
        if stereo_error:
            return {"ok": False, "error": stereo_error}
        from rdkit.Chem import AllChem

        reverse_map = {rd_idx: atom_id for atom_id, rd_idx in id_map.items()}
        initial_coordinates = request.get("coordinates", {})
        used_initial_coordinates = _apply_initial_coordinates(Chem, mol, id_map, initial_coordinates)
        original_atom_count = int(mol.GetNumAtoms())
        mol_h = Chem.AddHs(mol, addCoords=True)
        bond_targets, stereo_error = _apply_documented_double_bond_stereo(Chem, mol_h)
        if stereo_error:
            return {"ok": False, "error": stereo_error}
        stereo_targets = {"atoms": atom_targets, "bonds": bond_targets}
        if not used_initial_coordinates:
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            params.useRandomCoords = True
            embed_code = int(AllChem.EmbedMolecule(mol_h, params))
            if embed_code != 0:
                return {"ok": False, "error": "embed_failed", "detail": str(embed_code)}
        stereo_error = _validate_geometry_stereochemistry(
            Chem, mol_h, id_map, stereo_targets
        )
        if stereo_error:
            return {"ok": False, "error": stereo_error}

        method = requested
        ff = None
        if requested in {"MMFF94", "MMFF94S"}:
            variant = "MMFF94s" if requested == "MMFF94S" else "MMFF94"
            if not AllChem.MMFFHasAllMoleculeParams(mol_h):
                return {"ok": False, "error": "missing_mmff_parameters"}
            props = AllChem.MMFFGetMoleculeProperties(mol_h, mmffVariant=variant)
            ff = AllChem.MMFFGetMoleculeForceField(mol_h, props)
            method = variant
        elif requested == "UFF":
            ff = AllChem.UFFGetMoleculeForceField(mol_h)
            method = "UFF"
        else:
            return {"ok": False, "error": "unsupported_forcefield", "detail": requested}
        if ff is None:
            return {"ok": False, "error": "forcefield_unavailable"}
        raw_fixed_ids = request.get("fixed_atom_ids", [])
        fixed_atom_ids: set[int] = set()
        if isinstance(raw_fixed_ids, (list, tuple, set)):
            for raw_atom_id in raw_fixed_ids:
                try:
                    atom_id = int(raw_atom_id)
                    rdkit_index = id_map.get(atom_id)
                except (TypeError, ValueError):
                    continue
                if rdkit_index is not None and atom_id not in fixed_atom_ids:
                    ff.AddFixedPoint(int(rdkit_index))
                    fixed_atom_ids.add(atom_id)
        fixed_ids = sorted(fixed_atom_ids)
        constraint_summary = {
            "kind": "fixed_atom_positions",
            "atom_count": len(fixed_ids),
            "coordinate_unit": "Å",
        }
        if emit_frames:
            ff.Initialize()
        frames: list[dict[str, Any]] = []
        sequence = 0
        converged_code = 1
        completed = 0
        while completed < max_iters:
            block = min(steps_per_update, max_iters - completed)
            converged_code = int(ff.Minimize(maxIts=block))
            completed += block
            sequence += 1
            energy = float(ff.CalcEnergy())
            stereo_error = _validate_geometry_stereochemistry(
                Chem, mol_h, id_map, stereo_targets
            )
            if stereo_error:
                return {"ok": False, "error": stereo_error}
            conf = mol_h.GetConformer()
            generated = _generated_hydrogen_metadata(mol_h, conf, original_atom_count, reverse_map)
            frame = {
                "step": completed,
                "energy": energy,
                "energy_unit": "kcal/mol",
                "forcefield": method,
                "max_iters": max_iters,
                "steps_per_update": steps_per_update,
                "fixed_atom_ids": fixed_ids,
                "fixed_atom_count": len(fixed_ids),
                "constraints": constraint_summary,
                "converged": converged_code == 0,
                "positions": _positions_for_original_atoms(mol_h, reverse_map),
                "stereochemistry": _stereochemistry_metadata(atom_targets, bond_targets),
                **generated,
            }
            if emit_frames:
                _write_worker_event({"event": "frame", "sequence": sequence, **frame})
            else:
                frames.append(frame)
            if converged_code == 0:
                break
        energy = float(ff.CalcEnergy())
        conf = mol_h.GetConformer()
        positions: dict[str, list[float]] = {}
        for rd_idx, atom_id in reverse_map.items():
            pos = conf.GetAtomPosition(int(rd_idx))
            positions[str(atom_id)] = [float(pos.x), float(pos.y), float(pos.z)]
        generated = _generated_hydrogen_metadata(mol_h, conf, original_atom_count, reverse_map)
        return {
            "ok": True,
            "positions": positions,
            "metadata": {
                "source": "rdkit",
                "method": method,
                "energy": energy,
                "energy_unit": "kcal/mol",
                "forcefield": method,
                "converged": converged_code == 0,
                "max_iters": max_iters,
                "steps_per_update": steps_per_update,
                "fixed_atom_ids": fixed_ids,
                "fixed_atom_count": len(fixed_ids),
                "constraints": constraint_summary,
                "seed": seed,
                "used_initial_coordinates": used_initial_coordinates,
                "stereochemistry": _stereochemistry_metadata(atom_targets, bond_targets),
                **generated,
            },
            "frames": frames,
        }
    except Exception as exc:
        return {"ok": False, "error": "optimize3d_failed", "detail": str(exc)}


_EVENT_STREAM: TextIO | None = None


def _write_worker_event(event: dict[str, Any]) -> None:
    # En modo congelado el canal interno pasa un archivo de respuesta; los
    # eventos NDJSON deben escribirse ahí y no en stdout del proceso.
    stream = _EVENT_STREAM if _EVENT_STREAM is not None else sys.stdout
    stream.write(json.dumps(event, separators=(",", ":")) + "\n")
    stream.flush()


def _apply_initial_coordinates(Chem, mol, id_map: dict[int, int], coordinates: Any) -> bool:
    if not isinstance(coordinates, dict) or not coordinates:
        return False
    conf = Chem.Conformer(mol.GetNumAtoms())
    assigned = 0
    from rdkit.Geometry import Point3D

    for atom_id, rd_idx in id_map.items():
        coords = coordinates.get(str(atom_id), coordinates.get(atom_id))
        if not isinstance(coords, (list, tuple)) or len(coords) != 3:
            return False
        conf.SetAtomPosition(int(rd_idx), Point3D(float(coords[0]), float(coords[1]), float(coords[2])))
        assigned += 1
    if assigned != mol.GetNumAtoms():
        return False
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return True


def _handle_graph_clean2d_mode(Chem, request: dict[str, Any]) -> dict[str, Any]:
    """Genera coordenadas 2D optimizadas con Compute2DCoords.

    Preserva _chemuson_id como propiedad RDKit para que RemoveHs no
    pierda el mapeo a los IDs del grafo original.
    """
    mol, id_map, error = _build_mol_from_graph_payload(Chem, request)
    if mol is None:
        return {"ok": False, "error": error or "invalid_graph"}
    if mol.GetNumAtoms() == 0:
        return {"ok": False, "error": "empty_graph"}
    try:
        from rdkit.Chem import AllChem

        depict_method = _prepare_mol_for_2d_depiction(Chem, mol)

        # Asignar _chemuson_id a cada átomo RDKit para preservar el mapeo
        # incluso después de RemoveHs.
        for atom_id, rd_idx in id_map.items():
            mol.GetAtomWithIdx(rd_idx).SetIntProp("_chemuson_id", int(atom_id))

        reverse_map = {rd_idx: atom_id for atom_id, rd_idx in id_map.items()}

        positions: dict[str, list[float]] = {}

        mol_no_h = Chem.RemoveHs(Chem.Mol(mol), sanitize=False)
        if 0 < mol_no_h.GetNumAtoms() < mol.GetNumAtoms():
            AllChem.Compute2DCoords(mol_no_h)
            conf_no_h = mol_no_h.GetConformer()
            for atom in mol_no_h.GetAtoms():
                if not atom.HasProp("_chemuson_id"):
                    continue
                chemuson_id = int(atom.GetIntProp("_chemuson_id"))
                rd_idx = atom.GetIdx()
                pos = conf_no_h.GetAtomPosition(int(rd_idx))
                positions[str(chemuson_id)] = [float(pos.x), float(pos.y), 0.0]
        else:
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()
            for rd_idx, atom_id in reverse_map.items():
                pos = conf.GetAtomPosition(int(rd_idx))
                positions[str(atom_id)] = [float(pos.x), float(pos.y), 0.0]

        return {
            "ok": True,
            "positions": positions,
            "metadata": {"source": "rdkit", "method": f"compute2dcoords:{depict_method}"},
        }
    except Exception as exc:
        return {"ok": False, "error": "clean2d_failed", "detail": str(exc)}


def _prepare_mol_for_2d_depiction(Chem, mol) -> str:
    """Best-effort RDKit preparation for mixed aromatic drawings.

    Compute2DCoords can work with molecules that are not fully sanitized, but
    property caches and aromatic flags still need to be in a consistent enough
    state.  We try strict sanitization first, then controlled partial paths.
    """
    try:
        Chem.SanitizeMol(mol)
        return "sanitize"
    except Exception:
        pass
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    try:
        sanitize_ops = Chem.SanitizeFlags.SANITIZE_ALL
        if hasattr(Chem.SanitizeFlags, "SANITIZE_KEKULIZE"):
            sanitize_ops = sanitize_ops ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        Chem.SanitizeMol(mol, sanitizeOps=sanitize_ops)
        return "sanitize_without_kekulize"
    except Exception:
        pass
    try:
        Chem.SetAromaticity(mol)
        return "set_aromaticity_only"
    except Exception:
        return "unsanitized"


def _positions_for_original_atoms(mol_h, reverse_map: dict[int, int]) -> dict[str, list[float]]:
    conf = mol_h.GetConformer()
    positions: dict[str, list[float]] = {}
    for rd_idx, atom_id in reverse_map.items():
        pos = conf.GetAtomPosition(int(rd_idx))
        positions[str(atom_id)] = [float(pos.x), float(pos.y), float(pos.z)]
    return positions


def main(stdin: TextIO | None = None, stdout: TextIO | None = None) -> int:
    """Process one JSON request, optionally using explicit file streams when frozen."""
    input_stream = stdin if stdin is not None else sys.stdin
    output_stream = stdout if stdout is not None else sys.stdout
    global _EVENT_STREAM
    _EVENT_STREAM = output_stream
    if input_stream is None or output_stream is None:
        return _fail("worker_stdio_unavailable", output=output_stream)
    try:
        request = json.loads(input_stream.read() or "{}")
    except Exception:
        return _fail("invalid_json", output=output_stream)
    if not isinstance(request, dict):
        return _fail("invalid_request", output=output_stream)

    mode = str(request.get("mode", "graph") or "graph")
    if mode == "diagnostics":
        output_stream.write(json.dumps(_handle_diagnostics_mode()))
        return 0

    try:
        from rdkit import Chem
    except Exception as exc:  # pragma: no cover - entorno sin rdkit
        error = _rdkit_import_failure_code()
        return _fail(error, output=output_stream, detail=str(exc), python_executable=sys.executable, worker_file=__file__)

    if mode == "text":
        result = _handle_text_mode(Chem, request)
        output_stream.write(json.dumps(result))
        return 0
    if mode == "to_molblock":
        result = _handle_to_molblock_mode(Chem, request)
        output_stream.write(json.dumps(result))
        return 0
    if mode == "smiles_depict_candidates":
        result = _handle_smiles_depict_candidates_mode(Chem, request)
        output_stream.write(json.dumps(result))
        return 0
    if mode == "graph_to_molblock":
        result = _handle_graph_to_molblock_mode(Chem, request)
        output_stream.write(json.dumps(result))
        return 0
    if mode == "graph_to_smiles":
        result = _handle_graph_to_smiles_mode(Chem, request)
        output_stream.write(json.dumps(result))
        return 0
    if mode == "graph_to_inchi":
        result = _handle_graph_to_inchi_mode(Chem, request)
        output_stream.write(json.dumps(result))
        return 0
    if mode == "graph_descriptors":
        result = _handle_graph_descriptors_mode(Chem, request)
        output_stream.write(json.dumps(result))
        return 0
    if mode == "graph_conformer3d":
        result = _handle_graph_conformer3d_mode(Chem, request)
        output_stream.write(json.dumps(result))
        return 0
    if mode == "graph_conformer3d_stream":
        result = _handle_graph_conformer3d_mode(
            Chem, request, emit_progress=True
        )
        terminal_event = {
            "event": "done" if result.get("ok") else "error",
            **result,
        }
        _write_worker_event(terminal_event)
        return 0
    if mode == "graph_optimize3d":
        result = _handle_graph_optimize3d_mode(Chem, request)
        output_stream.write(json.dumps(result))
        return 0
    if mode == "graph_optimize3d_stream":
        result = _handle_graph_optimize3d_mode(Chem, request, emit_frames=True)
        if result.get("ok"):
            metadata = result.get("metadata", {})
            _write_worker_event(
                {
                    "event": "done",
                    "ok": True,
                    "method": metadata.get("method", "rdkit"),
                    "forcefield": metadata.get("forcefield", ""),
                    "energy": metadata.get("energy"),
                    "energy_unit": metadata.get("energy_unit", "kcal/mol"),
                    "converged": metadata.get("converged", False),
                    "positions": result.get("positions", {}),
                    "metadata": metadata,
                }
            )
        else:
            _write_worker_event({"event": "error", **result})
        return 0
    if mode == "graph_clean2d":
        result = _handle_graph_clean2d_mode(Chem, request)
        output_stream.write(json.dumps(result))
        return 0

    mol, id_map, error = _build_mol_from_graph_payload(Chem, request)
    if mol is None:
        return _fail(error or "invalid_graph", output=output_stream)
    chain = request.get("chain", [])
    if not isinstance(chain, list) or not chain:
        return _fail("missing_chain", output=output_stream)
    try:
        normalized_chain = [int(x) for x in chain]
        if mode == "advanced_graph":
            descriptors = _advanced_stereo_descriptors_for_chain(
                Chem,
                mol,
                id_map,
                normalized_chain,
                request,
            )
        else:
            descriptors = _stereo_descriptors_for_chain(Chem, mol, id_map, normalized_chain)
    except Exception as exc:
        return _fail("stereo_failed", output=output_stream, detail=str(exc))
    output_stream.write(json.dumps({"ok": True, "descriptors": descriptors}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
