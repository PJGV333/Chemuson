"""Worker aislado para tareas RDKit (invocado por `rdkit_safe.py`)."""

from __future__ import annotations

import json
import sys
from typing import Any

try:
    from pathlib import Path

    _SRC_ROOT = Path(__file__).resolve().parents[2]
    if str(_SRC_ROOT) not in sys.path:
        sys.path.insert(0, str(_SRC_ROOT))
    from chemuson.core.model import ATOMIC_NUMBERS as _ATOMIC_NUMBERS
except Exception:  # pragma: no cover - fallback defensivo en worker aislado
    _ATOMIC_NUMBERS = {}


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
            if bond_data.get("stereo_axial"):
                rd_bond.SetProp("_ChemusonStereoAxial", str(bond_data.get("stereo_axial")))
            if bond_data.get("stereo_endo_exo"):
                rd_bond.SetProp("_ChemusonStereoEndoExo", str(bond_data.get("stereo_endo_exo")))
            if bond_data.get("stereo_helical"):
                rd_bond.SetProp("_ChemusonStereoHelical", str(bond_data.get("stereo_helical")))
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
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    except Exception:
        centers = []
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
    return {
        "chiral_center_count": len(centers),
        "wedge_bond_count": wedge_count,
        "has_wedged_bonds": wedge_count > 0,
    }


def _handle_diagnostics_mode() -> dict[str, Any]:
    """Diagnóstico aislado de disponibilidad RDKit en el worker."""
    try:
        try:
            import rdkit

            rdkit_version = str(getattr(rdkit, "__version__", ""))
        except Exception:
            rdkit_version = ""
        payload = {
            "ok": True,
            "rdkit_version": rdkit_version,
            "python_executable": sys.executable,
            "sys_path_head": [str(item) for item in sys.path[:6]],
            "worker_file": __file__,
        }
        return payload
    except Exception as exc:
        return {
            "ok": False,
            "error": "rdkit_unavailable",
            "detail": str(exc),
            "python_executable": sys.executable,
            "worker_file": __file__,
        }


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
        Chem.SanitizeMol(mol)
        from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

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


def _handle_graph_conformer3d_mode(Chem, request: dict[str, Any]) -> dict[str, Any]:
    """Genera coordenadas 3D reales con EmbedMolecule + MMFF/UFF."""
    mol, id_map, error = _build_mol_from_graph_payload(Chem, request)
    if mol is None:
        return {"ok": False, "error": error or "invalid_graph"}
    if mol.GetNumAtoms() == 0:
        return {"ok": False, "error": "empty_graph"}
    try:
        Chem.SanitizeMol(mol)
        from rdkit.Chem import AllChem

        reverse_map = {rd_idx: atom_id for atom_id, rd_idx in id_map.items()}
        mol_h = Chem.AddHs(mol, addCoords=True)
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xC0FFEE
        params.useRandomCoords = True
        embed_code = int(AllChem.EmbedMolecule(mol_h, params))
        if embed_code != 0:
            params.useRandomCoords = True
            embed_code = int(AllChem.EmbedMolecule(mol_h, params))
        if embed_code != 0:
            return {"ok": False, "error": "embed_failed", "detail": str(embed_code)}

        method = "uff"
        energy = None
        try:
            if AllChem.MMFFHasAllMoleculeParams(mol_h):
                props = AllChem.MMFFGetMoleculeProperties(mol_h)
                ff = AllChem.MMFFGetMoleculeForceField(mol_h, props)
                method = "mmff"
            else:
                ff = AllChem.UFFGetMoleculeForceField(mol_h)
            if ff is not None:
                ff.Minimize(maxIts=200)
                energy = float(ff.CalcEnergy())
        except Exception:
            method = "embed-only"

        conf = mol_h.GetConformer()
        positions: dict[str, list[float]] = {}
        for rd_idx, atom_id in reverse_map.items():
            pos = conf.GetAtomPosition(int(rd_idx))
            positions[str(atom_id)] = [float(pos.x), float(pos.y), float(pos.z)]
        return {
            "ok": True,
            "positions": positions,
            "metadata": {"source": "rdkit", "method": method, "energy": energy},
        }
    except Exception as exc:
        return {"ok": False, "error": "conformer3d_failed", "detail": str(exc)}


def _handle_graph_optimize3d_mode(Chem, request: dict[str, Any]) -> dict[str, Any]:
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
        from rdkit.Chem import AllChem

        reverse_map = {rd_idx: atom_id for atom_id, rd_idx in id_map.items()}
        initial_coordinates = request.get("coordinates", {})
        used_initial_coordinates = _apply_initial_coordinates(Chem, mol, id_map, initial_coordinates)
        mol_h = Chem.AddHs(mol, addCoords=True)
        if not used_initial_coordinates:
            params = AllChem.ETKDGv3()
            params.randomSeed = seed
            params.useRandomCoords = True
            embed_code = int(AllChem.EmbedMolecule(mol_h, params))
            if embed_code != 0:
                return {"ok": False, "error": "embed_failed", "detail": str(embed_code)}

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
        frames: list[dict[str, Any]] = []
        converged_code = 1
        completed = 0
        while completed < max_iters:
            block = min(steps_per_update, max_iters - completed)
            converged_code = int(ff.Minimize(maxIts=block))
            completed += block
            energy = float(ff.CalcEnergy())
            frames.append(
                {
                    "step": completed,
                    "energy": energy,
                    "converged": converged_code == 0,
                    "positions": _positions_for_original_atoms(mol_h, reverse_map),
                }
            )
            if converged_code == 0:
                break
        energy = float(ff.CalcEnergy())
        conf = mol_h.GetConformer()
        positions: dict[str, list[float]] = {}
        for rd_idx, atom_id in reverse_map.items():
            pos = conf.GetAtomPosition(int(rd_idx))
            positions[str(atom_id)] = [float(pos.x), float(pos.y), float(pos.z)]
        return {
            "ok": True,
            "positions": positions,
            "metadata": {
                "source": "rdkit",
                "method": method,
                "energy": energy,
                "converged": converged_code == 0,
                "max_iters": max_iters,
                "steps_per_update": steps_per_update,
                "seed": seed,
                "used_initial_coordinates": used_initial_coordinates,
            },
            "frames": frames,
        }
    except Exception as exc:
        return {"ok": False, "error": "optimize3d_failed", "detail": str(exc)}


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


def main() -> int:
    try:
        request = json.loads(sys.stdin.read() or "{}")
    except Exception:
        return _fail("invalid_json")
    if not isinstance(request, dict):
        return _fail("invalid_request")

    mode = str(request.get("mode", "graph") or "graph")
    if mode == "diagnostics":
        sys.stdout.write(json.dumps(_handle_diagnostics_mode()))
        return 0

    try:
        from rdkit import Chem
    except Exception as exc:  # pragma: no cover - entorno sin rdkit
        return _fail("rdkit_unavailable", detail=str(exc), python_executable=sys.executable, worker_file=__file__)

    if mode == "text":
        result = _handle_text_mode(Chem, request)
        sys.stdout.write(json.dumps(result))
        return 0
    if mode == "to_molblock":
        result = _handle_to_molblock_mode(Chem, request)
        sys.stdout.write(json.dumps(result))
        return 0
    if mode == "smiles_depict_candidates":
        result = _handle_smiles_depict_candidates_mode(Chem, request)
        sys.stdout.write(json.dumps(result))
        return 0
    if mode == "graph_to_smiles":
        result = _handle_graph_to_smiles_mode(Chem, request)
        sys.stdout.write(json.dumps(result))
        return 0
    if mode == "graph_to_inchi":
        result = _handle_graph_to_inchi_mode(Chem, request)
        sys.stdout.write(json.dumps(result))
        return 0
    if mode == "graph_descriptors":
        result = _handle_graph_descriptors_mode(Chem, request)
        sys.stdout.write(json.dumps(result))
        return 0
    if mode == "graph_conformer3d":
        result = _handle_graph_conformer3d_mode(Chem, request)
        sys.stdout.write(json.dumps(result))
        return 0
    if mode == "graph_optimize3d":
        result = _handle_graph_optimize3d_mode(Chem, request)
        sys.stdout.write(json.dumps(result))
        return 0
    if mode == "graph_clean2d":
        result = _handle_graph_clean2d_mode(Chem, request)
        sys.stdout.write(json.dumps(result))
        return 0

    mol, id_map, error = _build_mol_from_graph_payload(Chem, request)
    if mol is None:
        return _fail(error or "invalid_graph")
    chain = request.get("chain", [])
    if not isinstance(chain, list) or not chain:
        return _fail("missing_chain")
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
        return _fail("stereo_failed", detail=str(exc))
    sys.stdout.write(json.dumps({"ok": True, "descriptors": descriptors}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
