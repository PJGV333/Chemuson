"""Extracción best-effort de estereoquímica avanzada.

Incluye soporte tolerante para:
- Axialidad (`R_a`/`S_a`) en enlaces/ejes marcados.
- Helicidad (`M`/`P`).
- Endo/exo y caras carbonílicas si/re cuando se dispone de anotaciones.

Si RDKit aislado está disponible, se consulta primero. Ante cualquier fallo,
se degrada de forma segura a anotaciones del grafo sin lanzar excepciones.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from .molview import MolView

if TYPE_CHECKING:  # pragma: no cover - typing only
    from .options import NameOptions

try:
    from chemuson.chemio.rdkit_safe import (
        advanced_stereo_descriptors_for_chain as safe_advanced_stereo_descriptors_for_chain,
    )
except Exception:  # pragma: no cover - dependencia opcional
    safe_advanced_stereo_descriptors_for_chain = None


def extract_advanced_stereo(
    view: MolView,
    *,
    locant_by_atom: dict[int, int] | None = None,
    chain: list[int] | None = None,
    opts: "NameOptions | None" = None,
) -> list[str]:
    """Devuelve descriptores avanzados normalizados para prefijo estereo.

    Args:
        view: Vista molecular.
        locant_by_atom: Mapeo opcional `atom_id -> locante`.
        chain: Cadena orientada (si existe), usada para mapeo/worker.
        opts: Opciones del motor para decidir aislamiento RDKit.

    Returns:
        Lista de tokens estereo avanzados (`M`, `2R_a`, `3endo`, `2si`, ...).
    """
    resolved_chain = [int(atom_id) for atom_id in (chain or []) if atom_id is not None]
    if locant_by_atom is None:
        if resolved_chain:
            locant_by_atom = {atom_id: idx + 1 for idx, atom_id in enumerate(resolved_chain)}
        else:
            ordered_atoms = sorted(view.atoms())
            locant_by_atom = {atom_id: idx + 1 for idx, atom_id in enumerate(ordered_atoms)}
            resolved_chain = list(ordered_atoms)

    descriptors: list[str] = []
    if (
        opts is not None
        and bool(getattr(opts, "rdkit_isolated", False))
        and callable(safe_advanced_stereo_descriptors_for_chain)
        and resolved_chain
    ):
        try:
            isolated = safe_advanced_stereo_descriptors_for_chain(view.graph, resolved_chain)
            descriptors.extend(_normalize_tokens(isolated))
        except Exception:
            # Degradación segura: se continúa con anotaciones locales.
            pass

    descriptors.extend(_advanced_from_annotations(view, locant_by_atom))
    return _dedupe_tokens(descriptors)


def _advanced_from_annotations(view: MolView, locant_by_atom: dict[int, int]) -> list[str]:
    """Recupera etiquetas avanzadas desde átomos/enlaces anotados."""
    graph = getattr(view, "graph", None)
    find_bond_between = getattr(graph, "find_bond_between", None)
    tokens: list[str] = []

    for atom_id, loc in sorted(locant_by_atom.items(), key=lambda item: item[1]):
        atom = view._get_atom(atom_id)  # noqa: SLF001 - acceso interno controlado.
        if atom is None:
            continue
        helicity = _normalize_helical(getattr(atom, "stereo_helical", None))
        if helicity:
            tokens.append(helicity)
        axial = _normalize_axial(getattr(atom, "stereo_axial", None))
        if axial:
            tokens.append(f"{loc}{axial}")
        face = _normalize_face(getattr(atom, "stereo_si_re", None))
        if face in {"si", "re"}:
            tokens.append(f"{loc}{face}")

    if callable(find_bond_between):
        atom_ids = list(locant_by_atom.keys())
        for i in range(len(atom_ids)):
            for j in range(i + 1, len(atom_ids)):
                a_atom = atom_ids[i]
                b_atom = atom_ids[j]
                bond = find_bond_between(a_atom, b_atom)
                if bond is None:
                    continue
                loc = min(locant_by_atom[a_atom], locant_by_atom[b_atom])
                axial = _normalize_axial(getattr(bond, "stereo_axial", None))
                if axial:
                    tokens.append(f"{loc}{axial}")
                endo_exo = _normalize_face(getattr(bond, "stereo_endo_exo", None))
                if endo_exo in {"endo", "exo"}:
                    tokens.append(f"{loc}{endo_exo}")
                helicity = _normalize_helical(getattr(bond, "stereo_helical", None))
                if helicity:
                    tokens.append(helicity)
    return _normalize_tokens(tokens)


def _normalize_tokens(tokens: list[str] | tuple[str, ...]) -> list[str]:
    """Normaliza variantes equivalentes (`Ra` -> `R_a`, etc.)."""
    normalized: list[str] = []
    for raw in tokens:
        token = str(raw or "").strip()
        if not token:
            continue
        if token in {"M", "P"}:
            normalized.append(token)
            continue

        match = re.fullmatch(r"(\d+)(.+)", token)
        if match:
            loc = match.group(1)
            tail = match.group(2)
            axial = _normalize_axial(tail)
            if axial:
                normalized.append(f"{loc}{axial}")
                continue
            face = _normalize_face(tail)
            if face in {"endo", "exo", "si", "re"}:
                normalized.append(f"{loc}{face}")
                continue

        axial_standalone = _normalize_axial(token)
        if axial_standalone:
            normalized.append(axial_standalone)
            continue
        face_standalone = _normalize_face(token)
        if face_standalone in {"endo", "exo", "si", "re"}:
            normalized.append(face_standalone)
            continue
        helicity = _normalize_helical(token)
        if helicity:
            normalized.append(helicity)
    return normalized


def _normalize_axial(value) -> str | None:
    text = str(value or "").strip()
    lowered = text.lower().replace("-", "").replace("_", "")
    if lowered == "ra":
        return "R_a"
    if lowered == "sa":
        return "S_a"
    return text if text in {"R_a", "S_a"} else None


def _normalize_helical(value) -> str | None:
    text = str(value or "").strip().upper()
    return text if text in {"M", "P"} else None


def _normalize_face(value) -> str | None:
    text = str(value or "").strip().lower()
    return text if text in {"endo", "exo", "si", "re"} else None


def _dedupe_tokens(tokens: list[str]) -> list[str]:
    """Elimina duplicados preservando orden de aparición."""
    result: list[str] = []
    seen: set[str] = set()
    for token in tokens:
        text = str(token or "").strip()
        if not text or text in seen:
            continue
        seen.add(text)
        result.append(text)
    return result

