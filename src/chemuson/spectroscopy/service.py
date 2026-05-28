from __future__ import annotations

"""Predicción espectral base con arquitectura de plugins."""

from dataclasses import dataclass, field
from typing import Protocol

from chemuson.core.model import MolGraph, bond_affects_valence


MONOISOTOPIC_MASSES = {
    "H": 1.007825032,
    "C": 12.0,
    "N": 14.003074005,
    "O": 15.99491462,
    "S": 31.972071,
    "P": 30.973761998,
    "F": 18.998403163,
    "Cl": 34.96885268,
    "Br": 78.9183376,
    "I": 126.904468,
    "Si": 27.9769265,
    "B": 11.009305,
}

_TYPICAL_VALENCE = {
    "C": 4,
    "N": 3,
    "O": 2,
    "S": 2,
    "P": 3,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
    "Si": 4,
    "B": 3,
}

_ELECTRONEGATIVE = {"N", "O", "F", "Cl", "Br", "I"}


@dataclass(frozen=True)
class ProtonNmrPeak:
    atom_id: int
    shift_ppm: float
    hydrogens: int
    environment: str
    confidence: float = 0.45


@dataclass(frozen=True)
class CarbonNmrPeak:
    atom_id: int
    shift_ppm: float
    environment: str
    confidence: float = 0.45


@dataclass(frozen=True)
class MassPeak:
    mz: float
    intensity: float
    label: str
    confidence: float = 0.70


@dataclass(frozen=True)
class SpectralPrediction:
    proton_nmr: list[ProtonNmrPeak] = field(default_factory=list)
    carbon_nmr: list[CarbonNmrPeak] = field(default_factory=list)
    mass_spectrum: list[MassPeak] = field(default_factory=list)
    source: str = "heuristic-v1"
    confidence: float = 0.45
    message: str = ""


class SpectrumPredictor(Protocol):
    name: str

    def predict(self, graph: MolGraph) -> SpectralPrediction:
        ...


_PREDICTORS: dict[str, SpectrumPredictor] = {}


def register_predictor(predictor: SpectrumPredictor) -> None:
    """Registra un predictor externo para reemplazar/extender las heurísticas."""
    _PREDICTORS[str(predictor.name)] = predictor


def predict_spectra(graph: MolGraph, *, predictor: str = "heuristic-v1") -> SpectralPrediction:
    """Predice ^1H/^13C NMR y picos MS básicos para un ``MolGraph``."""
    if predictor != "heuristic-v1":
        plugin = _PREDICTORS.get(predictor)
        if plugin is None:
            return SpectralPrediction(message=f"predictor_not_found:{predictor}")
        return plugin.predict(graph)

    try:
        from chemuson.chemio.rdkit_io import expand_abbreviations_for_calculation

        working = expand_abbreviations_for_calculation(graph)
    except Exception:
        working = graph

    proton = _predict_proton_nmr(working)
    carbon = _predict_carbon_nmr(working)
    mass = _predict_mass_spectrum(working)
    confidence = _prediction_confidence(proton, carbon, mass)
    return SpectralPrediction(proton, carbon, mass, confidence=confidence)


def _predict_proton_nmr(graph: MolGraph) -> list[ProtonNmrPeak]:
    peaks: list[ProtonNmrPeak] = []
    for atom_id, atom in sorted(graph.atoms.items()):
        hydrogens = _estimated_hydrogen_count(graph, atom_id)
        if hydrogens <= 0:
            continue
        element = str(atom.element)
        if element == "C":
            shift, env, confidence = _carbon_bound_proton_shift(graph, atom_id, hydrogens)
        elif element in {"O", "N", "S"}:
            shift, env, confidence = _hetero_bound_proton_shift(graph, atom_id, element)
        else:
            continue
        peaks.append(ProtonNmrPeak(atom_id, round(shift, 2), hydrogens, env, confidence))
    return peaks


def _predict_carbon_nmr(graph: MolGraph) -> list[CarbonNmrPeak]:
    peaks: list[CarbonNmrPeak] = []
    for atom_id, atom in sorted(graph.atoms.items()):
        if str(atom.element) != "C":
            continue
        shift, env, confidence = _carbon_shift(graph, atom_id)
        peaks.append(CarbonNmrPeak(atom_id, round(shift, 1), env, confidence))
    return peaks


def _predict_mass_spectrum(graph: MolGraph) -> list[MassPeak]:
    mass = 0.0
    for atom_id, atom in graph.atoms.items():
        element = str(atom.element)
        mass += MONOISOTOPIC_MASSES.get(element, 0.0)
        mass += _estimated_hydrogen_count(graph, atom_id) * MONOISOTOPIC_MASSES["H"]
    if mass <= 0.0:
        return []
    proton = MONOISOTOPIC_MASSES["H"]
    return [
        MassPeak(round(mass, 4), 100.0, "M+"),
        MassPeak(round(mass + proton, 4), 35.0, "[M+H]+"),
    ]


def _carbon_bound_proton_shift(graph: MolGraph, atom_id: int, hydrogens: int) -> tuple[float, str, float]:
    atom = graph.atoms[atom_id]
    if bool(getattr(atom, "is_aromatic", False)) or _has_aromatic_bond(graph, atom_id):
        return 7.2, "aromático", 0.58
    if _has_bond_order(graph, atom_id, 2):
        if _has_hetero_neighbor(graph, atom_id, {"O"}):
            return 9.8, "aldehídico", 0.56
        return 5.4, "vinílico", 0.50
    if _has_hetero_neighbor(graph, atom_id, _ELECTRONEGATIVE):
        return 3.5, "alquilo unido a heteroátomo", 0.55
    if _near_carbonyl(graph, atom_id):
        return 2.2, "alfa carbonilo", 0.53
    if hydrogens == 3:
        return 0.9, "metilo alifático", 0.48
    if hydrogens == 2:
        return 1.3, "metileno alifático", 0.48
    return 1.5, "metino alifático", 0.46


def _hetero_bound_proton_shift(graph: MolGraph, atom_id: int, element: str) -> tuple[float, str, float]:
    if element == "O":
        if _near_carbonyl_hetero(graph, atom_id):
            return 10.5, "ácido carboxílico O-H", 0.45
        return 2.5, "alcohol/fenol O-H intercambiable", 0.42
    if element == "N":
        return 3.0, "N-H intercambiable", 0.40
    return 2.2, f"{element}-H intercambiable", 0.38


def _carbon_shift(graph: MolGraph, atom_id: int) -> tuple[float, str, float]:
    if _has_bond_order(graph, atom_id, 2) and _has_hetero_neighbor(graph, atom_id, {"O", "N", "S"}):
        return 175.0, "carbonilo", 0.58
    if bool(getattr(graph.atoms[atom_id], "is_aromatic", False)) or _has_aromatic_bond(graph, atom_id):
        return 128.0, "aromático", 0.58
    if _has_bond_order(graph, atom_id, 2):
        return 115.0, "alqueno", 0.50
    if _has_hetero_neighbor(graph, atom_id, _ELECTRONEGATIVE):
        return 62.0, "C unido a heteroátomo", 0.55
    if _near_carbonyl(graph, atom_id):
        return 32.0, "alfa carbonilo", 0.52
    hydrogens = _estimated_hydrogen_count(graph, atom_id)
    if hydrogens == 3:
        return 15.0, "metilo alifático", 0.48
    if hydrogens == 2:
        return 28.0, "metileno alifático", 0.48
    return 35.0, "metino/cuaternario alifático", 0.44


def _estimated_hydrogen_count(graph: MolGraph, atom_id: int) -> int:
    atom = graph.atoms.get(atom_id)
    if atom is None:
        return 0
    element = str(atom.element)
    explicit_h = getattr(atom, "explicit_h", None)
    if explicit_h is not None:
        return max(0, int(explicit_h))
    assigned_h = graph.assigned_hydrogen_count(atom_id)
    if assigned_h:
        return assigned_h
    if element == "H":
        return 0
    valence = _TYPICAL_VALENCE.get(element)
    if valence is None:
        return 0
    bond_sum = _bond_order_sum(graph, atom_id)
    charge = int(getattr(atom, "charge", 0) or 0)
    estimated = int(round(float(valence) - bond_sum + charge))
    return max(0, estimated)


def _bond_order_sum(graph: MolGraph, atom_id: int) -> float:
    total = 0.0
    for bond in graph.bonds.values():
        if not bond_affects_valence(bond):
            continue
        if bond.a1_id != atom_id and bond.a2_id != atom_id:
            continue
        if bool(getattr(bond, "is_aromatic", False)):
            total += 1.5
        else:
            total += float(getattr(bond, "order", 1) or 1)
    return total


def _neighbor_ids(graph: MolGraph, atom_id: int) -> list[int]:
    ids: list[int] = []
    for bond in graph.bonds.values():
        if not bond_affects_valence(bond):
            continue
        if bond.a1_id == atom_id:
            ids.append(bond.a2_id)
        elif bond.a2_id == atom_id:
            ids.append(bond.a1_id)
    return ids


def _has_hetero_neighbor(graph: MolGraph, atom_id: int, elements: set[str]) -> bool:
    return any(str(graph.atoms.get(neigh_id).element) in elements for neigh_id in _neighbor_ids(graph, atom_id) if neigh_id in graph.atoms)


def _has_bond_order(graph: MolGraph, atom_id: int, order: int) -> bool:
    for bond in graph.bonds.values():
        if bond.a1_id == atom_id or bond.a2_id == atom_id:
            if int(getattr(bond, "order", 1) or 1) == int(order):
                return True
    return False


def _has_aromatic_bond(graph: MolGraph, atom_id: int) -> bool:
    return any(
        bool(getattr(bond, "is_aromatic", False))
        for bond in graph.bonds.values()
        if bond.a1_id == atom_id or bond.a2_id == atom_id
    )


def _near_carbonyl(graph: MolGraph, atom_id: int) -> bool:
    for neigh_id in _neighbor_ids(graph, atom_id):
        if str(graph.atoms.get(neigh_id).element) != "C":
            continue
        if _has_bond_order(graph, neigh_id, 2) and _has_hetero_neighbor(graph, neigh_id, {"O", "N", "S"}):
            return True
    return False


def _near_carbonyl_hetero(graph: MolGraph, atom_id: int) -> bool:
    for neigh_id in _neighbor_ids(graph, atom_id):
        if str(graph.atoms.get(neigh_id).element) != "C":
            continue
        if _has_bond_order(graph, neigh_id, 2) and _has_hetero_neighbor(graph, neigh_id, {"O"}):
            return True
    return False


def _prediction_confidence(
    proton: list[ProtonNmrPeak],
    carbon: list[CarbonNmrPeak],
    mass: list[MassPeak],
) -> float:
    values = [peak.confidence for peak in proton]
    values.extend(peak.confidence for peak in carbon)
    values.extend(peak.confidence for peak in mass)
    if not values:
        return 0.0
    return round(sum(values) / len(values), 2)
