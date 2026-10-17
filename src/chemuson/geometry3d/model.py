from __future__ import annotations

"""Modelos puros para conformaciones 3D y optimización molecular."""

from dataclasses import dataclass, field
import math
from enum import Enum
from typing import Any, Callable, Mapping


class ForceField(str, Enum):
    """Force fields soportados por backends externos."""

    UFF = "UFF"
    MMFF94 = "MMFF94"
    MMFF94S = "MMFF94s"
    GAFF = "GAFF"
    GHEMICAL = "Ghemical"


@dataclass(frozen=True)
class CoordinateSet3D:
    """Coordenadas 3D asociadas a IDs de átomos de un ``MolGraph``."""

    positions: Mapping[int, tuple[float, float, float]]
    source: str
    method: str
    energy: float | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def normalized_positions(self) -> dict[int, tuple[float, float, float]]:
        """Devuelve posiciones normalizadas a tipos primitivos."""
        return {
            int(atom_id): (float(coords[0]), float(coords[1]), float(coords[2]))
            for atom_id, coords in self.positions.items()
        }


def translated_generated_hydrogen_metadata(
    previous: CoordinateSet3D | None,
    new_positions: Mapping[Any, Any],
) -> dict[str, object]:
    """Conserva los H derivados y los traslada con su átomo documental padre."""
    if previous is None:
        return {}
    metadata = previous.metadata
    generated_atoms = metadata.get("generated_atoms", ())
    generated_bonds = metadata.get("generated_bonds", ())
    if not isinstance(generated_atoms, list) or not generated_atoms:
        return {}
    if not isinstance(generated_bonds, list):
        generated_bonds = []

    old_positions = previous.normalized_positions()
    normalized_new: dict[int, tuple[float, float, float]] = {}
    for atom_id, raw_position in new_positions.items():
        if not isinstance(raw_position, (list, tuple)) or len(raw_position) != 3:
            continue
        try:
            normalized_id = int(atom_id)
            position = (
                float(raw_position[0]),
                float(raw_position[1]),
                float(raw_position[2]),
            )
        except (TypeError, ValueError, OverflowError):
            continue
        if all(math.isfinite(component) for component in position):
            normalized_new[normalized_id] = position

    generated_ids = {
        str(item.get("id"))
        for item in generated_atoms
        if isinstance(item, Mapping) and item.get("id") is not None
    }
    parent_by_hydrogen: dict[str, int] = {}
    for bond in generated_bonds:
        if not isinstance(bond, Mapping):
            continue
        first, second = bond.get("a1_id"), bond.get("a2_id")
        if str(first) in generated_ids:
            hydrogen_id, raw_parent_id = str(first), second
        elif str(second) in generated_ids:
            hydrogen_id, raw_parent_id = str(second), first
        else:
            continue
        if raw_parent_id is None:
            continue
        try:
            parent_by_hydrogen[hydrogen_id] = int(raw_parent_id)
        except (TypeError, ValueError, OverflowError):
            continue

    translated_atoms: list[dict[str, object]] = []
    for item in generated_atoms:
        if not isinstance(item, Mapping):
            continue
        raw_coordinates = item.get("coords") or item.get("position")
        if not isinstance(raw_coordinates, (list, tuple)) or len(raw_coordinates) != 3:
            translated_atoms.append(dict(item))
            continue
        try:
            coordinates = (
                float(raw_coordinates[0]),
                float(raw_coordinates[1]),
                float(raw_coordinates[2]),
            )
        except (TypeError, ValueError, OverflowError):
            translated_atoms.append(dict(item))
            continue
        if not all(math.isfinite(value) for value in coordinates):
            translated_atoms.append(dict(item))
            continue
        parent_id = parent_by_hydrogen.get(str(item.get("id")))
        if parent_id in old_positions and parent_id in normalized_new:
            old_parent = old_positions[parent_id]
            new_parent = normalized_new[parent_id]
            coordinates = tuple(
                value + new_parent[index] - old_parent[index]
                for index, value in enumerate(coordinates)
            )
        translated_item = dict(item)
        translated_item["coords"] = list(coordinates)
        translated_atoms.append(translated_item)

    return {
        "generated_atoms": translated_atoms,
        "generated_bonds": [dict(item) for item in generated_bonds if isinstance(item, Mapping)],
    }


def transformed_generated_hydrogen_metadata(
    previous: CoordinateSet3D | None,
    transformed_parent_ids: set[int] | frozenset[int],
    transform: Callable[[tuple[float, float, float]], tuple[float, float, float]],
) -> dict[str, object]:
    """Aplica una transformación rígida a H derivados con su padre movido."""
    if previous is None:
        return {}
    metadata = previous.metadata
    generated_atoms = metadata.get("generated_atoms", ())
    generated_bonds = metadata.get("generated_bonds", ())
    if not isinstance(generated_atoms, list) or not generated_atoms:
        return {}
    if not isinstance(generated_bonds, list):
        generated_bonds = []

    generated_ids = {
        str(item.get("id"))
        for item in generated_atoms
        if isinstance(item, Mapping) and item.get("id") is not None
    }
    parent_by_hydrogen: dict[str, int] = {}
    for bond in generated_bonds:
        if not isinstance(bond, Mapping):
            continue
        first, second = bond.get("a1_id"), bond.get("a2_id")
        if str(first) in generated_ids:
            hydrogen_id, raw_parent_id = str(first), second
        elif str(second) in generated_ids:
            hydrogen_id, raw_parent_id = str(second), first
        else:
            continue
        try:
            if raw_parent_id is not None:
                parent_by_hydrogen[hydrogen_id] = int(raw_parent_id)
        except (TypeError, ValueError, OverflowError):
            continue

    transformed: list[dict[str, object]] = []
    for item in generated_atoms:
        if not isinstance(item, Mapping):
            continue
        transformed_item = dict(item)
        hydrogen_id = str(item.get("id"))
        parent_id = parent_by_hydrogen.get(hydrogen_id)
        raw_position = item.get("coords") or item.get("position")
        if (
            parent_id not in transformed_parent_ids
            or not isinstance(raw_position, (list, tuple))
            or len(raw_position) != 3
        ):
            transformed.append(transformed_item)
            continue
        try:
            position = (
                float(raw_position[0]),
                float(raw_position[1]),
                float(raw_position[2]),
            )
            if not all(math.isfinite(value) for value in position):
                transformed.append(transformed_item)
                continue
            raw_new_position = transform(position)
            new_position = (
                float(raw_new_position[0]),
                float(raw_new_position[1]),
                float(raw_new_position[2]),
            )
        except (IndexError, TypeError, ValueError, OverflowError) as exc:
            raise ValueError("invalid_generated_hydrogen_transform") from exc
        if not all(math.isfinite(value) for value in new_position):
            raise ValueError("invalid_generated_hydrogen_transform")
        transformed_item["coords"] = list(new_position)
        transformed.append(transformed_item)

    return {
        "generated_atoms": transformed,
        "generated_bonds": [dict(item) for item in generated_bonds if isinstance(item, Mapping)],
    }


@dataclass(frozen=True)
class OptimizationSettings:
    """Parámetros comunes para generación/optimización 3D."""

    forcefield: ForceField = ForceField.MMFF94
    max_iters: int = 200
    steps_per_update: int = 25
    timeout_s: float = 20.0
    convergence: float | None = None
    seed: int = 0xC0FFEE


@dataclass(frozen=True)
class OptimizationFrame:
    """Snapshot intermedio de una optimización molecular."""

    step: int
    coordinates: CoordinateSet3D
    energy: float | None = None
    converged: bool = False
    message: str = ""


@dataclass(frozen=True)
class OptimizationResult:
    """Resultado final de generación u optimización 3D."""

    coordinates: CoordinateSet3D | None
    converged: bool = False
    energy: float | None = None
    method: str = ""
    message: str = ""
    frames: tuple[OptimizationFrame, ...] = ()
    metadata: Mapping[str, object] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.coordinates is not None and bool(self.coordinates.positions)


@dataclass(frozen=True)
class SceneAtom3D:
    """Átomo listo para una escena/visor 3D, incluyendo átomos generados."""

    id: int | str
    element: str
    x: float
    y: float
    z: float
    radius: float
    color: str
    generated: bool = False


@dataclass(frozen=True)
class SceneBond3D:
    """Enlace listo para render 3D."""

    id: int | str
    a1_id: int | str
    a2_id: int | str
    order: int = 1
    is_aromatic: bool = False


@dataclass(frozen=True)
class SceneMolecule3D:
    """Representación desacoplada para un visor 3D futuro.

    ``MolGraph`` sigue siendo el modelo 2D editable. Esta escena puede incluir
    hidrógenos generados por backends externos sin mutar el documento 2D.
    """

    atoms: tuple[SceneAtom3D, ...]
    bonds: tuple[SceneBond3D, ...]
    source: str = ""


_VDW_RADII = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
}

_ELEMENT_COLORS = {
    "H": "#FFFFFF",
    "C": "#909090",
    "N": "#3050F8",
    "O": "#FF0D0D",
    "F": "#90E050",
    "P": "#FF8000",
    "S": "#FFFF30",
    "Cl": "#1FF01F",
    "Br": "#A62929",
    "I": "#940094",
}


def scene_molecule_from_graph(graph, coordinates: CoordinateSet3D) -> SceneMolecule3D:
    """Construye una escena 3D desde ``MolGraph`` y coordenadas asociadas."""
    positions = coordinates.normalized_positions()
    atoms: list[SceneAtom3D] = []
    for atom in sorted(getattr(graph, "atoms", {}).values(), key=lambda item: item.id):
        pos = positions.get(int(atom.id))
        if pos is None:
            continue
        element = str(atom.element)
        atoms.append(
            SceneAtom3D(
                id=int(atom.id),
                element=element,
                x=pos[0],
                y=pos[1],
                z=pos[2],
                radius=float(_VDW_RADII.get(element, 1.70)),
                color=str(_ELEMENT_COLORS.get(element, "#B0B0B0")),
                generated=False,
            )
        )
    atom_ids = {atom.id for atom in atoms}
    generated_atoms = coordinates.metadata.get("generated_atoms", ()) if hasattr(coordinates, "metadata") else ()
    if isinstance(generated_atoms, list):
        for item in generated_atoms:
            if not isinstance(item, Mapping):
                continue
            coords = item.get("coords") or item.get("position")
            if not isinstance(coords, (list, tuple)) or len(coords) != 3:
                continue
            element = str(item.get("element", "H"))
            try:
                position = tuple(float(component) for component in coords)
                atom_id = str(item.get("id", f"gen{len(atoms) + 1}"))
            except (TypeError, ValueError):
                continue
            if atom_id in atom_ids or not all(math.isfinite(component) for component in position):
                continue
            atoms.append(
                SceneAtom3D(
                    id=atom_id,
                    element=element,
                    x=position[0],
                    y=position[1],
                    z=position[2],
                    radius=float(_VDW_RADII.get(element, 1.20)),
                    color=str(_ELEMENT_COLORS.get(element, "#B0B0B0")),
                    generated=True,
                )
            )
            atom_ids.add(atom_id)
    bonds = list(
        SceneBond3D(
            int(bond.id),
            int(bond.a1_id),
            int(bond.a2_id),
            int(getattr(bond, "order", 1) or 1),
            bool(getattr(bond, "is_aromatic", False)),
        )
        for bond in sorted(getattr(graph, "bonds", {}).values(), key=lambda item: item.id)
        if int(bond.a1_id) in positions and int(bond.a2_id) in positions
    )
    generated_bonds = coordinates.metadata.get("generated_bonds", ()) if hasattr(coordinates, "metadata") else ()
    if isinstance(generated_bonds, list):
        bond_ids = {bond.id for bond in bonds}
        for item in generated_bonds:
            if not isinstance(item, Mapping):
                continue
            a1_id = item.get("a1_id")
            a2_id = item.get("a2_id")
            if a1_id not in atom_ids or a2_id not in atom_ids:
                continue
            bond_id = str(item.get("id", f"genB{len(bonds) + 1}"))
            if bond_id in bond_ids:
                continue
            try:
                order = max(1, int(item.get("order", 1) or 1))
            except (TypeError, ValueError):
                continue
            bonds.append(
                SceneBond3D(
                    bond_id,
                    a1_id,
                    a2_id,
                    order,
                    bool(item.get("is_aromatic", False)),
                )
            )
            bond_ids.add(bond_id)
    return SceneMolecule3D(atoms=tuple(atoms), bonds=tuple(bonds), source=coordinates.source)
