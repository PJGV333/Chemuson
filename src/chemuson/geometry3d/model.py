from __future__ import annotations

"""Modelos puros para conformaciones 3D y optimización molecular."""

from collections import deque
from dataclasses import dataclass, field, replace
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
    parent_atom_id: int | None = None


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


def _generated_records(metadata: Mapping[str, object], key: str) -> list[Mapping[str, object]]:
    raw_items = metadata.get(key, ())
    if not isinstance(raw_items, list):
        return []
    return [item for item in raw_items if isinstance(item, Mapping)]


def coordinate_sets_interpolatable(
    start: CoordinateSet3D, end: CoordinateSet3D
) -> bool:
    """Return whether two coordinate snapshots describe the same 3D scene."""
    start_positions = start.normalized_positions()
    end_positions = end.normalized_positions()
    if not start_positions or set(start_positions) != set(end_positions):
        return False
    start_metadata, end_metadata = start.metadata, end.metadata
    start_forcefield = start_metadata.get("forcefield")
    end_forcefield = end_metadata.get("forcefield")
    if start_forcefield and end_forcefield and str(start_forcefield) != str(end_forcefield):
        return False
    for key in ("source_element_overrides", "stereochemistry"):
        if start_metadata.get(key) != end_metadata.get(key):
            return False
    start_generated = _generated_records(start_metadata, "generated_atoms")
    end_generated = _generated_records(end_metadata, "generated_atoms")
    start_hydrogen_ids = {
        str(item.get("id"))
        for item in start_generated
        if str(item.get("element", "")) == "H" and item.get("id") is not None
    }
    end_hydrogen_ids = {
        str(item.get("id"))
        for item in end_generated
        if str(item.get("element", "")) == "H" and item.get("id") is not None
    }
    start_heavy_identity = tuple(
        sorted(
            (str(item.get("id", "")), str(item.get("element", "")))
            for item in start_generated
            if str(item.get("element", "")) != "H"
        )
    )
    end_heavy_identity = tuple(
        sorted(
            (str(item.get("id", "")), str(item.get("element", "")))
            for item in end_generated
            if str(item.get("element", "")) != "H"
        )
    )
    if start_heavy_identity != end_heavy_identity:
        return False
    start_bonds = _generated_records(start_metadata, "generated_bonds")
    end_bonds = _generated_records(end_metadata, "generated_bonds")
    start_structural_bonds = tuple(
        sorted(
            tuple(str(item.get(field, "")) for field in ("id", "a1_id", "a2_id", "order", "is_aromatic"))
            for item in start_bonds
            if str(item.get("a1_id", "")) not in start_hydrogen_ids
            and str(item.get("a2_id", "")) not in start_hydrogen_ids
        )
    )
    end_structural_bonds = tuple(
        sorted(
            tuple(str(item.get(field, "")) for field in ("id", "a1_id", "a2_id", "order", "is_aromatic"))
            for item in end_bonds
            if str(item.get("a1_id", "")) not in end_hydrogen_ids
            and str(item.get("a2_id", "")) not in end_hydrogen_ids
        )
    )
    return start_structural_bonds == end_structural_bonds


def interpolate_coordinate_sets(
    start: CoordinateSet3D,
    end: CoordinateSet3D,
    fraction: float,
) -> CoordinateSet3D:
    """Create visual-only coordinates between compatible measured snapshots."""
    if not coordinate_sets_interpolatable(start, end):
        raise ValueError("coordinate_snapshots_incompatible")
    amount = float(fraction)
    if not math.isfinite(amount):
        raise ValueError("interpolation_fraction_not_finite")
    amount = min(1.0, max(0.0, amount))
    start_positions = start.normalized_positions()
    end_positions = end.normalized_positions()
    positions: dict[int, tuple[float, float, float]] = {}
    for atom_id, start_point in start_positions.items():
        end_point = end_positions[atom_id]
        positions[atom_id] = (
            start_point[0] + (end_point[0] - start_point[0]) * amount,
            start_point[1] + (end_point[1] - start_point[1]) * amount,
            start_point[2] + (end_point[2] - start_point[2]) * amount,
        )

    start_generated: dict[str, Mapping[str, object]] = {}
    raw_start_generated = start.metadata.get("generated_atoms", ())
    if isinstance(raw_start_generated, list):
        start_generated = {
            str(item.get("id")): item
            for item in raw_start_generated
            if isinstance(item, Mapping) and item.get("id") is not None
        }
    raw_end_generated = end.metadata.get("generated_atoms", ())
    generated_atoms: list[dict[str, object]] = []
    if isinstance(raw_end_generated, list):
        for item in raw_end_generated:
            if not isinstance(item, Mapping):
                continue
            atom_id = str(item.get("id", ""))
            old = start_generated.get(atom_id)
            before = (old or {}).get("coords") or (old or {}).get("position")
            after = item.get("coords") or item.get("position")
            if (
                isinstance(before, (list, tuple))
                and isinstance(after, (list, tuple))
                and len(before) == 3
                and len(after) == 3
            ):
                key = "coords" if "coords" in item else "position"
                copied = dict(item)
                copied[key] = [
                    float(before[axis])
                    + (float(after[axis]) - float(before[axis])) * amount
                    for axis in range(3)
                ]
                generated_atoms.append(copied)
            else:
                generated_atoms.append(dict(item))

    metadata = dict(end.metadata)
    if isinstance(raw_end_generated, list):
        metadata["generated_atoms"] = generated_atoms
    metadata["visual_interpolation"] = True
    metadata["visual_interpolation_from_step"] = start.metadata.get("step")
    metadata["visual_interpolation_to_step"] = end.metadata.get("step")
    return CoordinateSet3D(
        positions=positions,
        source="visual_interpolation",
        method=end.method,
        energy=None,
        metadata=metadata,
    )


def scene_with_coordinates(
    scene: SceneMolecule3D, coordinates: CoordinateSet3D
) -> SceneMolecule3D:
    """Update render-scene geometry without mutating the document graph."""
    positions = coordinates.normalized_positions()
    source_ids = {atom.id for atom in scene.atoms if isinstance(atom.id, int)}
    raw_overrides = coordinates.metadata.get("source_element_overrides", {})
    overrides: dict[int, str] = {}
    if isinstance(raw_overrides, Mapping):
        for atom_id, element in raw_overrides.items():
            try:
                overrides[int(atom_id)] = str(element)
            except (TypeError, ValueError, OverflowError):
                continue

    generated_items: dict[str, tuple[Mapping[str, object], tuple[float, float, float]]] = {}
    raw_generated_atoms = coordinates.metadata.get("generated_atoms", ())
    if isinstance(raw_generated_atoms, list):
        for item in raw_generated_atoms:
            if not isinstance(item, Mapping) or item.get("id") is None:
                continue
            raw_position = item.get("coords") or item.get("position")
            if not isinstance(raw_position, (list, tuple)) or len(raw_position) != 3:
                continue
            try:
                point = tuple(float(value) for value in raw_position)
            except (TypeError, ValueError, OverflowError):
                continue
            if all(math.isfinite(value) for value in point):
                generated_items[str(item["id"])] = (item, (point[0], point[1], point[2]))

    raw_generated_bonds = coordinates.metadata.get("generated_bonds", ())
    generated_bonds: list[Mapping[str, object]] = (
        [item for item in raw_generated_bonds if isinstance(item, Mapping)]
        if isinstance(raw_generated_bonds, list)
        else []
    )
    existing_generated_ids = {atom.id for atom in scene.atoms if atom.generated}
    atoms: list[SceneAtom3D] = []
    for atom in scene.atoms:
        if atom.generated:
            continue
        if isinstance(atom.id, int):
            point = positions.get(atom.id)
            element = overrides.get(atom.id, atom.element)
        else:
            point = None
            element = atom.element
        atoms.append(
            replace(
                atom,
                element=element,
                x=point[0] if point is not None else atom.x,
                y=point[1] if point is not None else atom.y,
                z=point[2] if point is not None else atom.z,
                radius=float(_VDW_RADII.get(element, atom.radius)),
                color=str(_ELEMENT_COLORS.get(element, atom.color)),
            )
        )

    generated_parent_ids: dict[str, int] = {}
    for bond in generated_bonds:
        first, second = bond.get("a1_id"), bond.get("a2_id")
        if str(first) in generated_items and isinstance(second, int) and second in source_ids:
            generated_parent_ids[str(first)] = second
        if str(second) in generated_items and isinstance(first, int) and first in source_ids:
            generated_parent_ids[str(second)] = first
    for atom_id, (item, point) in generated_items.items():
        element = str(item.get("element", "H"))
        raw_parent = item.get("parent_atom_id")
        if isinstance(raw_parent, int):
            parent_id: int | None = raw_parent
        elif isinstance(raw_parent, (float, str)):
            try:
                parent_id = int(raw_parent)
            except (TypeError, ValueError, OverflowError):
                parent_id = generated_parent_ids.get(atom_id)
        else:
            parent_id = generated_parent_ids.get(atom_id)
        if parent_id not in source_ids:
            parent_id = generated_parent_ids.get(atom_id)
        atoms.append(
            SceneAtom3D(
                id=atom_id,
                element=element,
                x=point[0],
                y=point[1],
                z=point[2],
                radius=float(_VDW_RADII.get(element, 1.20)),
                color=str(_ELEMENT_COLORS.get(element, "#B0B0B0")),
                generated=True,
                parent_atom_id=parent_id,
            )
        )

    bonds = [
        bond
        for bond in scene.bonds
        if bond.a1_id not in existing_generated_ids and bond.a2_id not in existing_generated_ids
    ]
    atom_ids = {atom.id for atom in atoms}
    for item in generated_bonds:
        first_raw, second_raw = item.get("a1_id"), item.get("a2_id")
        if not isinstance(first_raw, (int, str)) or not isinstance(second_raw, (int, str)):
            continue
        first, second = first_raw, second_raw
        if first not in atom_ids or second not in atom_ids:
            continue
        raw_order = item.get("order", 1)
        try:
            order = max(1, min(3, int(raw_order))) if isinstance(raw_order, (int, float, str)) else 1
        except (TypeError, ValueError, OverflowError):
            order = 1
        raw_bond_id = item.get("id")
        bond_id = raw_bond_id if isinstance(raw_bond_id, (int, str)) else f"visual:{len(bonds)}"
        bonds.append(
            SceneBond3D(
                id=bond_id,
                a1_id=first,
                a2_id=second,
                order=order,
                is_aromatic=bool(item.get("is_aromatic", False)),
            )
        )
    return SceneMolecule3D(tuple(atoms), tuple(bonds), scene.source)


class VisualFramePlayback:
    """Bounded, deterministic playback of measured frames for visual display only."""

    def __init__(self, *, transition_s: float = 0.08, max_pending: int = 128) -> None:
        duration = float(transition_s)
        if not math.isfinite(duration) or duration <= 0.0:
            raise ValueError("transition_s must be positive and finite")
        self.transition_s = duration
        self.max_pending = max(1, int(max_pending))
        self._pending: deque[CoordinateSet3D] = deque()
        self._current: CoordinateSet3D | None = None
        self._from: CoordinateSet3D | None = None
        self._to: CoordinateSet3D | None = None
        self._started_at: float | None = None
        self._enabled = False
        self._finished = False
        self.dropped_frames = 0

    @property
    def enabled(self) -> bool:
        return self._enabled

    @property
    def active(self) -> bool:
        return self._to is not None or bool(self._pending)

    @property
    def drained(self) -> bool:
        return self._enabled and self._finished and not self.active

    @property
    def pending_count(self) -> int:
        return len(self._pending) + int(self._to is not None)

    @property
    def current(self) -> CoordinateSet3D | None:
        return self._current

    def begin(self, initial: CoordinateSet3D | None) -> None:
        self.cancel()
        self._current = initial
        self._enabled = True
        self.dropped_frames = 0

    def enqueue(self, coordinates: CoordinateSet3D) -> bool:
        if not self._enabled:
            return False
        previous = self._pending[-1] if self._pending else self._to or self._current
        if previous is not None and not coordinate_sets_interpolatable(previous, coordinates):
            self._pending.clear()
            self._from = None
            self._to = None
            self._started_at = None
            self._current = coordinates
            return False
        if len(self._pending) >= self.max_pending:
            self._pending.popleft()
            self.dropped_frames += 1
        self._pending.append(coordinates)
        self._finished = False
        return True

    def finish(self, final: CoordinateSet3D | None = None) -> None:
        if not self._enabled:
            self.begin(final)
        if final is not None:
            previous = self._pending[-1] if self._pending else self._to or self._current
            if previous is None:
                self._current = final
            elif not coordinate_sets_interpolatable(previous, final):
                self._pending.clear()
                self._from = None
                self._to = None
                self._started_at = None
                self._current = final
            elif not _coordinate_sets_same_geometry(previous, final):
                self.enqueue(final)
        self._finished = True

    def cancel(self) -> None:
        self._pending.clear()
        self._current = None
        self._from = None
        self._to = None
        self._started_at = None
        self._enabled = False
        self._finished = False

    def sample(self, now: float) -> CoordinateSet3D | None:
        timestamp = float(now)
        if not math.isfinite(timestamp):
            raise ValueError("playback timestamp must be finite")
        if not self._enabled:
            return self._current
        if self._to is None:
            if not self._pending:
                return self._current
            target = self._pending.popleft()
            if self._current is None or not coordinate_sets_interpolatable(self._current, target):
                self._current = target
                return self._current
            self._from = self._current
            self._to = target
            self._started_at = timestamp
        if self._from is None or self._to is None or self._started_at is None:
            return self._current
        fraction = (timestamp - self._started_at) / self.transition_s
        if fraction >= 1.0:
            self._current = self._to
            self._from = None
            self._to = None
            self._started_at = None
            return self._current
        return interpolate_coordinate_sets(self._from, self._to, fraction)


def _coordinate_sets_same_geometry(
    first: CoordinateSet3D, second: CoordinateSet3D, *, tolerance: float = 1.0e-9
) -> bool:
    if not coordinate_sets_interpolatable(first, second):
        return False
    first_positions = first.normalized_positions()
    second_positions = second.normalized_positions()
    if any(
        math.dist(first_positions[atom_id], second_positions[atom_id]) > tolerance
        for atom_id in first_positions
    ):
        return False
    first_generated: dict[str, Mapping[str, object]] = {}
    second_generated: dict[str, Mapping[str, object]] = {}
    raw_first_generated = first.metadata.get("generated_atoms", ())
    raw_second_generated = second.metadata.get("generated_atoms", ())
    if isinstance(raw_first_generated, list):
        first_generated = {
            str(item.get("id")): item
            for item in raw_first_generated
            if isinstance(item, Mapping) and item.get("id") is not None
        }
    if isinstance(raw_second_generated, list):
        second_generated = {
            str(item.get("id")): item
            for item in raw_second_generated
            if isinstance(item, Mapping) and item.get("id") is not None
        }
    if first_generated.keys() != second_generated.keys():
        return False
    for atom_id in first_generated:
        first_point = first_generated[atom_id].get("coords") or first_generated[atom_id].get("position")
        second_point = second_generated[atom_id].get("coords") or second_generated[atom_id].get("position")
        if not isinstance(first_point, (list, tuple)) or not isinstance(second_point, (list, tuple)):
            continue
        if len(first_point) != 3 or len(second_point) != 3:
            return False
        if math.dist(tuple(float(value) for value in first_point), tuple(float(value) for value in second_point)) > tolerance:
            return False
    return True


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
    raw_element_overrides = coordinates.metadata.get("source_element_overrides", {})
    source_element_overrides: dict[int, str] = {}
    if isinstance(raw_element_overrides, Mapping):
        for atom_id, element in raw_element_overrides.items():
            try:
                source_element_overrides[int(atom_id)] = str(element)
            except (TypeError, ValueError, OverflowError):
                continue
    atoms: list[SceneAtom3D] = []
    for atom in sorted(getattr(graph, "atoms", {}).values(), key=lambda item: item.id):
        pos = positions.get(int(atom.id))
        if pos is None:
            continue
        element = source_element_overrides.get(int(atom.id), str(atom.element))
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
    generated_bonds = coordinates.metadata.get("generated_bonds", ()) if hasattr(coordinates, "metadata") else ()
    generated_ids = {
        str(item.get("id"))
        for item in generated_atoms
        if isinstance(item, Mapping) and item.get("id") is not None
    } if isinstance(generated_atoms, list) else set()
    generated_parent_ids: dict[str, int] = {}
    if isinstance(generated_bonds, list):
        for item in generated_bonds:
            if not isinstance(item, Mapping):
                continue
            first_id, second_id = item.get("a1_id"), item.get("a2_id")
            for generated_id, candidate_parent in (
                (str(first_id), second_id),
                (str(second_id), first_id),
            ):
                if generated_id not in generated_ids or candidate_parent is None:
                    continue
                try:
                    parent_id = int(candidate_parent)
                except (TypeError, ValueError, OverflowError):
                    continue
                if parent_id in getattr(graph, "atoms", {}):
                    generated_parent_ids[generated_id] = parent_id
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
                    parent_atom_id=generated_parent_ids.get(atom_id),
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
