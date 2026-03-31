"""Modelos semánticos para diagramas electrónicos."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Mapping


RepresentationKind = Literal["boxes", "bar"]
ConnectorStyle = Literal["solid", "dashed"]
SemanticDiagramKind = Literal["atomic", "molecular_orbital", "ligand_field"]


def _normalize_box_occupancies(values: object, degeneracy: int) -> list[int]:
    """Normaliza ocupaciones por caja a enteros 0/1/2."""
    target = max(0, int(degeneracy))
    if target == 0:
        return []

    if isinstance(values, str):
        tokens = [token.strip() for token in values.split(",")]
    elif isinstance(values, (list, tuple)):
        tokens = list(values)
    else:
        tokens = []

    aliases = {
        "": 0,
        "0": 0,
        "empty": 0,
        "none": 0,
        "up": 1,
        "down": 1,
        "u": 1,
        "d": 1,
        "1": 1,
        "pair": 2,
        "paired": 2,
        "upup": 2,
        "downdown": 2,
        "2": 2,
    }

    normalized: list[int] = []
    for token in tokens[:target]:
        value = token
        if isinstance(token, str):
            value = aliases.get(token.strip().lower(), token)
        try:
            occupancy = int(value)
        except Exception:
            occupancy = 0
        normalized.append(max(0, min(2, occupancy)))

    if len(normalized) < target:
        normalized.extend([0] * (target - len(normalized)))
    return normalized


@dataclass
class DiagramLane:
    id: str
    title: str
    x: float

    def to_json_dict(self) -> dict[str, Any]:
        return {
            "id": str(self.id),
            "title": str(self.title),
            "x": float(self.x),
        }

    @classmethod
    def from_json_dict(cls, data: Mapping[str, Any]) -> "DiagramLane":
        return cls(
            id=str(data.get("id", "")),
            title=str(data.get("title", "")),
            x=float(data.get("x", 0.0)),
        )


@dataclass
class DiagramLevel:
    id: str
    lane_id: str
    energy: float
    label: str
    representation: RepresentationKind
    degeneracy: int
    occupancies: list[int]
    width: float = 96.0
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.id = str(self.id)
        self.lane_id = str(self.lane_id)
        self.label = str(self.label)
        self.representation = "bar" if self.representation == "bar" else "boxes"
        self.degeneracy = max(1, int(self.degeneracy))
        self.energy = float(self.energy)
        self.width = max(12.0, float(self.width))
        self.occupancies = _normalize_box_occupancies(self.occupancies, self.degeneracy)
        self.metadata = dict(self.metadata or {})

    def to_json_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "lane_id": self.lane_id,
            "energy": self.energy,
            "label": self.label,
            "representation": self.representation,
            "degeneracy": self.degeneracy,
            "occupancies": list(self.occupancies),
            "width": self.width,
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_json_dict(cls, data: Mapping[str, Any]) -> "DiagramLevel":
        return cls(
            id=str(data.get("id", "")),
            lane_id=str(data.get("lane_id", "")),
            energy=float(data.get("energy", 0.0)),
            label=str(data.get("label", "")),
            representation="bar" if data.get("representation") == "bar" else "boxes",
            degeneracy=max(1, int(data.get("degeneracy", 1))),
            occupancies=_normalize_box_occupancies(
                data.get("occupancies", ()),
                int(data.get("degeneracy", 1)),
            ),
            width=float(data.get("width", 96.0)),
            metadata=dict(data.get("metadata", {}) or {}),
        )


@dataclass
class DiagramConnector:
    source_level_id: str
    target_level_id: str
    style: ConnectorStyle = "dashed"

    def __post_init__(self) -> None:
        self.source_level_id = str(self.source_level_id)
        self.target_level_id = str(self.target_level_id)
        self.style = "solid" if self.style == "solid" else "dashed"

    def to_json_dict(self) -> dict[str, Any]:
        return {
            "source_level_id": self.source_level_id,
            "target_level_id": self.target_level_id,
            "style": self.style,
        }

    @classmethod
    def from_json_dict(cls, data: Mapping[str, Any]) -> "DiagramConnector":
        return cls(
            source_level_id=str(data.get("source_level_id", "")),
            target_level_id=str(data.get("target_level_id", "")),
            style="solid" if data.get("style") == "solid" else "dashed",
        )


@dataclass
class SemanticDiagram:
    kind: SemanticDiagramKind
    title: str
    lanes: list[DiagramLane]
    levels: list[DiagramLevel]
    connectors: list[DiagramConnector]
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        valid_kinds = {"atomic", "molecular_orbital", "ligand_field"}
        self.kind = self.kind if self.kind in valid_kinds else "atomic"
        self.title = str(self.title)
        self.lanes = [lane if isinstance(lane, DiagramLane) else DiagramLane.from_json_dict(lane) for lane in self.lanes]
        self.levels = [
            level if isinstance(level, DiagramLevel) else DiagramLevel.from_json_dict(level)
            for level in self.levels
        ]
        self.connectors = [
            connector
            if isinstance(connector, DiagramConnector)
            else DiagramConnector.from_json_dict(connector)
            for connector in self.connectors
        ]
        self.metadata = dict(self.metadata or {})

    def to_json_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "title": self.title,
            "lanes": [lane.to_json_dict() for lane in self.lanes],
            "levels": [level.to_json_dict() for level in self.levels],
            "connectors": [connector.to_json_dict() for connector in self.connectors],
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_json_dict(cls, data: Mapping[str, Any]) -> "SemanticDiagram":
        kind = str(data.get("kind", "atomic")).strip().lower()
        if kind not in {"atomic", "molecular_orbital", "ligand_field"}:
            kind = "atomic"
        return cls(
            kind=kind,  # type: ignore[arg-type]
            title=str(data.get("title", "")),
            lanes=[DiagramLane.from_json_dict(item) for item in list(data.get("lanes", []) or [])],
            levels=[DiagramLevel.from_json_dict(item) for item in list(data.get("levels", []) or [])],
            connectors=[
                DiagramConnector.from_json_dict(item)
                for item in list(data.get("connectors", []) or [])
            ],
            metadata=dict(data.get("metadata", {}) or {}),
        )
