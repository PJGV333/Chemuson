from __future__ import annotations

"""Extracción ligera de anotaciones polímero/Markush."""

from dataclasses import dataclass, field
from typing import Iterable

from chemuson.core.model import MolGraph


@dataclass(frozen=True)
class PolymerRepeat:
    """Corchete de repetición polimérica persistente."""

    repeat_label: str
    kind: str
    rect: tuple[float, float, float, float]


@dataclass(frozen=True)
class RGroupAtom:
    """Átomo de consulta tipo R/R1/R2 usado en Markush."""

    atom_id: int
    label: str


@dataclass(frozen=True)
class MarkushSummary:
    """Resumen reutilizable por UI/exportadores de polímeros y Markush."""

    r_groups: list[RGroupAtom] = field(default_factory=list)
    polymer_repeats: list[PolymerRepeat] = field(default_factory=list)

    @property
    def has_markush(self) -> bool:
        return bool(self.r_groups or self.polymer_repeats)


def summarize_markush(
    graph: MolGraph,
    *,
    bracket_items: Iterable[object] = (),
) -> MarkushSummary:
    """Extrae R-groups y corchetes de repetición desde modelo/canvas."""
    r_groups = [
        RGroupAtom(int(atom.id), str(atom.element))
        for atom in sorted(graph.atoms.values(), key=lambda item: item.id)
        if _is_r_group_label(str(atom.element))
    ]
    repeats: list[PolymerRepeat] = []
    for item in bracket_items:
        label = ""
        if hasattr(item, "repeat_label"):
            try:
                label = str(item.repeat_label() or "").strip()
            except Exception:
                label = ""
        if not label:
            continue
        try:
            rect = item.base_rect()
            kind = str(getattr(item, "_kind", "[]"))
            repeats.append(
                PolymerRepeat(
                    repeat_label=label,
                    kind=kind,
                    rect=(float(rect.x()), float(rect.y()), float(rect.width()), float(rect.height())),
                )
            )
        except Exception:
            continue
    return MarkushSummary(r_groups=r_groups, polymer_repeats=repeats)


def _is_r_group_label(label: str) -> bool:
    text = str(label or "").strip()
    if text == "R":
        return True
    return len(text) > 1 and text[0] == "R" and text[1:].isdigit()
