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
    repeat_min: int | None = None
    repeat_max: int | None = None


@dataclass(frozen=True)
class RGroupAtom:
    """Átomo de consulta tipo R/R1/R2 usado en Markush."""

    atom_id: int
    label: str
    allowed_substituents: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class MarkushSummary:
    """Resumen reutilizable por UI/exportadores de polímeros y Markush."""

    r_groups: list[RGroupAtom] = field(default_factory=list)
    polymer_repeats: list[PolymerRepeat] = field(default_factory=list)

    @property
    def has_markush(self) -> bool:
        return bool(self.r_groups or self.polymer_repeats)

    def as_dict(self) -> dict[str, object]:
        """Representación estable para UI/exportadores."""
        return {
            "r_groups": [
                {
                    "atom_id": r_group.atom_id,
                    "label": r_group.label,
                    "allowed_substituents": list(r_group.allowed_substituents),
                }
                for r_group in self.r_groups
            ],
            "polymer_repeats": [
                {
                    "repeat_label": repeat.repeat_label,
                    "kind": repeat.kind,
                    "rect": list(repeat.rect),
                    "repeat_min": repeat.repeat_min,
                    "repeat_max": repeat.repeat_max,
                }
                for repeat in self.polymer_repeats
            ],
        }


def summarize_markush(
    graph: MolGraph,
    *,
    bracket_items: Iterable[object] = (),
) -> MarkushSummary:
    """Extrae R-groups y corchetes de repetición desde modelo/canvas."""
    r_groups = [
        RGroupAtom(
            int(atom.id),
            str(atom.element),
            tuple(getattr(atom, "r_group_substituents", ()) or ()),
        )
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
                    repeat_min=_repeat_bound(label, 0),
                    repeat_max=_repeat_bound(label, 1),
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


def sanitize_r_group_substituents(value: str | list[str] | tuple[str, ...] | None) -> tuple[str, ...]:
    """Normaliza una lista de sustituyentes Markush simples."""
    if value is None:
        return ()
    if isinstance(value, str):
        raw_items = value.replace(";", ",").split(",")
    else:
        raw_items = list(value)
    out: list[str] = []
    seen: set[str] = set()
    for item in raw_items:
        text = "".join(ch for ch in str(item).strip() if ch.isalnum() or ch in {"-", "+", "_"})
        if not text or text in seen:
            continue
        seen.add(text)
        out.append(text)
    return tuple(out)


def set_r_group_substituents(graph: MolGraph, atom_ids: list[int], substituents: str | list[str] | tuple[str, ...]) -> tuple[str, ...]:
    """Asigna sustituyentes permitidos a átomos R/Rn del grafo."""
    values = sanitize_r_group_substituents(substituents)
    for atom_id in atom_ids:
        atom = graph.atoms.get(int(atom_id))
        if atom is None or not _is_r_group_label(str(atom.element)):
            continue
        atom.r_group_substituents = values
    return values


def _repeat_bound(label: str, index: int) -> int | None:
    text = str(label or "").strip()
    if text.isdigit():
        return int(text)
    for sep in ("-", ".."):
        if sep in text:
            parts = text.split(sep, 1)
            try:
                return int(parts[index].strip())
            except Exception:
                return None
    return None
