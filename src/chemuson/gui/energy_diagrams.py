"""Presets y utilidades para diagramas de niveles de energia."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

from PyQt6.QtCore import QSizeF

from chemuson.gui.diagram_models import (
    DiagramConnector,
    DiagramLane,
    DiagramLevel,
    SemanticDiagram,
)


ENERGY_OCCUPANCY_VALUES = {"empty", "up", "down", "pair", "upup", "downdown"}
ENERGY_LABEL_SIDES = {"left", "right"}
ENERGY_DIAGRAM_FAMILIES = {"row", "levels", "mo"}


@dataclass(frozen=True)
class EnergyDiagramPreset:
    """Preset de insercion para un diagrama de cajas de orbitales."""

    kind: str
    display_name: str
    family: str
    label: str
    label_side: str
    slot_count: int


ENERGY_DIAGRAM_PRESETS: dict[str, EnergyDiagramPreset] = {
    "sublevel_s": EnergyDiagramPreset(
        kind="sublevel_s",
        display_name="Subnivel s",
        family="row",
        label="s",
        label_side="left",
        slot_count=1,
    ),
    "custom_level": EnergyDiagramPreset(
        kind="custom_level",
        display_name="Nivel personalizable",
        family="row",
        label="",
        label_side="left",
        slot_count=1,
    ),
    "sublevel_p": EnergyDiagramPreset(
        kind="sublevel_p",
        display_name="Subnivel p",
        family="row",
        label="p",
        label_side="left",
        slot_count=3,
    ),
    "sublevel_d": EnergyDiagramPreset(
        kind="sublevel_d",
        display_name="Subnivel d",
        family="row",
        label="d",
        label_side="left",
        slot_count=5,
    ),
    "sublevel_f": EnergyDiagramPreset(
        kind="sublevel_f",
        display_name="Subnivel f",
        family="row",
        label="f",
        label_side="left",
        slot_count=7,
    ),
    "hybrid_sp": EnergyDiagramPreset(
        kind="hybrid_sp",
        display_name="Hibrido sp",
        family="row",
        label="sp",
        label_side="right",
        slot_count=2,
    ),
    "hybrid_sp2": EnergyDiagramPreset(
        kind="hybrid_sp2",
        display_name="Hibrido sp2",
        family="row",
        label="sp2",
        label_side="right",
        slot_count=3,
    ),
    "hybrid_sp3": EnergyDiagramPreset(
        kind="hybrid_sp3",
        display_name="Hibrido sp3",
        family="row",
        label="sp3",
        label_side="right",
        slot_count=4,
    ),
    "levels_aufbau": EnergyDiagramPreset(
        kind="levels_aufbau",
        display_name="Niveles de energia",
        family="levels",
        label="",
        label_side="left",
        slot_count=56,
    ),
    "mo_2s_2p": EnergyDiagramPreset(
        kind="mo_2s_2p",
        display_name="Orbital molecular 2s/2p",
        family="mo",
        label="",
        label_side="left",
        slot_count=16,
    ),
}

ENERGY_DIAGRAM_MENU_ORDER: tuple[str, ...] = (
    "sublevel_s",
    "custom_level",
    "sublevel_p",
    "sublevel_d",
    "sublevel_f",
    "hybrid_sp",
    "hybrid_sp2",
    "hybrid_sp3",
)

DEFAULT_ENERGY_DIAGRAM_KIND = "sublevel_s"

SUBSHELL_DEGENERACY: dict[str, int] = {
    "s": 1,
    "p": 3,
    "d": 5,
    "f": 7,
}
MADELUNG_SEQUENCE: tuple[str, ...] = (
    "1s",
    "2s",
    "2p",
    "3s",
    "3p",
    "4s",
    "3d",
    "4p",
    "5s",
    "4d",
    "5p",
    "6s",
    "4f",
    "5d",
    "6p",
    "7s",
    "5f",
    "6d",
    "7p",
)


def energy_diagram_preset(kind: str) -> EnergyDiagramPreset:
    """Devuelve el preset asociado o un fallback seguro."""
    return ENERGY_DIAGRAM_PRESETS.get(kind, ENERGY_DIAGRAM_PRESETS[DEFAULT_ENERGY_DIAGRAM_KIND])


def energy_diagram_display_name(kind: str) -> str:
    """Etiqueta legible para UI."""
    return energy_diagram_preset(kind).display_name


def energy_diagram_family(kind: str) -> str:
    """Familia de layout para un preset."""
    family = str(energy_diagram_preset(kind).family or "").strip().lower()
    return family if family in ENERGY_DIAGRAM_FAMILIES else "row"


def energy_diagram_supports_free_label(kind: str) -> bool:
    """Indica si el preset usa la etiqueta libre del objeto."""
    return energy_diagram_family(kind) == "row"


def energy_diagram_default_style_payload(kind: str) -> dict[str, object]:
    """Overrides de estilo por preset."""
    if kind == "levels_aufbau":
        return {
            "connector_color": "#C75A2C",
            "box_stroke_visible": True,
            "box_base_visible": False,
            "fill_visible": True,
        }
    if kind == "mo_2s_2p":
        return {
            "connector_color": "#7A7A7A",
            "box_stroke_visible": False,
            "box_base_visible": False,
            "fill_visible": False,
        }
    return {
        "box_stroke_visible": True,
        "box_base_visible": False,
        "fill_visible": True,
    }


def energy_diagram_tool_id(kind: str) -> str:
    """Convierte un preset a un id de herramienta."""
    return f"tool_energy_diagram_{kind}"


def energy_diagram_kind_from_tool_id(tool_id: str) -> str | None:
    """Extrae el preset desde un id de herramienta."""
    prefix = "tool_energy_diagram_"
    if not str(tool_id or "").startswith(prefix):
        return None
    kind = str(tool_id)[len(prefix) :]
    return kind if kind in ENERGY_DIAGRAM_PRESETS else None


def normalize_energy_label_side(value: object, *, default: str = "left") -> str:
    """Normaliza la posicion de etiqueta soportada."""
    normalized = str(value or "").strip().lower()
    if normalized in ENERGY_LABEL_SIDES:
        return normalized
    fallback = str(default or "").strip().lower()
    return fallback if fallback in ENERGY_LABEL_SIDES else "left"


def normalize_energy_occupancy(value: object) -> str:
    """Normaliza un estado de ocupacion electronica."""
    if isinstance(value, bool):
        return "up" if value else "empty"
    if isinstance(value, (int, float)):
        try:
            numeric = int(value)
        except Exception:
            numeric = 0
        if numeric <= 0:
            return "empty"
        if numeric == 1:
            return "up"
        return "pair"
    raw = str(value or "").strip().lower()
    aliases = {
        "": "empty",
        "0": "empty",
        "-": "empty",
        "empty": "empty",
        "none": "empty",
        "vacio": "empty",
        "vacío": "empty",
        "u": "up",
        "up": "up",
        "arriba": "up",
        "↑": "up",
        "uu": "upup",
        "upup": "upup",
        "arribaarriba": "upup",
        "↑↑": "upup",
        "d": "down",
        "down": "down",
        "abajo": "down",
        "↓": "down",
        "dd": "downdown",
        "downdown": "downdown",
        "abajoabajo": "downdown",
        "↓↓": "downdown",
        "2": "pair",
        "pair": "pair",
        "paired": "pair",
        "par": "pair",
        "ud": "pair",
        "du": "pair",
        "↑↓": "pair",
        "↓↑": "pair",
    }
    normalized = aliases.get(raw, raw)
    return normalized if normalized in ENERGY_OCCUPANCY_VALUES else "empty"


def normalize_energy_occupancies(
    values: object,
    *,
    kind: str | None = None,
    box_count: int | None = None,
) -> tuple[str, ...]:
    """Normaliza una secuencia de ocupaciones con el largo esperado."""
    target_count = int(box_count or 0)
    if target_count <= 0:
        target_count = energy_diagram_preset(kind or DEFAULT_ENERGY_DIAGRAM_KIND).slot_count

    if isinstance(values, str):
        tokens = [token.strip() for token in values.split(",")]
    elif isinstance(values, (list, tuple)):
        tokens = list(values)
    else:
        tokens = []

    normalized = [normalize_energy_occupancy(token) for token in tokens[:target_count]]
    if len(normalized) < target_count:
        normalized.extend(["empty"] * (target_count - len(normalized)))
    return tuple(normalized)


def default_energy_label(kind: str) -> str:
    """Etiqueta inicial por preset."""
    return energy_diagram_preset(kind).label


def default_energy_label_side(kind: str) -> str:
    """Lado inicial de la etiqueta."""
    return energy_diagram_preset(kind).label_side


def default_energy_occupancies(kind: str) -> tuple[str, ...]:
    """Ocupacion inicial vacia para el preset."""
    preset = energy_diagram_preset(kind)
    return tuple("empty" for _ in range(preset.slot_count))


def energy_diagram_default_size(kind: str, bond_length: float) -> QSizeF:
    """Tamano visual sugerido para insertar un diagrama."""
    preset = energy_diagram_preset(kind)
    family = energy_diagram_family(kind)
    unit = max(26.0, float(bond_length) * 0.72)
    if family == "levels":
        return QSizeF(unit * 12.5, unit * 15.5)
    if family == "mo":
        return QSizeF(unit * 12.5, unit * 9.5)
    box_height = unit
    box_width = unit * 0.84
    gap = unit * 0.16
    label_space = unit * 1.25
    width = box_width * preset.slot_count + gap * max(0, preset.slot_count - 1)
    if preset.label:
        width += label_space
    height = box_height + unit * 0.34
    return QSizeF(width, height)


def _normalized_electron_count(value: int, *, capacity: int) -> int:
    return max(0, min(int(value), max(0, int(capacity))))


def _fill_level_hund(electron_count: int, degeneracy: int) -> list[int]:
    """Llena un conjunto degenerado siguiendo Hund dentro del mismo nivel."""
    box_count = max(1, int(degeneracy))
    remaining = max(0, min(int(electron_count), box_count * 2))
    occupancies = [0] * box_count
    for index in range(min(remaining, box_count)):
        occupancies[index] = 1
    remaining -= min(remaining, box_count)
    for index in range(box_count):
        if remaining <= 0:
            break
        occupancies[index] += 1
        remaining -= 1
    return occupancies


def _sequential_level_fill(
    level_specs: list[dict[str, Any]],
    electron_count: int,
) -> dict[str, list[int]]:
    """Llena niveles de menor a mayor energía."""
    occupancies: dict[str, list[int]] = {}
    remaining = max(0, int(electron_count))
    for spec in level_specs:
        capacity = int(spec["degeneracy"]) * 2
        filled = min(remaining, capacity)
        occupancies[str(spec["id"])] = _fill_level_hund(filled, int(spec["degeneracy"]))
        remaining -= filled
    return occupancies


def _high_spin_level_fill(
    level_specs: list[dict[str, Any]],
    electron_count: int,
) -> dict[str, list[int]]:
    """Maximiza desapareados entre niveles antes de aparear."""
    occupancies = {
        str(spec["id"]): [0] * max(1, int(spec["degeneracy"]))
        for spec in level_specs
    }
    orbital_order = [
        (str(spec["id"]), orbital_index)
        for spec in level_specs
        for orbital_index in range(max(1, int(spec["degeneracy"])))
    ]
    remaining = max(0, int(electron_count))
    for level_id, orbital_index in orbital_order:
        if remaining <= 0:
            break
        occupancies[level_id][orbital_index] = 1
        remaining -= 1
    for level_id, orbital_index in orbital_order:
        if remaining <= 0:
            break
        if occupancies[level_id][orbital_index] < 2:
            occupancies[level_id][orbital_index] += 1
            remaining -= 1
    return occupancies


def _count_level_electrons(occupancies: dict[str, list[int]]) -> int:
    return sum(sum(level) for level in occupancies.values())


def _count_unpaired_electrons(occupancies: dict[str, list[int]]) -> int:
    return sum(1 for level in occupancies.values() for value in level if int(value) == 1)


SUBSHELL_ORDER: dict[str, int] = {
    "s": 0,
    "p": 1,
    "d": 2,
    "f": 3,
}
ATOMIC_NUMBERS: dict[str, int] = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
}
NEUTRAL_CONFIGURATION_EXCEPTIONS: dict[str, dict[str, int]] = {
    "Cr": {
        "1s": 2,
        "2s": 2,
        "2p": 6,
        "3s": 2,
        "3p": 6,
        "3d": 5,
        "4s": 1,
    },
    "Cu": {
        "1s": 2,
        "2s": 2,
        "2p": 6,
        "3s": 2,
        "3p": 6,
        "3d": 10,
        "4s": 1,
    },
}


def _level_electron_map(levels: list[DiagramLevel]) -> dict[str, list[int]]:
    return {str(level.id): list(level.occupancies) for level in levels}


def _magnetic_behavior(unpaired_electrons: int) -> str:
    return "diamagnetic" if int(unpaired_electrons) <= 0 else "paramagnetic"


def _format_metric_value(value: float | int) -> str:
    numeric = float(value)
    if abs(numeric - round(numeric)) <= 1e-6:
        return str(int(round(numeric)))
    return f"{numeric:.1f}".rstrip("0").rstrip(".")


def _configuration_sort_key(level: DiagramLevel) -> tuple[int, int, str]:
    label = str(level.id or level.label or "").strip()
    principal_text = []
    family = ""
    for char in label:
        if char.isdigit():
            principal_text.append(char)
            continue
        family = char.lower()
        break
    try:
        principal = int("".join(principal_text)) if principal_text else 999
    except Exception:
        principal = 999
    return principal, SUBSHELL_ORDER.get(family, 99), label


def configuration_string_from_levels(levels: list[DiagramLevel]) -> str:
    """Construye una configuración electrónica legible a partir de niveles semánticos."""
    parts: list[str] = []
    for level in sorted(levels, key=_configuration_sort_key):
        electrons = sum(int(value) for value in level.occupancies)
        if electrons <= 0:
            continue
        token = str(level.id or level.label or "").strip()
        if not token:
            continue
        parts.append(f"{token}{electrons}")
    return " ".join(parts)


def _builder_payload(name: str, params: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": str(name),
        "params": dict(params),
    }


def _cfse_string(levels: list[DiagramLevel], geometry: str) -> tuple[float | None, str]:
    totals = {str(level.id): sum(int(value) for value in level.occupancies) for level in levels}
    if geometry == "octahedral":
        value = totals.get("t2g", 0) * (-0.4) + totals.get("eg", 0) * 0.6
        return value, f"CFSE: {_format_metric_value(value)} Δo"
    if geometry == "tetrahedral":
        value = totals.get("e", 0) * (-0.6) + totals.get("t2", 0) * 0.4
        return value, f"CFSE: {_format_metric_value(value)} Δt"
    return None, "CFSE: n/a"


def refresh_semantic_diagram_metadata(diagram: SemanticDiagram) -> SemanticDiagram:
    """Recalcula líneas de resumen y métricas derivadas sin perder metadatos previos."""
    metadata = dict(diagram.metadata or {})
    occupancies = _level_electron_map(diagram.levels)
    unpaired_electrons = _count_unpaired_electrons(occupancies)
    magnetic_behavior = _magnetic_behavior(unpaired_electrons)
    metadata.setdefault("show_summary", True)

    if diagram.kind == "atomic":
        electron_count = int(
            metadata.get("electron_count", _count_level_electrons(occupancies)) or 0
        )
        configuration_string = configuration_string_from_levels(diagram.levels)
        derived_metrics = dict(metadata.get("derived_metrics", {}) or {})
        derived_metrics.update(
            {
                "configuration_string": configuration_string,
                "electron_count": electron_count,
                "unpaired_electrons": unpaired_electrons,
                "magnetic_behavior": magnetic_behavior,
            }
        )
        summary_lines = []
        if configuration_string:
            summary_lines.append(configuration_string)
        summary_lines.append(f"Electrons: {electron_count}")
        summary_lines.append(
            f"Unpaired electrons: {unpaired_electrons} ({magnetic_behavior})"
        )
        metadata.update(
            {
                "electron_count": electron_count,
                "unpaired_electrons": unpaired_electrons,
                "magnetic_behavior": magnetic_behavior,
                "configuration_string": configuration_string,
                "derived_metrics": derived_metrics,
                "summary_lines": summary_lines,
            }
        )
    elif diagram.kind == "molecular_orbital":
        molecular_levels = [
            level
            for level in diagram.levels
            if str(level.metadata.get("orbital_origin", "")) == "molecular"
        ]
        molecular_occupancies = _level_electron_map(molecular_levels)
        bonding_electrons = sum(
            sum(int(value) for value in level.occupancies)
            for level in molecular_levels
            if bool(level.metadata.get("bonding"))
        )
        antibonding_electrons = sum(
            sum(int(value) for value in level.occupancies)
            for level in molecular_levels
            if not bool(level.metadata.get("bonding"))
        )
        bond_order = (bonding_electrons - antibonding_electrons) / 2.0
        unpaired_electrons = _count_unpaired_electrons(molecular_occupancies)
        magnetic_behavior = _magnetic_behavior(unpaired_electrons)
        derived_metrics = dict(metadata.get("derived_metrics", {}) or {})
        derived_metrics.update(
            {
                "bond_order": bond_order,
                "unpaired_electrons": unpaired_electrons,
                "magnetic_behavior": magnetic_behavior,
            }
        )
        metadata.update(
            {
                "bond_order": bond_order,
                "unpaired_electrons": unpaired_electrons,
                "magnetic_behavior": magnetic_behavior,
                "derived_metrics": derived_metrics,
                "summary_lines": [
                    f"Bond order: {_format_metric_value(bond_order)}",
                    f"Unpaired electrons: {unpaired_electrons} ({magnetic_behavior})",
                ],
            }
        )
    else:
        geometry = str(metadata.get("geometry", "octahedral") or "octahedral")
        spin_mode = str(metadata.get("spin_mode", "high") or "high")
        d_electrons = int(metadata.get("d_electrons", _count_level_electrons(occupancies)) or 0)
        cfse_value, cfse_string = _cfse_string(diagram.levels, geometry)
        derived_metrics = dict(metadata.get("derived_metrics", {}) or {})
        derived_metrics.update(
            {
                "d_electrons": d_electrons,
                "geometry": geometry,
                "spin_mode": spin_mode,
                "unpaired_electrons": unpaired_electrons,
                "magnetic_behavior": magnetic_behavior,
                "cfse": cfse_value,
            }
        )
        summary_lines = [
            f"d-count: d{d_electrons}",
            f"Geometry: {geometry.replace('_', ' ')}",
            f"Spin mode: {spin_mode}",
            f"Unpaired electrons: {unpaired_electrons} ({magnetic_behavior})",
        ]
        if geometry in {"octahedral", "tetrahedral"}:
            summary_lines.append(cfse_string)
        elif geometry == "square_planar":
            summary_lines.append("CFSE: n/a")
        metadata.update(
            {
                "d_electrons": d_electrons,
                "geometry": geometry,
                "spin_mode": spin_mode,
                "unpaired_electrons": unpaired_electrons,
                "magnetic_behavior": magnetic_behavior,
                "derived_metrics": derived_metrics,
                "summary_lines": summary_lines,
            }
        )

    diagram.metadata = metadata
    return diagram


def _atomic_diagram_levels_from_electron_map(
    electron_map: dict[str, int],
    *,
    expanded_subshells: bool,
    max_n: int,
) -> tuple[list[DiagramLane], list[DiagramLevel]]:
    sequence = [
        label
        for label in MADELUNG_SEQUENCE
        if int(label[0]) <= max_n and int(electron_map.get(label, 0) or 0) > 0
    ]
    lane_positions = {"s": 0.0, "p": 132.0, "d": 276.0, "f": 444.0}
    used_families = [
        family
        for family in ("s", "p", "d", "f")
        if any(label.endswith(family) for label in sequence)
    ]
    lanes = [
        DiagramLane(id=f"{family}_lane", title=family, x=lane_positions[family])
        for family in used_families
    ]
    levels: list[DiagramLevel] = []
    for energy_index, label in enumerate(sequence):
        family = label[1]
        degeneracy = SUBSHELL_DEGENERACY[family]
        capacity_here = degeneracy * 2
        occupied_electrons = _normalized_electron_count(
            int(electron_map.get(label, 0) or 0),
            capacity=capacity_here,
        )
        levels.append(
            DiagramLevel(
                id=label,
                lane_id=f"{family}_lane",
                energy=float(energy_index),
                label=label,
                representation="bar" if (not expanded_subshells and degeneracy > 1) else "boxes",
                degeneracy=degeneracy,
                occupancies=_fill_level_hund(occupied_electrons, degeneracy),
                metadata={
                    "principal_quantum_number": int(label[0]),
                    "subshell": family,
                    "electron_capacity": capacity_here,
                },
            )
        )
    return lanes, levels


def _atomic_title(symbol: str | None, charge: int | None, electron_count: int, title: str | None) -> str:
    if title:
        return title
    if symbol:
        charge_value = int(charge or 0)
        if charge_value == 0:
            return f"{symbol} Atomic Diagram"
        sign = "+" if charge_value > 0 else "-"
        magnitude = abs(charge_value)
        charge_label = sign if magnitude == 1 else f"{magnitude}{sign}"
        return f"{symbol}{charge_label} Atomic Diagram"
    return f"Atomic Diagram ({electron_count} e-)"


def build_atomic_subshell_diagram(
    electron_count: int,
    title: str | None = None,
    expanded_subshells: bool = True,
    max_n: int = 7,
) -> SemanticDiagram:
    """Construye un diagrama semántico de subniveles atómicos."""
    max_n_value = max(1, min(7, int(max_n)))
    sequence = [
        label
        for label in MADELUNG_SEQUENCE
        if int(label[0]) <= max_n_value
    ]
    capacity = sum(SUBSHELL_DEGENERACY[label[1]] * 2 for label in sequence)
    total_electrons = _normalized_electron_count(electron_count, capacity=capacity)
    electron_map: dict[str, int] = {}
    remaining = total_electrons
    for label in sequence:
        if remaining <= 0:
            break
        family = label[1]
        capacity_here = SUBSHELL_DEGENERACY[family] * 2
        occupied_electrons = min(remaining, capacity_here)
        electron_map[label] = occupied_electrons
        remaining -= occupied_electrons

    lanes, levels = _atomic_diagram_levels_from_electron_map(
        electron_map,
        expanded_subshells=expanded_subshells,
        max_n=max_n_value,
    )
    diagram = SemanticDiagram(
        kind="atomic",
        title=title or f"Atomic Diagram ({total_electrons} e-)",
        lanes=lanes,
        levels=levels,
        connectors=[],
        metadata={
            "electron_count": total_electrons,
            "filling_rule": "aufbau_hund",
            "builder": _builder_payload(
                "build_atomic_subshell_diagram",
                {
                    "electron_count": total_electrons,
                    "title": title,
                    "expanded_subshells": bool(expanded_subshells),
                    "max_n": max_n_value,
                },
            ),
        },
    )
    return refresh_semantic_diagram_metadata(diagram)


def build_atomic_species_diagram(
    symbol: str,
    charge: int = 0,
    expanded_subshells: bool = True,
    title: str | None = None,
    use_known_exceptions: bool = True,
) -> SemanticDiagram:
    """Construye un diagrama atómico a partir de especie química y carga."""
    normalized_symbol = str(symbol or "").strip().capitalize()
    atomic_number = ATOMIC_NUMBERS.get(normalized_symbol)
    if atomic_number is None:
        raise ValueError(f"Unsupported atomic symbol: {symbol}")
    electron_count = max(0, atomic_number - int(charge))
    max_n_value = 7
    capacity = sum(SUBSHELL_DEGENERACY[label[1]] * 2 for label in MADELUNG_SEQUENCE)
    normalized_electrons = _normalized_electron_count(electron_count, capacity=capacity)

    exception_map = None
    if bool(use_known_exceptions) and int(charge) == 0:
        exception_map = NEUTRAL_CONFIGURATION_EXCEPTIONS.get(normalized_symbol)

    if exception_map is None or sum(exception_map.values()) != normalized_electrons:
        diagram = build_atomic_subshell_diagram(
            normalized_electrons,
            title=_atomic_title(normalized_symbol, charge, normalized_electrons, title),
            expanded_subshells=expanded_subshells,
            max_n=max_n_value,
        )
        metadata = dict(diagram.metadata or {})
    else:
        lanes, levels = _atomic_diagram_levels_from_electron_map(
            exception_map,
            expanded_subshells=expanded_subshells,
            max_n=max_n_value,
        )
        diagram = SemanticDiagram(
            kind="atomic",
            title=_atomic_title(normalized_symbol, charge, normalized_electrons, title),
            lanes=lanes,
            levels=levels,
            connectors=[],
            metadata={
                "electron_count": normalized_electrons,
                "filling_rule": "known_exception",
            },
        )
        metadata = dict(diagram.metadata or {})

    configuration_string = configuration_string_from_levels(diagram.levels)
    metadata.update(
        {
            "symbol": normalized_symbol,
            "charge": int(charge),
            "electron_count": normalized_electrons,
            "configuration_string": configuration_string,
            "builder": _builder_payload(
                "build_atomic_species_diagram",
                {
                    "symbol": normalized_symbol,
                    "charge": int(charge),
                    "expanded_subshells": bool(expanded_subshells),
                    "title": title,
                    "use_known_exceptions": bool(use_known_exceptions),
                },
            ),
        }
    )
    diagram.metadata = metadata
    return refresh_semantic_diagram_metadata(diagram)


def build_diatomic_mo_diagram(
    left_label: str,
    right_label: str,
    total_electrons: int,
    ordering: Literal["light_2p", "heavy_2p"] = "heavy_2p",
    include_core_1s: bool = False,
    title: str | None = None,
) -> SemanticDiagram:
    """Construye un diagrama MO diatómico simple."""
    ordering_value = "light_2p" if ordering == "light_2p" else "heavy_2p"
    mo_specs: list[dict[str, Any]] = []
    if include_core_1s:
        mo_specs.extend(
            [
                {"id": "sigma_1s", "label": "σ1s", "energy": 0.0, "degeneracy": 1, "bonding": True},
                {"id": "sigma_star_1s", "label": "σ*1s", "energy": 1.1, "degeneracy": 1, "bonding": False},
            ]
        )
    mo_specs.extend(
        [
            {"id": "sigma_2s", "label": "σ2s", "energy": 2.2, "degeneracy": 1, "bonding": True},
            {"id": "sigma_star_2s", "label": "σ*2s", "energy": 3.3, "degeneracy": 1, "bonding": False},
        ]
    )
    if ordering_value == "light_2p":
        mo_specs.extend(
            [
                {"id": "pi_2p", "label": "π2p", "energy": 4.2, "degeneracy": 2, "bonding": True},
                {"id": "sigma_2p", "label": "σ2p", "energy": 4.8, "degeneracy": 1, "bonding": True},
            ]
        )
    else:
        mo_specs.extend(
            [
                {"id": "sigma_2p", "label": "σ2p", "energy": 4.2, "degeneracy": 1, "bonding": True},
                {"id": "pi_2p", "label": "π2p", "energy": 4.8, "degeneracy": 2, "bonding": True},
            ]
        )
    mo_specs.extend(
        [
            {"id": "pi_star_2p", "label": "π*2p", "energy": 5.9, "degeneracy": 2, "bonding": False},
            {"id": "sigma_star_2p", "label": "σ*2p", "energy": 6.6, "degeneracy": 1, "bonding": False},
        ]
    )

    ao_specs: list[dict[str, Any]] = []
    if include_core_1s:
        ao_specs.append({"id": "1s", "label": "1s", "energy": 0.6, "degeneracy": 1})
    ao_specs.extend(
        [
            {"id": "2s", "label": "2s", "energy": 2.8, "degeneracy": 1},
            {"id": "2p", "label": "2p", "energy": 5.2, "degeneracy": 3},
        ]
    )
    atomic_capacity = sum(int(spec["degeneracy"]) * 2 for spec in ao_specs) * 2
    normalized_total = _normalized_electron_count(total_electrons, capacity=atomic_capacity)
    left_electrons = normalized_total // 2
    right_electrons = normalized_total - left_electrons

    left_fill = _sequential_level_fill(ao_specs, left_electrons)
    right_fill = _sequential_level_fill(ao_specs, right_electrons)
    mo_fill = _sequential_level_fill(mo_specs, normalized_total)

    lanes = [
        DiagramLane(id="left_atom", title=str(left_label or "Left"), x=-190.0),
        DiagramLane(id="molecule", title="Molecule", x=0.0),
        DiagramLane(id="right_atom", title=str(right_label or "Right"), x=190.0),
    ]
    levels: list[DiagramLevel] = []
    connectors: list[DiagramConnector] = []

    for spec in ao_specs:
        level_id = str(spec["id"])
        for lane_id, fill in (("left_atom", left_fill), ("right_atom", right_fill)):
            levels.append(
                DiagramLevel(
                    id=f"{lane_id}_{level_id}",
                    lane_id=lane_id,
                    energy=float(spec["energy"]),
                    label=str(spec["label"]),
                    representation="boxes",
                    degeneracy=int(spec["degeneracy"]),
                    occupancies=list(fill[level_id]),
                    metadata={
                        "orbital_family": level_id,
                        "orbital_origin": "atomic",
                    },
                )
            )

    for spec in mo_specs:
        levels.append(
            DiagramLevel(
                id=str(spec["id"]),
                lane_id="molecule",
                energy=float(spec["energy"]),
                label=str(spec["label"]),
                representation="boxes",
                degeneracy=int(spec["degeneracy"]),
                occupancies=list(mo_fill[str(spec["id"])]),
                metadata={
                    "orbital_origin": "molecular",
                    "bonding": bool(spec["bonding"]),
                },
            )
        )

    if include_core_1s:
        for lane_id in ("left_atom", "right_atom"):
            connectors.extend(
                [
                    DiagramConnector(f"{lane_id}_1s", "sigma_1s", "dashed"),
                    DiagramConnector(f"{lane_id}_1s", "sigma_star_1s", "dashed"),
                ]
            )
    for lane_id in ("left_atom", "right_atom"):
        connectors.extend(
            [
                DiagramConnector(f"{lane_id}_2s", "sigma_2s", "dashed"),
                DiagramConnector(f"{lane_id}_2s", "sigma_star_2s", "dashed"),
                DiagramConnector(f"{lane_id}_2p", "sigma_2p", "dashed"),
                DiagramConnector(f"{lane_id}_2p", "pi_2p", "dashed"),
                DiagramConnector(f"{lane_id}_2p", "pi_star_2p", "dashed"),
                DiagramConnector(f"{lane_id}_2p", "sigma_star_2p", "dashed"),
            ]
        )

    diagram = SemanticDiagram(
        kind="molecular_orbital",
        title=title or f"{left_label}-{right_label} MO Diagram",
        lanes=lanes,
        levels=levels,
        connectors=connectors,
        metadata={
            "total_electrons": normalized_total,
            "ordering": ordering_value,
            "builder": _builder_payload(
                "build_diatomic_mo_diagram",
                {
                    "left_label": str(left_label),
                    "right_label": str(right_label),
                    "total_electrons": normalized_total,
                    "ordering": ordering_value,
                    "include_core_1s": bool(include_core_1s),
                    "title": title,
                },
            ),
        },
    )
    return refresh_semantic_diagram_metadata(diagram)


def build_ligand_field_diagram(
    d_electrons: int,
    geometry: Literal["octahedral", "tetrahedral", "square_planar"] = "octahedral",
    spin_mode: Literal["high", "low"] = "high",
    title: str | None = None,
) -> SemanticDiagram:
    """Construye un diagrama simple de desdoblamiento de campo ligando."""
    geometry_value = (
        geometry
        if geometry in {"octahedral", "tetrahedral", "square_planar"}
        else "octahedral"
    )
    spin_value = "low" if spin_mode == "low" else "high"

    if geometry_value == "octahedral":
        level_specs = [
            {"id": "t2g", "label": "t2g", "energy": 0.0, "degeneracy": 3},
            {"id": "eg", "label": "eg", "energy": 1.6, "degeneracy": 2},
        ]
    elif geometry_value == "tetrahedral":
        level_specs = [
            {"id": "e", "label": "e", "energy": 0.0, "degeneracy": 2},
            {"id": "t2", "label": "t2", "energy": 1.2, "degeneracy": 3},
        ]
    else:
        # Orden razonable para complejos cuadrado-planares:
        # d_xz/d_yz < d_z2 < d_xy < d_x2-y2.
        level_specs = [
            {"id": "d_xz_d_yz", "label": "d_xz / d_yz", "energy": 0.0, "degeneracy": 2},
            {"id": "d_z2", "label": "d_z2", "energy": 1.0, "degeneracy": 1},
            {"id": "d_xy", "label": "d_xy", "energy": 2.0, "degeneracy": 1},
            {"id": "d_x2_y2", "label": "d_x2-y2", "energy": 3.0, "degeneracy": 1},
        ]

    total_electrons = _normalized_electron_count(d_electrons, capacity=10)
    if spin_value == "high":
        occupancies = _high_spin_level_fill(level_specs, total_electrons)
    else:
        occupancies = _sequential_level_fill(level_specs, total_electrons)

    levels = [
        DiagramLevel(
            id=str(spec["id"]),
            lane_id="ligand_field",
            energy=float(spec["energy"]),
            label=str(spec["label"]),
            representation="boxes",
            degeneracy=int(spec["degeneracy"]),
            occupancies=list(occupancies[str(spec["id"])]),
            metadata={"geometry": geometry_value},
        )
        for spec in level_specs
    ]

    diagram = SemanticDiagram(
        kind="ligand_field",
        title=title or f"{geometry_value.replace('_', ' ').title()} d{total_electrons}",
        lanes=[DiagramLane(id="ligand_field", title="", x=0.0)],
        levels=levels,
        connectors=[],
        metadata={
            "d_electrons": total_electrons,
            "geometry": geometry_value,
            "spin_mode": spin_value,
            "builder": _builder_payload(
                "build_ligand_field_diagram",
                {
                    "d_electrons": total_electrons,
                    "geometry": geometry_value,
                    "spin_mode": spin_value,
                    "title": title,
                },
            ),
        },
    )
    return refresh_semantic_diagram_metadata(diagram)
