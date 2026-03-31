"""Presets y utilidades para diagramas de niveles de energia."""

from __future__ import annotations

from dataclasses import dataclass

from PyQt6.QtCore import QSizeF


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
            "fill_visible": True,
        }
    if kind == "mo_2s_2p":
        return {
            "connector_color": "#7A7A7A",
            "box_stroke_visible": False,
            "fill_visible": False,
        }
    return {
        "box_stroke_visible": True,
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
