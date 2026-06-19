"""Modelos de datos base del editor molecular Chemuson.

Este módulo concentra las estructuras que representan el grafo molecular
(átomos y enlaces) y el estado químico de la interfaz. El resto de la
aplicación (GUI, persistencia y nomenclatura) interactúa con estas clases
para añadir, modificar y validar la química dibujada.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass, field, replace
from enum import Enum
from typing import Dict, Iterable, List, Optional, Set

# Marcadores internos para distinguir "no se especificó" de "se desea borrar".
_STROKE_UNSET = object()
_COLOR_UNSET = object()
_DONOR_UNSET = object()
_FLEX_CURVE_UNSET = object()
_PI_OFFSET_UNSET = object()
_OPACITY_UNSET = object()
_INTERACTION_KIND_UNSET = object()


class BondStyle(str, Enum):
    """Estilos de representación de enlaces en el lienzo."""
    PLAIN = "plain"
    BOLD = "bold"
    WEDGE = "wedge"
    HASHED = "hashed"
    WAVY = "wavy"
    FLEX = "flex"
    INTERACTION = "interaction"
    COORDINATION = "coordination"


class BondStereo(str, Enum):
    """Categorías de estereoquímica para enlaces dibujados."""
    NONE = "none"
    UP = "up"
    DOWN = "down"
    EITHER = "either"


VALENCE_NEUTRAL_BOND_STYLES = frozenset({BondStyle.COORDINATION, BondStyle.INTERACTION})
NON_STRUCTURAL_BOND_STYLES = frozenset({BondStyle.INTERACTION})


def normalize_bond_style(style: BondStyle | str | None) -> BondStyle:
    """Normaliza estilos de enlace con fallback seguro a `PLAIN`."""
    if isinstance(style, BondStyle):
        return style
    try:
        return BondStyle(style)
    except Exception:
        return BondStyle.PLAIN


def bond_style_affects_valence(style: BondStyle | str | None) -> bool:
    """Indica si un estilo debe contar para valencia/H implícitos."""
    return normalize_bond_style(style) not in VALENCE_NEUTRAL_BOND_STYLES


def bond_style_is_structural(style: BondStyle | str | None) -> bool:
    """Indica si un estilo representa conectividad química estructural."""
    return normalize_bond_style(style) not in NON_STRUCTURAL_BOND_STYLES


def bond_affects_valence(bond: object) -> bool:
    """Indica si un enlace concreto debe aportar a la valencia."""
    return bond_style_affects_valence(getattr(bond, "style", BondStyle.PLAIN))


def bond_is_structural(bond: object) -> bool:
    """Indica si un enlace concreto forma parte de la estructura química."""
    return bond_style_is_structural(getattr(bond, "style", BondStyle.PLAIN))


def normalize_opacity(value: object, default: float = 1.0) -> float:
    """Normaliza una opacidad al rango [0.0, 1.0]."""
    try:
        opacity = float(value)
    except Exception:
        opacity = float(default)
    return max(0.0, min(1.0, opacity))


def normalize_optional_opacity(value: object) -> Optional[float]:
    """Normaliza una opacidad opcional; `None` significa heredar del canvas."""
    if value is None:
        return None
    return normalize_opacity(value)


# Valencias "típicas" usadas por la UI para estimar hidrógenos implícitos.
# Nota: no es un validador químico estricto; se admiten estados cargados
# e hipervalentes en otras partes del flujo.
VALENCE_MAP = {
    "H": 1,
    "C": 4,
    "N": 3,
    "O": 2,
    "S": 2,
    "P": 3,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
}

# Elementos para los que se sugiere autocompletar H implícitos.
IMPLICIT_H_DEFAULT_ELEMENTS = {"B", "C", "N", "O", "Si", "P", "S"}

# Números atómicos usados para cálculo iso-electrónico.
ATOMIC_NUMBERS: Dict[str, int] = {
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
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm": 100,
    "Md": 101,
    "No": 102,
    "Lr": 103,
}

# Valencias permitidas por número atómico iso-electrónico.
# `-1` indica "sin límite" (usado para metales de transición/complejos).
ISO_VALENCE_MAP: Dict[int, tuple[int, ...]] = {
    1: (1,),
    2: (0,),
    3: (1,),
    4: (2,),
    5: (3,),
    6: (4,),
    7: (3, 5),
    8: (2,),
    9: (1,),
    10: (0,),
    11: (1,),
    12: (2,),
    13: (3,),
    14: (4,),
    15: (3, 5),
    16: (2, 4, 6),
    17: (1, 3, 5, 7),
    18: (0,),
    19: (1,),
    20: (2,),
    33: (3, 5),
    34: (2, 4, 6),
    35: (1, 3, 5, 7),
    51: (3, 5),
    52: (2, 4, 6),
    53: (1, 3, 5, 7),
    54: (2, 4, 6, 8),
}

UNLIMITED_VALENCE_ELEMENTS = {
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Lu",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
}

# Excepciones manuales útiles para estados iónicos comunes.
VALENCE_EXCEPTIONS: Dict[tuple[str, int], tuple[int, ...]] = {
    ("P", -2): (3, 5),
    ("S", -1): (2, 6),
}

SIMPLE_HYDROGEN_GROUP_LABELS: Dict[str, tuple[str, int]] = {
    "NH2": ("N", 2),
    "H2N": ("N", 2),
    "NH": ("N", 1),
    "HN": ("N", 1),
    "OH": ("O", 1),
    "HO": ("O", 1),
    "SH": ("S", 1),
    "HS": ("S", 1),
}


@dataclass
class Atom:
    """Representa un átomo en el grafo molecular."""
    id: int
    element: str
    x: float
    y: float
    charge: int = 0
    isotope: Optional[int] = None
    radical_electrons: int = 0
    oxidation_state: Optional[int] = None
    stereo_cip: Optional[str] = None
    stereo_axial: Optional[str] = None
    stereo_helical: Optional[str] = None
    stereo_si_re: Optional[str] = None
    explicit_h: Optional[int] = None
    group_h_cap: Optional[int] = None
    mapping: Optional[int] = None
    is_query: bool = False
    is_explicit: bool = False
    no_implicit: bool = False
    implicit_h: int = 0
    has_valence_error: bool = False
    label_scale: Optional[float] = None
    r_group_substituents: tuple[str, ...] = field(default_factory=tuple)
    is_coordination_center: bool = False
    sphere_radius: Optional[float] = None
    sphere_color: Optional[str] = None
    sphere_filled: bool = True
    sphere_transparent: bool = False
    opacity: Optional[float] = None

    @property
    def formal_charge(self) -> int:
        """Alias semántico para la carga formal del átomo."""
        return int(self.charge)

    @formal_charge.setter
    def formal_charge(self, value: int) -> None:
        self.charge = int(value)


@dataclass(frozen=True)
class ValidationCorrectionAction:
    """Acción sugerida de corrección expuesta como dato para la UI."""

    action_id: str
    label: str
    description: str = ""

    def as_dict(self) -> dict[str, str]:
        return {
            "id": self.action_id,
            "label": self.label,
            "description": self.description,
        }


@dataclass(frozen=True)
class ValidationIssue:
    """Detalle de validación química para un átomo específico."""
    atom_id: int
    element: str
    code: str
    message: str
    severity: str = "error"
    target_type: str = "atom"
    hint: Optional[str] = None
    suggestion: Optional[str] = None
    allowed_valences: tuple[int, ...] = field(default_factory=tuple)
    observed_valence: float = 0.0
    bond_order_sum: float = 0.0
    assigned_h: int = 0
    implicit_h: int = 0
    correction_actions: tuple[ValidationCorrectionAction, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        resolved_hint = self.hint if self.hint is not None else self.suggestion
        resolved_suggestion = self.suggestion if self.suggestion is not None else resolved_hint
        object.__setattr__(self, "hint", resolved_hint)
        object.__setattr__(self, "suggestion", resolved_suggestion)

    def tooltip_text(self) -> str:
        """Texto compacto para tooltips y paneles de inspección."""
        parts = [self.message]
        if self.hint:
            parts.append(f"Sugerencia: {self.hint}")
        return "\n".join(part for part in parts if part)

    def suggested_actions(self) -> tuple[dict[str, str], ...]:
        """Acciones sugeridas en formato de datos, fuera de `as_dict()`."""
        return tuple(action.as_dict() for action in self.correction_actions)

    def with_correction_actions(
        self,
        actions: Iterable[ValidationCorrectionAction],
    ) -> "ValidationIssue":
        """Devuelve una copia con acciones calculadas por un controlador."""
        return replace(self, correction_actions=tuple(actions))

    def as_dict(self) -> dict[str, object]:
        """Representación estable para UI, exportes y pruebas."""
        return {
            "target_type": self.target_type,
            "atom_id": self.atom_id,
            "element": self.element,
            "code": self.code,
            "severity": self.severity,
            "message": self.message,
            "hint": self.hint or "",
            "allowed_valences": list(self.allowed_valences),
            "observed_valence": self.observed_valence,
            "bond_order_sum": self.bond_order_sum,
            "assigned_h": self.assigned_h,
            "implicit_h": self.implicit_h,
        }


@dataclass
class Bond:
    """Representa un enlace químico entre dos átomos."""
    id: int
    a1_id: int
    a2_id: int
    order: int = 1
    style: BondStyle = BondStyle.PLAIN
    stereo: BondStereo = BondStereo.NONE
    stereo_ez: Optional[str] = None
    stereo_axial: Optional[str] = None
    stereo_endo_exo: Optional[str] = None
    stereo_helical: Optional[str] = None
    is_aromatic: bool = False
    display_order: Optional[int] = None
    is_query: bool = False
    ring_id: Optional[int] = None
    length_px: Optional[float] = None
    stroke_px: Optional[float] = None
    color: Optional[str] = None
    donor_atom_id: Optional[int] = None
    flex_curve_1: Optional[float] = None
    flex_curve_2: Optional[float] = None
    # Orientación manual de la línea pi en dobles enlaces: None (auto), +1 o -1.
    pi_offset_sign: Optional[int] = None
    opacity: Optional[float] = None
    interaction_kind: Optional[str] = None


@dataclass
class ChemState:
    """Estado químico y de visualización activo en la interfaz."""
    active_tool: str = "tool_select"
    active_orbital_kind: str = "p_shaded"
    active_energy_diagram_kind: str = "sublevel_s"
    active_bond_order: int = 1
    active_bond_style: BondStyle = BondStyle.PLAIN
    active_bond_stereo: BondStereo = BondStereo.NONE
    active_bond_mode: str = "increment"
    active_bond_aromatic: bool = False
    fixed_angles: bool = True
    fixed_lengths: bool = True
    angle_step_deg: int = 30
    bond_length: float = 40.0
    active_ring_size: int = 6
    active_ring_aromatic: bool = False
    active_ring_template: Optional[str] = None
    active_ring_anomeric: Optional[str] = None
    active_bracket_type: str = "[]"
    default_element: str = "C"
    selected_atoms: Set[int] = field(default_factory=set)
    selected_bonds: Set[int] = field(default_factory=set)
    # Display preferences
    show_implicit_carbons: bool = False  # Show C labels (False = hide like ChemDoodle)
    show_implicit_hydrogens: bool = False  # Show H labels
    use_aromatic_circles: bool = False  # Draw circles in aromatic rings
    label_font_family: str = "Arial"
    label_font_size: float = 11.0
    label_font_bold: bool = False
    label_font_italic: bool = False
    label_font_underline: bool = False
    canvas_opacity: float = 1.0
    use_element_colors: bool = True
    numbering_enabled: bool = False
    numbering_mode: str = "atoms"  # atoms | structures | both
    numbering_font_size: float = 10.0
    numbering_offset_x: float = 8.0
    numbering_offset_y: float = -8.0
    numbering_circle: bool = False
    numbering_color: str = "#1F2933"
    numbering_background: str = "#FFFFFF"
    numbering_include_in_export: bool = True


class MolGraph:
    """Grafo molecular mutable con operaciones de edición básicas."""

    def __init__(self) -> None:
        """Inicializa el grafo vacío y contadores internos de IDs."""
        self.atoms: Dict[int, Atom] = {}
        self.bonds: Dict[int, Bond] = {}
        self._next_atom_id = 1
        self._next_bond_id = 1

    def add_atom(
        self,
        element: str,
        x: float,
        y: float,
        atom_id: Optional[int] = None,
        charge: int = 0,
        formal_charge: Optional[int] = None,
        isotope: Optional[int] = None,
        radical_electrons: int = 0,
        oxidation_state: Optional[int] = None,
        stereo_cip: Optional[str] = None,
        stereo_axial: Optional[str] = None,
        stereo_helical: Optional[str] = None,
        stereo_si_re: Optional[str] = None,
        explicit_h: Optional[int] = None,
        group_h_cap: Optional[int] = None,
        mapping: Optional[int] = None,
        is_query: bool = False,
        is_explicit: bool = False,
        no_implicit: bool = False,
        is_coordination_center: bool = False,
        sphere_radius: Optional[float] = None,
        sphere_color: Optional[str] = None,
        sphere_filled: bool = True,
        sphere_transparent: bool = False,
        label_scale: Optional[float] = None,
        r_group_substituents: Optional[Iterable[str]] = None,
        opacity: Optional[float] = None,
    ) -> Atom:
        """Crea y registra un átomo en el grafo.

        Args:
            element: Símbolo del elemento químico (p. ej., "C", "O").
            x: Posición X en coordenadas del lienzo.
            y: Posición Y en coordenadas del lienzo.
            atom_id: ID explícito si se desea restaurar desde un archivo.
            charge: Alias retrocompatible para la carga formal.
            formal_charge: Carga formal del átomo (prioritaria sobre `charge`).
            isotope: Número másico si se desea mostrar el isótopo.
            radical_electrons: Número de electrones desapareados.
            oxidation_state: Estado de oxidación estimado/asignado.
            stereo_cip: Descriptor CIP absoluto (R/S) cuando se conoce.
            stereo_axial: Descriptor axial (R_a/S_a), best-effort.
            stereo_helical: Descriptor helicoidal (M/P), best-effort.
            stereo_si_re: Descriptor de cara carbonílica (si/re), best-effort.
            explicit_h: Número de hidrógenos explícitos asociados.
            group_h_cap: Límite superior de H para grupos abreviados simples (`NH2`, `OH`, `SH`).
            mapping: Índice de mapeo (útil en exportaciones/reacciones).
            is_query: Marca de átomo de consulta (SMARTS-like).
            is_explicit: Si el símbolo debe mostrarse aunque sea implícito.
            no_implicit: Si se desactiva el ajuste automático de H implícitos.
            is_coordination_center: Si el átomo debe renderizarse como centro de coordinación.
            sphere_radius: Radio visual para el modo esfera (si aplica).
            sphere_color: Color base de la esfera (hex), o `None` para automático.
            sphere_filled: Si la esfera se dibuja con relleno (gradiente).
            sphere_transparent: Si la esfera se dibuja transparente (solo etiqueta).
            label_scale: Escala local de etiqueta o `None` para heredar la global.
            opacity: Opacidad local del átomo o `None` para heredar del canvas.

        Returns:
            El átomo creado y almacenado en el diccionario interno.

        Side Effects:
            Incrementa el contador de IDs y modifica `self.atoms`.
        """
        if atom_id is None:
            atom_id = self._next_atom_id
            self._next_atom_id += 1
        else:
            self._next_atom_id = max(self._next_atom_id, atom_id + 1)
        resolved_sphere_radius = sphere_radius
        if is_coordination_center and resolved_sphere_radius is None:
            resolved_sphere_radius = 16.0
        resolved_sphere_color = sphere_color
        if is_coordination_center and resolved_sphere_color is None:
            resolved_sphere_color = "#D9DDE3"
        resolved_formal_charge = int(charge)
        if formal_charge is not None:
            resolved_formal_charge = int(formal_charge)
        resolved_element = element
        resolved_group_h_cap = group_h_cap
        resolved_explicit_h = explicit_h
        resolved_no_implicit = bool(no_implicit)
        simple_group = self._simple_hydrogen_group_spec(element)
        if simple_group is not None:
            resolved_element = simple_group[0]
            if resolved_group_h_cap is None:
                resolved_group_h_cap = simple_group[1]
            resolved_no_implicit = True
        atom = Atom(
            id=atom_id,
            element=resolved_element,
            x=x,
            y=y,
            charge=resolved_formal_charge,
            isotope=isotope,
            radical_electrons=int(radical_electrons or 0),
            oxidation_state=(
                int(oxidation_state)
                if oxidation_state is not None
                else None
            ),
            stereo_cip=(str(stereo_cip) if stereo_cip else None),
            stereo_axial=(str(stereo_axial) if stereo_axial else None),
            stereo_helical=(str(stereo_helical) if stereo_helical else None),
            stereo_si_re=(str(stereo_si_re) if stereo_si_re else None),
            explicit_h=resolved_explicit_h,
            group_h_cap=(
                None if resolved_group_h_cap is None else max(0, int(resolved_group_h_cap))
            ),
            mapping=mapping,
            is_query=is_query,
            is_explicit=is_explicit,
            no_implicit=resolved_no_implicit,
            label_scale=(None if label_scale is None else float(label_scale)),
            r_group_substituents=tuple(
                str(item).strip()
                for item in (r_group_substituents or ())
                if str(item).strip()
            ),
            is_coordination_center=is_coordination_center,
            sphere_radius=resolved_sphere_radius,
            sphere_color=resolved_sphere_color,
            sphere_filled=bool(sphere_filled),
            sphere_transparent=bool(sphere_transparent),
            opacity=normalize_optional_opacity(opacity),
        )
        self.atoms[atom_id] = atom
        return atom

    def remove_atom(self, atom_id: int) -> tuple[Atom, List[Bond]]:
        """Elimina un átomo y todos los enlaces conectados.

        Args:
            atom_id: Identificador del átomo a eliminar.

        Returns:
            Una tupla con el átomo eliminado y la lista de enlaces removidos.

        Side Effects:
            Modifica `self.atoms` y `self.bonds`, actualizando el grafo.
        """
        atom = self.atoms.pop(atom_id)
        removed_bonds: List[Bond] = []
        for bond_id, bond in list(self.bonds.items()):
            if bond.a1_id == atom_id or bond.a2_id == atom_id:
                removed_bonds.append(self.remove_bond(bond_id))
        return atom, removed_bonds

    def add_bond(
        self,
        a1_id: int,
        a2_id: int,
        order: int = 1,
        bond_id: Optional[int] = None,
        style: BondStyle = BondStyle.PLAIN,
        stereo: BondStereo = BondStereo.NONE,
        stereo_ez: Optional[str] = None,
        stereo_axial: Optional[str] = None,
        stereo_endo_exo: Optional[str] = None,
        stereo_helical: Optional[str] = None,
        is_aromatic: bool = False,
        display_order: Optional[int] = None,
        is_query: bool = False,
        ring_id: Optional[int] = None,
        length_px: Optional[float] = None,
        stroke_px: Optional[float] = None,
        color: Optional[str] = None,
        donor_atom_id: Optional[int] = None,
        flex_curve_1: Optional[float] = None,
        flex_curve_2: Optional[float] = None,
        pi_offset_sign: Optional[int] = None,
        opacity: Optional[float] = None,
        interaction_kind: Optional[str] = None,
    ) -> Bond:
        """Crea y registra un enlace entre dos átomos.

        Args:
            a1_id: ID del primer átomo.
            a2_id: ID del segundo átomo.
            order: Orden de enlace (1, 2, 3).
            bond_id: ID explícito si se restaura desde un archivo.
            style: Estilo visual del enlace.
            stereo: Estereoquímica dibujada (cuña, trazos, etc.).
            stereo_ez: Descriptor geométrico E/Z cuando aplique.
            stereo_axial: Descriptor axial de enlace (R_a/S_a), best-effort.
            stereo_endo_exo: Descriptor endo/exo en sistemas bicíclicos.
            stereo_helical: Descriptor helicoidal asociado al eje de enlace.
            is_aromatic: Marca si el enlace pertenece a un sistema aromático.
            display_order: Orden visual alternativo para dibujado.
            is_query: Indica si es enlace de consulta.
            ring_id: Identificador de anillo (si aplica).
            length_px: Longitud de dibujo fija (px).
            stroke_px: Grosor de línea (px).
            color: Color personalizado del enlace.
            donor_atom_id: ID del átomo donador para enlace coordinativo.
            flex_curve_1: Curvatura normalizada del primer control (estilo FLEX).
            flex_curve_2: Curvatura normalizada del segundo control (estilo FLEX).
            pi_offset_sign: Orientación manual de línea pi (+1/-1) o `None`.
            opacity: Opacidad local del enlace o `None` para heredar del canvas.
            interaction_kind: Semántica opcional para enlaces de interacción
                (`hydrogen_bond`, `pi_pi`, etc.).

        Returns:
            El enlace creado.

        Side Effects:
            Incrementa el contador de IDs y modifica `self.bonds`.
        """
        if bond_id is None:
            bond_id = self._next_bond_id
            self._next_bond_id += 1
        else:
            self._next_bond_id = max(self._next_bond_id, bond_id + 1)
        donor = donor_atom_id
        if style != BondStyle.COORDINATION:
            donor = None
        elif donor not in {a1_id, a2_id}:
            donor = None
        curve_1 = flex_curve_1
        curve_2 = flex_curve_2
        if style != BondStyle.FLEX:
            curve_1 = None
            curve_2 = None
        else:
            curve_1 = float(curve_1) if curve_1 is not None else None
            curve_2 = float(curve_2) if curve_2 is not None else None
        manual_pi: Optional[int] = None
        if pi_offset_sign in {-1, 1}:
            manual_pi = int(pi_offset_sign)
        resolved_interaction_kind = None
        if style in {BondStyle.INTERACTION, BondStyle.COORDINATION} and interaction_kind:
            resolved_interaction_kind = str(interaction_kind).strip().lower() or None

        bond = Bond(
            id=bond_id,
            a1_id=a1_id,
            a2_id=a2_id,
            order=order,
            style=style,
            stereo=stereo,
            stereo_ez=(str(stereo_ez) if stereo_ez else None),
            stereo_axial=(str(stereo_axial) if stereo_axial else None),
            stereo_endo_exo=(str(stereo_endo_exo) if stereo_endo_exo else None),
            stereo_helical=(str(stereo_helical) if stereo_helical else None),
            is_aromatic=is_aromatic,
            display_order=display_order,
            is_query=is_query,
            ring_id=ring_id,
            length_px=length_px,
            stroke_px=stroke_px,
            color=color,
            donor_atom_id=donor,
            flex_curve_1=curve_1,
            flex_curve_2=curve_2,
            pi_offset_sign=manual_pi,
            opacity=normalize_optional_opacity(opacity),
            interaction_kind=resolved_interaction_kind,
        )
        self.bonds[bond_id] = bond
        return bond

    def remove_bond(self, bond_id: int) -> Bond:
        """Elimina un enlace del grafo.

        Args:
            bond_id: Identificador del enlace.

        Returns:
            El enlace eliminado.

        Side Effects:
            Modifica el diccionario `self.bonds`.
        """
        return self.bonds.pop(bond_id)

    def get_atom(self, atom_id: int) -> Atom:
        """Obtiene un átomo por ID.

        Args:
            atom_id: Identificador del átomo.

        Returns:
            El átomo correspondiente.
        """
        return self.atoms[atom_id]

    def get_bond(self, bond_id: int) -> Bond:
        """Obtiene un enlace por ID.

        Args:
            bond_id: Identificador del enlace.

        Returns:
            El enlace correspondiente.
        """
        return self.bonds[bond_id]

    def find_bond_between(self, a1_id: int, a2_id: int) -> Optional[Bond]:
        """Busca un enlace existente entre dos átomos.

        Args:
            a1_id: ID del primer átomo.
            a2_id: ID del segundo átomo.

        Returns:
            El enlace si existe, o `None` en caso contrario.
        """
        for bond in self.bonds.values():
            if {bond.a1_id, bond.a2_id} == {a1_id, a2_id}:
                return bond
        return None

    def update_atom_position(self, atom_id: int, x: float, y: float) -> None:
        """Actualiza la posición de un átomo en el lienzo.

        Args:
            atom_id: Identificador del átomo.
            x: Nueva coordenada X.
            y: Nueva coordenada Y.

        Side Effects:
            Modifica el objeto `Atom` en `self.atoms`.
        """
        atom = self.atoms[atom_id]
        atom.x = x
        atom.y = y

    def update_atom_element(
        self,
        atom_id: int,
        element: str,
        is_explicit: Optional[bool] = None,
    ) -> None:
        """Cambia el elemento químico de un átomo.

        Args:
            atom_id: Identificador del átomo.
            element: Nuevo símbolo del elemento.
            is_explicit: Si se debe forzar la visibilidad del símbolo.

        Side Effects:
            Modifica `Atom.element` y opcionalmente `Atom.is_explicit`.
        """
        atom = self.atoms[atom_id]
        atom.element = element
        if is_explicit is not None:
            atom.is_explicit = is_explicit

    def update_atom_charge(self, atom_id: int, charge: int) -> None:
        """Actualiza la carga formal de un átomo.

        Args:
            atom_id: Identificador del átomo.
            charge: Nueva carga formal.

        Side Effects:
            Modifica `Atom.charge`.
        """
        atom = self.atoms[atom_id]
        atom.formal_charge = int(charge)

    def update_atom_no_implicit(self, atom_id: int, enabled: bool) -> None:
        """Activa/desactiva hidrógenos implícitos para un átomo."""
        atom = self.atoms[atom_id]
        atom.no_implicit = bool(enabled)

    def update_atom_label_scale(self, atom_id: int, label_scale: Optional[float]) -> None:
        """Actualiza la escala local de etiqueta de un átomo."""
        atom = self.atoms[atom_id]
        atom.label_scale = None if label_scale is None else float(label_scale)

    def update_atom_opacity(self, atom_id: int, opacity: Optional[float]) -> None:
        """Actualiza la opacidad local de un átomo."""
        atom = self.atoms[atom_id]
        atom.opacity = normalize_optional_opacity(opacity)

    def update_atom_group_h_cap(self, atom_id: int, group_h_cap: Optional[int]) -> None:
        """Actualiza el límite de H asociado a un grupo abreviado simple."""
        atom = self.atoms[atom_id]
        atom.group_h_cap = (
            None if group_h_cap is None else max(0, int(group_h_cap))
        )

    def atom_degree(self, atom_id: int) -> int:
        """Devuelve el número de enlaces estructurales conectados a un átomo."""
        return sum(
            1
            for bond in self.bonds.values()
            if bond_is_structural(bond)
            and (bond.a1_id == atom_id or bond.a2_id == atom_id)
        )

    def is_hidden_carbon_placeholder(self, atom_id: int) -> bool:
        """Indica si el átomo es un carbono implícito auto-generado y sin semántica propia."""
        atom = self.atoms.get(atom_id)
        if atom is None:
            return False
        if atom.element != "C" or atom.is_explicit:
            return False
        if atom.charge != 0 or atom.isotope is not None:
            return False
        if int(getattr(atom, "radical_electrons", 0) or 0) != 0:
            return False
        if atom.oxidation_state is not None or atom.explicit_h is not None:
            return False
        if getattr(atom, "group_h_cap", None) is not None:
            return False
        if atom.mapping is not None or atom.is_query or atom.no_implicit:
            return False
        if getattr(atom, "label_scale", None) is not None:
            return False
        if getattr(atom, "opacity", None) is not None:
            return False
        if getattr(atom, "is_coordination_center", False):
            return False
        if getattr(atom, "sphere_radius", None) is not None:
            return False
        if getattr(atom, "sphere_color", None) is not None:
            return False
        return True

    def is_disposable_orphan_atom(self, atom_id: int) -> bool:
        """Indica si el átomo puede eliminarse automáticamente al quedar aislado."""
        return self.is_hidden_carbon_placeholder(atom_id) and self.atom_degree(atom_id) <= 0

    def update_bond(
        self,
        bond_id: int,
        order: Optional[int] = None,
        style: Optional[BondStyle] = None,
        stereo: Optional[BondStereo] = None,
        stereo_ez: Optional[str] | object = _COLOR_UNSET,
        stereo_axial: Optional[str] | object = _COLOR_UNSET,
        stereo_endo_exo: Optional[str] | object = _COLOR_UNSET,
        stereo_helical: Optional[str] | object = _COLOR_UNSET,
        is_aromatic: Optional[bool] = None,
        display_order: Optional[int] = None,
        stroke_px: Optional[float] | object = _STROKE_UNSET,
        color: Optional[str] | object = _COLOR_UNSET,
        donor_atom_id: Optional[int] | object = _DONOR_UNSET,
        flex_curve_1: Optional[float] | object = _FLEX_CURVE_UNSET,
        flex_curve_2: Optional[float] | object = _FLEX_CURVE_UNSET,
        pi_offset_sign: Optional[int] | object = _PI_OFFSET_UNSET,
        opacity: Optional[float] | object = _OPACITY_UNSET,
        interaction_kind: Optional[str] | object = _INTERACTION_KIND_UNSET,
    ) -> Bond:
        """Actualiza propiedades de un enlace existente.

        Args:
            bond_id: Identificador del enlace a modificar.
            order: Nuevo orden de enlace.
            style: Estilo visual del enlace.
            stereo: Estereoquímica dibujada.
            stereo_ez: Descriptor E/Z; `None` limpia el valor.
            stereo_axial: Descriptor axial; `None` limpia el valor.
            stereo_endo_exo: Descriptor endo/exo; `None` limpia el valor.
            stereo_helical: Descriptor helicoidal; `None` limpia el valor.
            is_aromatic: Marca de aromaticidad.
            display_order: Orden visual alternativo.
            stroke_px: Grosor de línea; `None` limpia el valor.
            color: Color del enlace; `None` limpia el valor.
            donor_atom_id: ID del donador; `None` limpia el valor.
            flex_curve_1: Curvatura normalizada de control 1; `None` limpia.
            flex_curve_2: Curvatura normalizada de control 2; `None` limpia.
            pi_offset_sign: Orientación manual de línea pi (+1/-1), `None` limpia.
            opacity: Opacidad local del enlace; `None` hereda del canvas.
            interaction_kind: Semántica opcional para interacciones; `None` limpia.

        Returns:
            El enlace actualizado.

        Side Effects:
            Modifica el objeto `Bond` dentro de `self.bonds`.
        """
        bond = self.bonds[bond_id]
        if order is not None:
            bond.order = order
        if style is not None:
            bond.style = style
        if stereo is not None:
            bond.stereo = stereo
        if stereo_ez is not _COLOR_UNSET:
            bond.stereo_ez = None if stereo_ez is None else str(stereo_ez)
        if stereo_axial is not _COLOR_UNSET:
            bond.stereo_axial = None if stereo_axial is None else str(stereo_axial)
        if stereo_endo_exo is not _COLOR_UNSET:
            bond.stereo_endo_exo = None if stereo_endo_exo is None else str(stereo_endo_exo)
        if stereo_helical is not _COLOR_UNSET:
            bond.stereo_helical = None if stereo_helical is None else str(stereo_helical)
        if is_aromatic is not None:
            bond.is_aromatic = is_aromatic
        if display_order is not None:
            bond.display_order = display_order
        if stroke_px is not _STROKE_UNSET:
            bond.stroke_px = None if stroke_px is None else float(stroke_px)
        if color is not _COLOR_UNSET:
            bond.color = None if color is None else str(color)
        if donor_atom_id is not _DONOR_UNSET:
            bond.donor_atom_id = int(donor_atom_id) if donor_atom_id is not None else None
        if flex_curve_1 is not _FLEX_CURVE_UNSET:
            bond.flex_curve_1 = None if flex_curve_1 is None else float(flex_curve_1)
        if flex_curve_2 is not _FLEX_CURVE_UNSET:
            bond.flex_curve_2 = None if flex_curve_2 is None else float(flex_curve_2)
        if pi_offset_sign is not _PI_OFFSET_UNSET:
            if pi_offset_sign in {-1, 1}:
                bond.pi_offset_sign = int(pi_offset_sign)
            else:
                bond.pi_offset_sign = None
        if opacity is not _OPACITY_UNSET:
            bond.opacity = normalize_optional_opacity(opacity)
        if interaction_kind is not _INTERACTION_KIND_UNSET:
            bond.interaction_kind = None if interaction_kind is None else str(interaction_kind).strip().lower()
        if bond.style != BondStyle.COORDINATION:
            bond.donor_atom_id = None
        elif bond.donor_atom_id not in {bond.a1_id, bond.a2_id}:
            bond.donor_atom_id = None
        if bond.style != BondStyle.FLEX:
            bond.flex_curve_1 = None
            bond.flex_curve_2 = None
        if bond.style not in {BondStyle.INTERACTION, BondStyle.COORDINATION}:
            bond.interaction_kind = None
        return bond

    def update_bond_length(self, bond_id: int, length_px: Optional[float]) -> None:
        """Ajusta la longitud de dibujo del enlace.

        Args:
            bond_id: Identificador del enlace.
            length_px: Longitud en píxeles o `None` para usar la automática.

        Side Effects:
            Modifica el atributo `Bond.length_px`.
        """
        bond = self.bonds[bond_id]
        bond.length_px = length_px

    def clear(self) -> None:
        """Elimina todos los átomos y enlaces del grafo.

        Side Effects:
            Limpia `self.atoms`, `self.bonds` y reinicia contadores.
        """
        self.atoms.clear()
        self.bonds.clear()
        self._next_atom_id = 1
        self._next_bond_id = 1

    @property
    def formal_charge(self) -> int:
        """Carga formal total de la molécula."""
        return int(sum(int(atom.formal_charge) for atom in self.atoms.values()))

    def total_formal_charge(self) -> int:
        """Alias explícito para la carga formal total."""
        return self.formal_charge

    @staticmethod
    def _simple_hydrogen_group_spec(label: str) -> Optional[tuple[str, int]]:
        """Resuelve alias simples como `NH2`, `OH` o `SH` a elemento + cupo H."""
        cleaned = str(label or "").strip()
        if not cleaned:
            return None
        spec = SIMPLE_HYDROGEN_GROUP_LABELS.get(cleaned)
        if spec is not None:
            return spec
        return SIMPLE_HYDROGEN_GROUP_LABELS.get(cleaned.upper())

    def _effective_atom_element(self, atom: Atom) -> str:
        """Devuelve el elemento químico efectivo, incluyendo alias simples legacy."""
        spec = self._simple_hydrogen_group_spec(atom.element)
        if spec is not None:
            return spec[0]
        return atom.element

    def _effective_group_h_cap(self, atom: Atom) -> Optional[int]:
        """Devuelve el cupo H efectivo, incluso si el átomo viene de un alias legacy."""
        cap = getattr(atom, "group_h_cap", None)
        if cap is not None:
            return max(0, int(cap))
        spec = self._simple_hydrogen_group_spec(atom.element)
        if spec is None:
            return None
        return max(0, int(spec[1]))

    def _allowed_valences_for_atom(self, atom: Atom) -> List[int]:
        """Resuelve la lista de valencias permitidas para un átomo."""
        if bool(getattr(atom, "is_coordination_center", False)):
            return [-1]
        element = self._effective_atom_element(atom)
        if element in UNLIMITED_VALENCE_ELEMENTS:
            return [-1]

        charge = int(getattr(atom, "formal_charge", 0) or 0)
        explicit_override = VALENCE_EXCEPTIONS.get((element, charge))
        if explicit_override is not None:
            return list(explicit_override)

        atomic_number = ATOMIC_NUMBERS.get(element)
        if atomic_number is None:
            return []
        iso_z = atomic_number - charge
        if iso_z <= 0:
            return []

        valences = ISO_VALENCE_MAP.get(iso_z)
        if valences is not None:
            return list(valences)

        typical = VALENCE_MAP.get(element)
        if typical is not None:
            return [typical]
        return []

    def _bond_order_contribution(self, bond: Bond, aromatic_order: float = 1.5) -> float:
        """Contribución de un enlace a la suma de órdenes de valencia."""
        if not bond_affects_valence(bond):
            return 0.0
        if bond.is_aromatic:
            return float(aromatic_order)
        if bond.order <= 0:
            return 1.0
        return float(bond.order)

    def bond_order_sum(self, atom_id: int, aromatic_order: float = 1.5) -> float:
        """Suma de órdenes de enlace para un átomo."""
        total = 0.0
        for bond in self.bonds.values():
            if bond.a1_id != atom_id and bond.a2_id != atom_id:
                continue
            total += self._bond_order_contribution(bond, aromatic_order=aromatic_order)
        return total

    def _group_hydrogen_target(self, atom_id: int, atom: Atom) -> int:
        """Calcula H asignados dinámicamente para grupos simples como `NH2`."""
        cap = self._effective_group_h_cap(atom)
        if cap is None:
            return int(atom.explicit_h or 0)
        allowed = self._allowed_valences_for_atom(atom)
        if not allowed or -1 in allowed:
            return cap
        allowed_nonnegative = sorted(v for v in allowed if v >= 0)
        if not allowed_nonnegative:
            return 0
        preferred_valence = allowed_nonnegative[0]
        missing = float(preferred_valence) - float(self.bond_order_sum(atom_id, aromatic_order=1.5))
        if missing <= 1e-6:
            return 0
        rounded = int(round(missing))
        if rounded < 0:
            return 0
        if not math.isclose(float(rounded), missing, abs_tol=1e-4):
            return 0
        return min(cap, rounded)

    def assigned_hydrogen_count(self, atom_id: int) -> int:
        """Cuenta H asignados en la etiqueta del átomo, excluyendo átomos H vecinos."""
        atom = self.atoms.get(atom_id)
        if atom is None:
            return 0
        return max(0, int(self._group_hydrogen_target(atom_id, atom)))

    def explicit_hydrogen_count(self, atom_id: int) -> int:
        """Cuenta hidrógenos explícitos dibujados/asignados al átomo."""
        atom = self.atoms.get(atom_id)
        if atom is None:
            return 0
        count = self.assigned_hydrogen_count(atom_id)
        for bond in self.bonds.values():
            if not bond_affects_valence(bond):
                continue
            if bond.a1_id == atom_id:
                other_id = bond.a2_id
            elif bond.a2_id == atom_id:
                other_id = bond.a1_id
            else:
                continue
            other = self.atoms.get(other_id)
            if other is not None and other.element == "H":
                count += 1
        return max(0, int(count))

    def _aromatic_neighbors(self, atom_id: int) -> List[int]:
        """Devuelve vecinos conectados por enlaces aromáticos no coordinativos."""
        neighbors: List[int] = []
        for bond in self.bonds.values():
            if not bond_affects_valence(bond) or not bond.is_aromatic:
                continue
            if bond.a1_id == atom_id:
                neighbors.append(bond.a2_id)
            elif bond.a2_id == atom_id:
                neighbors.append(bond.a1_id)
        return neighbors

    def _shortest_aromatic_path_length(
        self,
        start_atom_id: int,
        end_atom_id: int,
        blocked_atom_id: int,
        max_edges: int,
    ) -> Optional[int]:
        """Busca la ruta aromática más corta entre dos vecinos excluyendo un átomo."""
        if start_atom_id == end_atom_id:
            return 0
        queue = deque([(start_atom_id, 0)])
        visited = {blocked_atom_id, start_atom_id}
        while queue:
            current_atom_id, distance = queue.popleft()
            if distance >= max_edges:
                continue
            for neighbor_id in self._aromatic_neighbors(current_atom_id):
                if neighbor_id == blocked_atom_id:
                    continue
                next_distance = distance + 1
                if neighbor_id == end_atom_id:
                    return next_distance
                if next_distance >= max_edges or neighbor_id in visited:
                    continue
                visited.add(neighbor_id)
                queue.append((neighbor_id, next_distance))
        return None

    def _is_pyrrolic_like_aromatic_n(self, atom_id: int, atom: Atom) -> bool:
        """Detecta N aromático neutro de cinco miembros que conserva un H implícito."""
        if self._effective_atom_element(atom) != "N" or int(getattr(atom, "formal_charge", 0) or 0) != 0:
            return False
        if bool(getattr(atom, "no_implicit", False)):
            return False
        if self.explicit_hydrogen_count(atom_id) > 0:
            return False

        aromatic_neighbors = self._aromatic_neighbors(atom_id)
        if len(aromatic_neighbors) != 2:
            return False

        valence_degree = 0
        for bond in self.bonds.values():
            if not bond_affects_valence(bond):
                continue
            if bond.a1_id == atom_id or bond.a2_id == atom_id:
                valence_degree += 1
        if valence_degree != 2:
            return False

        path_length = self._shortest_aromatic_path_length(
            aromatic_neighbors[0],
            aromatic_neighbors[1],
            blocked_atom_id=atom_id,
            max_edges=4,
        )
        return path_length == 3

    def _choose_implicit_h(self, atom: Atom, base_valence: float, allowed: List[int]) -> int:
        """Calcula H implícitos para llevar al estado de menor valencia permitida."""
        if self._effective_group_h_cap(atom) is not None:
            return 0
        if self._effective_atom_element(atom) not in IMPLICIT_H_DEFAULT_ELEMENTS:
            return 0
        if bool(getattr(atom, "no_implicit", False)):
            return 0
        allowed_nonnegative = sorted(v for v in allowed if v >= 0)
        if not allowed_nonnegative:
            return 0
        preferred_valence = allowed_nonnegative[0]
        missing = float(preferred_valence) - float(base_valence)
        if missing <= 1e-6:
            return 0
        rounded = int(round(missing))
        if rounded < 0:
            return 0
        if not math.isclose(float(rounded), missing, abs_tol=1e-4):
            return 0
        return rounded

    def implicit_h_count(self, atom_id: int) -> int:
        """Hidrógenos implícitos recomendados para un átomo."""
        atom = self.atoms.get(atom_id)
        if atom is None:
            return 0
        allowed = self._allowed_valences_for_atom(atom)
        if not allowed or -1 in allowed:
            atom.implicit_h = 0
            return 0
        if self._is_pyrrolic_like_aromatic_n(atom_id, atom):
            atom.implicit_h = 1
            return 1
        base_valence = (
            self.bond_order_sum(atom_id, aromatic_order=1.5)
            + self.assigned_hydrogen_count(atom_id)
        )
        implicit_h = self._choose_implicit_h(atom, base_valence, allowed)
        atom.implicit_h = int(max(0, implicit_h))
        return atom.implicit_h

    def validate(self) -> List[int]:
        """Valida valencias por átomo usando valencias iso-electrónicas.

        Returns:
            Lista de IDs de átomos con valencia no permitida.

        Side Effects:
            Actualiza `Atom.implicit_h` y `Atom.has_valence_error`.
        """
        errors: List[int] = []
        self.validate_detailed()
        for atom_id, atom in self.atoms.items():
            if atom.has_valence_error:
                errors.append(atom_id)
        return errors

    def validate_detailed(self) -> Dict[int, ValidationIssue]:
        """Valida valencias y devuelve detalles explicativos por átomo con error."""
        issues: Dict[int, ValidationIssue] = {}
        for atom_id, atom in self.atoms.items():
            allowed = self._allowed_valences_for_atom(atom)
            if not allowed:
                atom.implicit_h = 0
                atom.has_valence_error = False
                continue
            if -1 in allowed:
                atom.implicit_h = 0
                atom.has_valence_error = False
                continue

            bond_sum = self.bond_order_sum(atom_id)
            assigned_h = self.assigned_hydrogen_count(atom_id)
            base_valence = bond_sum + float(assigned_h)
            implicit_h = self._choose_implicit_h(atom, base_valence, allowed)
            total_valence = base_valence + float(implicit_h)
            is_valid = any(math.isclose(total_valence, float(v), abs_tol=1e-6) for v in allowed)

            atom.implicit_h = int(max(0, implicit_h))
            atom.has_valence_error = not is_valid
            if atom.has_valence_error:
                allowed_display = tuple(sorted(v for v in allowed if v >= 0))
                allowed_text = ", ".join(str(v) for v in allowed_display) or "N/D"
                min_allowed = min(allowed_display) if allowed_display else None
                max_allowed = max(allowed_display) if allowed_display else None
                preferred_allowed = allowed_display[0] if allowed_display else None
                if (
                    preferred_allowed is not None
                    and total_valence > float(preferred_allowed)
                ) or (max_allowed is not None and total_valence > float(max_allowed)):
                    code = "VALENCE_EXCEEDED"
                    hint = (
                        "Reduzca el orden de enlace, quite un H explícito/asignado "
                        "o aumente la carga formal si corresponde."
                    )
                elif min_allowed is not None and total_valence < float(min_allowed):
                    code = "VALENCE_UNDERFILLED"
                    hint = (
                        "Revise H implícitos, H explícitos o protonación para "
                        "alcanzar una valencia permitida."
                    )
                else:
                    code = "VALENCE_MISMATCH"
                    hint = (
                        "Revise aromaticidad, orden de enlace, carga formal o "
                        "protonación para ajustar la valencia."
                    )
                message = (
                    f"{atom.element} con valencia observada {total_valence:.2f} "
                    f"(enlaces {bond_sum:.2f} + H asignados {assigned_h} + "
                    f"H implícitos {int(max(0, implicit_h))}); "
                    f"valencias permitidas: {allowed_text}."
                )
                actions: list[ValidationCorrectionAction] = []
                if code == "VALENCE_EXCEEDED":
                    actions.extend(
                        [
                            ValidationCorrectionAction(
                                "reduce_selected_bond",
                                "Reducir enlace seleccionado",
                                "Reduce en 1 el orden del enlace seleccionado si mejora la valencia.",
                            ),
                            ValidationCorrectionAction(
                                "adjust_charge",
                                "Ajustar carga formal",
                                "Prueba +1/-1 y aplica solo una opción inequívoca.",
                            ),
                        ]
                    )
                if assigned_h > 0:
                    actions.append(
                        ValidationCorrectionAction(
                            "clear_assigned_h",
                            "Limpiar H asignados",
                            "Quita H asociados a la etiqueta si mejora la validación.",
                        )
                    )
                issues[atom_id] = ValidationIssue(
                    atom_id=atom_id,
                    element=atom.element,
                    code=code,
                    message=message,
                    severity="error",
                    hint=hint,
                    allowed_valences=allowed_display,
                    observed_valence=float(total_valence),
                    bond_order_sum=float(bond_sum),
                    assigned_h=int(assigned_h),
                    implicit_h=int(max(0, implicit_h)),
                    correction_actions=tuple(actions),
                )
        return issues
