"""Modelos de datos base del editor molecular Chemuson.

Este módulo concentra las estructuras que representan el grafo molecular
(átomos y enlaces) y el estado químico de la interfaz. El resto de la
aplicación (GUI, persistencia y nomenclatura) interactúa con estas clases
para añadir, modificar y validar la química dibujada.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional, Set

# Marcadores internos para distinguir "no se especificó" de "se desea borrar".
_STROKE_UNSET = object()
_COLOR_UNSET = object()


class BondStyle(str, Enum):
    """Estilos de representación de enlaces en el lienzo."""
    PLAIN = "plain"
    BOLD = "bold"
    WEDGE = "wedge"
    HASHED = "hashed"
    WAVY = "wavy"
    INTERACTION = "interaction"


class BondStereo(str, Enum):
    """Categorías de estereoquímica para enlaces dibujados."""
    NONE = "none"
    UP = "up"
    DOWN = "down"
    EITHER = "either"


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

# Valencias máximas (suma de órdenes de enlace) antes de marcar error.
# Se usa un umbral permisivo para patrones comunes hipervalentes:
# - P(V/VI): fosfatos, fosforanos, PF6-
# - S(IV/VI): sulfóxidos/sulfonas/sulfatos, SF6
# - Halógenos(III/V/VII): interhalógenos, oxoácidos (p. ej., IF7, ClO4-)
# - Xe(II/IV/VI/VIII): fluoruro de xenón y XeO4 en dibujos
MAX_VALENCE_MAP = {
    "H": 1,
    "C": 4,
    "N": 4,
    "O": 3,
    "F": 1,
    "Cl": 7,
    "Br": 7,
    "I": 7,
    "P": 6,
    "S": 6,
    "Xe": 8,
    "Se": 6,
    "Te": 6,
    "As": 6,
    "Sb": 6,
    "Bi": 6,
    "Si": 6,
    "Ge": 6,
    "Sn": 6,
    "Pb": 6,
    "B": 4,
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
    explicit_h: Optional[int] = None
    mapping: Optional[int] = None
    is_query: bool = False
    is_explicit: bool = False


@dataclass
class Bond:
    """Representa un enlace químico entre dos átomos."""
    id: int
    a1_id: int
    a2_id: int
    order: int = 1
    style: BondStyle = BondStyle.PLAIN
    stereo: BondStereo = BondStereo.NONE
    is_aromatic: bool = False
    display_order: Optional[int] = None
    is_query: bool = False
    ring_id: Optional[int] = None
    length_px: Optional[float] = None
    stroke_px: Optional[float] = None
    color: Optional[str] = None


@dataclass
class ChemState:
    """Estado químico y de visualización activo en la interfaz."""
    active_tool: str = "tool_atom"
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
    use_element_colors: bool = True


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
        isotope: Optional[int] = None,
        explicit_h: Optional[int] = None,
        mapping: Optional[int] = None,
        is_query: bool = False,
        is_explicit: bool = False,
    ) -> Atom:
        """Crea y registra un átomo en el grafo.

        Args:
            element: Símbolo del elemento químico (p. ej., "C", "O").
            x: Posición X en coordenadas del lienzo.
            y: Posición Y en coordenadas del lienzo.
            atom_id: ID explícito si se desea restaurar desde un archivo.
            charge: Carga formal del átomo.
            isotope: Número másico si se desea mostrar el isótopo.
            explicit_h: Número de hidrógenos explícitos asociados.
            mapping: Índice de mapeo (útil en exportaciones/reacciones).
            is_query: Marca de átomo de consulta (SMARTS-like).
            is_explicit: Si el símbolo debe mostrarse aunque sea implícito.

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
        atom = Atom(
            id=atom_id,
            element=element,
            x=x,
            y=y,
            charge=charge,
            isotope=isotope,
            explicit_h=explicit_h,
            mapping=mapping,
            is_query=is_query,
            is_explicit=is_explicit,
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
        is_aromatic: bool = False,
        display_order: Optional[int] = None,
        is_query: bool = False,
        ring_id: Optional[int] = None,
        length_px: Optional[float] = None,
        stroke_px: Optional[float] = None,
        color: Optional[str] = None,
    ) -> Bond:
        """Crea y registra un enlace entre dos átomos.

        Args:
            a1_id: ID del primer átomo.
            a2_id: ID del segundo átomo.
            order: Orden de enlace (1, 2, 3).
            bond_id: ID explícito si se restaura desde un archivo.
            style: Estilo visual del enlace.
            stereo: Estereoquímica dibujada (cuña, trazos, etc.).
            is_aromatic: Marca si el enlace pertenece a un sistema aromático.
            display_order: Orden visual alternativo para dibujado.
            is_query: Indica si es enlace de consulta.
            ring_id: Identificador de anillo (si aplica).
            length_px: Longitud de dibujo fija (px).
            stroke_px: Grosor de línea (px).
            color: Color personalizado del enlace.

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
        bond = Bond(
            id=bond_id,
            a1_id=a1_id,
            a2_id=a2_id,
            order=order,
            style=style,
            stereo=stereo,
            is_aromatic=is_aromatic,
            display_order=display_order,
            is_query=is_query,
            ring_id=ring_id,
            length_px=length_px,
            stroke_px=stroke_px,
            color=color,
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
        atom.charge = charge

    def update_bond(
        self,
        bond_id: int,
        order: Optional[int] = None,
        style: Optional[BondStyle] = None,
        stereo: Optional[BondStereo] = None,
        is_aromatic: Optional[bool] = None,
        display_order: Optional[int] = None,
        stroke_px: Optional[float] | object = _STROKE_UNSET,
        color: Optional[str] | object = _COLOR_UNSET,
    ) -> Bond:
        """Actualiza propiedades de un enlace existente.

        Args:
            bond_id: Identificador del enlace a modificar.
            order: Nuevo orden de enlace.
            style: Estilo visual del enlace.
            stereo: Estereoquímica dibujada.
            is_aromatic: Marca de aromaticidad.
            display_order: Orden visual alternativo.
            stroke_px: Grosor de línea; `None` limpia el valor.
            color: Color del enlace; `None` limpia el valor.

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
        if is_aromatic is not None:
            bond.is_aromatic = is_aromatic
        if display_order is not None:
            bond.display_order = display_order
        if stroke_px is not _STROKE_UNSET:
            bond.stroke_px = None if stroke_px is None else float(stroke_px)
        if color is not _COLOR_UNSET:
            bond.color = None if color is None else str(color)
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

    def validate(self) -> List[int]:
        """Valida valencias máximas según `MAX_VALENCE_MAP`.

        Calcula la suma de órdenes de enlace por átomo y reporta aquellos
        que superan la valencia máxima permitida.

        Returns:
            Lista de IDs de átomos que exceden la valencia permitida.

        Side Effects:
            No tiene efectos laterales; solo calcula y devuelve resultados.
        """
        bond_order_sum: Dict[int, int] = {atom_id: 0 for atom_id in self.atoms}
        for bond in self.bonds.values():
            if bond.a1_id in bond_order_sum:
                bond_order_sum[bond.a1_id] += bond.order
            if bond.a2_id in bond_order_sum:
                bond_order_sum[bond.a2_id] += bond.order

        errors: List[int] = []
        for atom_id, atom in self.atoms.items():
            expected = MAX_VALENCE_MAP.get(atom.element)
            if expected is None:
                continue
            if bond_order_sum.get(atom_id, 0) > expected:
                errors.append(atom_id)
        return errors
