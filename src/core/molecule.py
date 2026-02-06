"""Compatibilidad legada con un contenedor de grafo molecular.

Este módulo mantiene la clase `Molecule` como envoltorio obsoleto de
`MolGraph` para compatibilidad con código previo.
"""

from core.model import Atom, Bond, MolGraph


class Molecule(MolGraph):
    """
    Envoltorio obsoleto para compatibilidad.

    Use `core.model.MolGraph` en código nuevo.
    """

    def add_atom(self, symbol: str, x: float, y: float, is_explicit: bool = False) -> int:
        """Compatibilidad: delega en `MolGraph.add_atom` y retorna solo el ID.

        Args:
            symbol: Símbolo del elemento químico.
            x: Coordenada X del átomo.
            y: Coordenada Y del átomo.
            is_explicit: Si el símbolo debe mostrarse explícitamente.

        Returns:
            ID del átomo creado.

        Side Effects:
            Modifica el grafo interno mediante `MolGraph.add_atom`.
        """
        atom = super().add_atom(symbol, x, y, is_explicit=is_explicit)
        return atom.id

    def add_bond(self, atom1_idx: int, atom2_idx: int, order: int = 1) -> None:
        """Compatibilidad: delega en `MolGraph.add_bond`.

        Args:
            atom1_idx: ID del primer átomo.
            atom2_idx: ID del segundo átomo.
            order: Orden de enlace.

        Side Effects:
            Modifica el grafo interno mediante `MolGraph.add_bond`.
        """
        super().add_bond(atom1_idx, atom2_idx, order)


__all__ = ["Atom", "Bond", "Molecule"]
