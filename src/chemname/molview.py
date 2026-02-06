"""Vista de lectura sobre grafos moleculares heterogéneos.

`MolView` actúa como adaptador para distintos formatos de grafo, ofreciendo
una API uniforme a los módulos de nomenclatura.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Set, Tuple

from .errors import ChemNameNotSupported


class MolView:
    """Adaptador ligero sobre objetos tipo `MolGraph`."""

    def __init__(self, graph) -> None:
        """Inicializa la vista sobre un grafo arbitrario.

        Args:
            graph: Objeto con atributos compatibles (atoms/bonds/nodes/edges).

        Side Effects:
            Inicializa cachés internas para adyacencias y órdenes de enlace.
        """
        self.graph = graph
        self._adj: Optional[Dict[int, Set[int]]] = None
        self._bond_orders: Optional[Dict[frozenset[int], int]] = None
        self._bond_aromatic: Optional[Dict[frozenset[int], bool]] = None
        self._bond_order_list: Optional[List[int]] = None

    def atoms(self) -> List[int]:
        """Devuelve la lista de IDs atómicos disponibles.

        Returns:
            Lista de IDs de átomos detectados en el grafo.

        Raises:
            ChemNameNotSupported: Si el grafo no expone átomos accesibles.
        """
        graph = self.graph
        atoms_attr = getattr(graph, "atoms", None)
        if atoms_attr is not None:
            atoms_iter = atoms_attr() if callable(atoms_attr) else atoms_attr
            if isinstance(atoms_iter, dict):
                return list(atoms_iter.keys())
            return [self._atom_id(atom) for atom in atoms_iter]

        nodes_attr = getattr(graph, "nodes", None)
        if nodes_attr is not None:
            nodes_iter = nodes_attr() if callable(nodes_attr) else nodes_attr
            if isinstance(nodes_iter, dict):
                return list(nodes_iter.keys())
            return list(nodes_iter)

        get_atoms = getattr(graph, "get_atoms", None)
        if callable(get_atoms):
            return [self._atom_id(atom) for atom in get_atoms()]

        raise ChemNameNotSupported("Cannot read atoms from graph")

    def element(self, atom_id: int) -> str:
        """Obtiene el símbolo del elemento para un átomo.

        Args:
            atom_id: Identificador del átomo.

        Returns:
            Símbolo químico (por ejemplo, "C", "O").

        Raises:
            ChemNameNotSupported: Si no se puede resolver el elemento.
        """
        atom = self._get_atom(atom_id)
        if atom is None:
            raise ChemNameNotSupported("Cannot resolve atom element")
        if isinstance(atom, dict):
            element = atom.get("element") or atom.get("symbol")
            if element:
                return str(element)
        element = getattr(atom, "element", None)
        if element is not None:
            return str(element)
        symbol = getattr(atom, "symbol", None)
        if symbol is not None:
            return str(symbol)
        raise ChemNameNotSupported("Cannot resolve atom element")

    def neighbors(self, atom_id: int) -> List[int]:
        """Devuelve la lista de vecinos conectados al átomo.

        Args:
            atom_id: Identificador del átomo.

        Returns:
            Lista de IDs vecinos.
        """
        self._ensure_adjacency()
        return list(self._adj.get(atom_id, set()))

    def bond_order_between(self, atom_id_1: int, atom_id_2: int) -> int:
        """Obtiene el orden de enlace entre dos átomos.

        Args:
            atom_id_1: ID del primer átomo.
            atom_id_2: ID del segundo átomo.

        Returns:
            Orden de enlace (1 por defecto si no se conoce).
        """
        self._ensure_adjacency()
        return self._bond_orders.get(frozenset({atom_id_1, atom_id_2}), 1)

    def bond_is_aromatic(self, atom_id_1: int, atom_id_2: int) -> bool:
        """Indica si el enlace entre dos átomos es aromático.

        Args:
            atom_id_1: ID del primer átomo.
            atom_id_2: ID del segundo átomo.

        Returns:
            `True` si el enlace está marcado como aromático.
        """
        self._ensure_adjacency()
        return self._bond_aromatic.get(frozenset({atom_id_1, atom_id_2}), False)

    def bonds(self) -> List[Tuple[int, int, int]]:
        """Lista de enlaces como tuplas (a1, a2, orden)."""
        return list(self._iter_bonds())

    def bond_orders(self) -> List[int]:
        """Lista de órdenes de enlace presentes en el grafo."""
        self._ensure_adjacency()
        return list(self._bond_order_list or [])

    def is_acyclic(self) -> bool:
        """Determina si el grafo es acíclico (sin anillos).

        Returns:
            `True` si no se detectan ciclos, `False` en caso contrario.
        """
        atoms = self.atoms()
        visited: Set[int] = set()
        for start in atoms:
            if start in visited:
                continue
            stack: List[Tuple[int, Optional[int]]] = [(start, None)]
            while stack:
                node, parent = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                for nbr in self.neighbors(node):
                    if nbr == parent:
                        continue
                    if nbr in visited:
                        return False
                    stack.append((nbr, node))
        return True

    def _atom_id(self, atom) -> int:
        """Normaliza un objeto átomo a su ID entero.

        Args:
            atom: Objeto átomo o ID directo.

        Returns:
            Identificador entero del átomo.

        Raises:
            ChemNameNotSupported: Si no se puede extraer el ID.
        """
        if isinstance(atom, int):
            return atom
        atom_id = getattr(atom, "id", None)
        if atom_id is not None:
            return int(atom_id)
        atom_id = getattr(atom, "atom_id", None)
        if atom_id is not None:
            return int(atom_id)
        raise ChemNameNotSupported("Atom object has no id")

    def atom_charge(self, atom_id: int) -> int:
        """Obtiene la carga formal del átomo.

        Args:
            atom_id: Identificador del átomo.

        Returns:
            Carga formal (0 por defecto).
        """
        atom = self._get_atom(atom_id)
        if atom is None:
            return 0
        if isinstance(atom, dict):
            charge = atom.get("charge")
            return int(charge) if charge is not None else 0
        charge = getattr(atom, "charge", None)
        return int(charge) if charge is not None else 0

    def explicit_h(self, atom_id: int) -> int:
        """Obtiene el número de H explícitos asociados al átomo.

        Args:
            atom_id: Identificador del átomo.

        Returns:
            Número de hidrógenos explícitos (0 si no hay).
        """
        atom = self._get_atom(atom_id)
        if atom is None:
            return 0
        if isinstance(atom, dict):
            value = atom.get("explicit_h")
            return int(value) if value is not None else 0
        value = getattr(atom, "explicit_h", None)
        return int(value) if value is not None else 0

    def _get_atom(self, atom_id: int):
        """Resuelve un átomo desde el grafo usando múltiples convenciones.

        Args:
            atom_id: Identificador del átomo a recuperar.

        Returns:
            Objeto átomo o `None` si no se encuentra.
        """
        graph = self.graph
        get_atom = getattr(graph, "get_atom", None)
        if callable(get_atom):
            try:
                return get_atom(atom_id)
            except Exception:
                pass

        atoms_attr = getattr(graph, "atoms", None)
        if atoms_attr is not None:
            atoms_iter = atoms_attr() if callable(atoms_attr) else atoms_attr
            if isinstance(atoms_iter, dict):
                return atoms_iter.get(atom_id)
            for atom in atoms_iter:
                if self._atom_id(atom) == atom_id:
                    return atom

        nodes_attr = getattr(graph, "nodes", None)
        if nodes_attr is not None:
            try:
                return nodes_attr[atom_id]
            except Exception:
                pass

        return None

    def _iter_bonds(self) -> Iterable[Tuple[int, int, int]]:
        """Itera enlaces sin metadatos adicionales."""
        for a1, a2, order, _is_aromatic in self._iter_bonds_with_meta():
            yield a1, a2, order

    def _iter_bonds_with_meta(self) -> Iterable[Tuple[int, int, int, bool]]:
        """Itera enlaces con orden y aromaticidad si están disponibles."""
        graph = self.graph
        bonds_attr = getattr(graph, "bonds", None)
        if bonds_attr is not None and not callable(bonds_attr):
            bonds_iter = bonds_attr.values() if isinstance(bonds_attr, dict) else bonds_attr
            for bond in bonds_iter:
                a1 = getattr(bond, "a1_id", None)
                a2 = getattr(bond, "a2_id", None)
                if a1 is None or a2 is None:
                    a1 = getattr(bond, "a1", None)
                    a2 = getattr(bond, "a2", None)
                if a1 is None or a2 is None:
                    if isinstance(bond, (tuple, list)) and len(bond) >= 2:
                        a1, a2 = bond[0], bond[1]
                        order = bond[2] if len(bond) >= 3 else 1
                        yield int(a1), int(a2), int(order) if order is not None else 1, False
                    continue
                order = getattr(bond, "order", None)
                if order is None:
                    order = getattr(bond, "bond_order", None)
                if order is None:
                    order = 1
                is_aromatic = getattr(bond, "is_aromatic", False)
                yield int(a1), int(a2), int(order), bool(is_aromatic)
            return

        edges_attr = getattr(graph, "edges", None)
        if edges_attr is not None:
            edges_iter = edges_attr() if callable(edges_attr) else edges_attr
            for edge in edges_iter:
                if isinstance(edge, (tuple, list)) and len(edge) >= 2:
                    a1, a2 = edge[0], edge[1]
                    order = 1
                    is_aromatic = False
                    if len(edge) >= 3 and isinstance(edge[2], dict):
                        order = edge[2].get("order", 1)
                        is_aromatic = edge[2].get("is_aromatic", False)
                    yield int(a1), int(a2), int(order) if order is not None else 1, bool(is_aromatic)
            return

    def _ensure_adjacency(self) -> None:
        """Construye y cachea adyacencias y órdenes de enlace."""
        if self._adj is not None:
            return
        adjacency: Dict[int, Set[int]] = {atom_id: set() for atom_id in self.atoms()}
        bond_orders: Dict[frozenset[int], int] = {}
        bond_aromatic: Dict[frozenset[int], bool] = {}
        bond_order_list: List[int] = []

        bonds = list(self._iter_bonds_with_meta())
        if bonds:
            # La fuente principal es la lista de enlaces estructurados.
            for a1, a2, order, is_aromatic in bonds:
                adjacency.setdefault(a1, set()).add(a2)
                adjacency.setdefault(a2, set()).add(a1)
                bond_orders[frozenset({a1, a2})] = order
                bond_aromatic[frozenset({a1, a2})] = bool(is_aromatic)
                bond_order_list.append(order)
        else:
            # Ruta de respaldo si el grafo solo expone vecinos.
            neighbors_attr = getattr(self.graph, "neighbors", None)
            if neighbors_attr is not None:
                for atom_id in list(adjacency.keys()):
                    nbrs = neighbors_attr(atom_id) if callable(neighbors_attr) else neighbors_attr.get(atom_id, [])
                    for nbr in nbrs:
                        adjacency.setdefault(atom_id, set()).add(int(nbr))
                        adjacency.setdefault(int(nbr), set()).add(atom_id)

        self._adj = adjacency
        self._bond_orders = bond_orders
        self._bond_aromatic = bond_aromatic
        self._bond_order_list = bond_order_list
