"""Nombres de sustituyentes y cadenas padres para nomenclatura IUPAC-lite."""

from __future__ import annotations

from typing import Dict, Set

from chemcalc.valence import implicit_h_count

from .errors import ChemNameNotSupported
from .molview import MolView
from .rings import RingContext, ring_type_basic

# Mapeo de halógenos a prefijos de sustituyentes.
HALO_MAP: Dict[str, str] = {
    "F": "fluoro",
    "Cl": "chloro",
    "Br": "bromo",
    "I": "iodo",
}

# Nombres base de alcanos (cadena principal).
ALKANE_PARENT: Dict[int, str] = {
    1: "methane",
    2: "ethane",
    3: "propane",
    4: "butane",
    5: "pentane",
    6: "hexane",
    7: "heptane",
    8: "octane",
    9: "nonane",
    10: "decane",
    11: "undecane",
    12: "dodecane",
    13: "tridecane",
    14: "tetradecane",
    15: "pentadecane",
    16: "hexadecane",
    17: "heptadecane",
    18: "octadecane",
    19: "nonadecane",
    20: "eicosane",
}

# Nombres base para cicloalcanos simples.
CYCLO_PARENT: Dict[int, str] = {
    3: "cyclopropane",
    4: "cyclobutane",
    5: "cyclopentane",
    6: "cyclohexane",
    7: "cycloheptane",
    8: "cyclooctane",
}

# Nombres base de sustituyentes alquilo lineales.
ALKYL: Dict[int, str] = {
    1: "methyl",
    2: "ethyl",
    3: "propyl",
    4: "butyl",
    5: "pentyl",
    6: "hexyl",
    7: "heptyl",
    8: "octyl",
    9: "nonyl",
    10: "decyl",
    11: "undecyl",
    12: "dodecyl",
    13: "tridecyl",
    14: "tetradecyl",
    15: "pentadecyl",
    16: "hexadecyl",
    17: "heptadecyl",
    18: "octadecyl",
    19: "nonadecyl",
    20: "eicosyl",
}

# Nombres base de sustituyentes alcoxi lineales.
ALKOXY: Dict[int, str] = {
    1: "methoxy",
    2: "ethoxy",
    3: "propoxy",
}

# Sustituyentes alquilo ramificados soportados explícitamente.
BRANCHED_ALKYL = {
    "isopropyl",
    "sec-butyl",
    "tert-butyl",
    "neopentyl",
}


def alkyl_substituent_name(view: MolView, start_atom: int, chain_set: Set[int]) -> str:
    """Devuelve el nombre del sustituyente alquilo.

    Soporta patrones ramificados específicos y cadenas lineales.

    Args:
        view: Vista del grafo molecular.
        start_atom: Átomo inicial del sustituyente.
        chain_set: Conjunto de átomos de la cadena principal.

    Returns:
        Nombre del sustituyente (p. ej., "ethyl", "isopropyl").

    Raises:
        ChemNameNotSupported: Si el patrón no está soportado.
    """
    branch_nodes, adjacency = _collect_alkyl_branch(view, start_atom, chain_set)
    if not branch_nodes:
        raise ChemNameNotSupported("Unsupported alkyl branch")

    name = _classify_branched_alkyl(start_atom, branch_nodes, adjacency)
    if name is not None:
        return name

    if _is_linear_branch(adjacency):
        length = len(branch_nodes)
        name = ALKYL.get(length)
        if name is None:
            raise ChemNameNotSupported("Unsupported alkyl length")
        return name

    raise ChemNameNotSupported("Unsupported alkyl branch")


def ring_substituent_name(
    view: MolView,
    start_atom: int,
    parent_set: Set[int],
    ring_ctx: RingContext | None,
) -> str | None:
    """Detecta sustituyentes de anillo (fenilo, bencilo, ciclohexilo).

    Args:
        view: Vista del grafo molecular.
        start_atom: Átomo del sustituyente conectado a la cadena principal.
        parent_set: Conjunto de átomos de la cadena principal.
        ring_ctx: Contexto de anillos precomputado.

    Returns:
        Nombre del sustituyente de anillo o `None` si no aplica.
    """
    if ring_ctx is None:
        return None
    for ring in ring_ctx.atom_rings.get(start_atom, []):
        if not ring.isdisjoint(parent_set):
            continue
        rtype = ring_ctx.ring_types.get(ring)
        if rtype == "benzene":
            return "phenyl"
        if rtype == "cyclohexane":
            return "cyclohexyl"

    # bencilo: CH2 unido a un anillo bencénico y a la cadena principal.
    if view.element(start_atom) != "C":
        return None
    if start_atom in ring_ctx.atom_rings:
        return None
    heavy_neighbors = [nbr for nbr in view.neighbors(start_atom) if view.element(nbr) != "H"]
    if len(heavy_neighbors) != 2:
        return None
    parent_neighbor = next((nbr for nbr in heavy_neighbors if nbr in parent_set), None)
    ring_neighbor = next((nbr for nbr in heavy_neighbors if nbr not in parent_set), None)
    if parent_neighbor is None or ring_neighbor is None:
        return None
    for ring in ring_ctx.atom_rings.get(ring_neighbor, []):
        if ring.isdisjoint(parent_set) and ring_ctx.ring_types.get(ring) == "benzene":
            if view.bond_order_between(start_atom, ring_neighbor) == 1:
                return "benzyl"
    return None


def halomethyl_substituent_name(
    view: MolView, start_atom: int, parent_set: Set[int]
) -> str | None:
    """Detecta sustituyentes halometilo (p. ej., clorometilo).

    Args:
        view: Vista del grafo molecular.
        start_atom: Átomo inicial del sustituyente.
        parent_set: Conjunto de la cadena principal.

    Returns:
        Nombre halometilo o `None` si no aplica.
    """
    if view.element(start_atom) != "C":
        return None
    neighbors = [nbr for nbr in view.neighbors(start_atom) if nbr not in parent_set]
    halogens = [nbr for nbr in neighbors if view.element(nbr) in HALO_MAP]
    carbons = [nbr for nbr in neighbors if view.element(nbr) == "C"]
    if len(halogens) == 1 and not carbons:
        return f"{HALO_MAP[view.element(halogens[0])]}methyl"
    return None


def alkoxy_substituent_name(
    view: MolView,
    start_atom: int,
    parent_set: Set[int],
    max_len: int = 3,
) -> str | None:
    """Devuelve el nombre de un sustituyente alcoxi simple.

    Args:
        view: Vista del grafo molecular.
        start_atom: Átomo de oxígeno conectado al anillo/cadena.
        parent_set: Conjunto de la cadena principal.
        max_len: Longitud máxima del sustituyente alcoxi admitida.

    Returns:
        Nombre alcoxi o `None` si no aplica.

    Raises:
        ChemNameNotSupported: Si el patrón de alcoxi no es válido.
    """
    if view.element(start_atom) != "O":
        return None
    heavy_neighbors = [nbr for nbr in view.neighbors(start_atom) if view.element(nbr) != "H"]
    if len(heavy_neighbors) != 2:
        return None
    parent_neighbor = next((nbr for nbr in heavy_neighbors if nbr in parent_set), None)
    alkyl_neighbor = next((nbr for nbr in heavy_neighbors if nbr not in parent_set), None)
    if parent_neighbor is None or alkyl_neighbor is None:
        return None
    if view.bond_order_between(start_atom, parent_neighbor) != 1:
        raise ChemNameNotSupported("Unsupported alkoxy bond")
    if view.bond_order_between(start_atom, alkyl_neighbor) != 1:
        raise ChemNameNotSupported("Unsupported alkoxy bond")
    length = alkyl_length_linear(view, alkyl_neighbor, parent_set | {start_atom})
    if length > max_len:
        raise ChemNameNotSupported("Alkoxy substituent too large")
    name = ALKOXY.get(length)
    if name is None:
        raise ChemNameNotSupported("Unsupported alkoxy substituent")
    return name


def nitro_substituent_name(
    view: MolView, start_atom: int, parent_set: Set[int]
) -> str | None:
    """Detecta sustituyente nitro (-NO2) conectado a la cadena principal.

    Args:
        view: Vista del grafo molecular.
        start_atom: Átomo de nitrógeno del grupo nitro.
        parent_set: Conjunto de la cadena principal.

    Returns:
        "nitro" si coincide, o `None` si no aplica.
    """
    if view.element(start_atom) != "N":
        return None
    if implicit_h_count(view, start_atom) + view.explicit_h(start_atom) != 0:
        return None
    heavy_neighbors = [nbr for nbr in view.neighbors(start_atom) if view.element(nbr) != "H"]
    if len(heavy_neighbors) != 3:
        return None
    if not any(nbr in parent_set for nbr in heavy_neighbors):
        return None
    oxygen_neighbors = [nbr for nbr in heavy_neighbors if view.element(nbr) == "O"]
    if len(oxygen_neighbors) != 2:
        return None
    for o_id in oxygen_neighbors:
        order = view.bond_order_between(start_atom, o_id)
        if order not in {1, 2}:
            raise ChemNameNotSupported("Unsupported nitro substituent")
    return "nitro"


def amino_substituent_name(
    view: MolView, start_atom: int, parent_set: Set[int]
) -> str | None:
    """Detecta sustituyentes amino simples (-NH2, -NHR).

    Args:
        view: Vista del grafo molecular.
        start_atom: Átomo de nitrógeno.
        parent_set: Conjunto de la cadena principal.

    Returns:
        "amino" si coincide, o `None` si no aplica.
    """
    if view.element(start_atom) != "N":
        return None
    heavy_neighbors = [nbr for nbr in view.neighbors(start_atom) if view.element(nbr) != "H"]
    if len(heavy_neighbors) != 1:
        return None
    if heavy_neighbors[0] not in parent_set:
        return None
    if implicit_h_count(view, start_atom) + view.explicit_h(start_atom) < 1:
        return None
    if view.bond_order_between(start_atom, heavy_neighbors[0]) != 1:
        raise ChemNameNotSupported("Unsupported amino substituent")
    return "amino"


def _collect_alkyl_branch(
    view: MolView, start_atom: int, chain_set: Set[int], max_atoms: int = 5
) -> tuple[Set[int], Dict[int, Set[int]]]:
    """Recolecta nodos de un sustituyente alquilo conectado a la cadena.

    Args:
        view: Vista del grafo molecular.
        start_atom: Átomo inicial del sustituyente.
        chain_set: Conjunto de átomos de la cadena principal.
        max_atoms: Límite de tamaño del sustituyente para evitar combinatoria.

    Returns:
        Conjunto de nodos del sustituyente y su adyacencia interna.

    Raises:
        ChemNameNotSupported: Si el sustituyente no es válido o es demasiado grande.
    """
    if start_atom in chain_set:
        raise ChemNameNotSupported("Branch atom is part of the chain")

    branch_nodes: Set[int] = set()
    adjacency: Dict[int, Set[int]] = {}
    stack = [start_atom]
    parent: Dict[int, int | None] = {start_atom: None}
    while stack:
        node = stack.pop()
        if node in branch_nodes:
            continue
        if view.element(node) != "C":
            raise ChemNameNotSupported("Non-carbon in alkyl branch")
        branch_nodes.add(node)
        if len(branch_nodes) > max_atoms:
            raise ChemNameNotSupported("Alkyl branch too large")
        adjacency.setdefault(node, set())
        for nbr in view.neighbors(node):
            if nbr in chain_set:
                continue
            if view.element(nbr) == "H":
                continue
            if view.element(nbr) != "C":
                raise ChemNameNotSupported("Non-carbon in alkyl branch")
            if view.bond_order_between(node, nbr) != 1:
                raise ChemNameNotSupported("Unsaturated bond in alkyl branch")
            adjacency.setdefault(nbr, set()).add(node)
            adjacency[node].add(nbr)
            if nbr not in parent:
                parent[nbr] = node
                stack.append(nbr)
            elif parent[node] != nbr:
                raise ChemNameNotSupported("Alkyl branch cycle detected")
    return branch_nodes, adjacency


def _is_linear_branch(adjacency: Dict[int, Set[int]]) -> bool:
    """Indica si el sustituyente es lineal (sin ramificaciones)."""
    if not adjacency:
        return False
    degrees = [len(neighbors) for neighbors in adjacency.values()]
    if len(adjacency) == 1:
        return True
    return max(degrees) <= 2 and sum(1 for d in degrees if d == 1) == 2


def _classify_branched_alkyl(
    root: int, branch_nodes: Set[int], adjacency: Dict[int, Set[int]]
) -> str | None:
    """Clasifica sustituyentes alquilo ramificados conocidos.

    Args:
        root: Átomo de unión al esqueleto principal.
        branch_nodes: Nodos del sustituyente.
        adjacency: Adyacencia del sustituyente.

    Returns:
        Nombre reconocido (isopropyl, sec-butyl, etc.) o `None`.
    """
    if root not in branch_nodes:
        return None
    degrees = {node: len(neighbors) for node, neighbors in adjacency.items()}
    total = len(branch_nodes)

    # isopropil: carbono central con dos metilos.
    if total == 3 and degrees[root] == 2:
        if all(degrees[node] == 1 for node in branch_nodes if node != root):
            return "isopropyl"

    # tert-butil: carbono central con tres metilos.
    if total == 4 and degrees[root] == 3:
        if all(degrees[node] == 1 for node in branch_nodes if node != root):
            return "tert-butyl"

    # sec-butil: cadena de 4 con sustitución en carbono secundario.
    if total == 4 and degrees[root] == 2:
        neighbors = list(adjacency[root])
        if len(neighbors) != 2:
            return None
        degs = sorted(degrees[n] for n in neighbors)
        if degs == [1, 2]:
            chain_node = neighbors[0] if degrees[neighbors[0]] == 2 else neighbors[1]
            other = [n for n in adjacency[chain_node] if n != root]
            if len(other) == 1 and degrees[other[0]] == 1:
                return "sec-butyl"

    # neopentil: carbono terciario unido a cuatro metilos.
    if total == 5 and degrees[root] == 1:
        neighbor = next(iter(adjacency[root]), None)
        if neighbor is not None and degrees.get(neighbor) == 4:
            if all(
                degrees[node] == 1 for node in adjacency[neighbor] if node != root
            ):
                return "neopentyl"

    return None

def alkane_root(parent: str, use_a: bool = False) -> str:
    """Devuelve la raíz del nombre del alcano (sin 'ane').

    Args:
        parent: Nombre base (p. ej., "propane").
        use_a: Si se debe insertar la vocal 'a' en la raíz.

    Returns:
        Raíz modificada para combinar con insaturaciones.
    """
    if parent.endswith("ane"):
        root = parent[:-3]
    else:
        root = parent
    if use_a and not root.endswith("a"):
        root = f"{root}a"
    return root


def _format_unsaturations(
    unsaturations: list[tuple[int, int]], with_terminal_e: bool
) -> tuple[str, bool]:
    """Formatea insaturaciones y decide si insertar 'a' en la raíz.

    Args:
        unsaturations: Lista de (orden, locante) para dobles/triples.
        with_terminal_e: Si se conserva la 'e' final (p. ej., nitrilo).

    Returns:
        Tupla (descriptor, use_a) con el texto de insaturaciones y flag.
    """
    if not unsaturations:
        return "", False
    doubles = sorted(loc for order, loc in unsaturations if order == 2)
    triples = sorted(loc for order, loc in unsaturations if order == 3)
    if len(doubles) + len(triples) != len(unsaturations):
        raise ChemNameNotSupported("Unsupported bond order")

    # Caso mixto: un doble y un triple permitido.
    if doubles and triples:
        if len(doubles) != 1 or len(triples) != 1:
            raise ChemNameNotSupported("Unsupported unsaturation pattern")
        en_part = "en"
        yn_part = "yne" if with_terminal_e else "yn"
        descriptor = f"{doubles[0]}-{en_part}-{triples[0]}-{yn_part}"
        return descriptor, False

    # Solo dobles enlaces.
    if doubles:
        if len(doubles) > 3:
            raise ChemNameNotSupported("Too many double bonds")
        count = len(doubles)
        locants = ",".join(str(loc) for loc in doubles)
        if count == 1:
            suffix = "ene" if with_terminal_e else "en"
            return f"{locants}-{suffix}", False
        if count == 2:
            suffix = "diene" if with_terminal_e else "dien"
            return f"{locants}-{suffix}", True
        suffix = "triene" if with_terminal_e else "trien"
        return f"{locants}-{suffix}", True

    # Solo triples enlaces.
    if triples:
        if len(triples) > 2:
            raise ChemNameNotSupported("Too many triple bonds")
        count = len(triples)
        locants = ",".join(str(loc) for loc in triples)
        if count == 1:
            suffix = "yne" if with_terminal_e else "yn"
            return f"{locants}-{suffix}", False
        suffix = "diyne" if with_terminal_e else "diyn"
        return f"{locants}-{suffix}", True

    raise ChemNameNotSupported("Unsupported unsaturation pattern")


def parent_name(
    length: int,
    unsaturations: list[tuple[int, int]] | None = None,
    suffix: str | None = None,
    suffix_locant: int | None = None,
) -> str:
    """Construye el nombre del padre con insaturaciones y sufijo.

    Args:
        length: Longitud de la cadena principal.
        unsaturations: Lista de insaturaciones (orden, locante).
        suffix: Sufijo funcional (al, one, ol, etc.).
        suffix_locant: Locante del sufijo si corresponde.

    Returns:
        Nombre del padre con insaturaciones y sufijos aplicados.

    Raises:
        ChemNameNotSupported: Si la longitud o el patrón no están soportados.
    """
    parent = ALKANE_PARENT.get(length)
    if parent is None:
        raise ChemNameNotSupported("Unsupported parent length")

    unsaturations = unsaturations or []
    unsat_descriptor = ""
    use_a = False
    if unsaturations:
        with_terminal_e = suffix is None or suffix == "nitrile"
        unsat_descriptor, use_a = _format_unsaturations(unsaturations, with_terminal_e)

    if suffix is None:
        if not unsaturations:
            return parent
        root = alkane_root(parent, use_a=use_a)
        return f"{root}-{unsat_descriptor}"

    if suffix in {"al", "oic acid", "nitrile"}:
        if suffix_locant not in {None, 1}:
            raise ChemNameNotSupported("Unsupported suffix locant")
        if unsaturations:
            root = alkane_root(parent, use_a=use_a)
            base = f"{root}-{unsat_descriptor}"
        else:
            if suffix == "nitrile":
                base = parent
            else:
                base = parent[:-1] if parent.endswith("e") else parent
        if suffix == "oic acid":
            return f"{base}oic acid"
        return f"{base}{suffix}"

    if suffix_locant is None:
        raise ChemNameNotSupported("Missing suffix locant")

    if unsaturations:
        root = alkane_root(parent, use_a=use_a)
        base = f"{root}-{unsat_descriptor}"
    else:
        base = parent[:-1] if parent.endswith("e") else parent
    return f"{base}-{suffix_locant}-{suffix}"


def alkyl_length_linear(view: MolView, start_atom: int, chain_set: Set[int]) -> int:
    """Calcula la longitud de un sustituyente alquilo lineal.

    Args:
        view: Vista del grafo molecular.
        start_atom: Átomo inicial del sustituyente.
        chain_set: Conjunto de átomos de la cadena principal.

    Returns:
        Longitud del sustituyente lineal.

    Raises:
        ChemNameNotSupported: Si el sustituyente no es lineal o válido.
    """
    if start_atom in chain_set:
        raise ChemNameNotSupported("Branch atom is part of the chain")

    visited: Set[int] = set()
    current = start_atom
    prev = None
    length = 0

    while True:
        if current in visited:
            raise ChemNameNotSupported("Branch cycle detected")
        visited.add(current)

        if view.element(current) != "C":
            raise ChemNameNotSupported("Non-carbon in alkyl branch")

        length += 1
        neighbors = [nbr for nbr in view.neighbors(current) if nbr not in chain_set]
        next_candidates = []
        for nbr in neighbors:
            if nbr == prev:
                continue
            elem = view.element(nbr)
            if elem == "H":
                continue
            if elem != "C":
                raise ChemNameNotSupported("Non-carbon in alkyl branch")
            if view.bond_order_between(current, nbr) != 1:
                raise ChemNameNotSupported("Unsaturated bond in alkyl branch")
            next_candidates.append(nbr)

        if len(next_candidates) > 1:
            raise ChemNameNotSupported("Alkyl branch is branched")
        if not next_candidates:
            break

        prev, current = current, next_candidates[0]

    return length
