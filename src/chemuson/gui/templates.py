"""
Constructores de plantillas para proyecciones de carbohidratos.
"""
from __future__ import annotations

import math
from typing import List, Tuple

from chemuson.core.model import BondStyle, MolGraph


def _add_atom(
    graph: MolGraph,
    element: str,
    x: float,
    y: float,
    *,
    explicit: bool | None = None,
) -> int:
    """Añade un átomo al grafo y retorna su ID.

    Args:
        graph: Grafo molecular destino.
        element: Símbolo del elemento o etiqueta (p. ej., "OH").
        x: Coordenada X.
        y: Coordenada Y.
        explicit: Forzar visibilidad del símbolo.

    Returns:
        ID del átomo creado.
    """
    if explicit is None:
        explicit = element != "C"
    atom = graph.add_atom(element, x, y, is_explicit=explicit)
    return atom.id


def _add_bond(
    graph: MolGraph,
    a1_id: int,
    a2_id: int,
    order: int = 1,
    style: BondStyle = BondStyle.PLAIN,
) -> None:
    """Añade un enlace al grafo con estilo opcional."""
    graph.add_bond(a1_id, a2_id, order, style=style)


def _front_edge_indices(coords: List[Tuple[float, float]]) -> set[int]:
    """Devuelve índices de aristas "frontales" para resaltar en negrita.

    Se usa para dar profundidad a proyecciones (Haworth/silla).
    """
    if not coords:
        return set()
    ys = [y for _x, y in coords]
    center_y = sum(ys) / len(ys)
    height = max(ys) - min(ys)
    threshold = height * 0.05
    n = len(coords)
    bold: set[int] = set()
    for i in range(n):
        j = (i + 1) % n
        mid_y = (coords[i][1] + coords[j][1]) / 2.0
        if mid_y > center_y + threshold:
            bold.add(i)
    return bold


def build_linear_chain_template(bond_length: float) -> MolGraph:
    """Construye una plantilla de cadena lineal tipo Fischer (glucosa).

    Args:
        bond_length: Longitud base de enlace para el escalado.

    Returns:
        `MolGraph` con la proyección lineal y sustituyentes.
    """
    graph = MolGraph()
    length = float(bond_length)
    x0 = 0.0

    # Fischer Projection: Vertical Carbon backbone
    # Increased spacing to prevent label overlap
    y_positions = [
        -3.5 * length, # CHO
        -2.0 * length, # C2
        -0.7 * length, # C3
        0.7 * length,  # C4
        2.0 * length,  # C5
        3.5 * length,  # CH2OH
    ]
    chain_ids: List[int] = []

    # C1 (CHO)
    chain_ids.append(_add_atom(graph, "CHO", x0, y_positions[0], explicit=True))
    
    # C2-C5 Backbone
    for y in y_positions[1:5]:
        # Explicit Carbons ensure the "cross" is drawn
        chain_ids.append(_add_atom(graph, "C", x0, y, explicit=True))
        
    # C6 (CH2OH)
    chain_ids.append(_add_atom(graph, "CH2OH", x0, y_positions[5], explicit=True))

    # Connect backbone
    for i in range(len(chain_ids) - 1):
        _add_bond(graph, chain_ids[i], chain_ids[i + 1])

    substituent_offset = 1.2 * length
    # D-glucose Fischer: OH right, left, right, right
    # (OH direction, H direction) where 1=Right, -1=Left
    orientations = [(1, -1), (-1, 1), (1, -1), (1, -1)]
    
    # Add substituents to C2-C5
    for idx, (oh_dir, h_dir) in enumerate(orientations):
        base_id = chain_ids[idx + 1] # C2 is at index 1
        y = y_positions[idx + 1]
        
        oh_id = _add_atom(graph, "OH", x0 + oh_dir * substituent_offset, y, explicit=True)
        h_id = _add_atom(graph, "H", x0 + h_dir * substituent_offset, y, explicit=True)
        
        _add_bond(graph, base_id, oh_id)
        _add_bond(graph, base_id, h_id)

    return graph


def build_haworth_template(
    bond_length: float,
    *,
    anomeric_up: bool = True,
    bold_front: bool = True,
) -> MolGraph:
    """Construye una plantilla de anillo tipo Haworth (piranosa).

    Args:
        bond_length: Longitud base de enlace.
        anomeric_up: Si el OH anomérico apunta hacia arriba.
        bold_front: Si se resaltan los enlaces frontales.

    Returns:
        `MolGraph` con la proyección Haworth.
    """
    graph = MolGraph()
    L = float(bond_length)
    
    # Standard Haworth Geometry (Flattened Hexagon)
    # Coordinates normalized to approx L=1, centered
    # O is Top-Right (Back-Right)
    
    # Geometry Definitions:
    #      O(0)  -------  C1(1)
    #    /                 \
    #  C5(5)               C2(2)
    #    \                 /
    #     C4(4) ------- C3(3)
    
    # Vertices (relative to center)
    w = 0.8 * L # half-width top/bottom
    w_mid = 1.3 * L # half-width middle
    h_top = 0.5 * L # height top
    h_bot = 0.7 * L # height bottom
    
    # Note: Y increases downwards in screen coordinates
    # Top points have negative Y, Bottom points have positive Y
    
    points = [
        ( 0.4 * L, -h_top), # 0: O (Back-Right) - slightly indented? No, usually symmetric with C5.
        ( w_mid,    0.0),   # 1: C1 (Right)
        ( w,        h_bot), # 2: C2 (Front-Right)
        (-w,        h_bot), # 3: C3 (Front-Left)
        (-w_mid,    0.0),   # 4: C4 (Left)
        (-0.4 * L, -h_top), # 5: C5 (Back-Left)
    ]
    
    # Adjust O and C5 to form the "back" line properly
    # Actually Haworth often draws Furanose/Pyranose with straight horizontal back/front lines.
    # Pyranose (6-ring):
    #       O  ------ C1
    #     /            \
    #   C5              C2
    #     \            /
    #      C4 ------ C3
    # Let's try this specific "straight top/bottom" look which is very clean.
    
    pts = [
        ( 0.6 * L, -0.6 * L), # O
        ( 1.4 * L, -0.2 * L), # C1 (Right-Topish?) No, Haworth usually has points at sides.
    ]
    
    # Let's restart with the user-provided layout in the image (trapezoid-like)
    # Image 2 in Prompt: C1 right, O top-right segment.
    
    # Refined Coordinates for Haworth
    coords = [
        ( 0.5 * L, -0.7 * L), # 0: O
        ( 1.1 * L, -0.2 * L), # 1: C1
        ( 0.6 * L,  0.8 * L), # 2: C2
        (-0.6 * L,  0.8 * L), # 3: C3
        (-1.1 * L, -0.2 * L), # 4: C4
        (-0.5 * L, -0.7 * L), # 5: C5
    ] # Connect 0-1-2-3-4-5-0

    ring_ids: List[int] = []
    atoms = ["O", "C", "C", "C", "C", "C"]
    
    for i, (x, y) in enumerate(coords):
        chain_el = atoms[i]
        ring_ids.append(_add_atom(graph, chain_el, x, y, explicit=(chain_el!="C")))
        
    bold_edges = _front_edge_indices(coords) if bold_front else set()
    for i in range(6):
        style = BondStyle.BOLD if i in bold_edges else BondStyle.PLAIN
        _add_bond(graph, ring_ids[i], ring_ids[(i + 1) % 6], style=style)
        
    # Substituents
    # Vertical offsets
    v_up = -1.2 * L
    v_down = 1.2 * L
    
    # Beta-D-Glucose orientations
    # C1: OH Up (Beta)
    # C2: OH Down
    # C3: OH Up
    # C4: OH Down
    # C5: CH2OH Up
    
    anomeric_offset = v_up if anomeric_up else v_down
    subs = [
        # (RingAtomIdx, Element, OffsetY, OffsetX)
        (1, "OH", anomeric_offset,   0.0), # C1
        (2, "OH", v_down, 0.0), # C2
        (3, "OH", v_up,   0.0), # C3
        (4, "OH", v_down, 0.0), # C4
        (5, "CH2OH", v_up, 0.0), # C5
    ]
    
    for (idx, el, dy, dx) in subs:
        rx, ry = coords[idx]
        sid = _add_atom(graph, el, rx + dx, ry + dy, explicit=True)
        _add_bond(graph, ring_ids[idx], sid)
        
    return graph


def build_chair_template(
    bond_length: float,
    *,
    anomeric_up: bool = True,
    bold_front: bool = True,
) -> MolGraph:
    """Construye una plantilla de ciclohexano en conformación silla.

    Args:
        bond_length: Longitud base de enlace.
        anomeric_up: Parámetro legado; se conserva por compatibilidad.
        bold_front: Si se resaltan los enlaces frontales.

    Returns:
        `MolGraph` con un anillo de 6 carbonos en forma de silla.
    """
    graph = MolGraph()
    L = float(bond_length)

    # Compatibilidad con API previa (beta/alpha no aplica a una silla simple).
    _ = anomeric_up

    # Geometría "chair" de ciclohexano.
    # Coordenadas base normalizadas (longitud media de enlace ~= 1.0).
    base_coords = [
        (-1.35, 0.35),   # 0: extremo inferior izquierdo
        (-0.80, -0.55),  # 1: extremo superior izquierdo
        (0.15, -0.25),   # 2: valle superior central
        (1.10, -0.55),   # 3: extremo superior derecho
        (0.55, 0.35),    # 4: extremo inferior derecho
        (-0.35, 0.05),   # 5: valle inferior central
    ]
    coords = [(x * L, y * L) for (x, y) in base_coords]

    ring_ids: List[int] = []
    for x, y in coords:
        ring_ids.append(_add_atom(graph, "C", x, y, explicit=False))

    bold_edges = _front_edge_indices(coords) if bold_front else set()
    for i in range(6):
        style = BondStyle.BOLD if i in bold_edges else BondStyle.PLAIN
        _add_bond(graph, ring_ids[i], ring_ids[(i + 1) % 6], style=style)

    return graph
