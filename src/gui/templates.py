"""
Template builders for common carbohydrate projections.
"""
from __future__ import annotations

import math
from typing import List, Tuple

from core.model import MolGraph


def _add_atom(
    graph: MolGraph,
    element: str,
    x: float,
    y: float,
    *,
    explicit: bool | None = None,
) -> int:
    if explicit is None:
        explicit = element != "C"
    atom = graph.add_atom(element, x, y, is_explicit=explicit)
    return atom.id


def _add_bond(graph: MolGraph, a1_id: int, a2_id: int, order: int = 1) -> None:
    graph.add_bond(a1_id, a2_id, order)


def build_linear_chain_template(bond_length: float) -> MolGraph:
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


def build_haworth_template(bond_length: float) -> MolGraph:
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
        
    # Bonds
    for i in range(6):
        _add_bond(graph, ring_ids[i], ring_ids[(i + 1) % 6])
        
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
    
    subs = [
        # (RingAtomIdx, Element, OffsetY, OffsetX)
        (1, "OH", v_up,   0.0), # C1
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


def build_chair_template(bond_length: float) -> MolGraph:
    graph = MolGraph()
    L = float(bond_length)
    
    # 4C1 Chair Conformation (Zig-Zag)
    # We use a standard skew.
    
    # Ring Coordinates
    r_coords = [
        ( 0.9 * L, -0.9 * L), # 0: O (Top-Rightish)
        ( 1.7 * L,  0.0 * L), # 1: C1 (Right Point)
        ( 0.9 * L,  0.9 * L), # 2: C2 (Bottom-Rightish)
        (-0.9 * L,  0.9 * L), # 3: C3 (Bottom-Leftish)
        (-1.7 * L,  0.0 * L), # 4: C4 (Left Point)
        (-0.9 * L, -0.9 * L), # 5: C5 (Top-Leftish)
    ]
    # Connectivity 0-1-2-3-4-5-0
    
    ring_ids: List[int] = []
    atoms = ["O", "C", "C", "C", "C", "C"]
    
    for i, (x, y) in enumerate(r_coords):
        chain_el = atoms[i]
        # Explicit carbons usually not needed for ring, but helps debugging
        ring_ids.append(_add_atom(graph, chain_el, x, y, explicit=(chain_el!="C")))
        
    for i in range(6):
        _add_bond(graph, ring_ids[i], ring_ids[(i + 1) % 6])
        
    # Substituents (Beta-D-Glucose = All Equatorial)
    # Equatorial bonds point "out" away from the ring center
    
    # Vector logic:
    # C1 (Right): Eq is Right/Down? No, C1 Eq is Up/Right (Beta) or Down/Right (Alpha).
    # Wait, in 4C1 Beta-D-Glc:
    # C1 is Up-Right.
    # C2 is Down-Right.
    # C3 is Up-Left.
    # C4 is Down-Left.
    # C5 is Up-Left (CH2OH).
    
    off_weak = 0.5 * L
    off_strong = 1.0 * L
    
    # (dx, dy)
    sub_offsets = {
        1: ( off_strong, -off_weak), # C1: Right, Up
        2: ( off_strong,  off_weak), # C2: Right, Down
        3: (-off_strong,  off_weak), # C3: Left, Down? No, C3 is usually Left-Up?
                                     # Let's check: C1(R), C2(D), C3(L), C4(U), C5(D)... it alternates.
                                     # Actually Beta-D-Glucose, all equatorial are roughly in plane.
                                     # Let's use visual placement.
        # C3 is at bottom-left (-0.9, 0.9). Eq bond goes Left/Up? Or Left/Down?
        # Standard: C2 Eq is Down-Right. C3 Eq is Up-Left.
        # Wait, if C2 is Down-Right, C3 must be Up-Left to be trans-diaxial? No they are equatorial.
        # Equatorial bonds are roughly parallel to the ring bonds one bond away.
        # C2-Eq || O-C1 and C3-C4.
        # O-C1 is (0.8, 0.9). C3-C4 is (-0.8, -0.9).
        # So C2-Eq should be vector (0.8, 0.9) approx. -> Down Right. Correct.
        
        # C3-Eq || C1-C2 and C4-C5.
        # C1-C2 is (-0.8, -0.9). C4-C5 is (0.8, -0.9).
        # C1-C2 is (-0.8, 0.9). (1.7->0.9, 0.0->0.9).
        # So C3-Eq should be (-0.8, 0.9)? No parallelism is to C1-C2?
        # Actually simplest rule: "Equatorial goes OUT".
        
        3: (-1.0 * L,  0.2 * L), # C3: Left, slightly Down
        4: (-0.8 * L, -0.8 * L), # C4: Left, Up
        # Let's just hardcode what looks good.
    }
    
    # Re-defining visually good offsets
    subs = [
        (1, "OH",     1.0 * L, -0.4 * L), # C1: Right-Up
        (2, "OH",     0.6 * L,  1.0 * L), # C2: Right-Down
        (3, "OH",    -0.6 * L,  1.0 * L), # C3: Left-Down
        (4, "OH",    -1.0 * L, -0.4 * L), # C4: Left-Up (Back)
        (5, "CH2OH", -0.8 * L, -1.0 * L), # C5: Left-Top
    ]
    
    for (idx, el, dx, dy) in subs:
        rx, ry = r_coords[idx]
        sid = _add_atom(graph, el, rx + dx, ry + dy, explicit=True)
        _add_bond(graph, ring_ids[idx], sid)
        
    return graph
