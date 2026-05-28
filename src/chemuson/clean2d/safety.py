from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any


@dataclass
class Clean2DQualityReport:
    """Reporte de calidad comparando geometría antes/después de una limpieza."""
    atom_ids: set[int]
    before: dict[int, tuple[float, float]]
    after: dict[int, tuple[float, float]]
    target_bond_length: float
    mean_bond_length_before: float = 0.0
    mean_bond_length_after: float = 0.0
    bond_length_ratio: float = 1.0
    max_displacement: float = 0.0
    min_nonbonded_before: float = 0.0
    min_nonbonded_after: float = 0.0
    new_crossings: int = 0
    ring_degeneracy_before: float = 0.0
    ring_degeneracy_after: float = 0.0
    bounding_box_ratio: float = 1.0
    missing_coords: list[int] = field(default_factory=list)
    nan_coords: list[int] = field(default_factory=list)
    mean_displacement: float = 0.0
    is_cyclic: bool = False
    passed: bool = False
    rejection_reason: str = ""


def _bond_a1(bond: Any) -> int:
    return bond.a1_id if hasattr(bond, "a1_id") else bond.get("a1_id", 0)

def _bond_a2(bond: Any) -> int:
    return bond.a2_id if hasattr(bond, "a2_id") else bond.get("a2_id", 0)

def has_cycles(atom_ids: set[int], bonds: list[Any]) -> bool:
    if not atom_ids:
        return False
    adjacency: dict[int, list[int]] = {}
    for bond in bonds:
        a1 = _bond_a1(bond)
        a2 = _bond_a2(bond)
        if a1 not in atom_ids or a2 not in atom_ids:
            continue
        adjacency.setdefault(a1, []).append(a2)
        adjacency.setdefault(a2, []).append(a1)
    visited: set[int] = set()
    for start in atom_ids:
        if start in visited:
            continue
        stack: list[tuple[int, int]] = [(start, -1)]
        while stack:
            node, parent = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            for neigh in adjacency.get(node, []):
                if neigh == parent:
                    continue
                if neigh in visited:
                    return True
                stack.append((neigh, node))
    return False


def _find_cycles(atom_ids: set[int], bonds: list[Any]) -> list[set[int]]:
    adjacency: dict[int, set[int]] = {}
    for bond in bonds:
        a1 = _bond_a1(bond)
        a2 = _bond_a2(bond)
        if a1 not in atom_ids or a2 not in atom_ids:
            continue
        adjacency.setdefault(a1, set()).add(a2)
        adjacency.setdefault(a2, set()).add(a1)
    rings: list[set[int]] = []
    edges_processed: set[tuple[int, int]] = set()
    for a in sorted(adjacency):
        for b in list(adjacency.get(a, [])):
            key = (min(a, b), max(a, b))
            if key in edges_processed:
                continue
            edges_processed.add(key)
            path = _shortest_path(adjacency, a, b, key)
            if path and len(path) >= 3:
                ring = set(path)
                if not any(ring.issuperset(existing) for existing in rings):
                    rings.append(ring)
    return rings


def _shortest_path(
    adjacency: dict[int, set[int]],
    start: int,
    end: int,
    blocked: tuple[int, int],
) -> list[int] | None:
    visited: set[int] = set()
    queue: list[tuple[int, list[int]]] = [(start, [start])]
    while queue:
        node, path = queue.pop(0)
        if node == end and len(path) > 1:
            return path
        if node in visited:
            continue
        visited.add(node)
        for neigh in adjacency.get(node, []):
            edge = (min(node, neigh), max(node, neigh))
            if edge == blocked:
                continue
            if neigh in visited:
                continue
            queue.append((neigh, path + [neigh]))
    return None


def evaluate_clean2d_layout(
    atom_ids: set[int],
    bonds: list[Any],
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    target_len: float,
    *,
    is_cyclic: bool | None = None,
) -> Clean2DQualityReport:
    if is_cyclic is None:
        is_cyclic = has_cycles(atom_ids, bonds)

    report = Clean2DQualityReport(
        atom_ids=atom_ids,
        before=before,
        after=after,
        target_bond_length=max(1.0, target_len),
        is_cyclic=is_cyclic,
    )

    intersection = atom_ids & after.keys()
    missing = [aid for aid in atom_ids if aid not in after]
    report.missing_coords = missing
    nan_list = [aid for aid in intersection if _is_nan(after[aid])]
    report.nan_coords = nan_list
    if missing or nan_list:
        report.rejection_reason = "faltan_coordenadas"
        return report

    filtered_bonds = [
        b for b in bonds
        if _bond_a1(b) in intersection
        and _bond_a2(b) in intersection
    ]

    if filtered_bonds:
        report.mean_bond_length_before = _mean_bond_length(before, filtered_bonds)
        report.mean_bond_length_after = _mean_bond_length(after, filtered_bonds)
        if report.mean_bond_length_before > 1e-6:
            report.bond_length_ratio = (
                report.mean_bond_length_after / report.mean_bond_length_before
            )

    report.max_displacement = max_atom_displacement(before, after, intersection)
    displacements = [
        math.hypot(
            after[aid][0] - before[aid][0],
            after[aid][1] - before[aid][1],
        )
        for aid in intersection
    ]
    report.mean_displacement = sum(displacements) / len(displacements) if displacements else 0.0

    report.min_nonbonded_before = min_nonbonded_distance(before, filtered_bonds, intersection)
    report.min_nonbonded_after = min_nonbonded_distance(after, filtered_bonds, intersection)

    report.new_crossings = count_new_bond_crossings(before, after, filtered_bonds)

    if is_cyclic:
        cycles = _find_cycles(intersection, filtered_bonds)
        if cycles:
            report.ring_degeneracy_before = min(
                ring_degeneracy_score(before, ring) for ring in cycles
            )
            report.ring_degeneracy_after = min(
                ring_degeneracy_score(after, ring) for ring in cycles
            )

    report.bounding_box_ratio = _compute_bounding_box_ratio(before, after, intersection)

    return report


def is_clean2d_candidate_safe(report: Clean2DQualityReport, mode: str = "quick") -> bool:
    if report.missing_coords or report.nan_coords:
        report.passed = False
        report.rejection_reason = report.rejection_reason or "coordenadas_invalidas"
        return False

    t = max(1.0, report.target_bond_length)

    if report.mean_bond_length_after < t * 0.3 or report.mean_bond_length_after > t * 3.0:
        report.passed = False
        report.rejection_reason = "longitud_enlace_fuera_de_rango"
        return False

    max_disp = report.max_displacement
    disp_limit = t * 4.0 if mode == "publication" else t * 3.0
    if max_disp > disp_limit:
        report.passed = False
        report.rejection_reason = f"desplazamiento_maximo_excesivo:{max_disp:.1f}>{disp_limit:.1f}"
        return False

    if report.min_nonbonded_after < float("inf") and report.min_nonbonded_after < t * 0.25:
        report.passed = False
        report.rejection_reason = "colision_no_enlazada"
        return False

    if report.new_crossings > 0:
        report.passed = False
        report.rejection_reason = f"nuevos_cruces_enlaces:{report.new_crossings}"
        return False

    if report.is_cyclic and report.ring_degeneracy_after < 0.05 * t * t:
        report.passed = False
        report.rejection_reason = "anillo_colapsado"
        return False

    if report.bounding_box_ratio > 5.0 or report.bounding_box_ratio < 0.2:
        report.passed = False
        report.rejection_reason = f"cambio_caja_absurdo:{report.bounding_box_ratio:.2f}"
        return False

    if report.mean_displacement > t * 0.01 and report.mean_bond_length_before > 1e-6:
        bl_ratio = report.mean_bond_length_after / report.mean_bond_length_before
        if bl_ratio > 2.0 or bl_ratio < 0.5:
            report.passed = False
            report.rejection_reason = (
                f"cambio_longitud_enlace_demasiado_abrupto:{bl_ratio:.2f}"
            )
            return False

    report.passed = True
    report.rejection_reason = ""
    return True


def count_new_bond_crossings(
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    bonds: list[Any],
) -> int:
    before_crossings = _count_crossings(before, bonds)
    after_crossings = _count_crossings(after, bonds)
    return max(0, after_crossings - before_crossings)


def _count_crossings(
    positions: dict[int, tuple[float, float]],
    bonds: list[Any],
) -> int:
    crossings = 0
    for i, b1 in enumerate(bonds):
        a1 = _bond_a1(b1)
        a2 = _bond_a2(b1)
        if a1 not in positions or a2 not in positions:
            continue
        for b2 in bonds[i + 1:]:
            b3 = _bond_a1(b2)
            b4 = _bond_a2(b2)
            if b3 not in positions or b4 not in positions:
                continue
            shared = {a1, a2} & {b3, b4}
            if shared:
                continue
            if _segments_intersect(
                positions[a1], positions[a2],
                positions[b3], positions[b4],
            ):
                crossings += 1
    return crossings


def _segments_intersect(
    p1: tuple[float, float],
    p2: tuple[float, float],
    p3: tuple[float, float],
    p4: tuple[float, float],
) -> bool:
    def orient(a, b, c) -> float:
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    def on_segment(a, b, c) -> bool:
        return (
            min(a[0], b[0]) <= c[0] <= max(a[0], b[0])
            and min(a[1], b[1]) <= c[1] <= max(a[1], b[1])
        )

    o1 = orient(p1, p2, p3)
    o2 = orient(p1, p2, p4)
    o3 = orient(p3, p4, p1)
    o4 = orient(p3, p4, p2)

    if o1 == 0 and on_segment(p1, p2, p3):
        return False
    if o2 == 0 and on_segment(p1, p2, p4):
        return False
    if o3 == 0 and on_segment(p3, p4, p1):
        return False
    if o4 == 0 and on_segment(p3, p4, p2):
        return False

    return (o1 > 0) != (o2 > 0) and (o3 > 0) != (o4 > 0)


def min_nonbonded_distance(
    positions: dict[int, tuple[float, float]],
    bonds: list[Any],
    atom_ids: set[int],
) -> float:
    bonded_pairs: set[tuple[int, int]] = set()
    for bond in bonds:
        a1 = _bond_a1(bond)
        a2 = _bond_a2(bond)
        if a1 in atom_ids and a2 in atom_ids:
            bonded_pairs.add((min(a1, a2), max(a1, a2)))

    ids = sorted(atom_ids & set(positions.keys()))
    min_dist = float("inf")
    for i, a_id in enumerate(ids):
        ax, ay = positions[a_id]
        for b_id in ids[i + 1:]:
            if (min(a_id, b_id), max(a_id, b_id)) in bonded_pairs:
                continue
            bx, by = positions[b_id]
            dist = math.hypot(bx - ax, by - ay)
            if dist < min_dist:
                min_dist = dist
    return min_dist if min_dist < float("inf") else float("inf")


def ring_degeneracy_score(
    positions: dict[int, tuple[float, float]],
    ring_atoms: set[int],
) -> float:
    nodes = [n for n in ring_atoms if n in positions]
    if len(nodes) < 3:
        return 0.0
    area = _polygon_area([positions[n] for n in nodes])
    perimeter = sum(
        math.hypot(
            positions[nodes[i]][0] - positions[nodes[(i + 1) % len(nodes)]][0],
            positions[nodes[i]][1] - positions[nodes[(i + 1) % len(nodes)]][1],
        )
        for i in range(len(nodes))
    )
    if perimeter < 1e-6:
        return 0.0
    return area / (perimeter * perimeter / (4.0 * math.pi))


def _polygon_area(vertices: list[tuple[float, float]]) -> float:
    n = len(vertices)
    if n < 3:
        return 0.0
    area = 0.0
    for i in range(n):
        x1, y1 = vertices[i]
        x2, y2 = vertices[(i + 1) % n]
        area += x1 * y2 - x2 * y1
    return abs(area) * 0.5


def max_atom_displacement(
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    atom_ids: set[int],
) -> float:
    max_d = 0.0
    for aid in atom_ids:
        if aid not in before or aid not in after:
            continue
        d = math.hypot(
            after[aid][0] - before[aid][0],
            after[aid][1] - before[aid][1],
        )
        if d > max_d:
            max_d = d
    return max_d


def bond_length_stats(
    positions: dict[int, tuple[float, float]],
    bonds: list[Any],
) -> dict[str, float]:
    lengths = []
    for bond in bonds:
        a1 = _bond_a1(bond)
        a2 = _bond_a2(bond)
        if a1 not in positions or a2 not in positions:
            continue
        lengths.append(math.hypot(
            positions[a2][0] - positions[a1][0],
            positions[a2][1] - positions[a1][1],
        ))
    if not lengths:
        return {"mean": 0.0, "min": 0.0, "max": 0.0, "std": 0.0}
    mean = sum(lengths) / len(lengths)
    variance = sum((l - mean) ** 2 for l in lengths) / len(lengths)
    return {
        "mean": mean,
        "min": min(lengths),
        "max": max(lengths),
        "std": math.sqrt(variance),
    }


def _mean_bond_length(
    positions: dict[int, tuple[float, float]],
    bonds: list[Any],
) -> float:
    lengths = []
    for bond in bonds:
        a1 = _bond_a1(bond)
        a2 = _bond_a2(bond)
        if a1 in positions and a2 in positions:
            lengths.append(math.hypot(
                positions[a2][0] - positions[a1][0],
                positions[a2][1] - positions[a1][1],
            ))
    return sum(lengths) / len(lengths) if lengths else 0.0


def _compute_bounding_box_ratio(
    before: dict[int, tuple[float, float]],
    after: dict[int, tuple[float, float]],
    atom_ids: set[int],
) -> float:
    def bbox(coords: dict[int, tuple[float, float]]) -> tuple[float, float]:
        xs = [coords[a][0] for a in atom_ids if a in coords]
        ys = [coords[a][1] for a in atom_ids if a in coords]
        if not xs or not ys:
            return 1.0, 1.0
        return (max(xs) - min(xs)) or 1.0, (max(ys) - min(ys)) or 1.0

    bw, bh = bbox(before)
    aw, ah = bbox(after)
    before_diag = math.hypot(bw, bh)
    after_diag = math.hypot(aw, ah)
    if before_diag < 1e-6:
        return 1.0
    return after_diag / before_diag


def _is_nan(coord: tuple[float, float]) -> bool:
    return math.isnan(coord[0]) or math.isnan(coord[1])
