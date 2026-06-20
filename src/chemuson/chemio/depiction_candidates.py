from __future__ import annotations

from dataclasses import dataclass, field
import math

from chemuson.clean2d.geometry import count_crossings, cycle_basis, segments_intersect
from chemuson.clean2d.safety import has_cycles, min_nonbonded_distance, ring_degeneracy_score
from chemuson.core.layers import BlockKind, build_multilayer_chemical_graph
from chemuson.core.model import MolGraph, bond_affects_valence


@dataclass(frozen=True)
class DepictionCandidate:
    source: str
    graph: MolGraph
    score: float
    rejected: bool = False
    rejection_reason: str = ""
    metadata: dict[str, object] = field(default_factory=dict)


def score_imported_depiction(
    graph: MolGraph,
    *,
    target_bond_length: float = 40.0,
) -> tuple[float, dict[str, object]]:
    """Scores an imported 2D depiction; lower is better."""
    target = max(8.0, float(target_bond_length or 40.0))
    atom_ids = set(graph.atoms)
    bonds = [bond for bond in graph.bonds.values() if bond_affects_valence(bond)]
    coords = {atom_id: (float(atom.x), float(atom.y)) for atom_id, atom in graph.atoms.items()}
    metadata: dict[str, object] = {
        "atom_count": len(atom_ids),
        "bond_count": len(bonds),
    }
    score = 0.0
    if not atom_ids:
        metadata["reason"] = "missing_atoms"
        return 1_000_000.0, metadata
    if not bonds:
        metadata["reason"] = "missing_bonds"
        return 500_000.0 + len(atom_ids) * 100.0, metadata
    if any(not (math.isfinite(x) and math.isfinite(y)) for x, y in coords.values()):
        score += 300_000.0
        metadata["non_finite_coordinates"] = True

    lengths = [_distance(coords[bond.a1_id], coords[bond.a2_id]) for bond in bonds if bond.a1_id in coords and bond.a2_id in coords]
    near_zero = sum(1 for length in lengths if length <= max(1e-3, target * 0.03))
    mean_len = sum(lengths) / len(lengths) if lengths else 0.0
    variance = sum((length - mean_len) ** 2 for length in lengths) / len(lengths) if lengths else 0.0
    std_len = math.sqrt(variance)
    metadata.update({"mean_bond_length": mean_len, "bond_length_std": std_len, "near_zero_bonds": near_zero})
    if near_zero:
        score += near_zero * 20_000.0
    if mean_len <= 1e-6:
        score += 250_000.0
    else:
        score += abs(mean_len - target) / target * 40.0
        score += min(5000.0, (std_len / max(mean_len, 1e-6)) * 250.0)

    crossings = count_crossings(coords, bonds)
    metadata["bond_crossings"] = crossings
    score += crossings * 2500.0

    nonbonded = min_nonbonded_distance(coords, bonds, atom_ids)
    metadata["min_nonbonded_distance"] = nonbonded
    if nonbonded < target * 0.22:
        score += (target * 0.22 - nonbonded) / target * 5000.0
    elif nonbonded < target * 0.45:
        score += (target * 0.45 - nonbonded) / target * 500.0

    rings = cycle_basis(atom_ids, bonds, max_size=18)
    ring_scores = [ring_degeneracy_score(coords, set(ring)) for ring in rings]
    min_ring = min(ring_scores) if ring_scores else math.inf
    metadata.update({"ring_count": len(rings), "min_ring_degeneracy": min_ring})
    if min_ring < math.inf and min_ring < 0.10:
        score += (0.10 - min_ring) * 25_000.0
    elif min_ring < math.inf and min_ring < 0.25:
        score += (0.25 - min_ring) * 2500.0

    bbox = _bbox(coords)
    width, height = bbox[2] - bbox[0], bbox[3] - bbox[1]
    aspect = max(width, height) / max(1e-6, min(width, height)) if width > 0 and height > 0 else math.inf
    metadata.update({"bbox_width": width, "bbox_height": height, "aspect_ratio": aspect})
    if min(width, height) < target * 2.0:
        score += 15_000.0
    if aspect > 12.0 or aspect < 1.0:
        score += 800.0
    if len(atom_ids) >= 60 and 1.4 <= aspect <= 5.5:
        score -= 150.0

    block_score, block_metadata = _score_blocks(graph, target, coords, width, height, aspect, len(atom_ids))
    metadata.update(block_metadata)
    score += block_score
    donut_score, donut_metadata = block_donut_score(graph, target_bond_length=target)
    metadata.update(donut_metadata)
    metadata["donut_score"] = donut_score
    if int(metadata.get("block_count", 0) or 0) >= 4 or len(atom_ids) >= 35:
        score += donut_score * 900.0
    return score, metadata


def block_donut_score(graph: MolGraph, target_bond_length: float = 40.0) -> tuple[float, dict[str, object]]:
    """Scores square, radially uniform, block-rich cavity layouts."""
    target = max(8.0, float(target_bond_length or 40.0))
    coords = {atom_id: (float(atom.x), float(atom.y)) for atom_id, atom in graph.atoms.items()}
    if len(coords) < 4:
        return 0.0, {"donut_reason": "too_few_atoms"}
    try:
        layers = build_multilayer_chemical_graph(graph)
    except Exception as exc:
        return 0.0, {"donut_error": str(exc)}
    rigid_kinds = {BlockKind.AROMATIC_RING, BlockKind.FUSED_SYSTEM, BlockKind.MACROCYCLE, BlockKind.CYCLOPHANE}
    rigid_centroids: list[tuple[float, float]] = []
    counts: dict[str, int] = {}
    for block in layers.block_graph.blocks:
        counts[block.kind.value] = counts.get(block.kind.value, 0) + 1
        if block.kind in rigid_kinds:
            centroid = _centroid_for_atoms(coords, block.atom_ids)
            if centroid is not None:
                rigid_centroids.append(centroid)
    ring_count = sum(1 for motif in layers.motif_graph.motifs if getattr(motif, "kind", None).value == "ring")
    block_richness = (
        counts.get(BlockKind.AROMATIC_RING.value, 0)
        + counts.get(BlockKind.FUSED_SYSTEM.value, 0)
        + counts.get(BlockKind.MACROCYCLE.value, 0)
        + counts.get(BlockKind.CYCLOPHANE.value, 0)
        + counts.get(BlockKind.LINKER.value, 0)
    )
    metadata: dict[str, object] = {
        "donut_block_counts": counts,
        "donut_rigid_centroid_count": len(rigid_centroids),
        "donut_ring_count": ring_count,
        "donut_block_richness": block_richness,
    }
    if len(rigid_centroids) < 4 and not (len(coords) >= 35 and ring_count >= 4):
        return 0.0, metadata
    bbox = _bbox(coords)
    width, height = bbox[2] - bbox[0], bbox[3] - bbox[1]
    if width <= 1e-6 or height <= 1e-6:
        return 8.0, {**metadata, "donut_reason": "collapsed_bbox"}
    ratio = width / height
    square_factor = max(0.0, 1.0 - abs(math.log(max(ratio, 1e-6))) / math.log(1.55)) if 0.65 <= ratio <= 1.55 else 0.0
    cx = sum(x for x, _y in coords.values()) / len(coords)
    cy = sum(y for _x, y in coords.values()) / len(coords)
    nearest_atom = min(_distance((cx, cy), point) for point in coords.values())
    centroid_points = rigid_centroids or list(coords.values())
    nearest_rigid = min(_distance((cx, cy), point) for point in centroid_points)
    radial_uniformity = _radial_uniformity(centroid_points) if len(centroid_points) >= 3 else 1.0
    radial_factor = max(0.0, 1.0 - radial_uniformity / 0.38)
    hole_factor = min(2.0, max(nearest_atom, nearest_rigid) / max(target, 1e-6))
    area_ratio = _occupied_area_ratio(coords, width, height)
    sparse_factor = max(0.0, 0.70 - area_ratio)
    richness_factor = min(2.0, max(block_richness, ring_count) / 5.0)
    score = (square_factor * 2.0 + radial_factor * 2.0 + hole_factor * 1.4 + sparse_factor * 1.6) * richness_factor
    metadata.update(
        {
            "donut_bbox_width": width,
            "donut_bbox_height": height,
            "donut_bbox_ratio": ratio,
            "donut_square_factor": square_factor,
            "donut_radial_uniformity": radial_uniformity,
            "donut_nearest_atom_to_center": nearest_atom,
            "donut_nearest_rigid_to_center": nearest_rigid,
            "donut_occupied_area_ratio": area_ratio,
        }
    )
    return score, metadata


def _score_blocks(
    graph: MolGraph,
    target: float,
    coords: dict[int, tuple[float, float]],
    width: float,
    height: float,
    aspect: float,
    atom_count: int,
) -> tuple[float, dict[str, object]]:
    try:
        layers = build_multilayer_chemical_graph(graph)
    except Exception as exc:
        return 0.0, {"block_error": str(exc)}
    blocks = tuple(layers.block_graph.blocks)
    counts: dict[str, int] = {}
    centroids: list[tuple[float, float]] = []
    rigid_centroids: list[tuple[float, float]] = []
    for block in blocks:
        counts[block.kind.value] = counts.get(block.kind.value, 0) + 1
        centroid = _centroid_for_atoms(coords, block.atom_ids)
        if centroid is None:
            continue
        centroids.append(centroid)
        if block.kind in {BlockKind.AROMATIC_RING, BlockKind.FUSED_SYSTEM, BlockKind.MACROCYCLE, BlockKind.CYCLOPHANE}:
            rigid_centroids.append(centroid)
    score = 0.0
    metadata: dict[str, object] = {"block_counts": counts, "block_count": len(blocks)}
    complex_blocks = sum(counts.get(kind.value, 0) for kind in (BlockKind.AROMATIC_RING, BlockKind.FUSED_SYSTEM, BlockKind.MACROCYCLE, BlockKind.CYCLOPHANE, BlockKind.LINKER, BlockKind.STEREO_CENTER))
    if len(centroids) >= 3:
        nearest = _nearest_neighbor_distances(centroids)
        mean_nearest = sum(nearest) / len(nearest)
        metadata["mean_block_centroid_spacing"] = mean_nearest
        if mean_nearest < target * 1.15:
            score += (target * 1.15 - mean_nearest) / target * 2500.0
        else:
            score -= min(120.0, mean_nearest / target * 18.0)
    if len(rigid_centroids) >= 3:
        rigid_bbox = _bbox_from_points(rigid_centroids)
        rw, rh = rigid_bbox[2] - rigid_bbox[0], rigid_bbox[3] - rigid_bbox[1]
        rigid_aspect = max(rw, rh) / max(1e-6, min(rw, rh)) if rw > 0 and rh > 0 else math.inf
        metadata["rigid_block_aspect_ratio"] = rigid_aspect
        if rigid_aspect >= 1.6:
            score -= min(250.0, rigid_aspect * 35.0)
        if atom_count >= 60 and complex_blocks >= 6 and 0.75 <= width / max(height, 1e-6) <= 1.35:
            circularity = _radial_uniformity(rigid_centroids)
            metadata["rigid_block_radial_uniformity"] = circularity
            if circularity < 0.28:
                score += 1200.0
    if atom_count >= 60 and complex_blocks >= 6 and has_cycles(set(graph.atoms), [bond for bond in graph.bonds.values() if bond_affects_valence(bond)]):
        area_ratio = _occupied_area_ratio(coords, width, height)
        metadata["occupied_area_ratio"] = area_ratio
        if area_ratio < 0.18 and aspect < 1.6:
            score += 1800.0
    return score, metadata


def _distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(b[0] - a[0], b[1] - a[1])


def _bbox(coords: dict[int, tuple[float, float]]) -> tuple[float, float, float, float]:
    xs = [coord[0] for coord in coords.values()]
    ys = [coord[1] for coord in coords.values()]
    return min(xs), min(ys), max(xs), max(ys)


def _bbox_from_points(points: list[tuple[float, float]]) -> tuple[float, float, float, float]:
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    return min(xs), min(ys), max(xs), max(ys)


def _centroid_for_atoms(coords: dict[int, tuple[float, float]], atom_ids) -> tuple[float, float] | None:
    points = [coords[atom_id] for atom_id in atom_ids if atom_id in coords]
    if not points:
        return None
    return sum(x for x, _y in points) / len(points), sum(y for _x, y in points) / len(points)


def _nearest_neighbor_distances(points: list[tuple[float, float]]) -> list[float]:
    out: list[float] = []
    for idx, point in enumerate(points):
        distances = [_distance(point, other) for j, other in enumerate(points) if j != idx]
        if distances:
            out.append(min(distances))
    return out


def _radial_uniformity(points: list[tuple[float, float]]) -> float:
    cx = sum(x for x, _y in points) / len(points)
    cy = sum(y for _x, y in points) / len(points)
    radii = [math.hypot(x - cx, y - cy) for x, y in points]
    mean = sum(radii) / len(radii) if radii else 0.0
    if mean <= 1e-6:
        return 0.0
    variance = sum((radius - mean) ** 2 for radius in radii) / len(radii)
    return math.sqrt(variance) / mean


def _occupied_area_ratio(coords: dict[int, tuple[float, float]], width: float, height: float) -> float:
    if width <= 1e-6 or height <= 1e-6:
        return 0.0
    points = sorted(set(coords.values()))
    if len(points) < 3:
        return 0.0
    hull = _convex_hull(points)
    return _polygon_area(hull) / (width * height)


def _convex_hull(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lower: list[tuple[float, float]] = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[float, float]] = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def _polygon_area(points: list[tuple[float, float]]) -> float:
    if len(points) < 3:
        return 0.0
    area = 0.0
    for idx, (x1, y1) in enumerate(points):
        x2, y2 = points[(idx + 1) % len(points)]
        area += x1 * y2 - x2 * y1
    return abs(area) * 0.5


def _count_crossings(coords: dict[int, tuple[float, float]], bonds) -> int:
    return count_crossings(coords, list(bonds))


def _segments_intersect(p1, p2, p3, p4) -> bool:
    return segments_intersect(p1, p2, p3, p4)


def _cycle_basis(atom_ids: set[int], bonds, *, max_size: int) -> list[list[int]]:
    return cycle_basis(atom_ids, list(bonds), max_size=max_size)

