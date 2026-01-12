from __future__ import annotations

import math
from typing import Dict, Tuple

from core.model import BondStyle, BondStereo, MolGraph

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
except Exception:  # pragma: no cover - optional dependency at runtime
    Chem = None
    AllChem = None
    rdMolDraw2D = None


def _require_rdkit():
    if Chem is None or AllChem is None or rdMolDraw2D is None:
        raise RuntimeError("RDKit no disponible")


def _ensure_ring_info(mol) -> None:
    if hasattr(Chem, "FastFindRings"):
        Chem.FastFindRings(mol)
    else:
        Chem.GetSymmSSSR(mol)


def molgraph_to_rdkit_with_map(molgraph: MolGraph):
    _require_rdkit()
    rw = Chem.RWMol()
    id_map: Dict[int, int] = {}

    for atom in sorted(molgraph.atoms.values(), key=lambda a: a.id):
        rd_atom = Chem.Atom(atom.element)
        rd_atom.SetFormalCharge(atom.charge)
        if atom.isotope is not None:
            rd_atom.SetIsotope(atom.isotope)
        rd_idx = rw.AddAtom(rd_atom)
        id_map[atom.id] = rd_idx

    for bond in molgraph.bonds.values():
        if bond.is_aromatic:
            rw.GetAtomWithIdx(id_map[bond.a1_id]).SetIsAromatic(True)
            rw.GetAtomWithIdx(id_map[bond.a2_id]).SetIsAromatic(True)
            bond_type = Chem.BondType.AROMATIC
        elif bond.order == 2:
            bond_type = Chem.BondType.DOUBLE
        elif bond.order == 3:
            bond_type = Chem.BondType.TRIPLE
        else:
            bond_type = Chem.BondType.SINGLE
        
        # Robust check: don't add bond if it already exists in RDKit mol
        if rw.GetBondBetweenAtoms(id_map[bond.a1_id], id_map[bond.a2_id]) is None:
            rw.AddBond(id_map[bond.a1_id], id_map[bond.a2_id], bond_type)

    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for atom_id, idx in id_map.items():
        atom = molgraph.atoms[atom_id]
        conf.SetAtomPosition(idx, (atom.x, atom.y, 0.0))
    mol.AddConformer(conf, assignId=True)
    return mol, id_map


def molgraph_to_rdkit(molgraph: MolGraph):
    mol, _ = molgraph_to_rdkit_with_map(molgraph)
    return mol


def molgraph_to_smiles(molgraph: MolGraph) -> str:
    mol = molgraph_to_rdkit(molgraph)
    return Chem.MolToSmiles(mol, canonical=True)


def molgraph_to_molfile(molgraph: MolGraph) -> str:
    mol = molgraph_to_rdkit(molgraph)
    return Chem.MolToMolBlock(mol)


def molgraph_to_svg(molgraph: MolGraph, size: Tuple[int, int] = (300, 200)) -> str:
    _require_rdkit()
    mol = molgraph_to_rdkit(molgraph)
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def molfile_to_molgraph(molfile: str) -> MolGraph:
    _require_rdkit()
    mol = Chem.MolFromMolBlock(molfile, sanitize=True)
    return rdkit_to_molgraph(mol)


def smiles_to_molgraph(smiles: str) -> MolGraph:
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    return rdkit_to_molgraph(mol)


def rdkit_to_molgraph(mol) -> MolGraph:
    _require_rdkit()
    if mol is None:
        raise ValueError("Mol inválido")
    if mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()

    graph = MolGraph()
    idx_map: Dict[int, int] = {}

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        pos = conf.GetAtomPosition(idx)
        new_atom = graph.add_atom(atom.GetSymbol(), pos.x, pos.y)
        new_atom.charge = atom.GetFormalCharge()
        if atom.GetIsotope():
            new_atom.isotope = atom.GetIsotope()
        idx_map[idx] = new_atom.id

    for bond in mol.GetBonds():
        order = 1
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            order = 2
        elif bond.GetBondType() == Chem.BondType.TRIPLE:
            order = 3
        graph.add_bond(
            idx_map[bond.GetBeginAtomIdx()],
            idx_map[bond.GetEndAtomIdx()],
            order,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=bond.GetIsAromatic(),
        )

    _scale_to_default(graph)
    return graph


def _scale_to_default(graph: MolGraph, target: float = 40.0) -> None:
    if not graph.bonds:
        return
    lengths = []
    for bond in graph.bonds.values():
        a1 = graph.atoms[bond.a1_id]
        a2 = graph.atoms[bond.a2_id]
        dx = a2.x - a1.x
        dy = a2.y - a1.y
        lengths.append((dx * dx + dy * dy) ** 0.5)
    avg = sum(lengths) / len(lengths)
    if avg <= 0:
        return
    scale = target / avg
    for atom in graph.atoms.values():
        atom.x *= scale
        atom.y *= scale


def _fallback_aromatic_orders(
    molgraph: MolGraph,
    ring_counts: Dict[int, int] | None = None,
) -> Dict[int, int]:
    aromatic_bonds = [bond for bond in molgraph.bonds.values() if bond.is_aromatic]
    if not aromatic_bonds:
        return {}
    if ring_counts is None:
        ring_counts = {}

    adjacency: dict[int, list[int]] = {}
    bond_lookup: dict[frozenset[int], int] = {}
    for bond in aromatic_bonds:
        adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
        adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)
        bond_lookup[frozenset((bond.a1_id, bond.a2_id))] = bond.id

    color: dict[int, int] = {}
    for node in adjacency:
        if node in color:
            continue
        color[node] = 0
        stack = [node]
        while stack:
            current = stack.pop()
            for neighbor in adjacency.get(current, []):
                if neighbor in color:
                    if color[neighbor] == color[current]:
                        return {}
                else:
                    color[neighbor] = 1 - color[current]
                    stack.append(neighbor)

    u_nodes = [n for n, c in color.items() if c == 0]
    v_nodes = [n for n, c in color.items() if c == 1]
    u_index = {u: i for i, u in enumerate(u_nodes)}
    v_index = {v: i for i, v in enumerate(v_nodes)}

    degrees = {node: len(adjacency.get(node, [])) for node in adjacency}

    def edge_weight(u: int, v: int) -> int:
        du = degrees.get(u, 0)
        dv = degrees.get(v, 0)
        bond_id = bond_lookup.get(frozenset((u, v)))
        ring_count = ring_counts.get(bond_id, 1)
        outer_bonus = 4 if ring_count <= 1 else 0
        if du == 2 and dv == 2:
            return 3 + outer_bonus
        if du == 2 or dv == 2:
            return 2 + outer_bonus
        return 1 + outer_bonus

    class Edge:
        def __init__(self, to: int, rev: int, cap: int, cost: int) -> None:
            self.to = to
            self.rev = rev
            self.cap = cap
            self.cost = cost

    size = 2 + len(u_nodes) + len(v_nodes)
    src = 0
    sink = size - 1
    graph: list[list[Edge]] = [[] for _ in range(size)]

    def add_edge(frm: int, to: int, cap: int, cost: int) -> None:
        graph[frm].append(Edge(to, len(graph[to]), cap, cost))
        graph[to].append(Edge(frm, len(graph[frm]) - 1, 0, -cost))

    for u in u_nodes:
        add_edge(src, 1 + u_index[u], 1, 0)
    for v in v_nodes:
        add_edge(1 + len(u_nodes) + v_index[v], sink, 1, 0)

    for u in u_nodes:
        for v in adjacency.get(u, []):
            if v not in v_index:
                continue
            w = edge_weight(u, v)
            add_edge(1 + u_index[u], 1 + len(u_nodes) + v_index[v], 1, -w)

    while True:
        dist = [10 ** 9] * size
        in_queue = [False] * size
        prev_node = [-1] * size
        prev_edge = [-1] * size
        dist[src] = 0
        queue = [src]
        in_queue[src] = True
        while queue:
            node = queue.pop(0)
            in_queue[node] = False
            for i, edge in enumerate(graph[node]):
                if edge.cap <= 0:
                    continue
                nd = dist[node] + edge.cost
                if nd < dist[edge.to]:
                    dist[edge.to] = nd
                    prev_node[edge.to] = node
                    prev_edge[edge.to] = i
                    if not in_queue[edge.to]:
                        in_queue[edge.to] = True
                        queue.append(edge.to)
        if dist[sink] == 10 ** 9:
            break
        v = sink
        while v != src:
            edge = graph[prev_node[v]][prev_edge[v]]
            edge.cap -= 1
            graph[v][edge.rev].cap += 1
            v = prev_node[v]

    pair_u: dict[int, int | None] = {u: None for u in u_nodes}
    for u in u_nodes:
        node_u = 1 + u_index[u]
        for edge in graph[node_u]:
            if edge.to == src or edge.cap != 0:
                continue
            if edge.to == sink:
                continue
            if edge.to >= 1 + len(u_nodes) and edge.to < sink:
                v = v_nodes[edge.to - 1 - len(u_nodes)]
                pair_u[u] = v
                break

    orders: Dict[int, int] = {}
    for bond in aromatic_bonds:
        double = False
        if color.get(bond.a1_id) == 0:
            double = pair_u.get(bond.a1_id) == bond.a2_id
        else:
            double = pair_u.get(bond.a2_id) == bond.a1_id
        orders[bond.id] = 2 if double else 1
    return orders


def get_visual_properties(
    molgraph: MolGraph,
) -> Tuple[Dict[int, int], Dict[int, Tuple[float, float]]]:
    """
    Returns visual properties derived from RDKit:
    1. Visual bond orders (canonical Kekule form)
    2. Ring centers for bonds in rings (to resolve double bond direction)
    """
    _require_rdkit()
    try:
        mol, id_map = molgraph_to_rdkit_with_map(molgraph)
        rd_to_graph_id = {rd_idx: g_id for g_id, rd_idx in id_map.items()}
        bond_lookup = {
            frozenset((bond.a1_id, bond.a2_id)): bond.id
            for bond in molgraph.bonds.values()
        }

        _ensure_ring_info(mol)
        ring_info = mol.GetRingInfo()
        ring_bonds = {idx for ring in ring_info.BondRings() for idx in ring}
        ring_atoms = {idx for ring in ring_info.AtomRings() for idx in ring}
        ring_counts: Dict[int, int] = {}
        for bond_ring in ring_info.BondRings():
            for bond_idx in bond_ring:
                rd_bond = mol.GetBondWithIdx(bond_idx)
                a1_id = rd_to_graph_id[rd_bond.GetBeginAtomIdx()]
                a2_id = rd_to_graph_id[rd_bond.GetEndAtomIdx()]
                bond_id = bond_lookup.get(frozenset((a1_id, a2_id)))
                if bond_id is None:
                    continue
                ring_counts[bond_id] = ring_counts.get(bond_id, 0) + 1
        for bond in mol.GetBonds():
            if bond.GetIsAromatic() and bond.GetIdx() not in ring_bonds:
                bond.SetIsAromatic(False)
                bond.SetBondType(Chem.BondType.SINGLE)
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic() and atom.GetIdx() not in ring_atoms:
                atom.SetIsAromatic(False)
        for bond in mol.GetBonds():
            if bond.GetIsAromatic():
                bond.GetBeginAtom().SetIsAromatic(True)
                bond.GetEndAtom().SetIsAromatic(True)

        kekule_mol = Chem.Mol(mol)
        kekulized = False
        try:
            try:
                from rdkit import RDLogger
            except Exception:
                RDLogger = None
            if RDLogger is not None:
                RDLogger.DisableLog("rdApp.warning")
            try:
                Chem.Kekulize(kekule_mol, clearAromaticFlags=True)
                kekulized = True
            finally:
                if RDLogger is not None:
                    RDLogger.EnableLog("rdApp.warning")
        except Exception:
            kekulized = False

        bond_orders: Dict[int, int] = {}
        if kekulized:
            for bond in kekule_mol.GetBonds():
                a1_id = rd_to_graph_id[bond.GetBeginAtomIdx()]
                a2_id = rd_to_graph_id[bond.GetEndAtomIdx()]
                bond_id = bond_lookup.get(frozenset((a1_id, a2_id)))
                if bond_id is None:
                    continue
                order = 1
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    order = 2
                elif bond.GetBondType() == Chem.BondType.TRIPLE:
                    order = 3
                bond_orders[bond_id] = order
            for bond in molgraph.bonds.values():
                bond_orders.setdefault(bond.id, bond.order)
        else:
            bond_orders = {bond.id: bond.order for bond in molgraph.bonds.values()}
            bond_orders.update(_fallback_aromatic_orders(molgraph, ring_counts))

        bond_ring_centers: Dict[int, Tuple[float, float]] = {}
        atom_rings = ring_info.AtomRings()
        bond_rings = ring_info.BondRings()
        ring_candidates: dict[int, list[tuple[int, float, float, float]]] = {}

        for ring_index, bond_ring in enumerate(bond_rings):
            atom_ring = atom_rings[ring_index] if ring_index < len(atom_rings) else ()
            if not atom_ring:
                continue
            xs = []
            ys = []
            for rd_idx in atom_ring:
                atom_id = rd_to_graph_id[rd_idx]
                atom = molgraph.atoms[atom_id]
                xs.append(atom.x)
                ys.append(atom.y)
            cx = sum(xs) / len(xs)
            cy = sum(ys) / len(ys)
            ring_size = len(atom_ring)
            for bond_idx in bond_ring:
                rd_bond = mol.GetBondWithIdx(bond_idx)
                a1_id = rd_to_graph_id[rd_bond.GetBeginAtomIdx()]
                a2_id = rd_to_graph_id[rd_bond.GetEndAtomIdx()]
                bond_id = bond_lookup.get(frozenset((a1_id, a2_id)))
                if bond_id is None:
                    continue
                bond = molgraph.bonds[bond_id]
                midx = (molgraph.atoms[bond.a1_id].x + molgraph.atoms[bond.a2_id].x) / 2
                midy = (molgraph.atoms[bond.a1_id].y + molgraph.atoms[bond.a2_id].y) / 2
                dist = math.hypot(cx - midx, cy - midy)
                ring_candidates.setdefault(bond_id, []).append((ring_size, dist, cx, cy))

        for bond_id, candidates in ring_candidates.items():
            _, _, cx, cy = min(candidates, key=lambda item: (item[0], item[1]))
            bond_ring_centers[bond_id] = (cx, cy)

        return bond_orders, bond_ring_centers

    except Exception:
        return {}, {}


def get_aromatic_ring_geometries(molgraph: MolGraph) -> list[tuple[float, float, float]]:
    """Return ring centers/radii for aromatic rings (for circle rendering)."""
    _require_rdkit()
    try:
        mol, id_map = molgraph_to_rdkit_with_map(molgraph)
        rd_to_graph_id = {rd_idx: g_id for g_id, rd_idx in id_map.items()}
        bond_lookup = {
            frozenset((bond.a1_id, bond.a2_id)): bond
            for bond in molgraph.bonds.values()
        }

        _ensure_ring_info(mol)
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        bond_rings = ring_info.BondRings()
        geometries: list[tuple[float, float, float]] = []

        for ring_index, bond_ring in enumerate(bond_rings):
            atom_ring = atom_rings[ring_index] if ring_index < len(atom_rings) else ()
            if not atom_ring:
                continue
            is_aromatic = True
            for bond_idx in bond_ring:
                rd_bond = mol.GetBondWithIdx(bond_idx)
                a1_id = rd_to_graph_id[rd_bond.GetBeginAtomIdx()]
                a2_id = rd_to_graph_id[rd_bond.GetEndAtomIdx()]
                bond = bond_lookup.get(frozenset((a1_id, a2_id)))
                if bond is None or not bond.is_aromatic:
                    is_aromatic = False
                    break
            if not is_aromatic:
                continue
            xs = []
            ys = []
            for rd_idx in atom_ring:
                atom_id = rd_to_graph_id[rd_idx]
                atom = molgraph.atoms[atom_id]
                xs.append(atom.x)
                ys.append(atom.y)
            cx = sum(xs) / len(xs)
            cy = sum(ys) / len(ys)
            avg_dist = sum(
                math.hypot(x - cx, y - cy) for x, y in zip(xs, ys)
            ) / len(xs)
            radius = avg_dist * 0.55
            geometries.append((cx, cy, radius))

        return geometries
    except Exception:
        return []
