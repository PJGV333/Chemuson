from __future__ import annotations

from .errors import ChemNameInternalError, ChemNameNotSupported
from chemcalc.valence import implicit_h_count

from .locants import Sub, orientation_key, substituents_on_chain
from .molview import MolView
from .options import NameOptions
from .parent_chain import longest_carbon_chain, longest_chain_in_subset
from .render import render_name
from .ring_naming import (
    choose_hetero_ring_orientation,
    choose_ring_orientation,
    ring_substituents,
)
from .rings import (
    RingContext,
    classify_aromatic_ring,
    detect_naphthalene,
    build_ring_context,
    find_rings_simple,
    is_simple_ring,
    ring_order,
)
from .substituents import CYCLO_PARENT, HALO_MAP, alkyl_substituent_name, parent_name


def iupac_name(graph, opts: NameOptions = NameOptions()) -> str:
    """Public entry point for the lite naming engine."""
    try:
        return iupac_name_lite(graph, opts)
    except ChemNameNotSupported:
        if opts.return_nd_on_fail:
            return "N/D"
        raise
    except Exception as exc:  # pragma: no cover - defensive
        if opts.return_nd_on_fail:
            return "N/D"
        raise ChemNameInternalError(str(exc)) from exc


def iupac_name_lite(graph, opts: NameOptions) -> str:
    """Lite naming engine for acyclic hydrocarbons with simple substituents."""
    view = MolView(graph)
    ring_ctx = build_ring_context(view)
    rings = ring_ctx.rings
    if not rings:
        return _name_linear(view, opts)

    ring_atoms = set().union(*rings)
    chain = _longest_chain_excluding(view, ring_atoms)
    if chain:
        func = _find_functional_group(view, chain, allow_other_hetero=True)
        if func is not None:
            return _name_linear(
                view,
                opts,
                chain_override=chain,
                ring_ctx=ring_ctx,
                allow_rings=True,
            )

    # fused systems first
    try:
        return _name_cyclic(view, rings, opts, ring_ctx)
    except ChemNameNotSupported:
        pass

    chain_len = len(chain)
    benzene_rings = [r for r in rings if ring_ctx.ring_types.get(r) == "benzene"]
    cyclo_rings = [r for r in rings if ring_ctx.ring_types.get(r) == "cyclohexane"]

    if benzene_rings and chain_len <= 6:
        return _name_cyclic(view, benzene_rings, opts, ring_ctx)

    if cyclo_rings:
        ring_size = max(len(r) for r in cyclo_rings)
        if chain_len > ring_size:
            return _name_linear(
                view,
                opts,
                chain_override=chain,
                ring_ctx=ring_ctx,
                allow_rings=True,
            )
        return _name_cyclic(view, cyclo_rings, opts, ring_ctx)

    # default: chain if available, else fall back to cyclic
    if chain:
        return _name_linear(
            view,
            opts,
            chain_override=chain,
            ring_ctx=ring_ctx,
            allow_rings=True,
        )
    return _name_cyclic(view, rings, opts, ring_ctx)


def _name_linear(
    view: MolView,
    opts: NameOptions,
    chain_override: list[int] | None = None,
    ring_ctx: "RingContext | None" = None,
    allow_rings: bool = False,
) -> str:
    if not allow_rings and not view.is_acyclic():
        raise ChemNameNotSupported("Cyclic structures not supported")

    allowed_elements = {"C", "H", "F", "Cl", "Br", "I", "O", "N"}
    for atom_id in view.atoms():
        elem = view.element(atom_id)
        if elem not in allowed_elements:
            raise ChemNameNotSupported("Unsupported element")

    chain = chain_override or longest_carbon_chain(view)
    if not chain:
        raise ChemNameNotSupported("No carbon chain found")

    unsat_order, _ = _find_unsaturation(view, chain, strict_outside_chain=not allow_rings)
    func = _find_functional_group(view, chain)
    func_atom = func[1] if func else None
    func_suffix = func[0] if func else None
    ignore_atoms = {func[2]} if func else set()

    if allow_rings and ring_ctx is not None and func is None and unsat_order is None:
        chain = _trim_aromatic_linker(view, chain, ring_ctx)

    chain = _choose_oriented_chain(
        view,
        chain,
        opts,
        func_atom=func_atom,
        unsat_order=unsat_order,
        ignore_atoms=ignore_atoms,
        ring_ctx=ring_ctx,
    )

    substituents = substituents_on_chain(
        view, chain, ignore_atoms=ignore_atoms, ring_ctx=ring_ctx
    )
    func_locant = _locant_for_atom(chain, func_atom) if func_atom is not None else None
    unsat_locant = _unsaturation_locant_for_chain(view, chain, unsat_order)

    parent = parent_name(
        len(chain),
        unsat_order=unsat_order,
        unsat_locant=unsat_locant,
        suffix=func_suffix,
        suffix_locant=func_locant,
    )

    return render_name(substituents, parent)


def _longest_chain_excluding(view: MolView, excluded_atoms: set[int]) -> list[int]:
    allowed = {atom_id for atom_id in view.atoms() if atom_id not in excluded_atoms}
    return longest_chain_in_subset(view, allowed)


def _name_cyclic(
    view: MolView,
    rings: list[frozenset[int]],
    opts: NameOptions,
    ring_ctx: RingContext | None = None,
) -> str:
    if len(rings) == 1:
        ring_nodes = rings[0]
        if not is_simple_ring(view, ring_nodes):
            raise ChemNameNotSupported("Fused rings not supported")

        ring_atoms = ring_order(view, ring_nodes)
        if not ring_atoms:
            raise ChemNameNotSupported("Ring ordering failed")

        aromatic_info = classify_aromatic_ring(view, ring_nodes)
        if aromatic_info and aromatic_info.get("kind") == "benzene":
            return _name_benzene(view, ring_atoms, opts, ring_ctx)
        if aromatic_info and aromatic_info.get("kind") is not None:
            return _name_heteroaromatic(view, aromatic_info, opts, ring_ctx)

        return _name_cycloalkane(view, ring_atoms, opts, ring_ctx)

    fused_hetero = _detect_fused_hetero(view, rings)
    if fused_hetero is not None:
        return _name_fused_hetero(view, fused_hetero, opts)

    tricyclic = _detect_tricyclic(view, rings)
    if tricyclic is not None:
        return _name_tricyclic(view, tricyclic, opts)

    naph = detect_naphthalene(view, rings)
    if naph is not None:
        return _name_naphthalene(view, naph, opts)

    raise ChemNameNotSupported("Multiple rings not supported")


def _name_cycloalkane(
    view: MolView,
    ring_atoms: list[int],
    opts: NameOptions,
    ring_ctx: RingContext | None = None,
) -> str:
    for idx in range(len(ring_atoms)):
        if view.bond_order_between(ring_atoms[idx], ring_atoms[(idx + 1) % len(ring_atoms)]) != 1:
            raise ChemNameNotSupported("Unsaturated ring not supported")
    for atom_id in ring_atoms:
        if view.element(atom_id) != "C":
            raise ChemNameNotSupported("Non-carbon ring not supported")

    parent = CYCLO_PARENT.get(len(ring_atoms))
    if parent is None:
        raise ChemNameNotSupported("Unsupported ring size")

    oriented = choose_ring_orientation(view, ring_atoms, opts, ring_ctx=ring_ctx)
    substituents = ring_substituents(view, oriented, ring_ctx=ring_ctx)
    return render_name(substituents, parent, always_include_locant=False)


def _name_benzene(
    view: MolView,
    ring_atoms: list[int],
    opts: NameOptions,
    ring_ctx: RingContext | None = None,
) -> str:
    for atom_id in ring_atoms:
        if view.element(atom_id) != "C":
            raise ChemNameNotSupported("Unsupported aromatic ring")
    oriented = choose_ring_orientation(
        view,
        ring_atoms,
        opts,
        allow_hydroxy=True,
        allow_nitro=True,
        ring_ctx=ring_ctx,
    )
    substituents = ring_substituents(
        view,
        oriented,
        allow_hydroxy=True,
        allow_nitro=True,
        ring_ctx=ring_ctx,
    )
    return render_name(substituents, "benzene", always_include_locant=False)


def _name_heteroaromatic(
    view: MolView,
    info: dict,
    opts: NameOptions,
    ring_ctx: RingContext | None = None,
) -> str:
    kind = info.get("kind")
    ring_atoms = info.get("order") or []
    hetero_atoms = info.get("hetero_atoms") or []
    if not kind or not ring_atoms or not hetero_atoms:
        raise ChemNameNotSupported("Unsupported heteroaromatic ring")

    oriented = choose_hetero_ring_orientation(
        view,
        ring_atoms,
        hetero_atoms,
        opts,
        allow_hydroxy=True,
        allow_nitro=False,
        forbid_hetero_substituents=True,
        ring_ctx=ring_ctx,
    )
    substituents = ring_substituents(
        view,
        oriented,
        allow_hydroxy=True,
        allow_nitro=False,
        forbid_hetero_substituents=True,
        ring_ctx=ring_ctx,
    )
    return render_name(substituents, kind, always_include_locant=False)


def _name_naphthalene(view: MolView, info: dict, opts: NameOptions) -> str:
    ring_atoms = set(info.get("atoms", set()))
    fusion_atoms = set(info.get("fusion_atoms", ()))
    if len(ring_atoms) != 10 or len(fusion_atoms) != 2:
        raise ChemNameNotSupported("Unsupported fused ring")

    order = _perimeter_cycle(view, ring_atoms, fusion_atoms)
    if not order:
        raise ChemNameNotSupported("Naphthalene ordering failed")

    substituents: list[tuple[int, str]] = []
    for atom_id in ring_atoms:
        for nbr in view.neighbors(atom_id):
            if nbr in ring_atoms:
                continue
            if view.element(nbr) == "H":
                continue
            if atom_id in fusion_atoms:
                raise ChemNameNotSupported("Fusion-atom substitution not supported")
            name = _simple_substituent_name(view, atom_id, nbr, ring_atoms)
            substituents.append((atom_id, name))

    if len(substituents) > 2:
        raise ChemNameNotSupported("Multiple substitutions not supported")

    best_subs: list[Sub] = []
    best_key = None
    for locant_map in _naphthalene_locant_maps(order, fusion_atoms):
        subs = [Sub(name, locant_map[atom_id]) for atom_id, name in substituents]
        locants = sorted(sub.locant for sub in subs)
        key = orientation_key(subs, opts, primary_locants=locants)
        if best_key is None or key < best_key:
            best_key = key
            best_subs = subs

    return render_name(best_subs, "naphthalene", always_include_locant=True)


def _detect_tricyclic(view: MolView, rings: list[frozenset[int]]) -> dict | None:
    benzene_rings: list[frozenset[int]] = []
    ring_edges: dict[frozenset[int], set[frozenset[int]]] = {}
    for ring in rings:
        info = classify_aromatic_ring(view, ring)
        if not info or info.get("kind") != "benzene":
            continue
        benzene_rings.append(ring)
        edges = set()
        order = ring_order(view, ring)
        if not order:
            continue
        for idx in range(len(order)):
            a = order[idx]
            b = order[(idx + 1) % len(order)]
            edges.add(frozenset({a, b}))
        ring_edges[ring] = edges

    if len(benzene_rings) < 3:
        return None

    # build ring adjacency by shared edges
    adjacency: dict[frozenset[int], list[frozenset[int]]] = {r: [] for r in benzene_rings}
    shared_edges: dict[tuple[frozenset[int], frozenset[int]], set[frozenset[int]]] = {}
    for i in range(len(benzene_rings)):
        for j in range(i + 1, len(benzene_rings)):
            r1 = benzene_rings[i]
            r2 = benzene_rings[j]
            shared = ring_edges[r1] & ring_edges[r2]
            if shared:
                adjacency[r1].append(r2)
                adjacency[r2].append(r1)
                shared_edges[(r1, r2)] = shared

    # find chain of three rings
    for mid, neighbors in adjacency.items():
        if len(neighbors) != 2:
            continue
        r1, r2 = neighbors
        if r2 in adjacency.get(r1, []):
            # outer rings should not be directly connected
            continue
        union = set().union(mid, r1, r2)
        if len(union) != 14:
            continue
        if any(view.element(atom_id) != "C" for atom_id in union):
            continue
        fusion_edges = set()
        for pair, edges in shared_edges.items():
            if mid in pair and (r1 in pair or r2 in pair):
                fusion_edges |= edges
        if len(fusion_edges) != 2:
            continue
        triple_intersection = set(mid) & set(r1) & set(r2)
        kind = "phenanthrene" if triple_intersection else "anthracene"
        return {
            "kind": kind,
            "atoms": union,
            "fusion_edges": fusion_edges,
        }
    return None


def _name_tricyclic(view: MolView, info: dict, opts: NameOptions) -> str:
    ring_atoms = set(info.get("atoms", set()))
    fusion_edges = info.get("fusion_edges", set())
    kind = info.get("kind")
    if kind not in {"anthracene", "phenanthrene"}:
        raise ChemNameNotSupported("Unsupported fused ring system")
    if len(ring_atoms) != 14 or len(fusion_edges) != 2:
        raise ChemNameNotSupported("Unsupported fused ring system")

    fusion_atoms = set()
    for edge in fusion_edges:
        fusion_atoms |= set(edge)

    order = _perimeter_cycle_multi(view, ring_atoms, fusion_edges)
    if not order:
        raise ChemNameNotSupported("Fused ring ordering failed")

    substituents: list[tuple[int, str]] = []
    for atom_id in ring_atoms:
        for nbr in view.neighbors(atom_id):
            if nbr in ring_atoms:
                continue
            if view.element(nbr) == "H":
                continue
            if atom_id in fusion_atoms:
                raise ChemNameNotSupported("Fusion-atom substitution not supported")
            name = _simple_substituent_name(view, atom_id, nbr, ring_atoms)
            substituents.append((atom_id, name))
    if len(substituents) > 1:
        raise ChemNameNotSupported("Multiple substitutions not supported")

    best_subs: list[Sub] = []
    best_key = None
    for locant_map in _fused_locant_maps(order, fusion_atoms):
        subs = [Sub(name, locant_map[atom_id]) for atom_id, name in substituents]
        locants = sorted(sub.locant for sub in subs)
        key = orientation_key(subs, opts, primary_locants=locants)
        if best_key is None or key < best_key:
            best_key = key
            best_subs = subs

    return render_name(best_subs, kind, always_include_locant=True)


def _detect_fused_hetero(view: MolView, rings: list[frozenset[int]]) -> dict | None:
    ring_infos: dict[frozenset[int], dict] = {}
    for ring in rings:
        info = classify_aromatic_ring(view, ring)
        if info:
            ring_infos[ring] = info

    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            r1 = rings[i]
            r2 = rings[j]
            info1 = ring_infos.get(r1)
            info2 = ring_infos.get(r2)
            if not info1 or not info2:
                continue
            shared = r1 & r2
            if len(shared) != 2:
                continue
            a, b = tuple(shared)
            if view.bond_order_between(a, b) <= 0:
                continue
            kinds = {info1.get("kind"), info2.get("kind")}
            union = set(r1) | set(r2)
            if kinds == {"benzene", "pyridine"} and len(union) == 10:
                pyridine_ring = r1 if info1.get("kind") == "pyridine" else r2
                hetero_atom = next(
                    atom_id for atom_id in pyridine_ring if view.element(atom_id) == "N"
                )
                fusion_atoms = set(shared)
                alpha_atoms = set()
                for fusion in fusion_atoms:
                    for nbr in view.neighbors(fusion):
                        if nbr in union and nbr not in fusion_atoms:
                            alpha_atoms.add(nbr)
                kind = "quinoline" if hetero_atom in alpha_atoms else "isoquinoline"
                return {
                    "kind": kind,
                    "atoms": union,
                    "fusion_atoms": fusion_atoms,
                    "hetero_atom": hetero_atom,
                }
            if kinds == {"benzene", "pyrrole"} and len(union) == 9:
                pyrrole_ring = r1 if info1.get("kind") == "pyrrole" else r2
                hetero_atom = next(
                    atom_id for atom_id in pyrrole_ring if view.element(atom_id) == "N"
                )
                return {
                    "kind": "indole",
                    "atoms": union,
                    "fusion_atoms": set(shared),
                    "hetero_atom": hetero_atom,
                }
    return None


def _name_fused_hetero(view: MolView, info: dict, opts: NameOptions) -> str:
    kind = info.get("kind")
    ring_atoms = set(info.get("atoms", set()))
    fusion_atoms = set(info.get("fusion_atoms", set()))
    hetero_atom = info.get("hetero_atom")
    if not kind or hetero_atom is None:
        raise ChemNameNotSupported("Unsupported fused heteroaromatic")

    order = _perimeter_cycle(view, ring_atoms, fusion_atoms)
    if not order:
        raise ChemNameNotSupported("Fused hetero ordering failed")

    if hetero_atom in fusion_atoms:
        raise ChemNameNotSupported("Unsupported fused heteroaromatic")

    substituents: list[tuple[int, str]] = []
    for atom_id in ring_atoms:
        for nbr in view.neighbors(atom_id):
            if nbr in ring_atoms:
                continue
            if view.element(nbr) == "H":
                continue
            if atom_id in fusion_atoms or atom_id == hetero_atom:
                raise ChemNameNotSupported("Unsupported hetero substitution")
            name = _simple_substituent_name(view, atom_id, nbr, ring_atoms)
            substituents.append((atom_id, name))
    if len(substituents) > 1:
        raise ChemNameNotSupported("Multiple substitutions not supported")

    best_subs: list[Sub] = []
    best_key = None
    for locant_map in _hetero_fused_locant_maps(order, fusion_atoms, hetero_atom):
        subs = [Sub(name, locant_map[atom_id]) for atom_id, name in substituents]
        locants = sorted(sub.locant for sub in subs)
        key = orientation_key(subs, opts, primary_locants=locants)
        if best_key is None or key < best_key:
            best_key = key
            best_subs = subs

    return render_name(best_subs, kind, always_include_locant=True)


def _hetero_fused_locant_maps(
    order: list[int], fusion_atoms: set[int], hetero_atom: int
) -> list[dict[int, int]]:
    if hetero_atom not in order:
        return []
    idx0 = order.index(hetero_atom)
    base = order[idx0:] + order[:idx0]
    maps: list[dict[int, int]] = []
    for direction in (1, -1):
        if direction == 1:
            seq = list(base)
        else:
            rev = list(reversed(base))
            ridx = rev.index(hetero_atom)
            seq = rev[ridx:] + rev[:ridx]
        locant_map: dict[int, int] = {}
        loc = 1
        for atom_id in seq:
            if atom_id in fusion_atoms:
                continue
            locant_map[atom_id] = loc
            loc += 1
        maps.append(locant_map)
    return maps


def _perimeter_cycle_multi(
    view: MolView, ring_atoms: set[int], fusion_edges: set[frozenset[int]]
) -> list[int]:
    adjacency = {
        atom_id: sorted(nbr for nbr in view.neighbors(atom_id) if nbr in ring_atoms)
        for atom_id in ring_atoms
    }
    for edge in fusion_edges:
        a, b = tuple(edge)
        if b in adjacency.get(a, []):
            adjacency[a] = [nbr for nbr in adjacency[a] if nbr != b]
        if a in adjacency.get(b, []):
            adjacency[b] = [nbr for nbr in adjacency[b] if nbr != a]
    if any(len(nbrs) != 2 for nbrs in adjacency.values()):
        return []
    start = min(ring_atoms)
    return _cycle_order(adjacency, start)


def _fused_locant_maps(order: list[int], fusion_atoms: set[int]) -> list[dict[int, int]]:
    if not order:
        return []
    start_fusion = min(fusion_atoms) if fusion_atoms else order[0]
    if start_fusion not in order:
        start_fusion = order[0]
    idx0 = order.index(start_fusion)
    ordered = order[idx0:] + order[:idx0]
    maps: list[dict[int, int]] = []
    for direction in (1, -1):
        seq = ordered if direction == 1 else list(reversed(ordered))
        locant_map: dict[int, int] = {}
        loc = 1
        for atom_id in seq:
            if atom_id in fusion_atoms:
                continue
            locant_map[atom_id] = loc
            loc += 1
        maps.append(locant_map)
    return maps


def _perimeter_cycle(
    view: MolView, ring_atoms: set[int], fusion_atoms: set[int]
) -> list[int]:
    adjacency = {
        atom_id: sorted(nbr for nbr in view.neighbors(atom_id) if nbr in ring_atoms)
        for atom_id in ring_atoms
    }
    fusion = tuple(fusion_atoms)
    if len(fusion) != 2:
        return []
    a, b = fusion
    if b in adjacency.get(a, []):
        adjacency[a] = [nbr for nbr in adjacency[a] if nbr != b]
    if a in adjacency.get(b, []):
        adjacency[b] = [nbr for nbr in adjacency[b] if nbr != a]
    if any(len(nbrs) != 2 for nbrs in adjacency.values()):
        return []
    start = min(ring_atoms)
    return _cycle_order(adjacency, start)


def _cycle_order(adjacency: dict[int, list[int]], start: int) -> list[int]:
    order = [start]
    if not adjacency.get(start):
        return []
    next_node = adjacency[start][0]
    prev = start
    current = next_node
    while True:
        if current == start:
            break
        order.append(current)
        nbrs = adjacency[current]
        next_candidate = nbrs[0] if nbrs[1] == prev else nbrs[1]
        prev, current = current, next_candidate
        if len(order) > len(adjacency):
            return []
    if len(order) != len(adjacency):
        return []
    return order


def _naphthalene_locant_maps(
    order: list[int], fusion_atoms: set[int]
) -> list[dict[int, int]]:
    fusion_list = list(fusion_atoms)
    if len(fusion_list) != 2:
        return []
    maps: list[dict[int, int]] = []
    for start_fusion in fusion_list:
        idx0 = order.index(start_fusion)
        ordered = order[idx0:] + order[:idx0]
        for direction in (1, -1):
            seq = ordered if direction == 1 else list(reversed(ordered))
            try:
                idx_other = seq.index(fusion_list[1] if start_fusion == fusion_list[0] else fusion_list[0])
            except ValueError:
                continue
            if idx_other != 5:
                continue
            locant_map: dict[int, int] = {}
            loc = 1
            for i, atom_id in enumerate(seq):
                if atom_id in fusion_atoms:
                    continue
                locant_map[atom_id] = loc
                loc += 1
            if len(locant_map) == 8:
                maps.append(locant_map)
    return maps


def _simple_substituent_name(view: MolView, ring_atom: int, nbr: int, ring_set: set[int]) -> str:
    elem = view.element(nbr)
    if elem in HALO_MAP:
        return HALO_MAP[elem]
    if elem == "C":
        return alkyl_substituent_name(view, nbr, ring_set)
    raise ChemNameNotSupported("Unsupported substituent")


def _find_unsaturation(
    view: MolView, chain: list[int], strict_outside_chain: bool = True
) -> tuple[int | None, int | None]:
    chain_set = set(chain)
    index = {atom_id: idx for idx, atom_id in enumerate(chain)}
    unsat_order = None
    unsat_index = None

    for a1, a2, order in view.bonds():
        if order == 1:
            continue
        if order not in {2, 3}:
            raise ChemNameNotSupported("Unsupported bond order")
        if a1 in chain_set and a2 in chain_set and abs(index[a1] - index[a2]) == 1:
            idx = min(index[a1], index[a2])
            if unsat_order is None:
                unsat_order = order
                unsat_index = idx
            else:
                raise ChemNameNotSupported("Multiple unsaturations not supported")
        else:
            if strict_outside_chain:
                raise ChemNameNotSupported("Unsaturation outside parent chain")

    return unsat_order, unsat_index


def _find_functional_group(
    view: MolView, chain: list[int], allow_other_hetero: bool = False
) -> tuple[str, int, int] | None:
    chain_set = set(chain)
    candidates: list[tuple[str, int, int]] = []
    for atom_id in chain:
        for nbr in view.neighbors(atom_id):
            if nbr in chain_set:
                continue
            elem = view.element(nbr)
            if elem not in {"O", "N"}:
                continue
            if view.bond_order_between(atom_id, nbr) != 1:
                raise ChemNameNotSupported("Unsupported functional group bond order")
            h_total = implicit_h_count(view, nbr) + view.explicit_h(nbr)
            if h_total < 1:
                raise ChemNameNotSupported("Unsupported functional group")
            suffix = "ol" if elem == "O" else "amine"
            candidates.append((suffix, atom_id, nbr))

    if len(candidates) > 1:
        raise ChemNameNotSupported("Multiple functional groups not supported")
    if not candidates:
        if not allow_other_hetero:
            _assert_no_heteroatoms(view, chain_set)
        return None
    suffix, carbon_id, hetero_id = candidates[0]
    _assert_no_heteroatoms(view, chain_set, {hetero_id})
    return suffix, carbon_id, hetero_id


def _trim_aromatic_linker(
    view: MolView, chain: list[int], ring_ctx: RingContext
) -> list[int]:
    if len(chain) <= 1:
        return chain
    chain_set = set(chain)

    def should_trim(atom_id: int) -> bool:
        if view.element(atom_id) != "C":
            return False
        heavy_neighbors = [n for n in view.neighbors(atom_id) if view.element(n) != "H"]
        if len(heavy_neighbors) != 2:
            return False
        chain_neighbors = [n for n in heavy_neighbors if n in chain_set]
        ring_neighbors = [n for n in heavy_neighbors if n not in chain_set]
        if len(chain_neighbors) != 1 or len(ring_neighbors) != 1:
            return False
        for ring in ring_ctx.atom_rings.get(ring_neighbors[0], []):
            if ring_ctx.ring_types.get(ring) == "benzene":
                ring_size = len(ring)
                if len(chain) - 1 > ring_size:
                    return True
        return False

    start, end = chain[0], chain[-1]
    if should_trim(start):
        return chain[1:]
    if should_trim(end):
        return chain[:-1]
    return chain


def _assert_no_heteroatoms(
    view: MolView, chain_set: set[int], allowed: set[int] | None = None
) -> None:
    allowed = allowed or set()
    for atom_id in view.atoms():
        if atom_id in chain_set or atom_id in allowed:
            continue
        elem = view.element(atom_id)
        if elem in {"O", "N"}:
            raise ChemNameNotSupported("Unsupported heteroatom")


def _locant_for_atom(chain: list[int], atom_id: int | None) -> int | None:
    if atom_id is None:
        return None
    for idx, cid in enumerate(chain):
        if cid == atom_id:
            return idx + 1
    return None


def _unsaturation_locant_for_chain(
    view: MolView, chain: list[int], unsat_order: int | None
) -> int | None:
    if unsat_order is None:
        return None
    locants = []
    for idx in range(len(chain) - 1):
        if view.bond_order_between(chain[idx], chain[idx + 1]) == unsat_order:
            locants.append(idx + 1)
    if len(locants) != 1:
        raise ChemNameNotSupported("Multiple unsaturations not supported")
    return locants[0]


def _choose_oriented_chain(
    view: MolView,
    chain: list[int],
    opts: NameOptions,
    func_atom: int | None,
    unsat_order: int | None,
    ignore_atoms: set[int],
    ring_ctx: RingContext | None = None,
) -> list[int]:
    forward = list(chain)
    reverse = list(reversed(chain))

    f_primary = [_locant_for_atom(forward, func_atom)] if func_atom is not None else []
    r_primary = [_locant_for_atom(reverse, func_atom)] if func_atom is not None else []
    f_secondary = (
        [_unsaturation_locant_for_chain(view, forward, unsat_order)]
        if unsat_order is not None
        else []
    )
    r_secondary = (
        [_unsaturation_locant_for_chain(view, reverse, unsat_order)]
        if unsat_order is not None
        else []
    )

    f_subs = substituents_on_chain(
        view, forward, ignore_atoms=ignore_atoms, ring_ctx=ring_ctx
    )
    r_subs = substituents_on_chain(
        view, reverse, ignore_atoms=ignore_atoms, ring_ctx=ring_ctx
    )
    f_key = orientation_key(f_subs, opts, primary_locants=f_primary, secondary_locants=f_secondary)
    r_key = orientation_key(r_subs, opts, primary_locants=r_primary, secondary_locants=r_secondary)

    if r_key < f_key:
        return reverse
    return forward
