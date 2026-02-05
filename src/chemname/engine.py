from __future__ import annotations

from pathlib import Path

from .errors import ChemNameInternalError, ChemNameNotSupported
from chemcalc.valence import implicit_h_count

from .locants import Sub, orientation_key, substituents_on_chain
from .molview import MolView
from .options import NameOptions
from .parent_chain import longest_carbon_chain, longest_chain_in_subset
from .render import render_name
from .ring_naming import (
    enumerate_ring_numberings,
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
    _ring_aromatic_basic,
)
from .substituents import CYCLO_PARENT, HALO_MAP, alkyl_substituent_name, parent_name
from .template import TemplateMol, load_template
from .template_match import match_template_exact, select_template_mapping

_TEMPLATE_CACHE: dict[str, TemplateMol] = {}
_TEMPLATE_DIR = Path(__file__).resolve().parent / "templates"


def _load_template_cached(key: str, path: Path) -> TemplateMol:
    if key in _TEMPLATE_CACHE:
        return _TEMPLATE_CACHE[key]
    template = load_template(path)
    _TEMPLATE_CACHE[key] = template
    return template


def _pyrene_template(scheme: str) -> TemplateMol:
    scheme_key = scheme.lower()
    if scheme_key == "cas":
        filename = "pyrene_cas.mol"
    else:
        filename = "pyrene_iupac2004.mol"
    path = _TEMPLATE_DIR / "fused" / filename
    return _load_template_cached(f"pyrene:{scheme_key}", path)


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

    func = _find_functional_group(view, chain)
    func_atom = func[1] if func else None
    func_suffix = func[0] if func else None
    ignore_atoms = func[2] if func else set()

    unsaturations = _unsaturations_for_chain(view, chain)
    if allow_rings and ring_ctx is not None and func is None and not unsaturations:
        chain = _trim_aromatic_linker(view, chain, ring_ctx)

    chain = _choose_oriented_chain(
        view,
        chain,
        opts,
        func_atom=func_atom,
        ignore_atoms=ignore_atoms,
        ring_ctx=ring_ctx,
    )

    substituents = substituents_on_chain(
        view, chain, ignore_atoms=ignore_atoms, ring_ctx=ring_ctx
    )
    func_locant = _locant_for_atom(chain, func_atom) if func_atom is not None else None
    unsaturations = _unsaturations_for_chain(view, chain)

    parent = parent_name(
        len(chain),
        unsaturations=unsaturations,
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

    pyrene = _detect_pyrene_template(view, opts)
    if pyrene is not None:
        return _name_pyrene_template(view, pyrene, opts)

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
    ring_unsats: list[int] = []
    for idx in range(len(ring_atoms)):
        order = view.bond_order_between(ring_atoms[idx], ring_atoms[(idx + 1) % len(ring_atoms)])
        if order == 1:
            continue
        if order != 2:
            raise ChemNameNotSupported("Unsupported ring unsaturation")
        ring_unsats.append(idx + 1)
    for atom_id in ring_atoms:
        if view.element(atom_id) != "C":
            raise ChemNameNotSupported("Non-carbon ring not supported")

    parent = CYCLO_PARENT.get(len(ring_atoms))
    if parent is None:
        raise ChemNameNotSupported("Unsupported ring size")

    if ring_unsats:
        if len(ring_unsats) > 2:
            raise ChemNameNotSupported("Multiple unsaturations not supported")
        best = None
        best_key = None
        for numbering in enumerate_ring_numberings(ring_atoms):
            unsat_locants: list[int] = []
            for idx in range(len(numbering)):
                order = view.bond_order_between(
                    numbering[idx], numbering[(idx + 1) % len(numbering)]
                )
                if order == 2:
                    unsat_locants.append(idx + 1)
            subs = ring_substituents(view, numbering, ring_ctx=ring_ctx)
            key = orientation_key(subs, opts, primary_locants=sorted(unsat_locants))
            if best_key is None or key < best_key:
                best_key = key
                best = numbering
        oriented = list(best) if best is not None else list(ring_atoms)
        unsat_locants = []
        for idx in range(len(oriented)):
            if view.bond_order_between(oriented[idx], oriented[(idx + 1) % len(oriented)]) == 2:
                unsat_locants.append(idx + 1)
        if len(unsat_locants) == 1:
            parent_name_cyclo = parent[:-3] if parent.endswith("ane") else parent
            parent_name_cyclo = f"{parent_name_cyclo}-{unsat_locants[0]}-ene"
        else:
            locants = ",".join(str(loc) for loc in sorted(unsat_locants))
            parent_name_cyclo = parent[:-3] if parent.endswith("ane") else parent
            parent_name_cyclo = f"{parent_name_cyclo}-{locants}-diene"
        substituents = ring_substituents(view, oriented, ring_ctx=ring_ctx)
        return render_name(substituents, parent_name_cyclo, always_include_locant=False)

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
        allow_amino=True,
        allow_alkoxy=True,
        ring_ctx=ring_ctx,
    )
    substituents = ring_substituents(
        view,
        oriented,
        allow_hydroxy=True,
        allow_nitro=True,
        allow_amino=True,
        allow_alkoxy=True,
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
    preferred_start = info.get("preferred_start")
    if not kind or not ring_atoms or not hetero_atoms:
        raise ChemNameNotSupported("Unsupported heteroaromatic ring")

    oriented = choose_hetero_ring_orientation(
        view,
        ring_atoms,
        hetero_atoms,
        opts,
        allow_hydroxy=True,
        allow_nitro=False,
        allow_amino=True,
        allow_alkoxy=True,
        preferred_start_atoms=preferred_start,
        forbid_hetero_substituents=True,
        ring_ctx=ring_ctx,
    )
    substituents = ring_substituents(
        view,
        oriented,
        allow_hydroxy=True,
        allow_nitro=False,
        allow_amino=True,
        allow_alkoxy=True,
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


def _detect_pyrene(view: MolView, rings: list[frozenset[int]]) -> dict | None:
    benzene_rings: list[frozenset[int]] = []
    ring_edges: dict[frozenset[int], set[frozenset[int]]] = {}
    for ring in rings:
        info = classify_aromatic_ring(view, ring)
        if not info or info.get("kind") != "benzene":
            continue
        benzene_rings.append(ring)
        order = ring_order(view, ring)
        if not order:
            continue
        edges = set()
        for idx in range(len(order)):
            a = order[idx]
            b = order[(idx + 1) % len(order)]
            edges.add(frozenset({a, b}))
        ring_edges[ring] = edges

    if len(benzene_rings) < 4:
        return None

    adjacency: dict[frozenset[int], list[frozenset[int]]] = {r: [] for r in benzene_rings}
    for i in range(len(benzene_rings)):
        for j in range(i + 1, len(benzene_rings)):
            r1 = benzene_rings[i]
            r2 = benzene_rings[j]
            shared = ring_edges.get(r1, set()) & ring_edges.get(r2, set())
            if shared:
                adjacency[r1].append(r2)
                adjacency[r2].append(r1)

    visited: set[frozenset[int]] = set()
    for ring in benzene_rings:
        if ring in visited:
            continue
        stack = [ring]
        component: list[frozenset[int]] = []
        while stack:
            current = stack.pop()
            if current in visited:
                continue
            visited.add(current)
            component.append(current)
            for nbr in adjacency.get(current, []):
                if nbr not in visited:
                    stack.append(nbr)
        if len(component) != 4:
            continue
        if any(len([nbr for nbr in adjacency.get(r, []) if nbr in component]) != 2 for r in component):
            continue
        union = set().union(*component)
        if len(union) != 16:
            continue
        if any(view.element(atom_id) != "C" for atom_id in union):
            continue
        fusion_edges: set[frozenset[int]] = set()
        for i in range(len(component)):
            for j in range(i + 1, len(component)):
                shared = ring_edges.get(component[i], set()) & ring_edges.get(component[j], set())
                if shared:
                    fusion_edges |= shared
        if len(fusion_edges) != 4:
            continue
        return {
            "atoms": union,
            "fusion_edges": fusion_edges,
            "rings": component,
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


def _aromatic_components(view: MolView) -> list[set[int]]:
    rings = find_rings_simple(view)
    aromatic_atoms: set[int] = set()
    for ring in rings:
        if _ring_aromatic_basic(view, ring):
            aromatic_atoms |= set(ring)

    components: list[set[int]] = []
    visited: set[int] = set()
    for atom_id in aromatic_atoms:
        if atom_id in visited:
            continue
        stack = [atom_id]
        comp: set[int] = set()
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            comp.add(node)
            for nbr in view.neighbors(node):
                if nbr in aromatic_atoms and nbr not in visited:
                    stack.append(nbr)
        if comp:
            components.append(comp)
    return components


def _detect_pyrene_template(view: MolView, opts: NameOptions) -> dict | None:
    template = _pyrene_template(opts.fused_numbering_scheme)
    for comp in _aromatic_components(view):
        if len(comp) != len(template.atoms):
            continue
        mappings = match_template_exact(template, view, atom_ids=comp)
        if not mappings:
            continue
        mapping = select_template_mapping(template, view, mappings)
        if mapping is None:
            continue
        return {"template": template, "mapping": mapping}
    return None


def _name_pyrene_template(view: MolView, info: dict, opts: NameOptions) -> str:
    template: TemplateMol = info.get("template")
    mapping: dict[int, int] = info.get("mapping", {})
    if template is None or not mapping:
        raise ChemNameNotSupported("Pyrene template match failed")

    ring_atoms = set(mapping.values())
    inverse = {mol_id: t_idx for t_idx, mol_id in mapping.items()}
    substituents: list[tuple[int, str]] = []
    for atom_id in ring_atoms:
        for nbr in view.neighbors(atom_id):
            if nbr in ring_atoms:
                continue
            if view.element(nbr) == "H":
                continue
            t_idx = inverse.get(atom_id)
            if t_idx is None:
                raise ChemNameNotSupported("Invalid pyrene mapping")
            locant = template.locant_by_atom_idx.get(t_idx)
            if locant is None:
                raise ChemNameNotSupported("Unsupported pyrene substitution")
            name = _simple_substituent_name(view, atom_id, nbr, ring_atoms)
            substituents.append((locant, name))

    if len(substituents) > 2:
        raise ChemNameNotSupported("Multiple substitutions not supported")

    subs = [Sub(name, locant) for locant, name in substituents]
    return render_name(subs, "pyrene", always_include_locant=True)


def _name_pyrene(view: MolView, info: dict, opts: NameOptions) -> str:
    ring_atoms = set(info.get("atoms", set()))
    fusion_edges = info.get("fusion_edges", set())
    rings = info.get("rings", [])
    if len(ring_atoms) != 16 or len(fusion_edges) != 4:
        raise ChemNameNotSupported("Unsupported fused ring system")

    fusion_atoms = set()
    for edge in fusion_edges:
        fusion_atoms |= set(edge)

    order: list[int] = []
    if rings:
        ring_edges: list[set[frozenset[int]]] = []
        for ring in rings:
            order_ring = ring_order(view, ring)
            if not order_ring:
                continue
            edges = set()
            for idx in range(len(order_ring)):
                a = order_ring[idx]
                b = order_ring[(idx + 1) % len(order_ring)]
                edges.add(frozenset({a, b}))
            ring_edges.append(edges)
        edge_count: dict[frozenset[int], int] = {}
        for edges in ring_edges:
            for edge in edges:
                edge_count[edge] = edge_count.get(edge, 0) + 1
        perim_edges = [edge for edge, count in edge_count.items() if count == 1]
        adjacency: dict[int, list[int]] = {}
        for edge in perim_edges:
            a, b = tuple(edge)
            adjacency.setdefault(a, []).append(b)
            adjacency.setdefault(b, []).append(a)
        if adjacency and all(len(nbrs) == 2 for nbrs in adjacency.values()):
            start = min(adjacency)
            order = _cycle_order({k: sorted(v) for k, v in adjacency.items()}, start)

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
    locant_maps: list[dict[int, int]] = []
    if order:
        locant_maps = _fused_locant_maps(order, fusion_atoms)
    else:
        non_fusion = sorted(atom_id for atom_id in ring_atoms if atom_id not in fusion_atoms)
        locant_maps = [{atom_id: idx + 1 for idx, atom_id in enumerate(non_fusion)}]

    for locant_map in locant_maps:
        subs = [Sub(name, locant_map[atom_id]) for atom_id, name in substituents]
        locants = sorted(sub.locant for sub in subs)
        key = orientation_key(subs, opts, primary_locants=locants)
        if best_key is None or key < best_key:
            best_key = key
            best_subs = subs

    return render_name(best_subs, "pyrene", always_include_locant=True)


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
            if kinds == {"benzene", "furan"} and len(union) == 9:
                furan_ring = r1 if info1.get("kind") == "furan" else r2
                hetero_atom = next(
                    atom_id for atom_id in furan_ring if view.element(atom_id) == "O"
                )
                return {
                    "kind": "benzofuran",
                    "atoms": union,
                    "fusion_atoms": set(shared),
                    "hetero_atom": hetero_atom,
                }
            if kinds == {"benzene", "thiophene"} and len(union) == 9:
                thio_ring = r1 if info1.get("kind") == "thiophene" else r2
                hetero_atom = next(
                    atom_id for atom_id in thio_ring if view.element(atom_id) == "S"
                )
                return {
                    "kind": "benzothiophene",
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


def _unsaturations_for_chain(view: MolView, chain: list[int]) -> list[tuple[int, int]]:
    unsaturations: list[tuple[int, int]] = []
    for idx in range(len(chain) - 1):
        order = view.bond_order_between(chain[idx], chain[idx + 1])
        if order == 1:
            continue
        if order not in {2, 3}:
            raise ChemNameNotSupported("Unsupported bond order")
        unsaturations.append((order, idx + 1))
    return unsaturations


def _find_functional_group(
    view: MolView, chain: list[int], allow_other_hetero: bool = False
) -> tuple[str, int, set[int]] | None:
    chain_set = set(chain)
    acids: list[tuple[int, set[int]]] = []
    aldehydes: list[tuple[int, set[int]]] = []
    ketones: list[tuple[int, set[int]]] = []
    nitriles: list[tuple[int, set[int]]] = []
    alcohols: list[tuple[int, set[int]]] = []
    amines: list[tuple[int, set[int]]] = []

    for atom_id in chain:
        chain_neighbors = [nbr for nbr in view.neighbors(atom_id) if nbr in chain_set]
        outside_neighbors = [
            nbr
            for nbr in view.neighbors(atom_id)
            if nbr not in chain_set and view.element(nbr) != "H"
        ]

        carbonyl_oxygen = [
            nbr
            for nbr in outside_neighbors
            if view.element(nbr) == "O" and view.bond_order_between(atom_id, nbr) == 2
        ]
        single_oxygen = [
            nbr
            for nbr in outside_neighbors
            if view.element(nbr) == "O" and view.bond_order_between(atom_id, nbr) == 1
        ]
        nitrile_n = [
            nbr
            for nbr in outside_neighbors
            if view.element(nbr) == "N" and view.bond_order_between(atom_id, nbr) == 3
        ]

        acid_oxygen = None
        if len(carbonyl_oxygen) == 1 and len(single_oxygen) == 1:
            o_single = single_oxygen[0]
            h_total = implicit_h_count(view, o_single) + view.explicit_h(o_single)
            heavy_neighbors = [n for n in view.neighbors(o_single) if view.element(n) != "H"]
            if h_total >= 1 and len(heavy_neighbors) == 1 and len(chain_neighbors) == 1:
                acids.append((atom_id, {carbonyl_oxygen[0], o_single}))
                acid_oxygen = o_single

        if len(carbonyl_oxygen) == 1 and len(outside_neighbors) == 1:
            if len(chain_neighbors) == 1:
                aldehydes.append((atom_id, {carbonyl_oxygen[0]}))
            elif len(chain_neighbors) == 2:
                ketones.append((atom_id, {carbonyl_oxygen[0]}))

        if len(nitrile_n) == 1 and len(outside_neighbors) == 1 and len(chain_neighbors) == 1:
            nitriles.append((atom_id, {nitrile_n[0]}))

        for nbr in outside_neighbors:
            if acid_oxygen is not None and nbr == acid_oxygen:
                continue
            elem = view.element(nbr)
            if elem == "O" and view.bond_order_between(atom_id, nbr) == 1:
                h_total = implicit_h_count(view, nbr) + view.explicit_h(nbr)
                heavy_neighbors = [n for n in view.neighbors(nbr) if view.element(n) != "H"]
                if h_total >= 1 and len(heavy_neighbors) == 1:
                    alcohols.append((atom_id, {nbr}))
            if elem == "N" and view.bond_order_between(atom_id, nbr) == 1:
                h_total = implicit_h_count(view, nbr) + view.explicit_h(nbr)
                heavy_neighbors = [n for n in view.neighbors(nbr) if view.element(n) != "H"]
                if h_total >= 1 and len(heavy_neighbors) == 1:
                    amines.append((atom_id, {nbr}))

    def _select(
        group: list[tuple[int, set[int]]], suffix: str
    ) -> tuple[str, int, set[int]] | None:
        if not group:
            return None
        if len(group) > 1:
            raise ChemNameNotSupported("Multiple functional groups not supported")
        carbon_id, hetero_ids = group[0]
        return suffix, carbon_id, hetero_ids

    for group, suffix in (
        (acids, "oic acid"),
        (aldehydes, "al"),
        (ketones, "one"),
        (nitriles, "nitrile"),
        (alcohols, "ol"),
        (amines, "amine"),
    ):
        selected = _select(group, suffix)
        if selected is None:
            continue
        if any(other for other in (acids, aldehydes, ketones, nitriles, alcohols, amines) if other is not group and other):
            raise ChemNameNotSupported("Multiple functional groups not supported")
        suffix, carbon_id, hetero_ids = selected
        return suffix, carbon_id, hetero_ids

    return None


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


def _choose_oriented_chain(
    view: MolView,
    chain: list[int],
    opts: NameOptions,
    func_atom: int | None,
    ignore_atoms: set[int],
    ring_ctx: RingContext | None = None,
) -> list[int]:
    forward = list(chain)
    reverse = list(reversed(chain))

    f_primary = [_locant_for_atom(forward, func_atom)] if func_atom is not None else []
    r_primary = [_locant_for_atom(reverse, func_atom)] if func_atom is not None else []
    f_secondary = [loc for _order, loc in _unsaturations_for_chain(view, forward)]
    r_secondary = [loc for _order, loc in _unsaturations_for_chain(view, reverse)]

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
