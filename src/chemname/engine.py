from __future__ import annotations

from .errors import ChemNameInternalError, ChemNameNotSupported
from chemcalc.valence import implicit_h_count

from .locants import Sub, orientation_key, substituents_on_chain
from .molview import MolView
from .options import NameOptions
from .parent_chain import longest_carbon_chain
from .render import render_name
from .ring_naming import (
    choose_hetero_ring_orientation,
    choose_ring_orientation,
    ring_substituents,
)
from .rings import (
    classify_aromatic_ring,
    detect_naphthalene,
    find_rings_simple,
    is_simple_ring,
    ring_order,
)
from .substituents import ALKYL, CYCLO_PARENT, HALO_MAP, alkyl_length_linear, parent_name


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
    rings = find_rings_simple(view)
    if not rings:
        return _name_linear(view, opts)
    return _name_cyclic(view, rings, opts)


def _name_linear(view: MolView, opts: NameOptions) -> str:
    if not view.is_acyclic():
        raise ChemNameNotSupported("Cyclic structures not supported")

    allowed_elements = {"C", "H", "F", "Cl", "Br", "I", "O", "N"}
    for atom_id in view.atoms():
        elem = view.element(atom_id)
        if elem not in allowed_elements:
            raise ChemNameNotSupported("Unsupported element")

    chain = longest_carbon_chain(view)
    if not chain:
        raise ChemNameNotSupported("No carbon chain found")

    unsat_order, _ = _find_unsaturation(view, chain)
    func = _find_functional_group(view, chain)
    func_atom = func[1] if func else None
    func_suffix = func[0] if func else None
    ignore_atoms = {func[2]} if func else set()

    chain = _choose_oriented_chain(
        view,
        chain,
        opts,
        func_atom=func_atom,
        unsat_order=unsat_order,
        ignore_atoms=ignore_atoms,
    )

    substituents = substituents_on_chain(view, chain, ignore_atoms=ignore_atoms)
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


def _name_cyclic(view: MolView, rings: list[frozenset[int]], opts: NameOptions) -> str:
    if len(rings) == 1:
        ring_nodes = rings[0]
        if not is_simple_ring(view, ring_nodes):
            raise ChemNameNotSupported("Fused rings not supported")

        ring_atoms = ring_order(view, ring_nodes)
        if not ring_atoms:
            raise ChemNameNotSupported("Ring ordering failed")

        aromatic_info = classify_aromatic_ring(view, ring_nodes)
        if aromatic_info and aromatic_info.get("kind") == "benzene":
            return _name_benzene(view, ring_atoms, opts)
        if aromatic_info and aromatic_info.get("kind") is not None:
            return _name_heteroaromatic(view, aromatic_info, opts)

        return _name_cycloalkane(view, ring_atoms, opts)

    naph = detect_naphthalene(view, rings)
    if naph is not None:
        return _name_naphthalene(view, naph, opts)

    raise ChemNameNotSupported("Multiple rings not supported")


def _name_cycloalkane(view: MolView, ring_atoms: list[int], opts: NameOptions) -> str:
    for idx in range(len(ring_atoms)):
        if view.bond_order_between(ring_atoms[idx], ring_atoms[(idx + 1) % len(ring_atoms)]) != 1:
            raise ChemNameNotSupported("Unsaturated ring not supported")
    for atom_id in ring_atoms:
        if view.element(atom_id) != "C":
            raise ChemNameNotSupported("Non-carbon ring not supported")

    parent = CYCLO_PARENT.get(len(ring_atoms))
    if parent is None:
        raise ChemNameNotSupported("Unsupported ring size")

    oriented = choose_ring_orientation(view, ring_atoms, opts)
    substituents = ring_substituents(view, oriented)
    return render_name(substituents, parent, always_include_locant=False)


def _name_benzene(view: MolView, ring_atoms: list[int], opts: NameOptions) -> str:
    for atom_id in ring_atoms:
        if view.element(atom_id) != "C":
            raise ChemNameNotSupported("Unsupported aromatic ring")
    oriented = choose_ring_orientation(
        view,
        ring_atoms,
        opts,
        allow_hydroxy=True,
        allow_nitro=True,
    )
    substituents = ring_substituents(
        view,
        oriented,
        allow_hydroxy=True,
        allow_nitro=True,
    )
    return render_name(substituents, "benzene", always_include_locant=False)


def _name_heteroaromatic(view: MolView, info: dict, opts: NameOptions) -> str:
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
    )
    substituents = ring_substituents(
        view,
        oriented,
        allow_hydroxy=True,
        allow_nitro=False,
        forbid_hetero_substituents=True,
    )
    return render_name(substituents, kind, always_include_locant=False)


def _name_naphthalene(view: MolView, info: dict, opts: NameOptions) -> str:
    ring_atoms = set(info.get("atoms", set()))
    fusion_atoms = set(info.get("fusion_atoms", ()))
    if len(ring_atoms) != 10 or len(fusion_atoms) != 2:
        raise ChemNameNotSupported("Unsupported fused ring")

    # classify alpha/beta positions
    alpha_atoms: set[int] = set()
    for fusion in fusion_atoms:
        for nbr in view.neighbors(fusion):
            if nbr in ring_atoms and nbr not in fusion_atoms:
                alpha_atoms.add(nbr)
    beta_atoms = ring_atoms - fusion_atoms - alpha_atoms

    subs: list[Sub] = []
    for atom_id in ring_atoms:
        for nbr in view.neighbors(atom_id):
            if nbr in ring_atoms:
                continue
            if view.element(nbr) == "H":
                continue
            if atom_id in fusion_atoms:
                raise ChemNameNotSupported("Fusion-atom substitution not supported")
            name = _simple_substituent_name(view, atom_id, nbr, ring_atoms)
            locant = 1 if atom_id in alpha_atoms else 2
            subs.append(Sub(name, locant))

    if len(subs) > 1:
        raise ChemNameNotSupported("Multiple substitutions not supported")
    return render_name(subs, "naphthalene", always_include_locant=True)


def _simple_substituent_name(view: MolView, ring_atom: int, nbr: int, ring_set: set[int]) -> str:
    elem = view.element(nbr)
    if elem in HALO_MAP:
        return HALO_MAP[elem]
    if elem == "C":
        length = alkyl_length_linear(view, nbr, ring_set)
        name = ALKYL.get(length)
        if name is None:
            raise ChemNameNotSupported("Unsupported alkyl length")
        return name
    raise ChemNameNotSupported("Unsupported substituent")


def _find_unsaturation(view: MolView, chain: list[int]) -> tuple[int | None, int | None]:
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
            raise ChemNameNotSupported("Unsaturation outside parent chain")

    return unsat_order, unsat_index


def _find_functional_group(
    view: MolView, chain: list[int]
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
        _assert_no_heteroatoms(view, chain_set)
        return None
    suffix, carbon_id, hetero_id = candidates[0]
    _assert_no_heteroatoms(view, chain_set, {hetero_id})
    return suffix, carbon_id, hetero_id


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

    f_subs = substituents_on_chain(view, forward, ignore_atoms=ignore_atoms)
    r_subs = substituents_on_chain(view, reverse, ignore_atoms=ignore_atoms)
    f_key = orientation_key(f_subs, opts, primary_locants=f_primary, secondary_locants=f_secondary)
    r_key = orientation_key(r_subs, opts, primary_locants=r_primary, secondary_locants=r_secondary)

    if r_key < f_key:
        return reverse
    return forward
