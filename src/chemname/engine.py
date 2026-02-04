from __future__ import annotations

from .errors import ChemNameInternalError, ChemNameNotSupported
from chemcalc.valence import implicit_h_count

from .locants import orientation_key, substituents_on_chain
from .molview import MolView
from .options import NameOptions
from .parent_chain import longest_carbon_chain
from .render import render_name
from .substituents import parent_name


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

    unsat_order, unsat_index = _find_unsaturation(view, chain)
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
