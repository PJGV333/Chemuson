"""Motor principal de nomenclatura IUPAC-lite.

Orquesta la selección de cadena principal, detección de anillos y
composición final del nombre para las moléculas soportadas.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

from .errors import ChemNameInternalError, ChemNameNotSupported
from chemuson.chemcalc.valence import implicit_h_count
from chemuson.chemio.rdkit_io import molgraph_to_rdkit_with_map
from chemuson.utils.resources import open_resource_path

from .coordination import detect_coordination_name
from .functional_groups import detect_sulfonic_attachment
from .locants import Sub, orientation_key, substituents_on_chain
from .molview import MolView
from .options import NameOptions
from .parent_chain import longest_carbon_chain, longest_chain_in_subset
from .render import MULTIPLIER, render_name
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
from .substituents import (
    CYCLO_PARENT,
    HALO_MAP,
    alkyl_substituent_name,
    parent_name,
)
from .template import TemplateMol, load_template
from .template_match import match_template_exact, select_template_mapping

try:
    from rdkit import Chem
except Exception:  # pragma: no cover - runtime optional dependency
    Chem = None

try:
    from chemuson.chemio.rdkit_safe import stereo_descriptors_for_chain as safe_stereo_descriptors_for_chain
except Exception:  # pragma: no cover - runtime optional dependency
    safe_stereo_descriptors_for_chain = None

_TEMPLATE_CACHE: dict[str, TemplateMol] = {}


@dataclass
class FunctionalOccurrence:
    """Representa una ocurrencia de grupo funcional sobre la cadena."""

    kind: str
    atom_id: int
    aux_atom_ids: set[int]
    prefix_name: str
    suffix_name: str


@dataclass
class FunctionalSelection:
    """Grupo funcional seleccionado para sufijo y prefijos residuales."""

    suffix: str
    atom_id: int
    ignore_atoms: set[int]
    prefixes: list[tuple[str, int]]


def _load_template_cached(key: str, path: str | Path) -> TemplateMol:
    """Carga una plantilla con caché en memoria.

    Args:
        key: Clave única para identificar la plantilla.
        path: Ruta del archivo de plantilla.

    Returns:
        `TemplateMol` cargada o reutilizada desde caché.
    """
    if key in _TEMPLATE_CACHE:
        return _TEMPLATE_CACHE[key]
    template = load_template(path)
    _TEMPLATE_CACHE[key] = template
    return template


def _pyrene_template(scheme: str) -> TemplateMol:
    """Selecciona la plantilla de pireno según el esquema de numeración.

    Args:
        scheme: Nombre del esquema ("cas" o "iupac2004").

    Returns:
        Plantilla de pireno correspondiente.
    """
    scheme_key = scheme.lower()
    if scheme_key == "cas":
        filename = "pyrene_cas.mol"
    else:
        filename = "pyrene_iupac2004.mol"
    with open_resource_path("chemname", "templates", "fused", filename) as path:
        return _load_template_cached(f"pyrene:{scheme_key}", path)


def iupac_name(graph, opts: NameOptions = NameOptions()) -> str:
    """Punto de entrada público del motor de nombres.

    Args:
        graph: Grafo molecular a nombrar.
        opts: Opciones del motor.

    Returns:
        Nombre IUPAC-lite o "N/D" si `return_nd_on_fail` está habilitado.

    Raises:
        ChemNameNotSupported: Si la molécula no es soportada y no se pide "N/D".
        ChemNameInternalError: Si ocurre un error interno inesperado.
    """
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
    """Nomenclatura IUPAC-lite para moléculas sencillas.

    Decide si usar la ruta lineal o cíclica, priorizando sistemas fusionados
    y anillos aromáticos cuando corresponda.

    Args:
        graph: Grafo molecular a nombrar.
        opts: Opciones del motor.

    Returns:
        Nombre IUPAC-lite generado.
    """
    view = MolView(graph)
    if opts.allow_coordination:
        coordination_name = detect_coordination_name(view)
        if coordination_name:
            return coordination_name
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

    # Para anillos simples, priorizamos cadena cuando supera al anillo.
    if chain and len(rings) == 1:
        ring = rings[0]
        rtype = ring_ctx.ring_types.get(ring)
        if rtype == "benzene" and len(chain) > 6:
            return _name_linear(
                view,
                opts,
                chain_override=chain,
                ring_ctx=ring_ctx,
                allow_rings=True,
            )
        if rtype == "cyclohexane" and len(chain) > len(ring):
            return _name_linear(
                view,
                opts,
                chain_override=chain,
                ring_ctx=ring_ctx,
                allow_rings=True,
            )

    # Sistemas fusionados primero: si se reconoce, se retorna de inmediato.
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
    """Nombra una cadena lineal (o linealizada) con sustituyentes simples.

    Args:
        view: Vista del grafo molecular.
        opts: Opciones del motor.
        chain_override: Cadena principal preseleccionada (opcional).
        ring_ctx: Contexto de anillos (opcional, para sustituyentes).
        allow_rings: Permite enlazar anillos como sustituyentes.

    Returns:
        Nombre IUPAC-lite de la cadena.

    Raises:
        ChemNameNotSupported: Si hay elementos o estructuras no soportadas.
    """
    if not allow_rings and not view.is_acyclic():
        raise ChemNameNotSupported("Cyclic structures not supported")

    allowed_elements = {"C", "H", "F", "Cl", "Br", "I", "O", "N", "S", "P", "Si", "B"}
    for atom_id in view.atoms():
        elem = view.element(atom_id)
        if elem not in allowed_elements:
            raise ChemNameNotSupported("Unsupported element")

    # Selección de cadena principal (o la proporcionada por el llamador).
    chain = chain_override or longest_carbon_chain(view)
    if not chain:
        raise ChemNameNotSupported("No carbon chain found")

    # Identificar grupo funcional principal (ácido, aldehído, etc.).
    func = _find_functional_group(view, chain)
    func_atom = func.atom_id if func else None
    func_suffix = func.suffix if func else None
    ignore_atoms = set(func.ignore_atoms) if func else set()
    functional_prefixes = list(func.prefixes) if func else []

    # Determinar insaturaciones C=C/C#C dentro de la cadena.
    unsaturations = _unsaturations_for_chain(view, chain)
    if allow_rings and ring_ctx is not None and func is None and not unsaturations:
        # Recorta un carbono "puente" entre anillo aromático y cadena.
        chain = _trim_aromatic_linker(view, chain, ring_ctx)

    chain = _choose_oriented_chain(
        view,
        chain,
        opts,
        func_atom=func_atom,
        ignore_atoms=ignore_atoms,
        functional_prefixes=functional_prefixes,
        ring_ctx=ring_ctx,
    )

    # Recolectamos sustituyentes y calculamos el nombre del padre.
    substituents = substituents_on_chain(view, chain, ignore_atoms=ignore_atoms, ring_ctx=ring_ctx)
    substituents.extend(_isotope_substituents_on_chain(view, chain))
    for prefix_name, prefix_atom in functional_prefixes:
        locant = _locant_for_atom(chain, prefix_atom)
        if locant is None:
            continue
        substituents.append(Sub(prefix_name, locant))
    func_locant = _locant_for_atom(chain, func_atom) if func_atom is not None else None
    unsaturations = _unsaturations_for_chain(view, chain)

    parent = parent_name(
        len(chain),
        unsaturations=unsaturations,
        suffix=func_suffix,
        suffix_locant=func_locant,
    )

    stereo = _stereo_descriptors_for_linear(view, chain, opts)
    rendered = render_name(substituents, parent, stereo_descriptors=stereo)
    return _apply_radical_suffix_if_needed(view, chain, substituents, rendered, opts)


def _longest_chain_excluding(view: MolView, excluded_atoms: set[int]) -> list[int]:
    """Devuelve la cadena más larga excluyendo un conjunto de átomos.

    Args:
        view: Vista del grafo molecular.
        excluded_atoms: Átomos que no pueden pertenecer a la cadena.

    Returns:
        Lista de IDs de la cadena más larga encontrada.
    """
    allowed = {atom_id for atom_id in view.atoms() if atom_id not in excluded_atoms}
    return longest_chain_in_subset(view, allowed)


def _name_cyclic(
    view: MolView,
    rings: list[frozenset[int]],
    opts: NameOptions,
    ring_ctx: RingContext | None = None,
) -> str:
    """Nombra sistemas cíclicos simples o fusionados.

    Args:
        view: Vista del grafo molecular.
        rings: Lista de anillos detectados.
        opts: Opciones del motor.
        ring_ctx: Contexto de anillos para sustituyentes.

    Returns:
        Nombre IUPAC-lite de la estructura cíclica.

    Raises:
        ChemNameNotSupported: Si la topología no está soportada.
    """
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

    spiro_bicyclo = _detect_spiro_bicyclo(view, rings)
    if spiro_bicyclo is not None:
        return _name_spiro_bicyclo(view, spiro_bicyclo, opts)

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
    """Nombra cicloalcanos simples con sustituyentes.

    Args:
        view: Vista del grafo molecular.
        ring_atoms: Lista de átomos en orden circular.
        opts: Opciones del motor.
        ring_ctx: Contexto de anillos para sustituyentes.

    Returns:
        Nombre IUPAC-lite del cicloalcano.
    """
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
        best = None
        best_key = None
        # Probamos todas las numeraciones para minimizar locantes.
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
        count = len(unsat_locants)
        if count == 1:
            unsat_suffix = "ene"
        else:
            prefix = MULTIPLIER.get(count)
            if prefix is None:
                raise ChemNameNotSupported("Too many double bonds in ring")
            unsat_suffix = f"{prefix}ene"
        locants = ",".join(str(loc) for loc in sorted(unsat_locants))
        parent_name_cyclo = parent[:-3] if parent.endswith("ane") else parent
        parent_name_cyclo = f"{parent_name_cyclo}-{locants}-{unsat_suffix}"
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
    """Nombra benceno con sustituyentes simples.

    Args:
        view: Vista del grafo molecular.
        ring_atoms: Átomos del anillo en orden.
        opts: Opciones del motor.
        ring_ctx: Contexto de anillos.

    Returns:
        Nombre IUPAC-lite del benceno sustituido.
    """
    for atom_id in ring_atoms:
        if view.element(atom_id) != "C":
            raise ChemNameNotSupported("Unsupported aromatic ring")
    dione_name = _name_aromatic_dione(
        view,
        ring_atoms,
        parent_base="benzene",
        opts=opts,
        ring_ctx=ring_ctx,
    )
    if dione_name is not None:
        return dione_name
    sulfo_name = _name_benzene_sulfonic(view, ring_atoms, opts, ring_ctx)
    if sulfo_name is not None:
        return sulfo_name
    oriented = choose_ring_orientation(
        view,
        ring_atoms,
        opts,
        allow_hydroxy=True,
        allow_nitro=True,
        allow_amino=True,
        allow_alkoxy=True,
        allow_ester=True,
        allow_amide=True,
        allow_nitrile=True,
        ring_ctx=ring_ctx,
    )
    substituents = ring_substituents(
        view,
        oriented,
        allow_hydroxy=True,
        allow_nitro=True,
        allow_amino=True,
        allow_alkoxy=True,
        allow_ester=True,
        allow_amide=True,
        allow_nitrile=True,
        ring_ctx=ring_ctx,
    )
    return render_name(substituents, "benzene", always_include_locant=False)


def _name_benzene_sulfonic(
    view: MolView,
    ring_atoms: list[int],
    opts: NameOptions,
    ring_ctx: RingContext | None,
) -> str | None:
    """Nombra benzenos con grupo sulfónico/sulfonato como función principal."""
    ring_set = set(ring_atoms)
    matches: list[tuple[int, dict]] = []
    for atom_id in ring_atoms:
        for nbr in view.neighbors(atom_id):
            if nbr in ring_set or view.element(nbr) != "S":
                continue
            info = detect_sulfonic_attachment(view, nbr, ring_set)
            if info is None:
                continue
            if int(info.get("parent_atom", -1)) != atom_id:
                continue
            matches.append((nbr, info))
    if not matches:
        return None
    if len(matches) > 1:
        raise ChemNameNotSupported("Multiple sulfonic groups not supported")
    sulfur_atom, info = matches[0]
    kind = str(info.get("kind"))
    parent = "benzenesulfonic acid" if kind == "sulfonic_acid" else "benzenesulfonate"
    ignore_atoms = set(info.get("aux_atoms", set())) | {sulfur_atom}

    oriented = choose_ring_orientation(
        view,
        ring_atoms,
        opts,
        allow_hydroxy=True,
        allow_nitro=True,
        allow_amino=True,
        allow_alkoxy=True,
        allow_ester=True,
        allow_amide=True,
        allow_nitrile=True,
        ignore_atoms=ignore_atoms,
        ring_ctx=ring_ctx,
    )
    substituents = ring_substituents(
        view,
        oriented,
        allow_hydroxy=True,
        allow_nitro=True,
        allow_amino=True,
        allow_alkoxy=True,
        allow_ester=True,
        allow_amide=True,
        allow_nitrile=True,
        ignore_atoms=ignore_atoms,
        ring_ctx=ring_ctx,
    )
    return render_name(substituents, parent, always_include_locant=bool(substituents))


def _name_aromatic_dione(
    view: MolView,
    ring_order_atoms: list[int],
    parent_base: str,
    opts: NameOptions,
    ring_ctx: RingContext | None = None,
    fusion_atoms: set[int] | None = None,
) -> str | None:
    """Nombra anillos aromáticos con dos carbonilos exocíclicos como dionas."""
    ring_set = set(ring_order_atoms)
    carbonyl_atoms: list[int] = []
    carbonyl_oxygens: set[int] = set()
    for atom_id in ring_order_atoms:
        oxy = _exocyclic_carbonyl_oxygen(view, atom_id, ring_set)
        if oxy is None:
            continue
        carbonyl_atoms.append(atom_id)
        carbonyl_oxygens.add(oxy)
    carbonyl_pair = list(carbonyl_atoms)
    if len(carbonyl_pair) != 2:
        return None

    if parent_base == "benzene":
        best_name = None
        best_key = None
        for numbering in enumerate_ring_numberings(ring_order_atoms):
            locant_map = {atom_id: idx + 1 for idx, atom_id in enumerate(numbering)}
            dione_locants = sorted(locant_map[atom_id] for atom_id in carbonyl_pair)
            parent = f"{parent_base}-{dione_locants[0]},{dione_locants[1]}-dione"
            subs = ring_substituents(
                view,
                numbering,
                allow_hydroxy=True,
                allow_nitro=True,
                allow_amino=True,
                allow_alkoxy=True,
                allow_ester=True,
                allow_amide=True,
                allow_nitrile=True,
                ignore_atoms=set(carbonyl_oxygens),
                ring_ctx=ring_ctx,
            )
            key = orientation_key(subs, opts, primary_locants=dione_locants)
            if best_key is None or key < best_key:
                best_key = key
                best_name = render_name(subs, parent, always_include_locant=bool(subs))
        return best_name

    if parent_base == "naphthalene":
        if not fusion_atoms:
            return None
        locant_maps = _naphthalene_locant_maps(ring_order_atoms, set(fusion_atoms))
        if not locant_maps:
            return None

        substituents: list[tuple[int, str]] = []
        for atom_id in ring_set:
            for nbr in view.neighbors(atom_id):
                if nbr in ring_set or nbr in carbonyl_oxygens or view.element(nbr) == "H":
                    continue
                if atom_id in set(fusion_atoms):
                    raise ChemNameNotSupported("Fusion-atom substitution not supported")
                substituents.append((atom_id, _simple_substituent_name(view, atom_id, nbr, ring_set)))

        best_name = None
        best_key = None
        for locant_map in locant_maps:
            if any(atom_id not in locant_map for atom_id in carbonyl_pair):
                continue
            dione_locants = sorted(locant_map[atom_id] for atom_id in carbonyl_pair)
            parent = f"{parent_base}-{dione_locants[0]},{dione_locants[1]}-dione"
            subs = [Sub(name, locant_map[atom_id]) for atom_id, name in substituents]
            key = orientation_key(subs, opts, primary_locants=dione_locants)
            if best_key is None or key < best_key:
                best_key = key
                best_name = render_name(subs, parent, always_include_locant=bool(subs))
        return best_name

    return None


def _exocyclic_carbonyl_oxygen(view: MolView, atom_id: int, ring_set: set[int]) -> int | None:
    """Devuelve oxígeno de carbonilo exocíclico en un átomo del anillo, si existe."""
    oxygens = [
        nbr
        for nbr in view.neighbors(atom_id)
        if nbr not in ring_set and view.element(nbr) == "O" and view.bond_order_between(atom_id, nbr) == 2
    ]
    if len(oxygens) != 1:
        return None
    return oxygens[0]


def _name_heteroaromatic(
    view: MolView,
    info: dict,
    opts: NameOptions,
    ring_ctx: RingContext | None = None,
) -> str:
    """Nombra anillos heteroaromáticos simples (furano, piridina, etc.).

    Args:
        view: Vista del grafo molecular.
        info: Diccionario de clasificación del anillo.
        opts: Opciones del motor.
        ring_ctx: Contexto de anillos.

    Returns:
        Nombre IUPAC-lite del heteroaromático.
    """
    kind = info.get("kind")
    ring_atoms = info.get("order") or []
    hetero_atoms = info.get("hetero_atoms") or []
    preferred_start = info.get("preferred_start")
    hetero_priority = info.get("hetero_priority") or {}
    if not kind or not ring_atoms or not hetero_atoms:
        raise ChemNameNotSupported("Unsupported heteroaromatic ring")
    if kind in {"phosphabenzene", "silabenzene", "borabenzene"} and not opts.enable_exotic_hetero:
        raise ChemNameNotSupported("Experimental heteroaromatic disabled")

    if not preferred_start:
        priority_to_atoms: dict[int, list[int]] = {}
        for atom_id in hetero_atoms:
            priority = hetero_priority.get(atom_id, 9)
            priority_to_atoms.setdefault(priority, []).append(atom_id)
        if priority_to_atoms:
            preferred_start = priority_to_atoms[min(priority_to_atoms.keys())]

    oriented = choose_hetero_ring_orientation(
        view,
        ring_atoms,
        hetero_atoms,
        opts,
        allow_hydroxy=True,
        allow_nitro=False,
        allow_amino=True,
        allow_alkoxy=True,
        allow_ester=True,
        allow_amide=True,
        allow_nitrile=True,
        preferred_start_atoms=preferred_start,
        hetero_priority=hetero_priority,
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
        allow_ester=True,
        allow_amide=True,
        allow_nitrile=True,
        forbid_hetero_substituents=True,
        ring_ctx=ring_ctx,
    )
    return render_name(substituents, kind, always_include_locant=False)


def _name_naphthalene(view: MolView, info: dict, opts: NameOptions) -> str:
    """Nombra derivados simples de naftaleno.

    Args:
        view: Vista del grafo molecular.
        info: Información de naftaleno (átomos y fusiones).
        opts: Opciones del motor.

    Returns:
        Nombre IUPAC-lite con locantes según reglas de orientación.
    """
    ring_atoms = set(info.get("atoms", set()))
    fusion_atoms = set(info.get("fusion_atoms", ()))
    if len(ring_atoms) != 10 or len(fusion_atoms) != 2:
        raise ChemNameNotSupported("Unsupported fused ring")

    order = _perimeter_cycle(view, ring_atoms, fusion_atoms)
    if not order:
        raise ChemNameNotSupported("Naphthalene ordering failed")

    dione_name = _name_aromatic_dione(
        view,
        order,
        parent_base="naphthalene",
        opts=opts,
        ring_ctx=None,
        fusion_atoms=fusion_atoms,
    )
    if dione_name is not None:
        return dione_name

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
    """Detecta sistemas tricíclicos tipo antraceno/fenantreno.

    Args:
        view: Vista del grafo molecular.
        rings: Anillos detectados.

    Returns:
        Diccionario con metadatos del sistema o `None` si no coincide.
    """
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

    # Construimos adyacencia de anillos por aristas compartidas.
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

    # Buscamos una cadena de tres anillos con el patrón correcto.
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
    """Detecta sistemas tipo pireno por conectividad de anillos bencénicos.

    Args:
        view: Vista del grafo molecular.
        rings: Anillos detectados.

    Returns:
        Diccionario con metadatos del sistema o `None` si no coincide.
    """
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
    """Nombra antraceno o fenantreno con sustitución simple.

    Args:
        view: Vista del grafo molecular.
        info: Metadatos del sistema tricíclico.
        opts: Opciones del motor.

    Returns:
        Nombre IUPAC-lite del sistema tricíclico.
    """
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


def _detect_spiro_bicyclo(view: MolView, rings: list[frozenset[int]]) -> dict | None:
    """Detecta sistemas spiro y bicíclicos saturados.

    Args:
        view: Vista del grafo molecular.
        rings: Lista de anillos detectados.

    Returns:
        Diccionario descriptivo del sistema o `None` si no coincide.
    """
    if len(rings) < 2:
        return None

    # Spiro simple: dos anillos comparten un único átomo.
    if len(rings) == 2:
        r1, r2 = rings
        shared = r1 & r2
        if len(shared) == 1:
            union = set(r1) | set(r2)
            if all(view.element(atom_id) == "C" for atom_id in union) and _is_saturated_framework(
                view, union
            ):
                bridges = tuple(sorted((len(r1) - 1, len(r2) - 1)))
                return {
                    "kind": "spiro",
                    "bridges": bridges,
                    "parent_len": sum(bridges) + 1,
                    "atoms": union,
                }

    # Bicíclico: dos átomos puente y tres caminos internamente disjuntos.
    union = set().union(*rings)
    if not union:
        return None
    if any(view.element(atom_id) != "C" for atom_id in union):
        return None
    if not _is_saturated_framework(view, union):
        return None

    candidate_bridgeheads = [
        atom_id
        for atom_id in union
        if len([nbr for nbr in view.neighbors(atom_id) if nbr in union]) >= 3
    ]
    for a, b in combinations(candidate_bridgeheads, 2):
        paths = _all_simple_paths_between(view, a, b, union)
        if len(paths) < 3:
            continue
        for combo in combinations(paths, 3):
            if not _paths_are_internally_disjoint(combo, a, b):
                continue
            used_atoms = set().union(*(set(path) for path in combo))
            if used_atoms != union:
                continue
            bridges = tuple(sorted((len(path) - 2 for path in combo), reverse=True))
            if len(bridges) != 3:
                continue
            if sum(bridges) + 2 != len(union):
                continue
            return {
                "kind": "bicyclo",
                "bridges": bridges,
                "parent_len": len(union),
                "atoms": union,
            }
    return None


def _name_spiro_bicyclo(view: MolView, info: dict, opts: NameOptions) -> str:
    """Renderiza nombres base para sistemas spiro/bicíclicos.

    Args:
        view: Vista del grafo molecular.
        info: Metadatos devueltos por `_detect_spiro_bicyclo`.
        opts: Opciones del motor (reservado para ampliaciones).

    Returns:
        Nombre del sistema (`spiro[...]alkane` o `bicyclo[...]alkane`).
    """
    kind = info.get("kind")
    bridges = tuple(info.get("bridges", ()))
    parent_len = int(info.get("parent_len", 0))
    parent = parent_name(parent_len)
    if not parent.endswith("ane"):
        raise ChemNameNotSupported("Unsupported spiro/bicyclo parent")
    if kind == "spiro" and len(bridges) == 2:
        return f"spiro[{bridges[0]}.{bridges[1]}]{parent}"
    if kind == "bicyclo" and len(bridges) == 3:
        return f"bicyclo[{bridges[0]}.{bridges[1]}.{bridges[2]}]{parent}"
    raise ChemNameNotSupported("Unsupported spiro/bicyclo system")


def _is_saturated_framework(view: MolView, atom_set: set[int]) -> bool:
    """Valida que todos los enlaces internos sean sencillos no aromáticos."""
    for atom_id in atom_set:
        for nbr in view.neighbors(atom_id):
            if nbr not in atom_set or nbr < atom_id:
                continue
            if view.bond_is_aromatic(atom_id, nbr):
                return False
            if view.bond_order_between(atom_id, nbr) != 1:
                return False
    return True


def _all_simple_paths_between(
    view: MolView,
    start: int,
    end: int,
    allowed_atoms: set[int],
    max_paths: int = 80,
) -> list[list[int]]:
    """Enumera caminos simples acotados entre dos nodos."""
    stack: list[tuple[int, list[int]]] = [(start, [start])]
    paths: list[list[int]] = []
    max_len = len(allowed_atoms) + 1
    while stack and len(paths) < max_paths:
        current, path = stack.pop()
        for nbr in view.neighbors(current):
            if nbr not in allowed_atoms:
                continue
            if nbr == end:
                paths.append(path + [nbr])
                continue
            if nbr in path:
                continue
            if len(path) >= max_len:
                continue
            stack.append((nbr, path + [nbr]))
    unique: dict[tuple[int, ...], list[int]] = {}
    for path in paths:
        unique[tuple(path)] = path
    return sorted(unique.values(), key=len)


def _paths_are_internally_disjoint(paths: tuple[list[int], ...], a: int, b: int) -> bool:
    """Comprueba que los caminos solo compartan los nodos puente."""
    seen: set[int] = set()
    for path in paths:
        internal = set(path[1:-1])
        if a in internal or b in internal:
            return False
        if seen & internal:
            return False
        seen |= internal
    return True


def _aromatic_components(view: MolView) -> list[set[int]]:
    """Calcula componentes conexos aromáticos en el grafo.

    Args:
        view: Vista del grafo molecular.

    Returns:
        Lista de conjuntos de átomos aromáticos conectados.
    """
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
    """Intenta reconocer pireno usando una plantilla de numeración.

    Args:
        view: Vista del grafo molecular.
        opts: Opciones del motor.

    Returns:
        Diccionario con plantilla y mapeo o `None` si no coincide.
    """
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
    """Nombra pireno usando el mapeo de una plantilla.

    Args:
        view: Vista del grafo molecular.
        info: Diccionario con plantilla y mapeo.
        opts: Opciones del motor.

    Returns:
        Nombre IUPAC-lite del pireno sustituido.
    """
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

    subs = [Sub(name, locant) for locant, name in substituents]
    return render_name(subs, "pyrene", always_include_locant=True)


def _name_pyrene(view: MolView, info: dict, opts: NameOptions) -> str:
    """Nombra pireno sin plantilla (ruta alternativa).

    Args:
        view: Vista del grafo molecular.
        info: Metadatos del sistema tipo pireno.
        opts: Opciones del motor.

    Returns:
        Nombre IUPAC-lite del pireno sustituido.
    """
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
    """Detecta sistemas heteroaromáticos fusionados conocidos.

    Args:
        view: Vista del grafo molecular.
        rings: Anillos detectados.

    Returns:
        Diccionario con metadatos del sistema o `None` si no coincide.
    """
    ring_infos: dict[frozenset[int], dict] = {}
    for ring in rings:
        info = classify_aromatic_ring(view, ring)
        if info:
            ring_infos[ring] = info

    def _ordered_hetero_atoms(candidates: set[int]) -> list[int]:
        """Ordena heteroátomos por prioridad O > N > S y luego por id."""
        priorities = {"O": 0, "N": 1, "S": 2, "P": 3, "Si": 4, "B": 5}
        return sorted(candidates, key=lambda atom_id: (priorities.get(view.element(atom_id), 9), atom_id))

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
                    "hetero_atoms": [hetero_atom],
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
                    "hetero_atoms": [hetero_atom],
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
                    "hetero_atoms": [hetero_atom],
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
                    "hetero_atoms": [hetero_atom],
                }
            if kinds == {"benzene", "oxazole"} and len(union) == 9:
                hetero_atoms = _ordered_hetero_atoms(
                    {atom_id for atom_id in union if view.element(atom_id) in {"O", "N", "S"}}
                )
                return {
                    "kind": "benzoxazole",
                    "atoms": union,
                    "fusion_atoms": set(shared),
                    "hetero_atoms": hetero_atoms,
                }
            if kinds == {"benzene", "thiazole"} and len(union) == 9:
                hetero_atoms = _ordered_hetero_atoms(
                    {atom_id for atom_id in union if view.element(atom_id) in {"O", "N", "S"}}
                )
                return {
                    "kind": "benzothiazole",
                    "atoms": union,
                    "fusion_atoms": set(shared),
                    "hetero_atoms": hetero_atoms,
                }
            if kinds == {"benzene", "imidazole"} and len(union) == 9:
                hetero_atoms = _ordered_hetero_atoms(
                    {atom_id for atom_id in union if view.element(atom_id) in {"O", "N", "S"}}
                )
                return {
                    "kind": "benzimidazole",
                    "atoms": union,
                    "fusion_atoms": set(shared),
                    "hetero_atoms": hetero_atoms,
                }
            if kinds == {"benzene", "pyrazole"} and len(union) == 9:
                hetero_atoms = _ordered_hetero_atoms(
                    {atom_id for atom_id in union if view.element(atom_id) in {"O", "N", "S"}}
                )
                return {
                    "kind": "indazole",
                    "atoms": union,
                    "fusion_atoms": set(shared),
                    "hetero_atoms": hetero_atoms,
                }
            if kinds == {"benzene", "pyrazine"} and len(union) == 10:
                hetero_atoms = _ordered_hetero_atoms(
                    {atom_id for atom_id in union if view.element(atom_id) in {"O", "N", "S"}}
                )
                return {
                    "kind": "quinoxaline",
                    "atoms": union,
                    "fusion_atoms": set(shared),
                    "hetero_atoms": hetero_atoms,
                }
            triazole_kinds = {"triazole", "1,2,3-triazole", "1,2,4-triazole"}
            if len(union) == 9 and "benzene" in kinds and any(k in triazole_kinds for k in kinds):
                hetero_atoms = _ordered_hetero_atoms(
                    {atom_id for atom_id in union if view.element(atom_id) == "N"}
                )
                if len(hetero_atoms) == 3:
                    return {
                        "kind": "benzotriazole",
                        "atoms": union,
                        "fusion_atoms": set(shared),
                        "hetero_atoms": hetero_atoms,
                    }
    return None


def _name_fused_hetero(view: MolView, info: dict, opts: NameOptions) -> str:
    """Nombra sistemas heteroaromáticos fusionados simples.

    Args:
        view: Vista del grafo molecular.
        info: Metadatos del sistema fusionado.
        opts: Opciones del motor.

    Returns:
        Nombre IUPAC-lite del sistema fusionado.
    """
    kind = info.get("kind")
    ring_atoms = set(info.get("atoms", set()))
    fusion_atoms = set(info.get("fusion_atoms", set()))
    hetero_atoms = list(info.get("hetero_atoms") or [])
    legacy_hetero = info.get("hetero_atom")
    if legacy_hetero is not None and not hetero_atoms:
        hetero_atoms = [legacy_hetero]
    if not kind or not hetero_atoms:
        raise ChemNameNotSupported("Unsupported fused heteroaromatic")
    hetero_set = set(hetero_atoms)

    order = _perimeter_cycle(view, ring_atoms, fusion_atoms)
    if not order:
        raise ChemNameNotSupported("Fused hetero ordering failed")

    if hetero_set & fusion_atoms:
        raise ChemNameNotSupported("Unsupported fused heteroaromatic")

    substituents: list[tuple[int, str]] = []
    for atom_id in ring_atoms:
        for nbr in view.neighbors(atom_id):
            if nbr in ring_atoms:
                continue
            if view.element(nbr) == "H":
                continue
            if atom_id in fusion_atoms or atom_id in hetero_set:
                raise ChemNameNotSupported("Unsupported hetero substitution")
            name = _simple_substituent_name(view, atom_id, nbr, ring_atoms)
            substituents.append((atom_id, name))
    best_subs: list[Sub] = []
    best_key = None
    for locant_map in _hetero_fused_locant_maps(order, fusion_atoms, hetero_atoms):
        if any(atom_id not in locant_map for atom_id in hetero_atoms):
            continue
        subs = [Sub(name, locant_map[atom_id]) for atom_id, name in substituents]
        locants = sorted(sub.locant for sub in subs)
        hetero_locants = tuple(locant_map[atom_id] for atom_id in hetero_atoms)
        key = (hetero_locants, orientation_key(subs, opts, primary_locants=locants))
        if best_key is None or key < best_key:
            best_key = key
            best_subs = subs

    return render_name(best_subs, kind, always_include_locant=True)


def _hetero_fused_locant_maps(
    order: list[int], fusion_atoms: set[int], hetero_atoms: list[int]
) -> list[dict[int, int]]:
    """Genera mapas de locantes para anillos hetero fusionados.

    Args:
        order: Orden perimetral de átomos.
        fusion_atoms: Átomos de fusión entre anillos.
        hetero_atoms: Heteroátomos en orden de prioridad.

    Returns:
        Lista de mapas `atom_id -> locante`.
    """
    if not hetero_atoms:
        return []
    start_atom = hetero_atoms[0]
    if start_atom not in order:
        return []
    idx0 = order.index(start_atom)
    base = order[idx0:] + order[:idx0]
    maps: list[dict[int, int]] = []
    for direction in (1, -1):
        if direction == 1:
            seq = list(base)
        else:
            rev = list(reversed(base))
            ridx = rev.index(start_atom)
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
    """Obtiene el ciclo perimetral de un sistema fusionado múltiple.

    Args:
        view: Vista del grafo molecular.
        ring_atoms: Conjunto de átomos del sistema.
        fusion_edges: Enlaces internos de fusión que se eliminan del perímetro.

    Returns:
        Lista de átomos en orden perimetral o vacía si falla.
    """
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
    """Genera mapas de locantes para sistemas fusionados.

    Args:
        order: Orden perimetral de átomos.
        fusion_atoms: Átomos de fusión que no reciben locante.

    Returns:
        Lista de mapas `atom_id -> locante`.
    """
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
    """Obtiene el ciclo perimetral de un sistema fusionado simple.

    Args:
        view: Vista del grafo molecular.
        ring_atoms: Conjunto de átomos del sistema.
        fusion_atoms: Par de átomos de fusión.

    Returns:
        Lista de átomos en orden perimetral o vacía si falla.
    """
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
    """Devuelve el ciclo en orden dado una adyacencia 2-regular.

    Args:
        adjacency: Diccionario de vecinos (cada nodo con 2 vecinos).
        start: Nodo inicial.

    Returns:
        Orden del ciclo o lista vacía si la estructura no es válida.
    """
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
    """Genera mapas de locantes válidos para naftaleno.

    Args:
        order: Orden perimetral del naftaleno.
        fusion_atoms: Átomos de fusión.

    Returns:
        Lista de mapas `atom_id -> locante` válidos.
    """
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
    """Devuelve el nombre de un sustituyente simple en anillos fusionados.

    Args:
        view: Vista del grafo molecular.
        ring_atom: Átomo del anillo.
        nbr: Átomo sustituyente.
        ring_set: Conjunto de átomos del sistema.

    Returns:
        Nombre del sustituyente.
    """
    elem = view.element(nbr)
    if elem in HALO_MAP:
        return HALO_MAP[elem]
    if elem == "C":
        return alkyl_substituent_name(view, nbr, ring_set)
    raise ChemNameNotSupported("Unsupported substituent")


def _unsaturations_for_chain(view: MolView, chain: list[int]) -> list[tuple[int, int]]:
    """Detecta insaturaciones dentro de una cadena.

    Args:
        view: Vista del grafo molecular.
        chain: Lista de IDs de la cadena.

    Returns:
        Lista de (orden, locante) para dobles/triples enlaces.
    """
    unsaturations: list[tuple[int, int]] = []
    for idx in range(len(chain) - 1):
        order = view.bond_order_between(chain[idx], chain[idx + 1])
        if order == 1:
            continue
        if order not in {2, 3}:
            raise ChemNameNotSupported("Unsupported bond order")
        unsaturations.append((order, idx + 1))
    return unsaturations


def _stereo_descriptors_for_linear(
    view: MolView,
    chain: list[int],
    opts: NameOptions,
) -> list[str]:
    """Detecta designadores estereoquímicos para una cadena orientada.

    Usa RDKit cuando está disponible y, como respaldo, metadatos estéreo
    almacenados en el grafo.
    """
    locant_by_atom = {atom_id: idx + 1 for idx, atom_id in enumerate(chain)}
    descriptors: list[tuple[int, str]] = []

    if opts.rdkit_isolated and callable(safe_stereo_descriptors_for_chain):
        try:
            isolated = safe_stereo_descriptors_for_chain(view.graph, chain)
            if isolated:
                return list(isolated)
        except Exception:
            descriptors = []

    if Chem is not None:
        try:
            mol, id_map = molgraph_to_rdkit_with_map(view.graph)
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            Chem.FindPotentialStereoBonds(mol)
            idx_to_atom = {rd_idx: atom_id for atom_id, rd_idx in id_map.items()}

            centers = Chem.FindMolChiralCenters(
                mol,
                includeUnassigned=False,
                useLegacyImplementation=False,
            )
            for rd_idx, label in centers:
                atom_id = idx_to_atom.get(rd_idx)
                if atom_id not in locant_by_atom:
                    continue
                if label not in {"R", "S"}:
                    continue
                loc = locant_by_atom[atom_id]
                descriptors.append((loc, f"{loc}{label}"))

            for bond in mol.GetBonds():
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                stereo = bond.GetStereo()
                if stereo not in {Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ}:
                    continue
                a_atom = idx_to_atom.get(bond.GetBeginAtomIdx())
                b_atom = idx_to_atom.get(bond.GetEndAtomIdx())
                if a_atom not in locant_by_atom or b_atom not in locant_by_atom:
                    continue
                loc = min(locant_by_atom[a_atom], locant_by_atom[b_atom])
                label = "E" if stereo == Chem.BondStereo.STEREOE else "Z"
                descriptors.append((loc, f"{loc}{label}"))
        except Exception:
            descriptors = []

    if not descriptors:
        descriptors = _stereo_descriptors_from_annotations(view, locant_by_atom)

    descriptors.sort(key=lambda item: (item[0], item[1]))
    return [text for _loc, text in descriptors]


def _stereo_descriptors_from_annotations(
    view: MolView,
    locant_by_atom: dict[int, int],
) -> list[tuple[int, str]]:
    """Recupera descriptores R/S o E/Z almacenados como metadatos."""
    found: list[tuple[int, str]] = []
    for atom_id, loc in locant_by_atom.items():
        atom = view._get_atom(atom_id)  # noqa: SLF001 - lectura interna controlada.
        if atom is None:
            continue
        label = getattr(atom, "stereo_cip", None)
        if label in {"R", "S"}:
            found.append((loc, f"{loc}{label}"))

    graph = getattr(view, "graph", None)
    find_bond_between = getattr(graph, "find_bond_between", None)
    if callable(find_bond_between):
        atom_ids = list(locant_by_atom.keys())
        for i in range(len(atom_ids)):
            for j in range(i + 1, len(atom_ids)):
                a_atom = atom_ids[i]
                b_atom = atom_ids[j]
                bond = find_bond_between(a_atom, b_atom)
                if bond is None:
                    continue
                label = getattr(bond, "stereo_ez", None)
                if label not in {"E", "Z"}:
                    continue
                loc = min(locant_by_atom[a_atom], locant_by_atom[b_atom])
                found.append((loc, f"{loc}{label}"))
    return found


def _find_functional_group(
    view: MolView, chain: list[int], allow_other_hetero: bool = False
) -> FunctionalSelection | None:
    """Detecta grupos funcionales y selecciona sufijo por prioridad IUPAC.

    Args:
        view: Vista del grafo molecular.
        chain: Lista de IDs de la cadena.
        allow_other_hetero: Permite heteroátomos fuera de la cadena.

    Returns:
        `FunctionalSelection` o `None` si no se detectan grupos.
    """
    chain_set = set(chain)
    chain_index = {atom_id: idx for idx, atom_id in enumerate(chain)}
    occurrences: list[FunctionalOccurrence] = []

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
        single_nitrogen = [
            nbr
            for nbr in outside_neighbors
            if view.element(nbr) == "N" and view.bond_order_between(atom_id, nbr) == 1
        ]
        sulfur_neighbors = [
            nbr
            for nbr in outside_neighbors
            if view.element(nbr) == "S" and view.bond_order_between(atom_id, nbr) == 1
        ]

        acid_oxygen: int | None = None
        if len(carbonyl_oxygen) == 1 and len(single_oxygen) == 1:
            o_single = single_oxygen[0]
            h_total = implicit_h_count(view, o_single) + view.explicit_h(o_single)
            heavy_neighbors = [n for n in view.neighbors(o_single) if view.element(n) != "H"]
            if h_total >= 1 and len(heavy_neighbors) == 1 and len(chain_neighbors) == 1:
                occurrences.append(
                    FunctionalOccurrence(
                        kind="acid",
                        atom_id=atom_id,
                        aux_atom_ids={carbonyl_oxygen[0], o_single},
                        prefix_name="carboxy",
                        suffix_name="oic acid",
                    )
                )
                acid_oxygen = o_single
            elif (
                h_total == 0
                and len(heavy_neighbors) == 1
                and len(chain_neighbors) == 1
                and view.formal_charge(o_single) <= -1
            ):
                occurrences.append(
                    FunctionalOccurrence(
                        kind="carboxylate",
                        atom_id=atom_id,
                        aux_atom_ids={carbonyl_oxygen[0], o_single},
                        prefix_name="carboxylate",
                        suffix_name="oate",
                    )
                )
                acid_oxygen = o_single

        if len(carbonyl_oxygen) == 1 and len(chain_neighbors) == 1:
            for o_single in single_oxygen:
                if o_single == acid_oxygen:
                    continue
                h_total = implicit_h_count(view, o_single) + view.explicit_h(o_single)
                heavy_neighbors = [n for n in view.neighbors(o_single) if view.element(n) != "H"]
                if h_total == 0 and len(heavy_neighbors) == 2:
                    occurrences.append(
                        FunctionalOccurrence(
                            kind="ester",
                            atom_id=atom_id,
                            aux_atom_ids={carbonyl_oxygen[0], o_single},
                            prefix_name="alkoxycarbonyl",
                            suffix_name="oate",
                        )
                    )
                    break

        if len(carbonyl_oxygen) == 1 and len(chain_neighbors) == 1 and single_nitrogen:
            for n_atom in single_nitrogen:
                n_heavy = [n for n in view.neighbors(n_atom) if view.element(n) != "H"]
                if atom_id not in n_heavy:
                    continue
                extras = [n for n in n_heavy if n != atom_id]
                prefix = "amido" if not extras else "acylamido"
                occurrences.append(
                    FunctionalOccurrence(
                        kind="amide",
                        atom_id=atom_id,
                        aux_atom_ids={carbonyl_oxygen[0], n_atom},
                        prefix_name=prefix,
                        suffix_name="amide",
                    )
                )
                break

        for sulfur_atom in sulfur_neighbors:
            sulfo = detect_sulfonic_attachment(view, sulfur_atom, chain_set)
            if sulfo is None:
                continue
            kind = str(sulfo.get("kind"))
            aux = set(sulfo.get("aux_atoms", set())) | {sulfur_atom}
            if kind == "sulfonic_acid":
                occurrences.append(
                    FunctionalOccurrence(
                        kind="sulfonic_acid",
                        atom_id=atom_id,
                        aux_atom_ids=aux,
                        prefix_name="sulfo",
                        suffix_name="sulfonic acid",
                    )
                )
            elif kind == "sulfonate":
                occurrences.append(
                    FunctionalOccurrence(
                        kind="sulfonate",
                        atom_id=atom_id,
                        aux_atom_ids=aux,
                        prefix_name="sulfo",
                        suffix_name="sulfonate",
                    )
                )

        if len(carbonyl_oxygen) == 1 and len(outside_neighbors) == 1:
            if len(chain_neighbors) == 1:
                occurrences.append(
                    FunctionalOccurrence(
                        kind="aldehyde",
                        atom_id=atom_id,
                        aux_atom_ids={carbonyl_oxygen[0]},
                        prefix_name="oxo",
                        suffix_name="al",
                    )
                )
            elif len(chain_neighbors) == 2:
                occurrences.append(
                    FunctionalOccurrence(
                        kind="ketone",
                        atom_id=atom_id,
                        aux_atom_ids={carbonyl_oxygen[0]},
                        prefix_name="oxo",
                        suffix_name="one",
                    )
                )

        if len(nitrile_n) == 1 and len(outside_neighbors) == 1 and len(chain_neighbors) == 1:
            occurrences.append(
                FunctionalOccurrence(
                    kind="nitrile",
                    atom_id=atom_id,
                    aux_atom_ids={nitrile_n[0]},
                    prefix_name="cyano",
                    suffix_name="nitrile",
                )
            )

        for nbr in outside_neighbors:
            if acid_oxygen is not None and nbr == acid_oxygen:
                continue
            elem = view.element(nbr)
            if elem == "O" and view.bond_order_between(atom_id, nbr) == 1:
                if view.has_radical(nbr):
                    continue
                h_total = implicit_h_count(view, nbr) + view.explicit_h(nbr)
                heavy_neighbors = [n for n in view.neighbors(nbr) if view.element(n) != "H"]
                if h_total >= 1 and len(heavy_neighbors) == 1:
                    occurrences.append(
                        FunctionalOccurrence(
                            kind="alcohol",
                            atom_id=atom_id,
                            aux_atom_ids={nbr},
                            prefix_name="hydroxy",
                            suffix_name="ol",
                        )
                    )
            if elem == "N" and view.bond_order_between(atom_id, nbr) == 1:
                h_total = implicit_h_count(view, nbr) + view.explicit_h(nbr)
                heavy_neighbors = [n for n in view.neighbors(nbr) if view.element(n) != "H"]
                if h_total >= 1 and len(heavy_neighbors) == 1:
                    occurrences.append(
                        FunctionalOccurrence(
                            kind="amine",
                            atom_id=atom_id,
                            aux_atom_ids={nbr},
                            prefix_name="amino",
                            suffix_name="amine",
                        )
                    )
            if elem == "S" and view.bond_order_between(atom_id, nbr) == 1:
                h_total = implicit_h_count(view, nbr) + view.explicit_h(nbr)
                heavy_neighbors = [n for n in view.neighbors(nbr) if view.element(n) != "H"]
                if h_total >= 1 and len(heavy_neighbors) == 1:
                    occurrences.append(
                        FunctionalOccurrence(
                            kind="thiol",
                            atom_id=atom_id,
                            aux_atom_ids={nbr},
                            prefix_name="mercapto",
                            suffix_name="thiol",
                        )
                    )

    if not occurrences:
        return None

    # Desduplicar por (tipo, átomo, sufijo/prefijo).
    dedup: dict[tuple[str, int, str, str], FunctionalOccurrence] = {}
    for occ in occurrences:
        key = (occ.kind, occ.atom_id, occ.prefix_name, occ.suffix_name)
        if key not in dedup:
            dedup[key] = occ
    occurrences = list(dedup.values())

    priority = {
        "acid": 0,
        "sulfonic_acid": 0,
        "carboxylate": 1,
        "sulfonate": 1,
        "ester": 2,
        "amide": 3,
        "nitrile": 4,
        "aldehyde": 5,
        "ketone": 6,
        "alcohol": 7,
        "thiol": 7,
        "amine": 8,
    }
    primary = min(
        occurrences,
        key=lambda occ: (priority.get(occ.kind, 99), chain_index.get(occ.atom_id, 999)),
    )

    ignore_atoms: set[int] = set()
    prefixes: list[tuple[str, int]] = []
    for occ in occurrences:
        ignore_atoms |= occ.aux_atom_ids
        if occ is primary:
            continue
        prefixes.append((occ.prefix_name, occ.atom_id))
    prefixes.sort(key=lambda item: (chain_index.get(item[1], 999), item[0]))

    return FunctionalSelection(
        suffix=primary.suffix_name,
        atom_id=primary.atom_id,
        ignore_atoms=ignore_atoms,
        prefixes=prefixes,
    )


def _trim_aromatic_linker(
    view: MolView, chain: list[int], ring_ctx: RingContext
) -> list[int]:
    """Recorta un carbono puente entre cadena y anillo aromático.

    Esto evita elegir como cadena principal un átomo que solo conecta
    la cadena con un anillo aromático y no aporta longitud real.
    """
    if len(chain) <= 1:
        return chain
    chain_set = set(chain)

    def should_trim(atom_id: int) -> bool:
        """Decide si el átomo es un enlace puente hacia un anillo aromático."""
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
    """Lanza error si hay heteroátomos fuera de la cadena permitida."""
    allowed = allowed or set()
    for atom_id in view.atoms():
        if atom_id in chain_set or atom_id in allowed:
            continue
        elem = view.element(atom_id)
        if elem in {"O", "N"}:
            raise ChemNameNotSupported("Unsupported heteroatom")


def _locant_for_atom(chain: list[int], atom_id: int | None) -> int | None:
    """Obtiene el locante (1-based) de un átomo dentro de la cadena."""
    if atom_id is None:
        return None
    for idx, cid in enumerate(chain):
        if cid == atom_id:
            return idx + 1
    return None


def _isotope_substituents_on_chain(view: MolView, chain: list[int]) -> list[Sub]:
    """Genera prefijos isotópicos para átomos isotópicos de la cadena."""
    subs: list[Sub] = []
    for idx, atom_id in enumerate(chain, start=1):
        isotope = view.isotope(atom_id)
        if isotope is None:
            continue
        elem = view.element(atom_id)
        if elem == "H":
            if isotope == 2:
                name = "deuterio"
            elif isotope == 3:
                name = "tritio"
            else:
                name = f"({isotope}H)"
        else:
            name = f"({isotope}{elem})"
        subs.append(Sub(name, idx))
    return subs


def _apply_radical_suffix_if_needed(
    view: MolView,
    chain: list[int],
    substituents: list[Sub],
    rendered: str,
    opts: NameOptions,
) -> str:
    """Aplica degradación segura para radicales no capturados como sustituyente."""
    radical_atoms = [atom_id for atom_id in view.atoms() if view.has_radical(atom_id)]
    if not radical_atoms:
        return rendered
    if any(sub.name == "oxyl" for sub in substituents):
        return rendered
    chain_set = set(chain)
    if all(atom_id in chain_set and view.element(atom_id) == "C" for atom_id in radical_atoms):
        if rendered.endswith("e"):
            return f"{rendered[:-1]}yl"
        return f"{rendered}-yl"
    if opts.strict:
        raise ChemNameNotSupported("Unsupported radical placement")
    return "N/D"


def _choose_oriented_chain(
    view: MolView,
    chain: list[int],
    opts: NameOptions,
    func_atom: int | None,
    ignore_atoms: set[int],
    functional_prefixes: list[tuple[str, int]] | None = None,
    ring_ctx: RingContext | None = None,
) -> list[int]:
    """Selecciona orientación de la cadena minimizando locantes.

    Args:
        view: Vista del grafo molecular.
        chain: Cadena candidata.
        opts: Opciones del motor.
        func_atom: Átomo funcional principal (si aplica).
        ignore_atoms: Átomos a ignorar al calcular sustituyentes.
        functional_prefixes: Prefijos funcionales adicionales (nombre, atom_id).
        ring_ctx: Contexto de anillos.

    Returns:
        Cadena orientada (lista de IDs).
    """
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
    for name, atom_id in functional_prefixes or []:
        f_loc = _locant_for_atom(forward, atom_id)
        if f_loc is not None:
            f_subs.append(Sub(name, f_loc))
        r_loc = _locant_for_atom(reverse, atom_id)
        if r_loc is not None:
            r_subs.append(Sub(name, r_loc))
    f_key = orientation_key(f_subs, opts, primary_locants=f_primary, secondary_locants=f_secondary)
    r_key = orientation_key(r_subs, opts, primary_locants=r_primary, secondary_locants=r_secondary)

    if r_key < f_key:
        return reverse
    return forward
