"""Detección de plantillas especiales (carbohidratos y esteroides)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from chemuson.chemcalc.valence import implicit_h_count
from chemuson.chemio.rdkit_io import molgraph_to_rdkit_with_map
from chemuson.utils.resources import open_resource_path

from .locants import Sub
from .molview import MolView
from .render import render_name
from .template import TemplateMol, load_template

try:
    from rdkit import Chem
except Exception:  # pragma: no cover - dependencia opcional
    Chem = None


@dataclass(frozen=True)
class _SpecialTemplateSpec:
    """Metadatos de una plantilla especial."""

    key: str
    base_name: str
    filename: str
    locants: dict[int, int]


class SpecialMapping(dict):
    """Mapping `template_idx -> atom_id` enriquecido con metadatos."""

    def __init__(
        self,
        *args,
        template_key: str,
        base_name: str,
        locants: dict[int, int],
        **kwargs,
    ) -> None:
        super().__init__(*args, **kwargs)
        self.template_key = str(template_key)
        self.base_name = str(base_name)
        self.locants = {int(k): int(v) for k, v in locants.items()}


_SPECIAL_TEMPLATES: tuple[_SpecialTemplateSpec, ...] = (
    _SpecialTemplateSpec(
        key="alpha_d_glucopyranose",
        base_name="alpha-d-glucopyranose",
        filename="alpha_d_glucopyranose.mol",
        locants={idx: idx + 1 for idx in range(13)},
    ),
    _SpecialTemplateSpec(
        key="beta_d_glucopyranose",
        base_name="beta-d-glucopyranose",
        filename="beta_d_glucopyranose.mol",
        locants={idx: idx + 1 for idx in range(13)},
    ),
    _SpecialTemplateSpec(
        key="beta_d_fructofuranose",
        base_name="beta-d-fructofuranose",
        filename="beta_d_fructofuranose.mol",
        locants={idx: idx + 1 for idx in range(13)},
    ),
    _SpecialTemplateSpec(
        key="d_ribose",
        base_name="d-ribose",
        filename="d_ribose.mol",
        locants={idx: idx + 1 for idx in range(11)},
    ),
    _SpecialTemplateSpec(
        key="androstane",
        base_name="androstane",
        filename="androstane_core.mol",
        locants={idx: idx + 1 for idx in range(17)},
    ),
    _SpecialTemplateSpec(
        key="cholestane",
        base_name="cholestane",
        filename="cholestane_core.mol",
        locants={idx: idx + 1 for idx in range(25)},
    ),
)

_TEMPLATE_CACHE: dict[str, TemplateMol] = {}
_TEMPLATE_MOLBLOCK_CACHE: dict[str, str] = {}


def detect_special_template(view: MolView) -> tuple[str, SpecialMapping] | None:
    """Devuelve `(name, mapping)` si se detecta una plantilla especial."""
    loaded: dict[str, TemplateMol] = {}

    def _template_for(spec: _SpecialTemplateSpec) -> TemplateMol:
        template = loaded.get(spec.key)
        if template is None:
            template = _load_special_template(spec)
            loaded[spec.key] = template
        return template

    # Paso 1: priorizar matches exactos por plantilla completa.
    for spec in _SPECIAL_TEMPLATES:
        mapping = _match_template_exact_component(view, _template_for(spec), spec)
        if mapping is not None:
            return spec.base_name, mapping

    # Paso 2: para subestructura, priorizar núcleos más grandes primero.
    ordered_specs = sorted(
        _SPECIAL_TEMPLATES,
        key=lambda spec: len(_template_for(spec).atoms),
        reverse=True,
    )
    for spec in ordered_specs:
        mapping = _match_template_core(view, _template_for(spec), spec)
        if mapping is None:
            mapping = _match_template_core_rdkit(view, spec)
        if mapping is not None:
            return spec.base_name, mapping
    return None


def detect_special_templates(
    view: MolView,
) -> tuple[str, SpecialMapping, dict[str, Any]] | None:
    """Compatibilidad: retorna también metadata mínima."""
    detected = detect_special_template(view)
    if detected is None:
        return None
    name, mapping = detected
    return name, mapping, {"template_key": mapping.template_key}


def name_special_substituted(view: MolView, mapping: dict[int, int]) -> str:
    """Renderiza nombre de plantilla con sustitución simple si aplica."""
    special_map = _coerce_special_mapping(mapping)
    base_name = special_map.base_name
    core_atoms = set(special_map.values())
    mol_locants = {
        int(mol_atom): int(special_map.locants.get(t_idx, t_idx + 1))
        for t_idx, mol_atom in special_map.items()
    }

    substituents: list[Sub] = []
    for atom_id in sorted(core_atoms):
        locant = mol_locants.get(atom_id)
        if locant is None:
            continue
        for nbr in view.neighbors(atom_id):
            if nbr in core_atoms or view.element(nbr) == "H":
                continue
            name = _special_substituent_name(view, atom_id, nbr, core_atoms)
            if name is None:
                return f"{base_name} (substituted)"
            substituents.append(Sub(name, locant))

    if not substituents:
        return base_name
    return render_name(substituents, base_name, always_include_locant=True)


def _coerce_special_mapping(mapping: dict[int, int]) -> SpecialMapping:
    """Normaliza el mapping a `SpecialMapping`."""
    if isinstance(mapping, SpecialMapping):
        return mapping
    return SpecialMapping(
        {int(k): int(v) for k, v in mapping.items()},
        template_key="special",
        base_name="special",
        locants={int(k): int(k) + 1 for k in mapping.keys()},
    )


def _load_special_template(spec: _SpecialTemplateSpec) -> TemplateMol:
    """Carga plantilla especial con caché."""
    if spec.key in _TEMPLATE_CACHE:
        return _TEMPLATE_CACHE[spec.key]
    with open_resource_path("chemname", "templates", "special", spec.filename) as path:
        template = load_template(path)
    _TEMPLATE_CACHE[spec.key] = template
    return template


def _load_special_molblock(spec: _SpecialTemplateSpec) -> str:
    """Carga el contenido MOL de la plantilla para ruta RDKit."""
    if spec.key in _TEMPLATE_MOLBLOCK_CACHE:
        return _TEMPLATE_MOLBLOCK_CACHE[spec.key]
    with open_resource_path("chemname", "templates", "special", spec.filename) as path:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    _TEMPLATE_MOLBLOCK_CACHE[spec.key] = text
    return text


def _match_template_exact_component(
    view: MolView,
    template: TemplateMol,
    spec: _SpecialTemplateSpec,
) -> SpecialMapping | None:
    """Match exacto de la plantilla contra un componente del grafo."""
    for component in _connected_components(view):
        if len(component) != len(template.atoms):
            continue
        mappings = _match_template_with_backtracking(
            view=view,
            template=template,
            candidate_atoms=component,
            allow_extra_degree=False,
        )
        if not mappings:
            continue
        best = _select_best_mapping(view, mappings, allow_externals=False)
        if best is not None:
            return SpecialMapping(
                best,
                template_key=spec.key,
                base_name=spec.base_name,
                locants=spec.locants,
            )
    return None


def _match_template_core(
    view: MolView,
    template: TemplateMol,
    spec: _SpecialTemplateSpec,
) -> SpecialMapping | None:
    """Match de núcleo con sustituyentes permitidos fuera del mapping."""
    atom_ids = set(view.atoms())
    if len(atom_ids) < len(template.atoms):
        return None
    mappings = _match_template_with_backtracking(
        view=view,
        template=template,
        candidate_atoms=atom_ids,
        allow_extra_degree=True,
    )
    if not mappings:
        return None
    best = _select_best_mapping(view, mappings, allow_externals=True)
    if best is None:
        return None
    return SpecialMapping(
        best,
        template_key=spec.key,
        base_name=spec.base_name,
        locants=spec.locants,
    )


def _match_template_core_rdkit(
    view: MolView,
    spec: _SpecialTemplateSpec,
) -> SpecialMapping | None:
    """Ruta opcional RDKit para subestructura si está disponible."""
    if Chem is None:
        return None
    try:
        mol, id_map = molgraph_to_rdkit_with_map(view.graph)
    except Exception:
        return None
    try:
        template_mol = Chem.MolFromMolBlock(_load_special_molblock(spec), sanitize=False)
    except Exception:
        template_mol = None
    if template_mol is None:
        return None

    try:
        matches = mol.GetSubstructMatches(template_mol, uniquify=True, useChirality=False)
    except Exception:
        return None
    if not matches:
        return None

    idx_to_atom = {rd_idx: atom_id for atom_id, rd_idx in id_map.items()}
    rd_mapping = matches[0]
    mapped = {
        int(t_idx): int(idx_to_atom[int(rd_idx)])
        for t_idx, rd_idx in enumerate(rd_mapping)
        if int(rd_idx) in idx_to_atom
    }
    if len(mapped) != len(rd_mapping):
        return None
    return SpecialMapping(
        mapped,
        template_key=spec.key,
        base_name=spec.base_name,
        locants=spec.locants,
    )


def _match_template_with_backtracking(
    view: MolView,
    template: TemplateMol,
    candidate_atoms: set[int],
    allow_extra_degree: bool,
) -> list[dict[int, int]]:
    """Busca mapeos isomórficos por backtracking."""
    if len(candidate_atoms) < len(template.atoms):
        return []

    mol_degree = {
        atom_id: len([nbr for nbr in view.neighbors(atom_id) if nbr in candidate_atoms])
        for atom_id in candidate_atoms
    }
    mol_aromatic = _aromatic_atoms(view, candidate_atoms)
    mol_bonds = _bond_table(view, candidate_atoms)
    template_bonds = _template_bond_table(template)

    order = sorted(
        range(len(template.atoms)),
        key=lambda idx: (-template.atoms[idx].degree, template.atoms[idx].element, idx),
    )
    candidates: dict[int, list[int]] = {}
    for t_idx, t_atom in enumerate(template.atoms):
        options: list[int] = []
        for atom_id in candidate_atoms:
            if view.element(atom_id) != t_atom.element:
                continue
            if t_atom.aromatic and not mol_aromatic.get(atom_id, False):
                continue
            if allow_extra_degree:
                if mol_degree.get(atom_id, 0) < t_atom.degree:
                    continue
            else:
                if mol_degree.get(atom_id, 0) != t_atom.degree:
                    continue
            options.append(atom_id)
        if not options:
            return []
        candidates[t_idx] = options

    mapping: dict[int, int] = {}
    used: set[int] = set()
    results: list[dict[int, int]] = []

    def backtrack(pos: int) -> None:
        if pos >= len(order):
            if _mapping_respects_template(mapping, template_bonds, mol_bonds):
                results.append(dict(mapping))
            return
        t_idx = order[pos]
        for m_atom in candidates[t_idx]:
            if m_atom in used:
                continue
            if not _partial_mapping_is_compatible(
                t_idx,
                m_atom,
                mapping,
                template_bonds,
                mol_bonds,
            ):
                continue
            mapping[t_idx] = m_atom
            used.add(m_atom)
            backtrack(pos + 1)
            used.remove(m_atom)
            mapping.pop(t_idx, None)

    backtrack(0)
    return results


def _template_bond_table(template: TemplateMol) -> dict[tuple[int, int], tuple[int, bool]]:
    """Convierte enlaces de plantilla a tabla rápida."""
    table: dict[tuple[int, int], tuple[int, bool]] = {}
    for bond in template.bonds:
        key = (bond.a, bond.b) if bond.a < bond.b else (bond.b, bond.a)
        table[key] = (int(bond.order), bool(bond.aromatic))
    return table


def _bond_table(view: MolView, atom_ids: set[int]) -> dict[tuple[int, int], tuple[int, bool]]:
    """Tabla de enlaces `(a,b)->(order, aromatic)` sobre un subconjunto."""
    table: dict[tuple[int, int], tuple[int, bool]] = {}
    for atom_id in atom_ids:
        for nbr in view.neighbors(atom_id):
            if nbr not in atom_ids:
                continue
            key = (atom_id, nbr) if atom_id < nbr else (nbr, atom_id)
            if key in table:
                continue
            table[key] = (
                int(view.bond_order_between(atom_id, nbr)),
                bool(view.bond_is_aromatic(atom_id, nbr)),
            )
    return table


def _aromatic_atoms(view: MolView, atom_ids: set[int]) -> dict[int, bool]:
    """Marca átomos con al menos un enlace aromático."""
    aromatic = {atom_id: False for atom_id in atom_ids}
    for atom_id in atom_ids:
        for nbr in view.neighbors(atom_id):
            if nbr not in atom_ids:
                continue
            if view.bond_is_aromatic(atom_id, nbr):
                aromatic[atom_id] = True
                break
    return aromatic


def _partial_mapping_is_compatible(
    template_idx: int,
    mol_atom: int,
    mapping: dict[int, int],
    template_bonds: dict[tuple[int, int], tuple[int, bool]],
    mol_bonds: dict[tuple[int, int], tuple[int, bool]],
) -> bool:
    """Valida un nuevo par parcial `template_idx -> mol_atom`."""
    for t_other, m_other in mapping.items():
        t_key = (template_idx, t_other) if template_idx < t_other else (t_other, template_idx)
        m_key = (mol_atom, m_other) if mol_atom < m_other else (m_other, mol_atom)
        t_bond = template_bonds.get(t_key)
        m_bond = mol_bonds.get(m_key)
        if t_bond is None:
            if m_bond is not None:
                return False
            continue
        if m_bond is None:
            return False
        t_order, t_aromatic = t_bond
        m_order, m_aromatic = m_bond
        if t_aromatic:
            if not (m_aromatic or m_order in {1, 2}):
                return False
        elif int(m_order) != int(t_order):
            return False
    return True


def _mapping_respects_template(
    mapping: dict[int, int],
    template_bonds: dict[tuple[int, int], tuple[int, bool]],
    mol_bonds: dict[tuple[int, int], tuple[int, bool]],
) -> bool:
    """Comprueba que el mapeo final respete enlaces y no-enlaces."""
    t_indices = sorted(mapping.keys())
    for i in range(len(t_indices)):
        for j in range(i + 1, len(t_indices)):
            t_a = t_indices[i]
            t_b = t_indices[j]
            t_key = (t_a, t_b)
            m_a = mapping[t_a]
            m_b = mapping[t_b]
            m_key = (m_a, m_b) if m_a < m_b else (m_b, m_a)
            t_bond = template_bonds.get(t_key)
            m_bond = mol_bonds.get(m_key)
            if t_bond is None:
                if m_bond is not None:
                    return False
                continue
            if m_bond is None:
                return False
            t_order, t_aromatic = t_bond
            m_order, m_aromatic = m_bond
            if t_aromatic:
                if not (m_aromatic or m_order in {1, 2}):
                    return False
            elif int(m_order) != int(t_order):
                return False
    return True


def _select_best_mapping(
    view: MolView,
    mappings: list[dict[int, int]],
    allow_externals: bool,
) -> dict[int, int] | None:
    """Elige mapeo por menor sustitución externa y orden estable."""
    if not mappings:
        return None
    best = None
    best_key = None
    for mapping in mappings:
        mapped_atoms = set(mapping.values())
        external_count = 0
        if allow_externals:
            for atom_id in mapped_atoms:
                for nbr in view.neighbors(atom_id):
                    if nbr in mapped_atoms or view.element(nbr) == "H":
                        continue
                    external_count += 1
        map_key = tuple(mapping[idx] for idx in sorted(mapping.keys()))
        key = (external_count, map_key)
        if best_key is None or key < best_key:
            best_key = key
            best = mapping
    return best


def _connected_components(view: MolView) -> list[set[int]]:
    """Componentes conexos del grafo completo."""
    remaining = set(view.atoms())
    components: list[set[int]] = []
    while remaining:
        start = next(iter(remaining))
        stack = [start]
        component: set[int] = set()
        while stack:
            atom_id = stack.pop()
            if atom_id in component:
                continue
            component.add(atom_id)
            for nbr in view.neighbors(atom_id):
                if nbr not in component:
                    stack.append(nbr)
        remaining -= component
        components.append(component)
    return components


def _special_substituent_name(
    view: MolView,
    parent_atom: int,
    substituent_atom: int,
    core_atoms: set[int],
) -> str | None:
    """Detecta sustituyentes simples permitidos sobre plantillas especiales."""
    elem = view.element(substituent_atom)
    order = view.bond_order_between(parent_atom, substituent_atom)
    heavy_neighbors = [nbr for nbr in view.neighbors(substituent_atom) if view.element(nbr) != "H"]
    external_heavy = [nbr for nbr in heavy_neighbors if nbr not in core_atoms and nbr != parent_atom]

    if elem == "O":
        if order == 2 and not external_heavy:
            return "oxo"
        if order == 1 and not external_heavy:
            h_total = implicit_h_count(view, substituent_atom) + view.explicit_h(substituent_atom)
            if h_total >= 1:
                return "hydroxy"
    if elem == "N" and order == 1 and not external_heavy:
        h_total = implicit_h_count(view, substituent_atom) + view.explicit_h(substituent_atom)
        if h_total >= 1:
            return "amino"
    return None
