"""Conversión entre el grafo interno y formatos químicos externos.

Este módulo integra RDKit (si está disponible) para exportar/importar
SMILES, MOL y SVG. Incluye rutas de respaldo para exportar sin RDKit.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple

from chemuson.core.model import ATOMIC_NUMBERS, BondStyle, BondStereo, MolGraph

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except Exception:  # pragma: no cover - optional dependency at runtime
        rdMolStandardize = None
except Exception:  # pragma: no cover - optional dependency at runtime
    Chem = None
    AllChem = None
    rdMolDraw2D = None
    rdMolStandardize = None


def _rdkit_available() -> bool:
    """Indica si RDKit está disponible en tiempo de ejecución."""
    return Chem is not None and AllChem is not None and rdMolDraw2D is not None


def _require_rdkit():
    """Lanza un error si RDKit no está instalado."""
    if not _rdkit_available():
        raise RuntimeError("RDKit no disponible")


@dataclass
class StrictValidationResult:
    """Resultado de validación/normalización estricta basada en RDKit."""

    is_valid: bool
    errors: list[str]
    normalized_smiles: Optional[str] = None
    normalized_graph: Optional[MolGraph] = None


def molgraph_to_rdkit_with_map(molgraph: MolGraph):
    """Convierte un `MolGraph` a RDKit `Mol` y mapea IDs internos.

    Args:
        molgraph: Grafo molecular de Chemuson.

    Returns:
        Tupla `(mol, id_map)` donde `mol` es un RDKit Mol y `id_map` mapea
        IDs de átomos del grafo a índices RDKit.

    Raises:
        RuntimeError: Si RDKit no está disponible.

    Side Effects:
        No modifica el grafo; construye objetos RDKit en memoria.
    """
    _require_rdkit()
    rw = Chem.RWMol()
    id_map: Dict[int, int] = {}

    for atom in sorted(molgraph.atoms.values(), key=lambda a: a.id):
        element = atom.element
        try:
            atomic_number = Chem.GetPeriodicTable().GetAtomicNumber(element)
        except Exception:
            atomic_number = 0
        if atomic_number <= 0:
            rd_atom = Chem.Atom(0)
            rd_atom.SetProp("atomLabel", element)
            rd_atom.SetProp("dummyLabel", element)
        else:
            rd_atom = Chem.Atom(element)
        rd_atom.SetFormalCharge(atom.charge)
        if atom.isotope is not None:
            rd_atom.SetIsotope(atom.isotope)
        rd_idx = rw.AddAtom(rd_atom)
        id_map[atom.id] = rd_idx

    for bond in molgraph.bonds.values():
        # RDKit distingue enlaces aromáticos y órdenes discretos.
        begin_id = bond.a1_id
        end_id = bond.a2_id
        if bond.style == BondStyle.COORDINATION:
            bond_type = Chem.BondType.DATIVE
            donor = bond.donor_atom_id
            if donor == bond.a2_id:
                begin_id, end_id = bond.a2_id, bond.a1_id
            elif donor == bond.a1_id:
                begin_id, end_id = bond.a1_id, bond.a2_id
        elif bond.is_aromatic:
            rw.GetAtomWithIdx(id_map[bond.a1_id]).SetIsAromatic(True)
            rw.GetAtomWithIdx(id_map[bond.a2_id]).SetIsAromatic(True)
            bond_type = Chem.BondType.AROMATIC
        elif bond.order == 2:
            bond_type = Chem.BondType.DOUBLE
        elif bond.order == 3:
            bond_type = Chem.BondType.TRIPLE
        else:
            bond_type = Chem.BondType.SINGLE
        rw.AddBond(id_map[begin_id], id_map[end_id], bond_type)

    mol = rw.GetMol()
    # Preservar coordenadas 2D del editor.
    conf = Chem.Conformer(mol.GetNumAtoms())
    for atom_id, idx in id_map.items():
        atom = molgraph.atoms[atom_id]
        conf.SetAtomPosition(idx, (atom.x, atom.y, 0.0))
    mol.AddConformer(conf, assignId=True)
    return mol, id_map


def molgraph_to_rdkit(molgraph: MolGraph):
    """Convierte un `MolGraph` a RDKit `Mol` ignorando el mapeo de IDs.

    Args:
        molgraph: Grafo molecular de Chemuson.

    Returns:
        Mol de RDKit con conformación 2D.
    """
    mol, _ = molgraph_to_rdkit_with_map(molgraph)
    return mol


def molgraph_to_smiles(molgraph: MolGraph) -> str:
    """Genera SMILES desde un `MolGraph`.

    Usa RDKit si está disponible; si falla, utiliza un escritor interno.

    Args:
        molgraph: Grafo molecular de Chemuson.

    Returns:
        SMILES canónico o aproximado en el camino de respaldo.
    """
    if _graph_requires_rdkit_fallback(molgraph):
        return _molgraph_to_smiles_fallback(molgraph)
    if _rdkit_available():
        try:
            mol = molgraph_to_rdkit(molgraph)
            return Chem.MolToSmiles(mol, canonical=True)
        except Exception:
            # RDKit can reject some hypervalent depictions (e.g., interhalogens, noble gases).
            # Fall back to the internal writer so the editor can still export something useful.
            return _molgraph_to_smiles_fallback(molgraph)
    return _molgraph_to_smiles_fallback(molgraph)


def molgraph_to_molfile(molgraph: MolGraph) -> str:
    """Genera un bloque MOL (V2000) desde un `MolGraph`.

    Args:
        molgraph: Grafo molecular de Chemuson.

    Returns:
        Cadena MOL. Si RDKit falla, usa un exportador interno básico.
    """
    if _graph_requires_rdkit_fallback(molgraph):
        return _molgraph_to_molfile_fallback(molgraph)
    if _rdkit_available():
        try:
            mol = molgraph_to_rdkit(molgraph)
            return Chem.MolToMolBlock(mol)
        except Exception:
            return _molgraph_to_molfile_fallback(molgraph)
    return _molgraph_to_molfile_fallback(molgraph)


def molgraph_to_svg(molgraph: MolGraph, size: Tuple[int, int] = (300, 200)) -> str:
    """Renderiza un `MolGraph` como SVG con RDKit.

    Args:
        molgraph: Grafo molecular de Chemuson.
        size: Tamaño del lienzo (ancho, alto) en píxeles.

    Returns:
        SVG como cadena.
    """
    _require_rdkit()
    mol = molgraph_to_rdkit(molgraph)
    drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _looks_like_counts_line(line: str) -> bool:
    """Heurística para identificar la línea de conteo de un MOL V2000/V3000."""
    text = str(line).rstrip()
    if not text:
        return False
    if "V2000" in text or "V3000" in text:
        return True
    padded = f"{text:<6}"[:6]
    left = padded[:3].strip()
    right = padded[3:6].strip()
    return left.isdigit() and right.isdigit()


def _parse_counts(line: str) -> Optional[tuple[int, int]]:
    """Extrae número de átomos/enlaces desde la línea de conteo."""
    text = str(line or "")
    if not text:
        return None
    try:
        return int(text[0:3]), int(text[3:6])
    except Exception:
        pass
    ints = re.findall(r"-?\d+", text)
    if len(ints) < 2:
        return None
    try:
        return int(ints[0]), int(ints[1])
    except Exception:
        return None


def _find_counts_line_index(lines: list[str]) -> Optional[int]:
    """Localiza la línea de conteo CTAB en las primeras líneas del bloque."""
    limit = min(len(lines), 8)
    for idx in range(limit):
        if _looks_like_counts_line(lines[idx]):
            if _parse_counts(lines[idx]) is not None:
                return idx
    return None


def _is_float_token(token: str) -> bool:
    """Verifica si un token representa un número real."""
    try:
        float(token)
        return True
    except Exception:
        return False


def _parse_atom_line(line: str) -> tuple[float, float, str]:
    """Parsea coordenadas XY y símbolo desde una línea de átomo MOL."""
    tokens = line.split()
    if (
        len(tokens) >= 4
        and _is_float_token(tokens[0])
        and _is_float_token(tokens[1])
        and _is_float_token(tokens[2])
    ):
        x = float(tokens[0])
        y = float(tokens[1])
        symbol = tokens[3].strip() or "C"
        return x, y, symbol
    try:
        x = float(line[0:10])
        y = float(line[10:20])
    except Exception as exc:
        raise ValueError(f"Línea de átomo inválida: {line!r}") from exc
    symbol = line[31:34].strip() if len(line) >= 34 else ""
    if not symbol:
        symbol = "C"
    return x, y, symbol


def _parse_bond_line(line: str) -> tuple[int, int, int]:
    """Parsea índices de átomos y tipo de enlace desde una línea MOL."""
    try:
        return int(line[0:3]), int(line[3:6]), int(line[6:9])
    except Exception:
        tokens = line.split()
        if len(tokens) < 3:
            raise ValueError(f"Línea de enlace inválida: {line!r}")
        return int(tokens[0]), int(tokens[1]), int(tokens[2])


def _mol_symbol_requires_fallback(symbol: str) -> bool:
    """Indica si un símbolo no es elemento estándar y conviene parseo local."""
    clean = str(symbol or "").strip()
    if not clean:
        return True
    if clean in ATOMIC_NUMBERS:
        return False
    if clean in {"*", "A", "Q", "L", "LP", "Du", "R", "R#"}:
        return False
    return True


def _graph_requires_rdkit_fallback(molgraph: MolGraph) -> bool:
    """Determina si conviene evitar RDKit por símbolos no estándar."""
    for atom in molgraph.atoms.values():
        if _mol_symbol_requires_fallback(atom.element):
            return True
    return False


def _should_use_molfile_fallback(molfile: str) -> bool:
    """Detecta MOL con pseudoátomos/estructura no compatible con RDKit parser."""
    text = str(molfile or "")
    if not text:
        return True
    lines = text.splitlines()
    counts_idx = _find_counts_line_index(lines)
    if counts_idx is None:
        return True
    counts = _parse_counts(lines[counts_idx])
    if counts is None:
        return True
    atom_count, _ = counts
    atom_start = counts_idx + 1
    if atom_count < 0 or atom_start + atom_count > len(lines):
        return True
    for i in range(atom_count):
        _, _, symbol = _parse_atom_line(lines[atom_start + i])
        if _mol_symbol_requires_fallback(symbol):
            return True
    # Molblocks con alias `A  n` pueden venir de pseudoátomos serializados
    # como `*` + etiqueta. RDKit tiende a devolver `*` y se pierde la etiqueta;
    # forzamos parser interno para preservar el alias original.
    for line in lines:
        if line.startswith("A  "):
            return True
    return False


def _molfile_to_molgraph_fallback(molfile: str) -> MolGraph:
    """Importa MOL mediante parser interno tolerante (sin RDKit)."""
    text = str(molfile or "")
    if not text.strip():
        raise ValueError("Mol inválido")
    lines = text.splitlines()
    counts_idx = _find_counts_line_index(lines)
    if counts_idx is None:
        raise ValueError("Mol inválido")
    counts = _parse_counts(lines[counts_idx])
    if counts is None:
        raise ValueError("Mol inválido")
    atom_count, bond_count = counts
    if atom_count < 0 or bond_count < 0:
        raise ValueError("Mol inválido")

    atom_start = counts_idx + 1
    bond_start = atom_start + atom_count
    if bond_start + bond_count > len(lines):
        raise ValueError("Mol inválido")

    graph = MolGraph()
    atom_ids: list[int] = []
    for offset in range(atom_count):
        x, y, symbol = _parse_atom_line(lines[atom_start + offset])
        atom = graph.add_atom(symbol, x, y, is_explicit=(symbol != "C"))
        atom_ids.append(atom.id)

    for offset in range(bond_count):
        a1_idx, a2_idx, bond_type = _parse_bond_line(lines[bond_start + offset])
        if not (1 <= a1_idx <= len(atom_ids) and 1 <= a2_idx <= len(atom_ids)):
            continue
        order = 1
        aromatic = False
        if bond_type == 2:
            order = 2
        elif bond_type == 3:
            order = 3
        elif bond_type == 4:
            aromatic = True
            order = 1
        graph.add_bond(
            atom_ids[a1_idx - 1],
            atom_ids[a2_idx - 1],
            order=order,
            style=BondStyle.PLAIN,
            stereo=BondStereo.NONE,
            is_aromatic=aromatic,
        )

    aliases: dict[int, str] = {}
    idx = bond_start + bond_count
    while idx < len(lines):
        line = lines[idx]
        if line.startswith("M  END"):
            break
        if line.startswith("A  "):
            try:
                atom_idx = int(line[3:].strip())
            except Exception:
                parts = line.split()
                atom_idx = int(parts[1]) if len(parts) > 1 else -1
            if idx + 1 < len(lines):
                aliases[atom_idx] = lines[idx + 1].strip()
                idx += 2
                continue
        if line.startswith("M  CHG"):
            tokens = line.split()
            if len(tokens) >= 3:
                try:
                    count = int(tokens[2])
                except Exception:
                    count = 0
                cursor = 3
                for _ in range(count):
                    if cursor + 1 >= len(tokens):
                        break
                    try:
                        atom_idx = int(tokens[cursor])
                        charge = int(tokens[cursor + 1])
                    except Exception:
                        break
                    cursor += 2
                    if 1 <= atom_idx <= len(atom_ids):
                        graph.get_atom(atom_ids[atom_idx - 1]).charge = charge
        idx += 1

    for atom_idx, label in aliases.items():
        if 1 <= atom_idx <= len(atom_ids):
            atom = graph.get_atom(atom_ids[atom_idx - 1])
            atom.element = label
            atom.is_explicit = label != "C"

    _scale_to_default(graph)
    return graph


def normalize_molblock_header(molfile: str) -> str:
    """Normaliza cabecera CTAB para tolerar MOL con línea de comentario omitida.

    Algunos bloques guardados en bibliotecas pueden perder una línea de
    encabezado al aplicar `strip()`, desplazando la línea de conteo de CTAB.
    Esta función inserta una línea vacía cuando detecta ese patrón.
    """
    text = str(molfile or "")
    if not text:
        return ""
    normalized = text.replace("\r\n", "\n").replace("\r", "\n")
    lines = normalized.split("\n")
    if len(lines) >= 3 and _looks_like_counts_line(lines[2]):
        line4 = lines[3] if len(lines) > 3 else ""
        if not _looks_like_counts_line(line4):
            lines.insert(2, "")
    return "\n".join(lines)


def molfile_to_molgraph(molfile: str) -> MolGraph:
    """Importa un bloque MOL (V2000) a `MolGraph`.

    Args:
        molfile: Texto MOL.

    Returns:
        Grafo molecular equivalente.
    """
    normalized = normalize_molblock_header(molfile)
    if _should_use_molfile_fallback(normalized):
        return _molfile_to_molgraph_fallback(normalized)
    if not _rdkit_available():
        return _molfile_to_molgraph_fallback(normalized)
    try:
        mol = Chem.MolFromMolBlock(normalized, sanitize=True)
        if mol is not None:
            return rdkit_to_molgraph(mol)
    except Exception:
        pass
    try:
        mol = Chem.MolFromMolBlock(normalized, sanitize=False)
        if mol is not None:
            return rdkit_to_molgraph(mol)
    except Exception:
        pass
    return _molfile_to_molgraph_fallback(normalized)


def smiles_to_molgraph(smiles: str) -> MolGraph:
    """Importa un SMILES a `MolGraph` usando RDKit.

    Args:
        smiles: Cadena SMILES.

    Returns:
        Grafo molecular equivalente.
    """
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    return rdkit_to_molgraph(mol)


def rdkit_to_molgraph(mol) -> MolGraph:
    """Convierte un RDKit Mol a `MolGraph`.

    Args:
        mol: Instancia de RDKit Mol.

    Returns:
        Grafo molecular con coordenadas 2D.

    Raises:
        ValueError: Si el Mol es `None`.
    """
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
        symbol = atom.GetSymbol()
        new_atom = graph.add_atom(
            symbol,
            pos.x,
            pos.y,
            is_explicit=symbol != "C",
        )
        new_atom.charge = atom.GetFormalCharge()
        if atom.GetIsotope():
            new_atom.isotope = atom.GetIsotope()
        idx_map[idx] = new_atom.id

    for bond in mol.GetBonds():
        order = 1
        style = BondStyle.PLAIN
        donor_atom_id = None
        bond_type = bond.GetBondType()
        if bond_type == Chem.BondType.DOUBLE:
            order = 2
        elif bond_type == Chem.BondType.TRIPLE:
            order = 3
        elif bond_type == Chem.BondType.DATIVE:
            style = BondStyle.COORDINATION
            donor_atom_id = idx_map[bond.GetBeginAtomIdx()]
        graph.add_bond(
            idx_map[bond.GetBeginAtomIdx()],
            idx_map[bond.GetEndAtomIdx()],
            order,
            style=style,
            stereo=BondStereo.NONE,
            is_aromatic=bond.GetIsAromatic(),
            donor_atom_id=donor_atom_id,
        )

    _scale_to_default(graph)
    return graph


def strict_validate_and_normalize(molgraph: MolGraph) -> StrictValidationResult:
    """Valida y normaliza de forma estricta un `MolGraph` usando RDKit.

    La validación estricta ejecuta sanitización de RDKit y, cuando está
    disponible, aplica un pipeline de normalización (`Normalize`/`Reionize`).
    Está orientada a chequeos avanzados de especies iónicas y complejos.
    """
    _require_rdkit()
    mol, _ = molgraph_to_rdkit_with_map(molgraph)

    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:
        errors = [f"SanitizeMol: {exc}"]
        try:
            problems = Chem.DetectChemistryProblems(mol)
        except Exception:
            problems = []
        for problem in problems:
            message = ""
            if hasattr(problem, "Message"):
                try:
                    message = problem.Message()
                except Exception:
                    message = str(problem)
            else:
                message = str(problem)
            errors.append(message)
        return StrictValidationResult(is_valid=False, errors=errors)

    normalized = Chem.Mol(mol)
    if rdMolStandardize is not None:
        normalized = rdMolStandardize.Cleanup(normalized)
        normalized = rdMolStandardize.Normalize(normalized)
        normalized = rdMolStandardize.Reionize(normalized)
        Chem.SanitizeMol(normalized)

    smiles = Chem.MolToSmiles(normalized, canonical=True)
    normalized_graph = rdkit_to_molgraph(Chem.Mol(normalized))
    return StrictValidationResult(
        is_valid=True,
        errors=[],
        normalized_smiles=smiles,
        normalized_graph=normalized_graph,
    )


def kekulize_display_orders(
    molgraph: MolGraph, seed_atoms: Optional[Iterable[int]] = None
) -> Optional[Dict[int, int]]:
    """Calcula órdenes de dibujo para enlaces aromáticos.

    Usa la kekulización de RDKit y, opcionalmente, restringe el cálculo
    a los átomos conectados a una lista de semillas.

    Args:
        molgraph: Grafo molecular de Chemuson.
        seed_atoms: IDs de átomos desde los cuales expandir el cálculo.

    Returns:
        Diccionario de `bond_id -> display_order` o `None` si falla.
    """
    if Chem is None:
        return None
    try:
        mol, id_map = molgraph_to_rdkit_with_map(molgraph)
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        return None

    aromatic_atoms: Optional[set[int]] = None
    if seed_atoms is not None:
        # Expandimos las semillas a todo el componente aromático conectado.
        seeds = set(seed_atoms)
        if seeds:
            adjacency: dict[int, list[int]] = {}
            for bond in molgraph.bonds.values():
                if bond.is_aromatic:
                    adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
                    adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)
            aromatic_atoms = set()
            stack = [atom_id for atom_id in seeds if atom_id in adjacency]
            while stack:
                node = stack.pop()
                if node in aromatic_atoms:
                    continue
                aromatic_atoms.add(node)
                for neighbor in adjacency.get(node, []):
                    if neighbor not in aromatic_atoms:
                        stack.append(neighbor)

    display_orders: Dict[int, int] = {}
    for bond in molgraph.bonds.values():
        if not bond.is_aromatic:
            continue
        if aromatic_atoms is not None and bond.a1_id not in aromatic_atoms:
            continue
        rd_bond = mol.GetBondBetweenAtoms(id_map[bond.a1_id], id_map[bond.a2_id])
        if rd_bond is None:
            continue
        bond_type = rd_bond.GetBondType()
        display_orders[bond.id] = 2 if bond_type == Chem.BondType.DOUBLE else 1
    return display_orders


def _scale_to_default(graph: MolGraph, target: float = 40.0) -> None:
    """Escala el grafo para aproximar una longitud de enlace objetivo.

    Args:
        graph: Grafo molecular a escalar.
        target: Longitud promedio deseada (px).

    Side Effects:
        Modifica las coordenadas de los átomos del grafo.
    """
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


@dataclass
class _SmilesNode:
    """Nodo auxiliar para emitir SMILES en el escritor de respaldo."""
    atom_id: int
    symbol: str
    children: list[tuple[str, "_SmilesNode"]]
    ring_closures: list[tuple[str, int]]


def _molgraph_to_smiles_fallback(molgraph: MolGraph) -> str:
    """Emite un SMILES básico sin RDKit.

    Este exportador es deliberadamente simple: maneja aromaticidad básica,
    ramificaciones y cierres de anillo con numeración incremental.

    Args:
        molgraph: Grafo molecular de Chemuson.

    Returns:
        Cadena SMILES aproximada.
    """
    if not molgraph.atoms:
        return ""
    adjacency: dict[int, list[tuple[int, Bond]]] = {}
    for bond in molgraph.bonds.values():
        adjacency.setdefault(bond.a1_id, []).append((bond.a2_id, bond))
        adjacency.setdefault(bond.a2_id, []).append((bond.a1_id, bond))
    for neighbors in adjacency.values():
        neighbors.sort(key=lambda item: item[0])

    aromatic_atoms: set[int] = set()
    for bond in molgraph.bonds.values():
        if bond.is_aromatic:
            aromatic_atoms.add(bond.a1_id)
            aromatic_atoms.add(bond.a2_id)

    nodes_by_id: dict[int, _SmilesNode] = {}
    edge_handled: set[tuple[int, int]] = set()
    ring_counter = 1

    def bond_symbol(bond: Bond) -> str:
        """Convierte un enlace en su símbolo SMILES."""
        if bond.is_aromatic:
            return ":"
        if bond.order == 2:
            return "="
        if bond.order == 3:
            return "#"
        return ""

    def atom_symbol(atom_id: int) -> str:
        """Construye el símbolo SMILES del átomo con carga/isótopo."""
        atom = molgraph.atoms[atom_id]
        element = atom.element
        aromatic = atom_id in aromatic_atoms and element in {"B", "C", "N", "O", "P", "S"}
        symbol = element.lower() if aromatic else element
        charge = ""
        if atom.charge > 0:
            charge = f"+{atom.charge}" if atom.charge > 1 else "+"
        elif atom.charge < 0:
            magnitude = abs(atom.charge)
            charge = f"-{magnitude}" if magnitude > 1 else "-"
        isotope = f"{atom.isotope}" if atom.isotope is not None else ""
        explicit_h = ""
        if atom.explicit_h is not None and atom.explicit_h > 0:
            explicit_h = "H" if atom.explicit_h == 1 else f"H{atom.explicit_h}"
        if atom.mapping is not None:
            return f"[{isotope}{symbol}{explicit_h}{charge}:{atom.mapping}]"
        if isotope or explicit_h or charge or symbol not in {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"}:
            return f"[{isotope}{symbol}{explicit_h}{charge}]"
        return symbol

    def build(node_id: int, parent_id: Optional[int]) -> _SmilesNode:
        """Construye el árbol DFS y registra cierres de anillo."""
        nonlocal ring_counter
        node = nodes_by_id.get(node_id)
        if node is None:
            node = _SmilesNode(
                atom_id=node_id,
                symbol=atom_symbol(node_id),
                children=[],
                ring_closures=[],
            )
            nodes_by_id[node_id] = node
        for neighbor_id, bond in adjacency.get(node_id, []):
            if neighbor_id == parent_id:
                continue
            edge_key = (min(node_id, neighbor_id), max(node_id, neighbor_id))
            if neighbor_id not in nodes_by_id:
                child = build(neighbor_id, node_id)
                node.children.append((bond_symbol(bond), child))
            else:
                # Enlace ya visto: representamos un cierre de anillo.
                if edge_key in edge_handled:
                    continue
                edge_handled.add(edge_key)
                ring_id = ring_counter
                ring_counter += 1
                node.ring_closures.append((bond_symbol(bond), ring_id))
                other = nodes_by_id[neighbor_id]
                other.ring_closures.append((bond_symbol(bond), ring_id))
        return node

    def ring_token(symbol: str, ring_id: int) -> str:
        """Convierte ID de anillo en token SMILES, usando % para >=10."""
        if ring_id >= 10:
            return f"{symbol}%{ring_id}"
        return f"{symbol}{ring_id}"

    def emit(node: _SmilesNode) -> str:
        """Emite SMILES recorriendo el árbol con ramificaciones."""
        text = node.symbol
        for symbol, ring_id in node.ring_closures:
            text += ring_token(symbol, ring_id)
        for idx, (symbol, child) in enumerate(node.children):
            branch = emit(child)
            if idx == 0:
                text += f"{symbol}{branch}"
            else:
                text += f"({symbol}{branch})"
        return text

    seen: set[int] = set()
    parts: list[str] = []
    for atom_id in sorted(molgraph.atoms.keys()):
        if atom_id in seen:
            continue
        node = build(atom_id, None)
        for node_id in nodes_by_id:
            seen.add(node_id)
        parts.append(emit(node))
    return ".".join(parts)


def _molgraph_to_molfile_fallback(molgraph: MolGraph) -> str:
    """Serializa un `MolGraph` a un MOL V2000 mínimo sin RDKit.

    Args:
        molgraph: Grafo molecular de Chemuson.

    Returns:
        Cadena con el bloque MOL (V2000) básico.
    """
    atoms = [molgraph.atoms[atom_id] for atom_id in sorted(molgraph.atoms.keys())]
    bonds = [molgraph.bonds[bond_id] for bond_id in sorted(molgraph.bonds.keys())]
    atom_index = {atom.id: idx + 1 for idx, atom in enumerate(atoms)}
    aliases: list[tuple[int, str]] = []

    lines = [
        "Chemuson",
        "Chemuson",
        "",
        f"{len(atoms):>3}{len(bonds):>3}  0  0  0  0  0  0  0  0999 V2000",
    ]

    for idx, atom in enumerate(atoms, start=1):
        symbol = atom.element.strip() or "C"
        symbol_for_line = symbol
        if _mol_symbol_requires_fallback(symbol) or len(symbol) > 3:
            symbol_for_line = "*"
            aliases.append((idx, symbol))
        lines.append(
            f"{atom.x:>10.4f}{atom.y:>10.4f}{0.0:>10.4f} {symbol_for_line:<3} 0  0  0  0  0  0  0  0  0  0  0  0"
        )

    for bond in bonds:
        bond_type = 4 if bond.is_aromatic else max(1, min(3, bond.order))
        lines.append(
            f"{atom_index[bond.a1_id]:>3}{atom_index[bond.a2_id]:>3}{bond_type:>3}  0  0  0  0"
        )

    charges: list[tuple[int, int]] = []
    for atom in atoms:
        if atom.charge:
            charges.append((atom_index[atom.id], atom.charge))
    for i in range(0, len(charges), 8):
        chunk = charges[i : i + 8]
        parts = [f"{len(chunk):>3}"]
        for idx, charge in chunk:
            parts.append(f"{idx:>3}{charge:>3}")
        lines.append(f"M  CHG{''.join(parts)}")

    for idx, label in aliases:
        lines.append(f"A  {idx:>3}")
        lines.append(label)

    lines.append("M  END")
    return "\n".join(lines)
