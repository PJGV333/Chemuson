from __future__ import annotations

"""Servicio Name→Structure con fuentes offline y PubChem."""

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Protocol
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

from chemuson.core.model import MolGraph


COMMON_NAME_SMILES: dict[str, str] = {
    "agua": "O",
    "water": "O",
    "metano": "C",
    "methane": "C",
    "etano": "CC",
    "ethane": "CC",
    "etanol": "CCO",
    "ethanol": "CCO",
    "metanol": "CO",
    "methanol": "CO",
    "acetona": "CC(=O)C",
    "acetone": "CC(=O)C",
    "benceno": "c1ccccc1",
    "benzene": "c1ccccc1",
    "tolueno": "Cc1ccccc1",
    "toluene": "Cc1ccccc1",
    "fenol": "Oc1ccccc1",
    "phenol": "Oc1ccccc1",
    "acido acetico": "CC(=O)O",
    "ácido acético": "CC(=O)O",
    "acetic acid": "CC(=O)O",
}


@dataclass(frozen=True)
class NameToStructureResult:
    """Resultado trazable de Name→Structure."""

    query: str
    graph: MolGraph | None
    source: str
    confidence: float
    smiles: str = ""
    resolved_name: str = ""
    message: str = ""
    from_cache: bool = False

    @property
    def ok(self) -> bool:
        return self.graph is not None


class NameConnector(Protocol):
    source: str

    def resolve(self, name: str, timeout_s: float = 8.0) -> NameToStructureResult:
        ...


class StaticNameConnector:
    """Conector offline para nombres comunes frecuentes."""

    source = "offline-common"

    def __init__(self, entries: dict[str, str] | None = None) -> None:
        self._entries = {self._key(key): value for key, value in (entries or COMMON_NAME_SMILES).items()}

    def resolve(self, name: str, timeout_s: float = 8.0) -> NameToStructureResult:
        query = str(name or "").strip()
        smiles = self._entries.get(self._key(query), "")
        if not smiles:
            return NameToStructureResult(query, None, self.source, 0.0, message="not_found")
        graph, error = _smiles_to_graph(smiles, timeout_s=timeout_s)
        if graph is None:
            return NameToStructureResult(query, None, self.source, 0.0, smiles=smiles, message=error)
        return NameToStructureResult(
            query=query,
            graph=graph,
            source=self.source,
            confidence=0.72,
            smiles=smiles,
            resolved_name=query,
        )

    @staticmethod
    def _key(value: str) -> str:
        return " ".join(str(value or "").strip().lower().split())


class PubChemNameConnector:
    """Conector online PubChem PUG REST."""

    source = "pubchem"

    def __init__(self, cache_path: Path | None = None) -> None:
        self.cache_path = cache_path or Path.home() / ".chemuson" / "name2structure_cache.json"

    def resolve(self, name: str, timeout_s: float = 8.0) -> NameToStructureResult:
        query = str(name or "").strip()
        if not query:
            return NameToStructureResult(query, None, self.source, 0.0, message="empty_query")

        cached = self._read_cache().get(self._cache_key(query), {})
        if isinstance(cached, dict) and cached.get("smiles"):
            graph, error = _smiles_to_graph(str(cached["smiles"]), timeout_s=timeout_s)
            if graph is not None:
                return NameToStructureResult(
                    query=query,
                    graph=graph,
                    source=self.source,
                    confidence=float(cached.get("confidence", 0.86)),
                    smiles=str(cached["smiles"]),
                    resolved_name=str(cached.get("resolved_name", query)),
                    from_cache=True,
                )
            if error:
                return NameToStructureResult(query, None, self.source, 0.0, message=error)

        try:
            smiles, resolved_name = self._fetch_smiles(query, timeout_s=timeout_s)
        except (HTTPError, URLError, TimeoutError) as exc:
            return NameToStructureResult(query, None, self.source, 0.0, message=exc.__class__.__name__)
        except Exception as exc:
            return NameToStructureResult(query, None, self.source, 0.0, message=str(exc))

        graph, error = _smiles_to_graph(smiles, timeout_s=timeout_s)
        if graph is None:
            return NameToStructureResult(
                query,
                None,
                self.source,
                0.0,
                smiles=smiles,
                resolved_name=resolved_name,
                message=error,
            )
        self._write_cache_entry(query, smiles, resolved_name)
        return NameToStructureResult(
            query=query,
            graph=graph,
            source=self.source,
            confidence=0.86,
            smiles=smiles,
            resolved_name=resolved_name,
        )

    def _fetch_smiles(self, query: str, timeout_s: float) -> tuple[str, str]:
        encoded = quote(query, safe="")
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{encoded}/property/IsomericSMILES,CanonicalSMILES,IUPACName/JSON"
        )
        request = Request(url, headers={"User-Agent": "Chemuson/Name2Structure"})
        with urlopen(request, timeout=max(1.0, float(timeout_s))) as response:
            payload = json.loads(response.read().decode("utf-8"))
        props = payload.get("PropertyTable", {}).get("Properties", [])
        if not props:
            raise ValueError("not_found")
        first = props[0]
        smiles = str(first.get("IsomericSMILES") or first.get("CanonicalSMILES") or "").strip()
        if not smiles:
            raise ValueError("empty_smiles")
        resolved_name = str(first.get("IUPACName") or query).strip()
        return smiles, resolved_name

    def _read_cache(self) -> dict[str, object]:
        try:
            return json.loads(self.cache_path.read_text(encoding="utf-8"))
        except Exception:
            return {}

    def _write_cache_entry(self, query: str, smiles: str, resolved_name: str) -> None:
        cache = self._read_cache()
        cache[self._cache_key(query)] = {
            "smiles": smiles,
            "resolved_name": resolved_name,
            "confidence": 0.86,
        }
        try:
            self.cache_path.parent.mkdir(parents=True, exist_ok=True)
            self.cache_path.write_text(
                json.dumps(cache, ensure_ascii=False, indent=2, sort_keys=True),
                encoding="utf-8",
            )
        except Exception:
            return

    @staticmethod
    def _cache_key(query: str) -> str:
        return " ".join(str(query or "").strip().lower().split())


def resolve_name_to_structure(
    name: str,
    *,
    allow_network: bool = True,
    timeout_s: float = 8.0,
    connectors: list[NameConnector] | None = None,
) -> NameToStructureResult:
    """Resuelve un nombre común/sistemático a ``MolGraph``."""
    query = str(name or "").strip()
    if not query:
        return NameToStructureResult(query, None, "none", 0.0, message="empty_query")

    active_connectors = connectors or [StaticNameConnector()]
    if connectors is None and allow_network:
        active_connectors.append(PubChemNameConnector())

    last = NameToStructureResult(query, None, "none", 0.0, message="not_found")
    for connector in active_connectors:
        result = connector.resolve(query, timeout_s=timeout_s)
        if result.ok:
            return result
        last = result
    return last


def _smiles_to_graph(smiles: str, timeout_s: float) -> tuple[MolGraph | None, str]:
    from chemuson.chemio.rdkit_safe import smiles_to_molgraph_isolated

    graph, error = smiles_to_molgraph_isolated(smiles, timeout_s=timeout_s)
    if graph is None:
        fallback = _fallback_common_smiles_to_graph(smiles)
        if fallback is not None:
            return fallback, ""
        return None, str(error or "conversion_failed")
    return graph, ""


def _fallback_common_smiles_to_graph(smiles: str) -> MolGraph | None:
    """Construye grafos para el subconjunto offline sin depender de RDKit."""
    key = str(smiles or "").strip()
    builders = {
        "O": _graph_water,
        "C": lambda: _graph_chain(["C"]),
        "CC": lambda: _graph_chain(["C", "C"]),
        "CCO": lambda: _graph_chain(["C", "C", "O"]),
        "CO": lambda: _graph_chain(["C", "O"]),
        "CC(=O)C": _graph_acetone,
        "c1ccccc1": lambda: _graph_benzene(),
        "Cc1ccccc1": _graph_toluene,
        "Oc1ccccc1": _graph_phenol,
        "CC(=O)O": _graph_acetic_acid,
    }
    builder = builders.get(key)
    if builder is None:
        return None
    return builder()


def _graph_chain(elements: list[str]) -> MolGraph:
    graph = MolGraph()
    atoms = [graph.add_atom(element, float(idx) * 42.0, 0.0, is_explicit=element != "C") for idx, element in enumerate(elements)]
    for left, right in zip(atoms, atoms[1:]):
        graph.add_bond(left.id, right.id, order=1)
    return graph


def _graph_water() -> MolGraph:
    graph = MolGraph()
    graph.add_atom("O", 0.0, 0.0, is_explicit=True)
    return graph


def _graph_acetone() -> MolGraph:
    graph = MolGraph()
    c1 = graph.add_atom("C", 0.0, 0.0)
    c2 = graph.add_atom("C", 42.0, 0.0)
    o = graph.add_atom("O", 42.0, -42.0, is_explicit=True)
    c3 = graph.add_atom("C", 84.0, 0.0)
    graph.add_bond(c1.id, c2.id, order=1)
    graph.add_bond(c2.id, o.id, order=2)
    graph.add_bond(c2.id, c3.id, order=1)
    return graph


def _graph_acetic_acid() -> MolGraph:
    graph = MolGraph()
    c1 = graph.add_atom("C", 0.0, 0.0)
    c2 = graph.add_atom("C", 42.0, 0.0)
    o1 = graph.add_atom("O", 42.0, -42.0, is_explicit=True)
    o2 = graph.add_atom("O", 84.0, 0.0, is_explicit=True, explicit_h=1)
    graph.add_bond(c1.id, c2.id, order=1)
    graph.add_bond(c2.id, o1.id, order=2)
    graph.add_bond(c2.id, o2.id, order=1)
    return graph


def _graph_benzene() -> MolGraph:
    graph = MolGraph()
    import math

    radius = 42.0
    atoms = []
    for idx in range(6):
        angle = math.radians(60.0 * idx - 30.0)
        atom = graph.add_atom("C", math.cos(angle) * radius, math.sin(angle) * radius)
        atom.is_aromatic = True
        atoms.append(atom)
    for idx in range(6):
        graph.add_bond(atoms[idx].id, atoms[(idx + 1) % 6].id, order=1, is_aromatic=True)
    return graph


def _graph_toluene() -> MolGraph:
    graph = _graph_benzene()
    anchor = graph.get_atom(1)
    methyl = graph.add_atom("C", anchor.x + 42.0, anchor.y - 42.0)
    graph.add_bond(anchor.id, methyl.id, order=1)
    return graph


def _graph_phenol() -> MolGraph:
    graph = _graph_benzene()
    anchor = graph.get_atom(1)
    oxygen = graph.add_atom("O", anchor.x + 42.0, anchor.y - 42.0, is_explicit=True, explicit_h=1)
    graph.add_bond(anchor.id, oxygen.id, order=1)
    return graph
