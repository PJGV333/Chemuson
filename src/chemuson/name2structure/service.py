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
        return None, str(error or "conversion_failed")
    return graph, ""
