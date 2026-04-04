from __future__ import annotations

from typing import Callable

from chemuson.core.model import MolGraph


def safe_smiles(
    graph: MolGraph,
    *,
    molgraph_to_smiles_fn: Callable[[MolGraph], str],
) -> str:
    """Exporta SMILES en modo tolerante."""
    try:
        return molgraph_to_smiles_fn(graph)
    except Exception:
        return ""


def graph_from_template_payload(
    template: dict,
    *,
    molfile_to_molgraph_fn: Callable[[str], MolGraph],
    smiles_to_molgraph_fn: Callable[[str], MolGraph],
) -> MolGraph:
    """Convierte una plantilla almacenada a `MolGraph`."""
    molblock = str(template.get("molblock", ""))
    if molblock:
        try:
            return molfile_to_molgraph_fn(molblock)
        except Exception:
            smiles = str(template.get("smiles", "")).strip()
            if smiles:
                return smiles_to_molgraph_fn(smiles)
            raise
    smiles = str(template.get("smiles", "")).strip()
    if smiles:
        return smiles_to_molgraph_fn(smiles)
    raise ValueError("Plantilla sin contenido químico")
