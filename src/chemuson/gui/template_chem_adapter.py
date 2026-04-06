"""Adaptadores químicos para convertir entre plantillas y `MolGraph`."""

from __future__ import annotations

from typing import Callable

from chemuson.core.model import MolGraph
from chemuson.gui.template_conversion import graph_from_template_payload


class TemplateChemAdapter:
    """Aísla las conversiones químicas de la fachada pública TemplateLibrary."""

    def __init__(
        self,
        *,
        molgraph_to_molfile_fn: Callable[[MolGraph], str],
        safe_smiles_fn: Callable[[MolGraph], str],
        molfile_to_molgraph_fn: Callable[[str], MolGraph],
        smiles_to_molgraph_fn: Callable[[str], MolGraph],
    ) -> None:
        self._molgraph_to_molfile_fn = molgraph_to_molfile_fn
        self._safe_smiles_fn = safe_smiles_fn
        self._molfile_to_molgraph_fn = molfile_to_molgraph_fn
        self._smiles_to_molgraph_fn = smiles_to_molgraph_fn

    def add_template_from_graph(self, repository, name: str, category: str, graph: MolGraph) -> dict:
        if not graph.atoms:
            raise ValueError("La plantilla está vacía")
        molblock = self._molgraph_to_molfile_fn(graph)
        return repository.add_template(
            name,
            category,
            molblock,
            smiles=self._safe_smiles_fn(graph),
        )

    def graph_from_template(self, template: dict) -> MolGraph:
        return graph_from_template_payload(
            template,
            molfile_to_molgraph_fn=self._molfile_to_molgraph_fn,
            smiles_to_molgraph_fn=self._smiles_to_molgraph_fn,
        )
