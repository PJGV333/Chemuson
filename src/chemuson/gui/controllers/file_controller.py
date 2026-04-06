from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Callable

from PyQt6.QtWidgets import QApplication, QFileDialog, QMessageBox, QWidget

from chemuson.chemio.persistence import PersistenceManager
from chemuson.chemio.rdkit_io import (
    molfile_to_molgraph,
    molgraph_to_molfile,
)
from chemuson.gui.canvas import ChemusonCanvas


@dataclass(slots=True)
class FileWorkflowContext:
    """Dependencias mínimas para abrir/guardar documentos desde la shell principal."""

    parent: QWidget
    canvas: ChemusonCanvas
    tabs: object
    tab_manager: object
    create_document_tab: Callable[[bool], ChemusonCanvas]
    apply_toolbar_defaults: Callable[[ChemusonCanvas], None]
    set_active_canvas: Callable[[ChemusonCanvas], None]
    before_canvas_discard: Callable[[ChemusonCanvas], None]
    add_recent_file: Callable[[str], None]
    refresh_recent_menu: Callable[[], None]
    update_total_charge_indicator: Callable[[], None]
    show_status: Callable[[str], None]


class FileController:
    """Flujos de archivo desacoplados de `ChemusonWindow`."""

    def load_file_into_canvas(self, filepath: str, canvas: ChemusonCanvas) -> None:
        """Carga un archivo químico dentro de un canvas específico."""
        if filepath.lower().endswith(".cmsn"):
            canvas.clear_canvas()
            PersistenceManager.load_from_file(filepath, canvas)
            return

        with open(filepath, "r", encoding="utf-8") as fh:
            molfile = fh.read()
        graph = molfile_to_molgraph(molfile)
        canvas.clear_canvas()
        canvas._insert_molgraph(graph)

    def open_recent_file(self, context: FileWorkflowContext, filepath: str) -> None:
        """Abre una entrada reciente si el archivo sigue existiendo."""
        if not filepath or not os.path.exists(filepath):
            QMessageBox.warning(context.parent, "Archivo no encontrado", "El archivo no existe.")
            context.refresh_recent_menu()
            return
        self.open_file_path(context, filepath)

    def open_file_path(self, context: FileWorkflowContext, filepath: str) -> None:
        """Abre un archivo en una pestaña nueva manteniendo el documento actual intacto."""
        if not filepath:
            return
        canvas = context.create_document_tab(True)
        context.apply_toolbar_defaults(canvas)
        try:
            self.load_file_into_canvas(filepath, canvas)
            canvas.undo_stack.setClean()
            context.tab_manager.set_canvas_file_path(canvas, filepath)
            context.tabs.setCurrentWidget(canvas)
            context.set_active_canvas(canvas)
            context.add_recent_file(filepath)
            context.update_total_charge_indicator()
            context.show_status(f"Abierto: {filepath}")
        except Exception as exc:
            context.before_canvas_discard(canvas)
            context.tab_manager.discard_canvas(canvas)
            if context.tabs.count() == 0:
                replacement = context.create_document_tab(True)
                context.apply_toolbar_defaults(replacement)
                context.set_active_canvas(replacement)
            QMessageBox.critical(context.parent, "Error", f"No se pudo abrir el archivo:\n{exc}")

    def save_file(self, context: FileWorkflowContext) -> None:
        """Guarda el documento activo en formato nativo o MOL."""
        filepath = context.tab_manager.file_path_for(context.canvas)
        selected_filter = ""
        if not filepath:
            filepath, selected_filter = QFileDialog.getSaveFileName(
                context.parent,
                "Guardar archivo",
                "",
                "Archivo de Chemuson (*.cmsn);;Archivo MOL (*.mol);;Todos los archivos (*.*)",
            )
        if not filepath:
            return
        try:
            autosave_manager = context.tab_manager.autosave_manager_for(context.canvas)
            if filepath.lower().endswith(".mol") or filepath.lower().endswith(".sdf") or "MOL" in selected_filter:
                molfile = molgraph_to_molfile(context.canvas.graph)
                with open(filepath, "w", encoding="utf-8") as fh:
                    fh.write(molfile)
            else:
                if not filepath.lower().endswith(".cmsn"):
                    filepath += ".cmsn"
                PersistenceManager.save_to_file(filepath, context.canvas)
            context.tab_manager.set_canvas_file_path(context.canvas, filepath)
            context.canvas.undo_stack.setClean()
            if autosave_manager is not None:
                autosave_manager.cleanup_after_save()
            context.add_recent_file(filepath)
            context.show_status(f"Guardado: {filepath}")
        except Exception as exc:
            QMessageBox.critical(context.parent, "Error", f"No se pudo guardar:\n{exc}")

    def copy_as(self, context: FileWorkflowContext, format_name: str) -> None:
        """Copia la estructura actual en el formato textual solicitado."""
        try:
            clipboard = QApplication.clipboard()
            graph, _bbox = context.canvas._analysis_graph_and_bbox()
            target_graph = graph if graph is not None else context.canvas.graph

            if format_name == "smiles":
                from chemuson.chemio import rdkit_io

                text = rdkit_io.molgraph_to_smiles(target_graph) if target_graph.atoms else ""
            elif format_name == "molfile":
                text = molgraph_to_molfile(target_graph) if target_graph.atoms else ""
            elif format_name == "inchi":
                if target_graph.atoms:
                    from rdkit.Chem.inchi import MolToInchi
                    from chemuson.chemio import rdkit_io

                    text = MolToInchi(rdkit_io.molgraph_to_rdkit(target_graph))
                else:
                    text = ""
            else:
                text = ""

            clipboard.setText(text)
            context.show_status(f"Copiado como {format_name.upper()}")
        except Exception as exc:
            context.show_status(f"Error: {exc}")
