from __future__ import annotations

import copy
from dataclasses import dataclass
from typing import Callable, Optional

from PyQt6.QtCore import QObject, QThread, pyqtSignal, pyqtSlot
from PyQt6.QtWidgets import QApplication, QFileDialog, QInputDialog, QMessageBox, QWidget

from chemuson.core.model import MolGraph
from chemuson.gui.template_library import DEFAULT_CATEGORY_USER, TemplateLibrary
from chemuson.gui.templates import build_linear_chain_template


@dataclass(slots=True)
class TemplateControllerContext:
    """Dependencias mínimas para operaciones de plantillas."""

    parent: QWidget
    canvas: object
    template_library: TemplateLibrary
    show_status: Callable[[str], None]
    refresh_template_views: Callable[[], None]
    insert_template: Callable[[str, MolGraph], None]


class _SmilesExportWorker(QObject):
    """Worker para exportar SMILES sin bloquear la interfaz."""

    finished = pyqtSignal(int, str, str)

    def __init__(self, job_id: int, graph: MolGraph) -> None:
        super().__init__()
        self._job_id = int(job_id)
        self._graph = graph

    @pyqtSlot()
    def run(self) -> None:
        try:
            from chemuson.chemio.rdkit_io import molgraph_to_smiles_isolated_or_error

            smiles = molgraph_to_smiles_isolated_or_error(self._graph)
            self.finished.emit(self._job_id, smiles, "")
        except Exception as exc:
            self.finished.emit(self._job_id, "", str(exc))


@dataclass(slots=True)
class _SmilesExportJob:
    """Estado mínimo de una exportación SMILES en segundo plano."""

    thread: QThread
    worker: _SmilesExportWorker
    parent: QWidget
    show_status: Callable[[str], None]


class TemplateController(QObject):
    """Gestiona operaciones de biblioteca de plantillas y SMILES."""

    def __init__(self) -> None:
        super().__init__()
        self._next_smiles_export_job_id = 1
        self._smiles_export_jobs: dict[int, _SmilesExportJob] = {}

    def start_template_insert_by_id(
        self,
        context: TemplateControllerContext,
        template_id: str,
        *,
        place_now: bool = False,
    ) -> None:
        try:
            template = context.template_library.get_template(template_id)
            graph = context.template_library.graph_from_template(template_id)
            label = str(template.get("name", "Plantilla")).strip() or "Plantilla"
            if place_now:
                target = context.canvas._last_scene_pos
                if target is None:
                    target = context.canvas.mapToScene(context.canvas.viewport().rect().center())
                context.canvas._insert_molgraph_at(graph, target)
                context.canvas.cancel_template_insert_mode()
                context.show_status(f"Plantilla '{label}' insertada")
            else:
                context.insert_template(label, graph)
        except Exception as exc:
            QMessageBox.critical(context.parent, "Error", f"No se pudo cargar la plantilla:\n{exc}")

    def prompt_template_destination(
        self,
        context: TemplateControllerContext,
    ) -> Optional[tuple[str, str]]:
        name, ok = QInputDialog.getText(context.parent, "Guardar plantilla", "Nombre de la plantilla:")
        if not ok:
            return None
        clean_name = " ".join(name.strip().split())
        if not clean_name:
            return None
        categories = context.template_library.categories()
        if not categories:
            categories = [DEFAULT_CATEGORY_USER]
        if DEFAULT_CATEGORY_USER not in categories:
            categories.append(DEFAULT_CATEGORY_USER)
        choices = categories + ["+ Nueva categoría..."]
        default_idx = choices.index(DEFAULT_CATEGORY_USER) if DEFAULT_CATEGORY_USER in choices else 0
        category_choice, ok = QInputDialog.getItem(
            context.parent,
            "Guardar plantilla",
            "Categoría:",
            choices,
            default_idx,
            False,
        )
        if not ok:
            return None
        if category_choice == "+ Nueva categoría...":
            new_category, ok = QInputDialog.getText(
                context.parent,
                "Nueva categoría",
                "Nombre de la categoría:",
            )
            if not ok:
                return None
            category = context.template_library.add_category(new_category)
        else:
            category = context.template_library.add_category(category_choice)
        return clean_name, category

    def on_save_template(self, context: TemplateControllerContext) -> None:
        try:
            atom_ids, bonds = context.canvas._selected_structure_ids()
            graph_to_save = (
                context.canvas._build_selection_graph(atom_ids, bonds)
                if atom_ids
                else context.canvas.graph
            )
            if not graph_to_save.atoms:
                QMessageBox.warning(context.parent, "Aviso", "No hay nada para guardar.")
                return
            target = self.prompt_template_destination(context)
            if target is None:
                return
            name, category = target
            context.template_library.add_template_from_graph(name, category, graph_to_save)
            context.refresh_template_views()
            context.show_status(f"Plantilla guardada: {name} ({category})")
        except Exception as exc:
            QMessageBox.critical(context.parent, "Error", f"Error al guardar plantilla:\n{exc}")

    def on_new_template_category(self, context: TemplateControllerContext) -> None:
        name, ok = QInputDialog.getText(context.parent, "Nueva categoría", "Nombre de la categoría:")
        if not ok:
            return
        category = " ".join(name.strip().split())
        if not category:
            return
        context.template_library.add_category(category)
        context.refresh_template_views()
        context.show_status(f"Categoría creada: {category}")

    def on_rename_template_category(
        self,
        context: TemplateControllerContext,
        old_name: str,
    ) -> None:
        if not old_name:
            return
        new_name, ok = QInputDialog.getText(
            context.parent,
            "Renombrar categoría",
            "Nuevo nombre:",
            text=old_name,
        )
        if not ok:
            return
        clean_new = " ".join(new_name.strip().split())
        if not clean_new:
            return
        try:
            context.template_library.rename_category(old_name, clean_new)
            context.refresh_template_views()
            context.show_status(f"Categoría renombrada: {old_name} -> {clean_new}")
        except Exception as exc:
            QMessageBox.critical(
                context.parent,
                "Error",
                f"No se pudo renombrar la categoría:\n{exc}",
            )

    def on_delete_template_category(
        self,
        context: TemplateControllerContext,
        name: str,
    ) -> None:
        if not name:
            return
        reply = QMessageBox.question(
            context.parent,
            "Eliminar categoría",
            (
                f"¿Eliminar la categoría '{name}'?\n"
                f"Las plantillas se moverán a '{DEFAULT_CATEGORY_USER}'."
            ),
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        try:
            context.template_library.delete_category(
                name,
                fallback_category=DEFAULT_CATEGORY_USER,
            )
            context.refresh_template_views()
            context.show_status(f"Categoría eliminada: {name}")
        except Exception as exc:
            QMessageBox.critical(
                context.parent,
                "Error",
                f"No se pudo eliminar la categoría:\n{exc}",
            )

    def on_rename_template(
        self,
        context: TemplateControllerContext,
        template_id: str,
    ) -> None:
        if not template_id:
            return
        try:
            template = context.template_library.get_template(template_id)
        except Exception:
            return
        current_name = str(template.get("name", "Plantilla"))
        new_name, ok = QInputDialog.getText(
            context.parent,
            "Renombrar plantilla",
            "Nuevo nombre:",
            text=current_name,
        )
        if not ok:
            return
        clean_new = " ".join(new_name.strip().split())
        if not clean_new:
            return
        try:
            context.template_library.rename_template(template_id, clean_new)
            context.refresh_template_views()
            context.show_status(f"Plantilla renombrada: {clean_new}")
        except Exception as exc:
            QMessageBox.critical(
                context.parent,
                "Error",
                f"No se pudo renombrar la plantilla:\n{exc}",
            )

    def on_delete_template(
        self,
        context: TemplateControllerContext,
        template_id: str,
    ) -> None:
        if not template_id:
            return
        try:
            template = context.template_library.get_template(template_id)
        except Exception:
            return
        name = str(template.get("name", "Plantilla"))
        reply = QMessageBox.question(
            context.parent,
            "Eliminar plantilla",
            f"¿Eliminar la plantilla '{name}'?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        try:
            context.template_library.delete_template(template_id)
            context.refresh_template_views()
            context.show_status(f"Plantilla eliminada: {name}")
        except Exception as exc:
            QMessageBox.critical(
                context.parent,
                "Error",
                f"No se pudo eliminar la plantilla:\n{exc}",
            )

    def on_import_template_library(self, context: TemplateControllerContext) -> None:
        filepath, _ = QFileDialog.getOpenFileName(
            context.parent,
            "Importar biblioteca de plantillas",
            "",
            "Template Library (*.json)",
        )
        if not filepath:
            return
        choice = QMessageBox.question(
            context.parent,
            "Importar biblioteca",
            "¿Combinar con la biblioteca actual? (No = reemplazar)",
            QMessageBox.StandardButton.Yes
            | QMessageBox.StandardButton.No
            | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Yes,
        )
        if choice == QMessageBox.StandardButton.Cancel:
            return
        merge = choice == QMessageBox.StandardButton.Yes
        try:
            added = context.template_library.import_from_file(filepath, merge=merge)
            context.refresh_template_views()
            if merge:
                context.show_status(f"Biblioteca importada: {added} plantilla(s) agregadas.")
            else:
                context.show_status("Biblioteca importada (reemplazo completo).")
        except Exception as exc:
            QMessageBox.critical(
                context.parent,
                "Error",
                f"No se pudo importar la biblioteca:\n{exc}",
            )

    def on_export_template_library(self, context: TemplateControllerContext) -> None:
        filepath, _ = QFileDialog.getSaveFileName(
            context.parent,
            "Exportar biblioteca de plantillas",
            "chemuson_templates.json",
            "Template Library (*.json)",
        )
        if not filepath:
            return
        try:
            context.template_library.export_to_file(filepath)
            context.show_status("Biblioteca de plantillas exportada.")
        except Exception as exc:
            QMessageBox.critical(
                context.parent,
                "Error",
                f"No se pudo exportar la biblioteca:\n{exc}",
            )

    def on_insert_linear_chain(self, context: TemplateControllerContext) -> None:
        graph = build_linear_chain_template(context.canvas.state.bond_length)
        context.insert_template("Cadena lineal", graph)

    def on_import_smiles(self, context: TemplateControllerContext) -> None:
        smiles, ok = QInputDialog.getText(context.parent, "Importar SMILES", "SMILES:")
        if not ok or not smiles.strip():
            return
        try:
            from chemuson.chemio.rdkit_io import smiles_to_molgraph

            graph = smiles_to_molgraph(smiles.strip())
            context.canvas._insert_molgraph(graph)
            context.show_status("SMILES insertado")
        except Exception as exc:
            QMessageBox.critical(
                context.parent,
                "Error",
                f"No se pudo importar SMILES:\n{exc}",
            )

    def on_export_smiles(self, context: TemplateControllerContext) -> None:
        try:
            atom_ids, bonds = context.canvas._selected_structure_ids()
            target_graph = (
                context.canvas._build_selection_graph(atom_ids, bonds)
                if atom_ids
                else context.canvas.graph
            )
            if not target_graph.atoms:
                QMessageBox.warning(context.parent, "Aviso", "No hay estructura para exportar.")
                return
            self._start_smiles_export_job(context, target_graph)
        except Exception as exc:
            QMessageBox.critical(
                context.parent,
                "Error",
                f"No se pudo exportar SMILES:\n{exc}",
            )

    def _start_smiles_export_job(
        self,
        context: TemplateControllerContext,
        graph: MolGraph,
    ) -> None:
        """Inicia la exportación SMILES fuera del hilo de la UI."""
        if self._smiles_export_jobs:
            context.show_status("Exportación SMILES en curso...")
            return

        job_id = self._next_smiles_export_job_id
        self._next_smiles_export_job_id += 1
        thread = QThread(context.parent)
        worker = _SmilesExportWorker(job_id, copy.deepcopy(graph))
        worker.moveToThread(thread)

        self._smiles_export_jobs[job_id] = _SmilesExportJob(
            thread=thread,
            worker=worker,
            parent=context.parent,
            show_status=context.show_status,
        )
        thread.started.connect(worker.run)
        worker.finished.connect(self._on_smiles_export_finished)
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        context.show_status("Exportando SMILES...")
        thread.start()

    @pyqtSlot(int, str, str)
    def _on_smiles_export_finished(self, job_id: int, smiles: str, error: str) -> None:
        """Muestra el resultado de una exportación SMILES terminada."""
        job = self._smiles_export_jobs.pop(int(job_id), None)
        if job is None:
            return
        if error:
            QMessageBox.critical(
                job.parent,
                "Error",
                f"No se pudo exportar SMILES:\n{error}",
            )
            job.show_status("Error al exportar SMILES")
            return
        QApplication.clipboard().setText(smiles)
        QMessageBox.information(job.parent, "SMILES", smiles)
        job.show_status("SMILES exportado y copiado al portapapeles")
