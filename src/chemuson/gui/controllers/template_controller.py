from __future__ import annotations

import os
from typing import Optional

from PyQt6.QtWidgets import QFileDialog, QInputDialog, QMessageBox

from chemuson.gui.template_library import DEFAULT_CATEGORY_USER
from chemuson.gui.templates import build_linear_chain_template


class TemplateController:
    """Gestiona operaciones de biblioteca de plantillas y SMILES."""

    def start_template_insert_by_id(self, window, template_id: str, *, place_now: bool = False) -> None:
        try:
            template = window.template_library.get_template(template_id)
            graph = window.template_library.graph_from_template(template_id)
            label = str(template.get("name", "Plantilla")).strip() or "Plantilla"
            if place_now:
                target = window.canvas._last_scene_pos
                if target is None:
                    target = window.canvas.mapToScene(window.canvas.viewport().rect().center())
                window.canvas._insert_molgraph_at(graph, target)
                window.canvas.cancel_template_insert_mode()
                window.statusBar().showMessage(f"Plantilla '{label}' insertada")
            else:
                window._insert_template(label, graph)
        except Exception as exc:
            QMessageBox.critical(window, "Error", f"No se pudo cargar la plantilla:\n{exc}")

    def prompt_template_destination(self, window) -> Optional[tuple[str, str]]:
        name, ok = QInputDialog.getText(window, "Guardar plantilla", "Nombre de la plantilla:")
        if not ok:
            return None
        clean_name = " ".join(name.strip().split())
        if not clean_name:
            return None
        categories = window.template_library.categories()
        if not categories:
            categories = [DEFAULT_CATEGORY_USER]
        if DEFAULT_CATEGORY_USER not in categories:
            categories.append(DEFAULT_CATEGORY_USER)
        choices = categories + ["+ Nueva categoría..."]
        default_idx = choices.index(DEFAULT_CATEGORY_USER) if DEFAULT_CATEGORY_USER in choices else 0
        category_choice, ok = QInputDialog.getItem(
            window,
            "Guardar plantilla",
            "Categoría:",
            choices,
            default_idx,
            False,
        )
        if not ok:
            return None
        category = category_choice
        if category_choice == "+ Nueva categoría...":
            new_category, ok = QInputDialog.getText(window, "Nueva categoría", "Nombre de la categoría:")
            if not ok:
                return None
            category = window.template_library.add_category(new_category)
        else:
            category = window.template_library.add_category(category_choice)
        return clean_name, category

    def on_save_template(self, window) -> None:
        try:
            atom_ids, bonds = window.canvas._selected_structure_ids()
            graph_to_save = (
                window.canvas._build_selection_graph(atom_ids, bonds)
                if atom_ids
                else window.canvas.graph
            )
            if not graph_to_save.atoms:
                QMessageBox.warning(window, "Aviso", "No hay nada para guardar.")
                return
            target = self.prompt_template_destination(window)
            if target is None:
                return
            name, category = target
            window.template_library.add_template_from_graph(name, category, graph_to_save)
            window._refresh_template_views()
            window.statusBar().showMessage(f"Plantilla guardada: {name} ({category})")
        except Exception as exc:
            QMessageBox.critical(window, "Error", f"Error al guardar plantilla:\n{exc}")

    def on_new_template_category(self, window) -> None:
        name, ok = QInputDialog.getText(window, "Nueva categoría", "Nombre de la categoría:")
        if not ok:
            return
        category = " ".join(name.strip().split())
        if not category:
            return
        window.template_library.add_category(category)
        window._refresh_template_views()
        window.statusBar().showMessage(f"Categoría creada: {category}")

    def on_rename_template_category(self, window, old_name: str) -> None:
        if not old_name:
            return
        new_name, ok = QInputDialog.getText(
            window,
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
            window.template_library.rename_category(old_name, clean_new)
            window._refresh_template_views()
            window.statusBar().showMessage(f"Categoría renombrada: {old_name} -> {clean_new}")
        except Exception as exc:
            QMessageBox.critical(window, "Error", f"No se pudo renombrar la categoría:\n{exc}")

    def on_delete_template_category(self, window, name: str) -> None:
        if not name:
            return
        reply = QMessageBox.question(
            window,
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
            window.template_library.delete_category(name, fallback_category=DEFAULT_CATEGORY_USER)
            window._refresh_template_views()
            window.statusBar().showMessage(f"Categoría eliminada: {name}")
        except Exception as exc:
            QMessageBox.critical(window, "Error", f"No se pudo eliminar la categoría:\n{exc}")

    def on_rename_template(self, window, template_id: str) -> None:
        if not template_id:
            return
        try:
            template = window.template_library.get_template(template_id)
        except Exception:
            return
        current_name = str(template.get("name", "Plantilla"))
        new_name, ok = QInputDialog.getText(
            window,
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
            window.template_library.rename_template(template_id, clean_new)
            window._refresh_template_views()
            window.statusBar().showMessage(f"Plantilla renombrada: {clean_new}")
        except Exception as exc:
            QMessageBox.critical(window, "Error", f"No se pudo renombrar la plantilla:\n{exc}")

    def on_delete_template(self, window, template_id: str) -> None:
        if not template_id:
            return
        try:
            template = window.template_library.get_template(template_id)
        except Exception:
            return
        name = str(template.get("name", "Plantilla"))
        reply = QMessageBox.question(
            window,
            "Eliminar plantilla",
            f"¿Eliminar la plantilla '{name}'?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        try:
            window.template_library.delete_template(template_id)
            window._refresh_template_views()
            window.statusBar().showMessage(f"Plantilla eliminada: {name}")
        except Exception as exc:
            QMessageBox.critical(window, "Error", f"No se pudo eliminar la plantilla:\n{exc}")

    def on_import_template_library(self, window) -> None:
        filepath, _ = QFileDialog.getOpenFileName(
            window,
            "Importar biblioteca de plantillas",
            "",
            "Template Library (*.json)",
        )
        if not filepath:
            return
        choice = QMessageBox.question(
            window,
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
            added = window.template_library.import_from_file(filepath, merge=merge)
            window._refresh_template_views()
            if merge:
                window.statusBar().showMessage(f"Biblioteca importada: {added} plantilla(s) agregadas.")
            else:
                window.statusBar().showMessage("Biblioteca importada (reemplazo completo).")
        except Exception as exc:
            QMessageBox.critical(window, "Error", f"No se pudo importar la biblioteca:\n{exc}")

    def on_export_template_library(self, window) -> None:
        filepath, _ = QFileDialog.getSaveFileName(
            window,
            "Exportar biblioteca de plantillas",
            "chemuson_templates.json",
            "Template Library (*.json)",
        )
        if not filepath:
            return
        try:
            window.template_library.export_to_file(filepath)
            window.statusBar().showMessage("Biblioteca de plantillas exportada.")
        except Exception as exc:
            QMessageBox.critical(window, "Error", f"No se pudo exportar la biblioteca:\n{exc}")

    def on_insert_linear_chain(self, window) -> None:
        graph = build_linear_chain_template(window.canvas.state.bond_length)
        window._insert_template("Cadena lineal", graph)

    def on_import_smiles(self, window) -> None:
        smiles, ok = QInputDialog.getText(window, "Importar SMILES", "SMILES:")
        if not ok or not smiles.strip():
            return
        try:
            from chemuson.chemio.rdkit_io import smiles_to_molgraph

            graph = smiles_to_molgraph(smiles.strip())
            window.canvas._insert_molgraph(graph)
            window.statusBar().showMessage("SMILES insertado")
        except Exception as exc:
            QMessageBox.critical(window, "Error", f"No se pudo importar SMILES:\n{exc}")

    def on_export_smiles(self, window) -> None:
        try:
            from chemuson.chemio.rdkit_io import molgraph_to_smiles

            smiles = molgraph_to_smiles(window.canvas.graph)
            QMessageBox.information(window, "SMILES", smiles)
        except Exception as exc:
            QMessageBox.critical(window, "Error", f"No se pudo exportar SMILES:\n{exc}")
