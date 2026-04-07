from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Callable

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction, QIcon, QPainter, QPen, QPixmap, QColor
from PyQt6.QtWidgets import QMessageBox, QWidget


@dataclass(slots=True)
class TemplateBrowserContext:
    parent: QWidget
    template_library: object
    templates_menu: object
    templates_dock: object
    action_save_template: QAction
    action_template_new_category: QAction
    action_template_import_library: QAction
    action_template_export_library: QAction
    preview_cache: dict[str, QIcon]
    show_status: Callable[[str], None]
    start_template_insert_by_id: Callable[[str], None]
    insert_template: Callable[[str, object], None]
    default_category_user: str


class TemplateBrowserService:
    """Presentación y migración de biblioteca de plantillas."""

    def migrate_legacy_templates(self, context: TemplateBrowserContext) -> None:
        try:
            templates_dir = self._get_templates_dir()
            files = [
                path
                for path in sorted(os.listdir(templates_dir))
                if path.lower().endswith(".mol")
            ]
            if not files:
                return
            existing = {
                (
                    str(template.get("name", "")).strip().lower(),
                    str(template.get("molblock", "")).strip(),
                )
                for template in context.template_library.list_templates()
            }
            changed = False
            for filename in files:
                filepath = os.path.join(templates_dir, filename)
                try:
                    with open(filepath, "r", encoding="utf-8") as handle:
                        molblock = handle.read().strip()
                except Exception:
                    continue
                if not molblock:
                    continue
                name = os.path.splitext(filename)[0].replace("_", " ").strip() or "Plantilla"
                signature = (name.lower(), molblock)
                if signature in existing:
                    continue
                context.template_library.add_template(
                    name=name,
                    category=context.default_category_user,
                    molblock=molblock,
                )
                existing.add(signature)
                changed = True
            if changed:
                context.show_status("Plantillas legadas importadas a la biblioteca.")
        except Exception:
            return

    def template_preview_icon(
        self,
        context: TemplateBrowserContext,
        template_id: str,
    ) -> QIcon:
        cache = context.preview_cache.get(template_id)
        if cache is not None:
            return cache
        icon = QIcon()
        try:
            graph = context.template_library.graph_from_template(template_id)
            if not graph.atoms:
                context.preview_cache[template_id] = icon
                return icon
            pixmap = QPixmap(88, 56)
            pixmap.fill(Qt.GlobalColor.transparent)
            painter = QPainter(pixmap)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)

            xs = [atom.x for atom in graph.atoms.values()]
            ys = [atom.y for atom in graph.atoms.values()]
            min_x, max_x = min(xs), max(xs)
            min_y, max_y = min(ys), max(ys)
            width = max(1.0, max_x - min_x)
            height = max(1.0, max_y - min_y)
            margin = 8.0
            scale = min(
                (pixmap.width() - 2.0 * margin) / width,
                (pixmap.height() - 2.0 * margin) / height,
            )

            def map_point(x: float, y: float) -> tuple[float, float]:
                sx = margin + (x - min_x) * scale
                sy = margin + (max_y - y) * scale
                return sx, sy

            bond_pen = QPen(QColor("#222222"))
            bond_pen.setWidth(2)
            bond_pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            painter.setPen(bond_pen)
            for bond in graph.bonds.values():
                a1 = graph.get_atom(bond.a1_id)
                a2 = graph.get_atom(bond.a2_id)
                x1, y1 = map_point(a1.x, a1.y)
                x2, y2 = map_point(a2.x, a2.y)
                painter.drawLine(int(round(x1)), int(round(y1)), int(round(x2)), int(round(y2)))

            color_map = {
                "N": QColor("#2A56D1"),
                "O": QColor("#D11E1E"),
                "S": QColor("#D48400"),
                "P": QColor("#E66A00"),
                "Cl": QColor("#18B81F"),
                "Br": QColor("#A04300"),
                "I": QColor("#5E2A88"),
            }
            font = painter.font()
            font.setPointSize(8)
            painter.setFont(font)
            for atom in graph.atoms.values():
                if atom.element == "C":
                    continue
                x, y = map_point(atom.x, atom.y)
                painter.setPen(color_map.get(atom.element, QColor("#1A1A1A")))
                painter.drawText(int(round(x)) - 6, int(round(y)) + 4, atom.element)
            painter.end()
            icon = QIcon(pixmap)
        except Exception:
            icon = QIcon()
        context.preview_cache[template_id] = icon
        return icon

    def refresh_template_views(self, context: TemplateBrowserContext) -> None:
        context.preview_cache.clear()
        grouped = context.template_library.grouped_templates()
        grouped_with_icons: list[dict] = []
        for group in grouped:
            templates = []
            for template in group.get("templates", []):
                entry = dict(template)
                template_id = str(entry.get("id", "")).strip()
                if template_id:
                    entry["icon"] = self.template_preview_icon(context, template_id)
                templates.append(entry)
            grouped_with_icons.append({"name": group.get("name", ""), "templates": templates})
        context.templates_dock.set_templates(grouped_with_icons)
        self.refresh_templates_menu(context)

    def refresh_templates_menu(self, context: TemplateBrowserContext) -> None:
        context.templates_menu.clear()
        grouped = context.template_library.grouped_templates()
        total_templates = 0
        for group in grouped:
            category = str(group.get("name", "")).strip()
            if not category:
                continue
            submenu = context.templates_menu.addMenu(category)
            templates = list(group.get("templates", []))
            if not templates:
                empty = QAction("(Vacío)", context.parent)
                empty.setEnabled(False)
                submenu.addAction(empty)
                continue
            for template in templates:
                template_id = str(template.get("id", "")).strip()
                label = str(template.get("name", "")).strip() or "Plantilla"
                if not template_id:
                    continue
                action = QAction(label, context.parent)
                action.triggered.connect(
                    lambda checked=False, tid=template_id: context.start_template_insert_by_id(tid)
                )
                submenu.addAction(action)
                total_templates += 1
        if total_templates == 0:
            empty = QAction("(Sin plantillas)", context.parent)
            empty.setEnabled(False)
            context.templates_menu.addAction(empty)
        context.templates_menu.addSeparator()
        context.templates_menu.addAction(context.action_save_template)
        context.templates_menu.addAction(context.action_template_new_category)
        context.templates_menu.addAction(context.action_template_import_library)
        context.templates_menu.addAction(context.action_template_export_library)

    def insert_template_from_file(self, context: TemplateBrowserContext, filepath: str) -> None:
        try:
            from chemuson.chemio.rdkit_io import molfile_to_molgraph

            with open(filepath, "r", encoding="utf-8") as handle:
                molblock = handle.read()
            graph = molfile_to_molgraph(molblock)
            name = os.path.splitext(os.path.basename(filepath))[0]
            context.insert_template(name, graph)
        except Exception as exc:
            QMessageBox.critical(context.parent, "Error", f"Error al cargar plantilla:\n{exc}")

    @staticmethod
    def _get_templates_dir() -> str:
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        templates_dir = os.path.join(base_dir, "templates")
        if not os.path.exists(templates_dir):
            os.makedirs(templates_dir)
        return templates_dir
