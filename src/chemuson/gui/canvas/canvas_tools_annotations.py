from __future__ import annotations

from copy import deepcopy
from typing import Optional

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QTextDocument
from PyQt6.QtWidgets import QColorDialog, QGraphicsItem, QInputDialog

from chemuson.gui.commands import (
    ChangeCanvasOpacityCommand,
    ConfigureEnergyDiagramItemsCommand,
    EditSemanticDiagramCommand,
    StyleOrbitalItemsCommand,
)
from chemuson.gui.composite_diagram_item import CompositeDiagramItem
from chemuson.gui.energy_diagrams import normalize_energy_occupancies
from chemuson.gui.items import EnergyDiagramItem, OrbitalAnnotationItem

class CanvasToolsAnnotationsMixin:
    @staticmethod
    def _orbital_part_display_name(name: str) -> str:
        labels = {
            "sphere": "Esfera",
            "top": "Lóbulo superior",
            "bottom": "Lóbulo inferior",
            "left": "Lóbulo izquierdo",
            "right": "Lóbulo derecho",
            "major": "Lóbulo principal",
            "minor": "Lóbulo secundario",
            "minor_left": "Lóbulo menor izquierdo",
            "minor_right": "Lóbulo menor derecho",
            "bond": "Nube central",
            "upper": "Nube superior",
            "lower": "Nube inferior",
            "ring": "Toroide",
            "ring_upper": "Toroide superior",
            "ring_lower": "Toroide inferior",
        }
        if name.startswith("lobe_"):
            suffix = name.split("_", 1)[1]
            if suffix.isdigit():
                return f"Lóbulo {int(suffix) + 1}"
        return labels.get(name, name.replace("_", " "))

    def _single_orbital_target(
        self,
        clicked_item: Optional[QGraphicsItem] = None,
    ) -> Optional[OrbitalAnnotationItem]:
        if isinstance(clicked_item, OrbitalAnnotationItem):
            return clicked_item
        selected = self._selected_orbital_items()
        return selected[0] if len(selected) == 1 else None

    def _choose_orbital_part(
        self,
        item: OrbitalAnnotationItem,
        *,
        title: str = "Seleccionar lóbulo",
    ) -> Optional[str]:
        part_names = item.part_names()
        if not part_names:
            return None
        display_to_name = {
            self._orbital_part_display_name(name): name
            for name in part_names
        }
        value, ok = QInputDialog.getItem(
            self,
            title,
            "Lóbulo:",
            list(display_to_name.keys()),
            0,
            False,
        )
        if not ok:
            return None
        return display_to_name.get(value)

    def _push_orbital_style_change(
        self,
        updates: dict[OrbitalAnnotationItem, dict[str, dict[str, object]]],
        *,
        text: str,
    ) -> bool:
        before: dict[OrbitalAnnotationItem, dict[str, dict[str, object]]] = {}
        after: dict[OrbitalAnnotationItem, dict[str, dict[str, object]]] = {}
        for item, payload in updates.items():
            if item.scene() is not self.scene:
                continue
            current = item.part_styles()
            if current == payload:
                continue
            before[item] = current
            after[item] = payload
        if not after:
            return False
        self.undo_stack.push(StyleOrbitalItemsCommand(self, before, after, text))
        return True

    def _prompt_selected_orbital_color(self) -> None:
        orbitals = self._selected_orbital_items()
        if not orbitals:
            return
        initial = QColor("#111111")
        for item in orbitals:
            for name in item.part_names():
                color_value = item.effective_part_style(name).get("color")
                if color_value:
                    candidate = QColor(str(color_value))
                    if candidate.isValid():
                        initial = candidate
                        break
            else:
                continue
            break
        color = QColorDialog.getColor(initial, self, "Seleccionar color de orbital")
        if not color.isValid():
            return
        updates: dict[OrbitalAnnotationItem, dict[str, dict[str, object]]] = {}
        for item in orbitals:
            payload = item.part_styles()
            for name in item.part_names():
                part_payload = dict(payload.get(name, {}))
                part_payload["color"] = color.name(QColor.NameFormat.HexRgb)
                payload[name] = part_payload
            updates[item] = payload
        self._push_orbital_style_change(updates, text="Change orbital color")

    def _reset_selected_orbital_color(self) -> None:
        orbitals = self._selected_orbital_items()
        if not orbitals:
            return
        updates: dict[OrbitalAnnotationItem, dict[str, dict[str, object]]] = {}
        for item in orbitals:
            payload = item.part_styles()
            for name in list(payload.keys()):
                part_payload = dict(payload.get(name, {}))
                part_payload.pop("color", None)
                if part_payload:
                    payload[name] = part_payload
                else:
                    payload.pop(name, None)
            updates[item] = payload
        self._push_orbital_style_change(updates, text="Reset orbital color")

    def _prompt_selected_orbital_opacity(self) -> None:
        orbitals = self._selected_orbital_items()
        if not orbitals:
            return
        current = 100
        values = [self.effective_item_opacity(item) for item in orbitals if item.scene() is self.scene]
        if values:
            current = int(round(sum(values) / float(len(values)) * 100.0))
        value, ok = QInputDialog.getInt(
            self,
            "Opacidad de orbital",
            "Opacidad (%):",
            current,
            0,
            100,
            1,
        )
        if not ok:
            return
        target = max(0.0, min(1.0, float(value) / 100.0))
        item_values = {item: target for item in orbitals if item.scene() is self.scene}
        if not item_values:
            return
        self.undo_stack.push(
            ChangeCanvasOpacityCommand(
                self.model,
                self,
                item_values=item_values,
                text="Change orbital opacity",
            )
        )

    def _prompt_orbital_part_color(self, item: OrbitalAnnotationItem) -> None:
        part_name = self._choose_orbital_part(item, title="Color por lóbulo")
        if not part_name:
            return
        current = item.effective_part_style(part_name).get("color")
        initial = QColor(str(current)) if current else QColor("#111111")
        color = QColorDialog.getColor(initial, self, f"Color de {self._orbital_part_display_name(part_name)}")
        if not color.isValid():
            return
        payload = item.part_styles()
        part_payload = dict(payload.get(part_name, {}))
        part_payload["color"] = color.name(QColor.NameFormat.HexRgb)
        payload[part_name] = part_payload
        self._push_orbital_style_change({item: payload}, text="Change orbital lobe color")

    def _prompt_orbital_part_opacity(self, item: OrbitalAnnotationItem) -> None:
        part_name = self._choose_orbital_part(item, title="Opacidad por lóbulo")
        if not part_name:
            return
        current = item.effective_part_style(part_name).get("opacity")
        percent = 100 if current is None else int(round(float(current) * 100.0))
        value, ok = QInputDialog.getInt(
            self,
            "Opacidad por lóbulo",
            f"{self._orbital_part_display_name(part_name)} (%):",
            percent,
            0,
            100,
            1,
        )
        if not ok:
            return
        payload = item.part_styles()
        part_payload = dict(payload.get(part_name, {}))
        opacity = max(0.0, min(1.0, float(value) / 100.0))
        if abs(opacity - 1.0) <= 1e-6:
            part_payload.pop("opacity", None)
        else:
            part_payload["opacity"] = opacity
        if part_payload:
            payload[part_name] = part_payload
        else:
            payload.pop(part_name, None)
        self._push_orbital_style_change({item: payload}, text="Change orbital lobe opacity")

    def _prompt_orbital_part_offset(self, item: OrbitalAnnotationItem) -> None:
        part_name = self._choose_orbital_part(item, title="Mover lóbulo")
        if not part_name:
            return
        current = item.effective_part_style(part_name)
        current_x = float(current.get("offset_x", 0.0) or 0.0)
        current_y = float(current.get("offset_y", 0.0) or 0.0)
        offset_x, ok = QInputDialog.getDouble(
            self,
            "Mover lóbulo",
            f"{self._orbital_part_display_name(part_name)} desplazamiento X:",
            current_x,
            -80.0,
            80.0,
            1,
        )
        if not ok:
            return
        offset_y, ok = QInputDialog.getDouble(
            self,
            "Mover lóbulo",
            f"{self._orbital_part_display_name(part_name)} desplazamiento Y:",
            current_y,
            -80.0,
            80.0,
            1,
        )
        if not ok:
            return
        payload = item.part_styles()
        part_payload = dict(payload.get(part_name, {}))
        if abs(offset_x) <= 1e-6:
            part_payload.pop("offset_x", None)
        else:
            part_payload["offset_x"] = float(offset_x)
        if abs(offset_y) <= 1e-6:
            part_payload.pop("offset_y", None)
        else:
            part_payload["offset_y"] = float(offset_y)
        if part_payload:
            payload[part_name] = part_payload
        else:
            payload.pop(part_name, None)
        self._push_orbital_style_change({item: payload}, text="Move orbital lobe")

    def _reset_orbital_part_style(self, item: OrbitalAnnotationItem) -> None:
        part_name = self._choose_orbital_part(item, title="Restablecer lóbulo")
        if not part_name:
            return
        payload = item.part_styles()
        payload.pop(part_name, None)
        self._push_orbital_style_change({item: payload}, text="Reset orbital lobe style")

    def _single_energy_diagram_target(
        self,
        clicked_item: Optional[QGraphicsItem] = None,
    ) -> Optional[EnergyDiagramItem]:
        """Devuelve un diagrama único para acciones de edición puntual."""
        if isinstance(clicked_item, EnergyDiagramItem):
            return clicked_item
        selected = self._selected_energy_diagram_items()
        return selected[0] if len(selected) == 1 else None

    def _single_semantic_diagram_target(
        self,
        clicked_item: Optional[QGraphicsItem] = None,
    ) -> Optional[CompositeDiagramItem]:
        """Devuelve un diagrama semántico único para acciones de edición."""
        if isinstance(clicked_item, CompositeDiagramItem):
            return clicked_item
        selected = self._selected_semantic_diagram_items()
        return selected[0] if len(selected) == 1 else None

    def _set_semantic_diagram_summary_visible(
        self,
        item: CompositeDiagramItem,
        visible: bool,
    ) -> bool:
        """Alterna la visibilidad del resumen textual inferior con undo/redo."""
        if item is None or item.scene() is not self.scene:
            return False
        before_payload = item.to_json()
        before_payload["opacity"] = self.item_raw_opacity(item)
        current_visible = bool(
            item.semantic_diagram.metadata.get("show_summary", True)
        )
        if current_visible == bool(visible):
            return False
        after_payload = deepcopy(before_payload)
        semantic_payload = dict(after_payload.get("semantic_diagram", {}) or {})
        metadata = dict(semantic_payload.get("metadata", {}) or {})
        metadata["show_summary"] = bool(visible)
        semantic_payload["metadata"] = metadata
        after_payload["semantic_diagram"] = semantic_payload
        self.undo_stack.push(
            EditSemanticDiagramCommand(
                self,
                item,
                before_payload,
                after_payload,
                text="Toggle semantic diagram summary",
            )
        )
        return True

    @staticmethod
    def _plain_text_from_markup(value: object) -> str:
        raw = str(value or "")
        if not raw:
            return ""
        if not Qt.mightBeRichText(raw):
            return raw
        document = QTextDocument()
        document.setHtml(raw)
        return document.toPlainText()

    def _prompt_rich_text_value(
        self,
        title: str,
        label: str,
        initial_text: str,
    ) -> tuple[str, bool]:
        """Abre el editor enriquecido del main window o usa fallback plano."""
        window = self.window()
        prompt = getattr(window, "_open_rich_text_value_dialog", None)
        if callable(prompt):
            return prompt(title=title, label=label, initial_text=initial_text)
        value, ok = QInputDialog.getText(self, title, label, text=str(initial_text or ""))
        return str(value or ""), bool(ok)

    def _prompt_semantic_diagram_title(self, item: CompositeDiagramItem) -> None:
        """Solicita un nuevo título para el diagrama semántico."""
        value, ok = self._prompt_rich_text_value(
            "Título del diagrama",
            "Título:",
            str(item.semantic_diagram.title or ""),
        )
        if not ok:
            return
        if item.set_diagram_title(str(value or "")):
            self._sync_selection_from_scene()

    def _prompt_semantic_diagram_lane_title(
        self,
        item: CompositeDiagramItem,
        lane_id: Optional[str] = None,
    ) -> None:
        """Edita la etiqueta de un carril del diagrama semántico."""
        target_lane_id = str(lane_id or "")
        lane = None
        if target_lane_id:
            lane = next(
                (candidate for candidate in item.semantic_diagram.lanes if candidate.id == target_lane_id),
                None,
            )
        if lane is None:
            lanes = [candidate for candidate in item.semantic_diagram.lanes]
            if not lanes:
                return
            lane_labels = [
                self._plain_text_from_markup(candidate.title or candidate.id or f"Carril {index + 1}")
                for index, candidate in enumerate(lanes)
            ]
            selected_label, ok = QInputDialog.getItem(
                self,
                "Etiqueta de carril",
                "Carril:",
                lane_labels,
                0,
                False,
            )
            if not ok:
                return
            chosen_index = lane_labels.index(selected_label)
            lane = lanes[chosen_index]
        value, ok = self._prompt_rich_text_value(
            "Etiqueta de carril",
            "Etiqueta:",
            str(lane.title or ""),
        )
        if not ok:
            return
        if item.set_lane_title(lane.id, str(value or "")):
            self._sync_selection_from_scene()

    def _prompt_semantic_diagram_level_label(
        self,
        item: CompositeDiagramItem,
        level_id: Optional[str] = None,
    ) -> None:
        """Edita la etiqueta de un nivel del diagrama semántico."""
        target_level_id = str(level_id or "")
        level = None
        if target_level_id:
            level = next(
                (candidate for candidate in item.semantic_diagram.levels if candidate.id == target_level_id),
                None,
            )
        if level is None:
            levels = [candidate for candidate in item.semantic_diagram.levels]
            if not levels:
                return
            level_labels = [
                self._plain_text_from_markup(candidate.label or candidate.id or f"Nivel {index + 1}")
                for index, candidate in enumerate(levels)
            ]
            selected_label, ok = QInputDialog.getItem(
                self,
                "Etiqueta de nivel",
                "Nivel:",
                level_labels,
                0,
                False,
            )
            if not ok:
                return
            chosen_index = level_labels.index(selected_label)
            level = levels[chosen_index]
        value, ok = self._prompt_rich_text_value(
            "Etiqueta de nivel",
            "Etiqueta:",
            str(level.label or ""),
        )
        if not ok:
            return
        if item.set_level_label(level.id, str(value or "")):
            self._sync_selection_from_scene()

    def _push_energy_diagram_config_change(
        self,
        updates: dict[EnergyDiagramItem, dict[str, object]],
        *,
        text: str,
    ) -> bool:
        """Empaqueta cambios de configuración de diagramas en undo/redo."""
        before: dict[EnergyDiagramItem, dict[str, object]] = {}
        after: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item, payload in updates.items():
            if item.scene() is not self.scene:
                continue
            current = item.config_payload()
            if current == payload:
                continue
            before[item] = current
            after[item] = payload
        if not after:
            return False
        self.undo_stack.push(ConfigureEnergyDiagramItemsCommand(self, before, after, text))
        return True

    def _prompt_energy_diagram_label(self, item: EnergyDiagramItem) -> None:
        """Solicita una nueva etiqueta para un diagrama."""
        if not item.supports_free_label():
            return
        value, ok = self._prompt_rich_text_value(
            "Etiqueta de diagrama",
            "Etiqueta:",
            item.label(),
        )
        if not ok:
            return
        payload = item.config_payload()
        payload["label"] = str(value or "")
        self._push_energy_diagram_config_change(
            {item: payload},
            text="Change energy diagram label",
        )

    def _prompt_energy_diagram_occupancies(self, item: EnergyDiagramItem) -> None:
        """Edita la ocupación electrónica caja por caja."""
        value, ok = QInputDialog.getText(
            self,
            "Ocupación electrónica",
            "Ocupaciones separadas por comas (empty, up, down, pair, upup, downdown):",
            text=", ".join(item.occupancies()),
        )
        if not ok:
            return
        payload = item.config_payload()
        payload["occupancies"] = list(
            normalize_energy_occupancies(
                value,
                kind=item.kind(),
                box_count=item.slot_count(),
            )
        )
        self._push_energy_diagram_config_change(
            {item: payload},
            text="Change energy diagram occupancies",
        )

    def _prompt_selected_energy_diagram_box_count(self) -> None:
        """Solicita un nuevo numero de cajas para filas de orbitales."""
        items = [item for item in self._selected_energy_diagram_items() if item.family() == "row"]
        if not items:
            return
        current = int(items[0].slot_count())
        value, ok = QInputDialog.getInt(
            self,
            "Numero de cajas",
            "Numero de cajas:",
            current,
            1,
            20,
        )
        if not ok:
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            payload = item.config_payload()
            payload["slot_count"] = int(value)
            updates[item] = payload
        self._push_energy_diagram_config_change(
            updates,
            text="Change energy diagram box count",
        )

    def _set_selected_energy_diagram_label_side(self, side: str) -> None:
        """Aplica la posición de etiqueta a los diagramas seleccionados."""
        items = self._selected_energy_diagram_items()
        if not items:
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            if not item.supports_free_label():
                continue
            payload = item.config_payload()
            payload["label_side"] = side
            updates[item] = payload
        self._push_energy_diagram_config_change(
            updates,
            text="Change energy diagram label side",
        )

    def _prompt_selected_energy_diagram_color(self, key: str, title: str) -> None:
        """Solicita un color para un canal visual de diagramas seleccionados."""
        items = self._selected_energy_diagram_items()
        if not items:
            return
        initial = QColor("#222222")
        for item in items:
            candidate = QColor(str(item.effective_style().get(key, "")))
            if candidate.isValid():
                initial = candidate
                break
        color = QColorDialog.getColor(initial, self, title)
        if not color.isValid():
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            payload = item.config_payload()
            style_payload = dict(payload.get("style_payload", {}) or {})
            style_payload[key] = color.name(QColor.NameFormat.HexRgb)
            payload["style_payload"] = style_payload
            updates[item] = payload
        self._push_energy_diagram_config_change(updates, text=f"Change {key}")

    def _set_selected_energy_diagram_fill_visible(self, visible: bool) -> None:
        """Muestra u oculta el fondo de los diagramas seleccionados."""
        items = self._selected_energy_diagram_items()
        if not items:
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            payload = item.config_payload()
            style_payload = dict(payload.get("style_payload", {}) or {})
            style_payload["fill_visible"] = bool(visible)
            payload["style_payload"] = style_payload
            updates[item] = payload
        self._push_energy_diagram_config_change(
            updates,
            text="Toggle energy diagram fill",
        )

    def _set_selected_energy_diagram_box_stroke_visible(self, visible: bool) -> None:
        """Muestra u oculta el contorno de cajas/niveles editables."""
        items = self._selected_energy_diagram_items()
        if not items:
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            payload = item.config_payload()
            style_payload = dict(payload.get("style_payload", {}) or {})
            style_payload["box_stroke_visible"] = bool(visible)
            if visible:
                style_payload["box_base_visible"] = False
            payload["style_payload"] = style_payload
            updates[item] = payload
        self._push_energy_diagram_config_change(
            updates,
            text="Toggle energy diagram box outline",
        )

    def _set_selected_energy_diagram_box_base_visible(self, visible: bool) -> None:
        """Muestra u oculta la base horizontal de las cajas."""
        items = self._selected_energy_diagram_items()
        if not items:
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            payload = item.config_payload()
            style_payload = dict(payload.get("style_payload", {}) or {})
            style_payload["box_base_visible"] = bool(visible)
            if visible:
                style_payload["box_stroke_visible"] = False
            payload["style_payload"] = style_payload
            updates[item] = payload
        self._push_energy_diagram_config_change(
            updates,
            text="Toggle energy diagram box base",
        )

    def _reset_selected_energy_diagram_style(self) -> None:
        """Elimina overrides de color/fondo de los diagramas seleccionados."""
        items = self._selected_energy_diagram_items()
        if not items:
            return
        updates: dict[EnergyDiagramItem, dict[str, object]] = {}
        for item in items:
            payload = item.config_payload()
            payload["style_payload"] = {}
            updates[item] = payload
        self._push_energy_diagram_config_change(
            updates,
            text="Reset energy diagram style",
        )
