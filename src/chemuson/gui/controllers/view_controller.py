from __future__ import annotations

from dataclasses import dataclass
from dataclasses import replace

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QVBoxLayout,
)


@dataclass(frozen=True, slots=True)
class CanvasSizePreset:
    label: str
    width: int
    height: int


class ViewController:
    """Coordina estado de vista, numeración y tamaño del lienzo."""

    LETTER_PORTRAIT = CanvasSizePreset("Carta (vertical)", 816, 1056)
    LETTER_LANDSCAPE = CanvasSizePreset("Carta (horizontal)", 1056, 816)
    A4_PORTRAIT = CanvasSizePreset("A4 (vertical)", 794, 1123)
    A4_LANDSCAPE = CanvasSizePreset("A4 (horizontal)", 1123, 794)
    A3_PORTRAIT = CanvasSizePreset("A3 (vertical)", 1123, 1587)
    A3_LANDSCAPE = CanvasSizePreset("A3 (horizontal)", 1587, 1123)

    def load_numbering_preferences(self, window) -> None:
        mode = str(window._settings.value("numbering/mode", "atoms") or "atoms").strip().lower()
        if mode not in {"atoms", "structures", "both"}:
            mode = "atoms"
        include_export = window._setting_bool(
            window._settings.value("numbering/include_export", True),
            True,
        )
        window._numbering_default_mode = mode
        window._numbering_default_include_export = bool(include_export)
        self.apply_default_numbering_to_canvas(window, window.canvas)

    def save_numbering_preferences(self, window) -> None:
        window._settings.remove("numbering/enabled")
        window._numbering_default_mode = str(window.canvas.state.numbering_mode)
        window._numbering_default_include_export = bool(
            window.canvas.state.numbering_include_in_export
        )
        window._settings.setValue("numbering/mode", str(window._numbering_default_mode))
        window._settings.setValue(
            "numbering/include_export",
            bool(window._numbering_default_include_export),
        )

    def sync_numbering_actions(self, window) -> None:
        mode = str(window.canvas.state.numbering_mode or "").strip().lower()
        if mode not in {"atoms", "structures", "both"}:
            mode = "atoms"
        window.canvas.state.numbering_mode = mode
        state_to_action = {
            "atoms": window.action_numbering_mode_atoms,
            "structures": window.action_numbering_mode_structures,
            "both": window.action_numbering_mode_both,
        }
        actions = [
            window.action_numbering_enabled,
            window.action_numbering_mode_atoms,
            window.action_numbering_mode_structures,
            window.action_numbering_mode_both,
            window.action_numbering_export,
        ]
        previous_blocks = []
        for action in actions:
            previous_blocks.append((action, action.blockSignals(True)))
        try:
            window.action_numbering_enabled.setChecked(bool(window.canvas.state.numbering_enabled))
            for key, action in state_to_action.items():
                action.setChecked(key == mode)
            window.action_numbering_export.setChecked(
                bool(window.canvas.state.numbering_include_in_export)
            )
            window.action_numbering_recalculate.setEnabled(
                bool(window.canvas.state.numbering_enabled)
            )
        finally:
            for action, was_blocked in previous_blocks:
                action.blockSignals(was_blocked)

    def apply_default_numbering_to_canvas(self, window, canvas) -> None:
        canvas.set_numbering_mode(window._numbering_default_mode)
        canvas.set_numbering_enabled(False)
        canvas.set_numbering_include_in_export(window._numbering_default_include_export)

    def set_canvas_size(self, window, width: int, height: int) -> None:
        window.canvas.set_paper_size(width, height)
        window.statusBar().showMessage(f"Lienzo: {width} x {height} px")

    def on_toggle_numbering(self, window, checked: bool) -> None:
        window.canvas.set_numbering_enabled(checked)
        self.sync_numbering_actions(window)
        self.save_numbering_preferences(window)
        window.statusBar().showMessage("Numeración: visible" if checked else "Numeración: oculta")

    def on_set_numbering_mode(self, window, mode: str) -> None:
        window.canvas.set_numbering_mode(mode)
        self.sync_numbering_actions(window)
        self.save_numbering_preferences(window)
        labels = {
            "atoms": "átomos",
            "structures": "estructuras",
            "both": "átomos y estructuras",
        }
        active_mode = window.canvas.state.numbering_mode
        window.statusBar().showMessage(
            f"Numeración: {labels.get(active_mode, active_mode)}"
        )

    def on_recalculate_numbering(self, window) -> None:
        window.canvas.recompute_numbering(force_reset=True)
        window.statusBar().showMessage("Numeración recalculada")

    def on_toggle_numbering_export(self, window, checked: bool) -> None:
        window.canvas.set_numbering_include_in_export(checked)
        self.sync_numbering_actions(window)
        self.save_numbering_preferences(window)
        window.statusBar().showMessage(
            "Exportación: numeración incluida"
            if checked
            else "Exportación: numeración excluida"
        )

    def on_toggle_carbons(self, window, checked: bool) -> None:
        window.canvas.state.show_implicit_carbons = checked
        window.canvas.refresh_atom_visibility()
        window.statusBar().showMessage("Carbonos: visibles" if checked else "Carbonos: ocultos")

    def on_toggle_hydrogens(self, window, checked: bool) -> None:
        window.canvas.state.show_implicit_hydrogens = checked
        window.canvas.refresh_atom_visibility()
        window.statusBar().showMessage(
            "Hidrógenos: visibles" if checked else "Hidrógenos: ocultos"
        )

    def on_toggle_aromatic_circles(self, window, checked: bool) -> None:
        window.canvas.state.use_aromatic_circles = checked
        window.canvas.refresh_aromatic_circles()
        window.statusBar().showMessage(
            "Aromáticos: círculos" if checked else "Aromáticos: Kekulé"
        )

    def on_toggle_rules(self, window, checked: bool) -> None:
        window.canvas.set_show_rulers(checked)
        window.statusBar().showMessage("Reglas: visibles" if checked else "Reglas: ocultas")

    def on_toggle_crosshair(self, window, checked: bool) -> None:
        window.canvas.set_show_grid(checked)
        window.statusBar().showMessage(
            "Cuadrícula: visible" if checked else "Cuadrícula: oculta"
        )

    def on_zoom_in(self, window) -> None:
        window.canvas.zoom_in()

    def on_zoom_out(self, window) -> None:
        window.canvas.zoom_out()

    def on_zoom_reset(self, window) -> None:
        window.canvas.resetTransform()
        window.canvas.center_on_paper()
        window.statusBar().showMessage("Zoom: 100%")

    def open_canvas_custom_size_dialog(self, window) -> None:
        dialog = QDialog(window)
        dialog.setWindowTitle("Tamaño de lienzo")

        px_per_in = 96.0
        px_per_cm = px_per_in / 2.54

        width_spin = QDoubleSpinBox()
        height_spin = QDoubleSpinBox()
        width_spin.setDecimals(2)
        height_spin.setDecimals(2)

        unit_combo = QComboBox()
        unit_combo.addItems(["cm", "px", "in"])
        unit_combo.setCurrentText("cm")

        def apply_unit_settings(unit: str) -> None:
            if unit == "px":
                width_spin.setRange(200.0, 20000.0)
                height_spin.setRange(200.0, 20000.0)
                width_spin.setDecimals(0)
                height_spin.setDecimals(0)
            elif unit == "in":
                width_spin.setRange(1.0, 200.0)
                height_spin.setRange(1.0, 200.0)
                width_spin.setDecimals(2)
                height_spin.setDecimals(2)
            else:
                width_spin.setRange(1.0, 500.0)
                height_spin.setRange(1.0, 500.0)
                width_spin.setDecimals(2)
                height_spin.setDecimals(2)

        def to_px(value: float, unit: str) -> float:
            if unit == "px":
                return value
            if unit == "in":
                return value * px_per_in
            return value * px_per_cm

        def from_px(value: float, unit: str) -> float:
            if unit == "px":
                return value
            if unit == "in":
                return value / px_per_in
            return value / px_per_cm

        def on_unit_changed(text: str) -> None:
            old_unit = str(unit_combo.property("last_unit") or "cm")
            unit_combo.setProperty("last_unit", text)
            current_px_w = to_px(width_spin.value(), old_unit)
            current_px_h = to_px(height_spin.value(), old_unit)
            apply_unit_settings(text)
            width_spin.setValue(from_px(current_px_w, text))
            height_spin.setValue(from_px(current_px_h, text))

        apply_unit_settings("cm")
        width_spin.setValue(window.canvas.paper_width / px_per_cm)
        height_spin.setValue(window.canvas.paper_height / px_per_cm)
        unit_combo.setProperty("last_unit", "cm")
        unit_combo.currentTextChanged.connect(on_unit_changed)

        form = QFormLayout()
        form.addRow("Unidad:", unit_combo)
        form.addRow("Ancho:", width_spin)
        form.addRow("Alto:", height_spin)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)

        layout = QVBoxLayout(dialog)
        layout.addLayout(form)
        layout.addWidget(buttons)

        if dialog.exec() != QDialog.DialogCode.Accepted:
            return

        unit = unit_combo.currentText()
        if unit == "px":
            width_px = int(width_spin.value())
            height_px = int(height_spin.value())
        elif unit == "in":
            width_px = int(round(width_spin.value() * px_per_in))
            height_px = int(round(height_spin.value() * px_per_in))
        else:
            width_px = int(round(width_spin.value() * px_per_cm))
            height_px = int(round(height_spin.value() * px_per_cm))

        self.set_canvas_size(window, width_px, height_px)

    def apply_bond_caps(self, window, bond_caps: str) -> None:
        if bond_caps == "round":
            cap_style = Qt.PenCapStyle.RoundCap
            join_style = Qt.PenJoinStyle.RoundJoin
        else:
            cap_style = Qt.PenCapStyle.FlatCap
            join_style = Qt.PenJoinStyle.MiterJoin
        if (
            window.canvas.drawing_style.cap_style != cap_style
            or window.canvas.drawing_style.join_style != join_style
        ):
            style = replace(
                window.canvas.drawing_style,
                cap_style=cap_style,
                join_style=join_style,
            )
            window.canvas.apply_drawing_style(style)

    def apply_appearance_settings(self, window, prefs: dict) -> None:
        bond_caps = prefs.get("bond_caps")
        if bond_caps:
            self.apply_bond_caps(window, bond_caps)

    def update_total_charge_indicator(self, window) -> None:
        charge = int(window.canvas.model.total_formal_charge())
        charge_text = f"+{charge}" if charge > 0 else str(charge)
        window._total_charge_label.setText(f"Carga total: {charge_text}")
