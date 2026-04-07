from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QLineEdit,
    QMessageBox,
    QSpinBox,
    QWidget,
)

from chemuson.gui.commands import EditSemanticDiagramCommand
from chemuson.gui.composite_diagram_item import CompositeDiagramItem
from chemuson.gui.energy_diagrams import (
    build_atomic_species_diagram,
    build_atomic_subshell_diagram,
    build_diatomic_mo_diagram,
    build_ligand_field_diagram,
)


@dataclass(slots=True)
class SemanticDiagramWorkflowContext:
    parent: QWidget
    canvas: object
    visible_center: Callable[[], object]


class SemanticDiagramWorkflow:
    """Inserción y edición reconstruible de diagramas semánticos."""

    def builder_spec(self, item: CompositeDiagramItem | None) -> tuple[str, dict]:
        if item is None:
            return "", {}
        builder = dict(item.semantic_diagram.metadata.get("builder", {}) or {})
        return str(builder.get("name", "") or ""), dict(builder.get("params", {}) or {})

    def item_payload(self, context: SemanticDiagramWorkflowContext, item: CompositeDiagramItem) -> dict:
        payload = item.to_json()
        opacity_getter = getattr(context.canvas, "item_raw_opacity", None)
        if callable(opacity_getter):
            payload["opacity"] = opacity_getter(item)
        return payload

    def apply_result(
        self,
        context: SemanticDiagramWorkflowContext,
        diagram,
        *,
        existing_item: CompositeDiagramItem | None = None,
    ) -> None:
        if existing_item is None:
            context.canvas.insert_semantic_diagram(diagram, context.visible_center())
            return
        before_payload = self.item_payload(context, existing_item)
        after_payload = dict(before_payload)
        after_payload["semantic_diagram"] = diagram.to_json_dict()
        if before_payload.get("semantic_diagram") == after_payload.get("semantic_diagram"):
            return
        context.canvas.undo_stack.push(
            EditSemanticDiagramCommand(
                context.canvas,
                existing_item,
                before_payload,
                after_payload,
            )
        )

    def insert_preset(self, context: SemanticDiagramWorkflowContext, preset_name: str) -> None:
        context.canvas.insert_semantic_preset(preset_name, context.visible_center())

    def open_atomic_diagram_dialog(
        self,
        context: SemanticDiagramWorkflowContext,
        existing_item: CompositeDiagramItem | None = None,
    ) -> None:
        _builder_name, builder_params = self.builder_spec(existing_item)
        dialog = QDialog(context.parent)
        dialog.setWindowTitle("Atomic diagram")
        layout = QFormLayout(dialog)
        electron_count = QSpinBox(dialog)
        electron_count.setRange(0, 118)
        electron_count.setValue(int(builder_params.get("electron_count", 8) or 8))
        expanded_subshells = QCheckBox("Expanded subshells", dialog)
        expanded_subshells.setChecked(bool(builder_params.get("expanded_subshells", True)))
        layout.addRow("Electron count:", electron_count)
        layout.addRow("", expanded_subshells)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=dialog,
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addRow(buttons)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        diagram = build_atomic_subshell_diagram(
            electron_count.value(),
            title=builder_params.get("title"),
            expanded_subshells=expanded_subshells.isChecked(),
            max_n=int(builder_params.get("max_n", 7) or 7),
        )
        self.apply_result(context, diagram, existing_item=existing_item)

    def open_atomic_species_diagram_dialog(
        self,
        context: SemanticDiagramWorkflowContext,
        existing_item: CompositeDiagramItem | None = None,
    ) -> None:
        _builder_name, builder_params = self.builder_spec(existing_item)
        dialog = QDialog(context.parent)
        dialog.setWindowTitle("Atomic species diagram")
        layout = QFormLayout(dialog)
        symbol = QLineEdit(str(builder_params.get("symbol", "O") or "O"), dialog)
        charge = QSpinBox(dialog)
        charge.setRange(-6, 6)
        charge.setValue(int(builder_params.get("charge", 0) or 0))
        expanded_subshells = QCheckBox("Expanded subshells", dialog)
        expanded_subshells.setChecked(bool(builder_params.get("expanded_subshells", True)))
        use_known_exceptions = QCheckBox("Use known exceptions", dialog)
        use_known_exceptions.setChecked(bool(builder_params.get("use_known_exceptions", True)))
        layout.addRow("Symbol:", symbol)
        layout.addRow("Charge:", charge)
        layout.addRow("", expanded_subshells)
        layout.addRow("", use_known_exceptions)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=dialog,
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addRow(buttons)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            diagram = build_atomic_species_diagram(
                symbol=symbol.text().strip(),
                charge=charge.value(),
                expanded_subshells=expanded_subshells.isChecked(),
                title=builder_params.get("title"),
                use_known_exceptions=use_known_exceptions.isChecked(),
            )
        except ValueError as exc:
            QMessageBox.information(context.parent, "Atomic species diagram", str(exc))
            return
        self.apply_result(context, diagram, existing_item=existing_item)

    def open_diatomic_mo_diagram_dialog(
        self,
        context: SemanticDiagramWorkflowContext,
        existing_item: CompositeDiagramItem | None = None,
    ) -> None:
        _builder_name, builder_params = self.builder_spec(existing_item)
        dialog = QDialog(context.parent)
        dialog.setWindowTitle("Diatomic MO diagram")
        layout = QFormLayout(dialog)
        left_label = QLineEdit(str(builder_params.get("left_label", "A") or "A"), dialog)
        right_label = QLineEdit(str(builder_params.get("right_label", "B") or "B"), dialog)
        total_electrons = QSpinBox(dialog)
        total_electrons.setRange(0, 20)
        total_electrons.setValue(int(builder_params.get("total_electrons", 10) or 10))
        ordering = QComboBox(dialog)
        ordering.addItem("light_2p", "light_2p")
        ordering.addItem("heavy_2p", "heavy_2p")
        ordering_value = str(builder_params.get("ordering", "heavy_2p") or "heavy_2p")
        ordering.setCurrentIndex(0 if ordering_value == "light_2p" else 1)
        include_core_1s = QCheckBox("Include core 1s", dialog)
        include_core_1s.setChecked(bool(builder_params.get("include_core_1s", False)))
        layout.addRow("Left label:", left_label)
        layout.addRow("Right label:", right_label)
        layout.addRow("Total electrons:", total_electrons)
        layout.addRow("Ordering:", ordering)
        layout.addRow("", include_core_1s)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=dialog,
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addRow(buttons)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        diagram = build_diatomic_mo_diagram(
            left_label=left_label.text().strip() or "A",
            right_label=right_label.text().strip() or "B",
            total_electrons=total_electrons.value(),
            ordering=ordering.currentData(),
            include_core_1s=include_core_1s.isChecked(),
            title=builder_params.get("title"),
        )
        self.apply_result(context, diagram, existing_item=existing_item)

    def open_ligand_field_diagram_dialog(
        self,
        context: SemanticDiagramWorkflowContext,
        existing_item: CompositeDiagramItem | None = None,
    ) -> None:
        _builder_name, builder_params = self.builder_spec(existing_item)
        dialog = QDialog(context.parent)
        dialog.setWindowTitle("Ligand field diagram")
        layout = QFormLayout(dialog)
        d_electrons = QSpinBox(dialog)
        d_electrons.setRange(0, 10)
        d_electrons.setValue(int(builder_params.get("d_electrons", 6) or 6))
        geometry = QComboBox(dialog)
        geometry.addItem("octahedral", "octahedral")
        geometry.addItem("tetrahedral", "tetrahedral")
        geometry.addItem("square_planar", "square_planar")
        geometry_value = str(builder_params.get("geometry", "octahedral") or "octahedral")
        geometry.setCurrentIndex(max(0, geometry.findData(geometry_value)))
        spin_mode = QComboBox(dialog)
        spin_mode.addItem("high", "high")
        spin_mode.addItem("low", "low")
        spin_value = str(builder_params.get("spin_mode", "high") or "high")
        spin_mode.setCurrentIndex(0 if spin_value == "high" else 1)
        layout.addRow("d electron count:", d_electrons)
        layout.addRow("Geometry:", geometry)
        layout.addRow("Spin mode:", spin_mode)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=dialog,
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addRow(buttons)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        diagram = build_ligand_field_diagram(
            d_electrons=d_electrons.value(),
            geometry=geometry.currentData(),
            spin_mode=spin_mode.currentData(),
            title=builder_params.get("title"),
        )
        self.apply_result(context, diagram, existing_item=existing_item)

    def edit_selected_semantic_diagram(self, context: SemanticDiagramWorkflowContext) -> None:
        item = context.canvas.selected_semantic_diagram_item()
        if item is None:
            return
        builder_name, _builder_params = self.builder_spec(item)
        if not builder_name:
            QMessageBox.information(
                context.parent,
                "Edit Electronic Diagram",
                "This diagram can be edited only at label/occupancy level because its original builder parameters are not available.",
            )
            return
        if builder_name == "build_atomic_subshell_diagram":
            self.open_atomic_diagram_dialog(context, existing_item=item)
            return
        if builder_name == "build_atomic_species_diagram":
            self.open_atomic_species_diagram_dialog(context, existing_item=item)
            return
        if builder_name == "build_diatomic_mo_diagram":
            self.open_diatomic_mo_diagram_dialog(context, existing_item=item)
            return
        if builder_name == "build_ligand_field_diagram":
            self.open_ligand_field_diagram_dialog(context, existing_item=item)
            return
        QMessageBox.information(
            context.parent,
            "Edit Electronic Diagram",
            "This electronic diagram builder is not available for full reconstruction.",
        )
