"""Qt interface for the isolated Molecular Depiction Visual Lab."""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

from PyQt6.QtCore import QSignalBlocker, Qt
from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSlider,
    QSpinBox,
    QSplitter,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from chemuson.gui.theme import apply_theme

from .diagnostics import GeometryDiagnostics
from .exporting import export_case_bundle, export_gallery
from .parameters import (
    CAP_STYLES,
    JOIN_STYLES,
    RANGES,
    DepictionStyleParameters,
    ParameterValidationError,
)
from .named_presets import EDITORIAL_AB3, get_named_preset
from .profiles import (
    ProfileError,
    compare_profiles,
    load_profile,
    profile_directory,
    safe_profile_filename,
    save_profile,
)
from .renderer import ProductionCanvasAdapter
from .scenes import CASES, case_graph, graph_snapshot


@dataclass(frozen=True)
class NumericControlSpec:
    key: str
    label: str
    tooltip: str
    minimum: float
    maximum: float
    optional: bool = False
    integer: bool = False


NUMERIC_CONTROLS = (
    NumericControlSpec("bond_stroke_ratio", "Grosor de trazo", "Grosor base de enlaces simples y trazos de la cuña hash, como fracción de la longitud de enlace de referencia.", *RANGES["bond_stroke_ratio"]),
    NumericControlSpec("double_bond_offset_ratio", "Separación del doble", "Offset del trazo pi respecto del sigma para enlaces dobles plain.", *RANGES["double_bond_offset_ratio"]),
    NumericControlSpec("double_bond_inner_trim_ratio", "Recorte interno doble", "Longitud recortada de los extremos de la línea pi doble no simétrica.", *RANGES["double_bond_inner_trim_ratio"]),
    NumericControlSpec("bold_double_bond_offset_ratio", "Separación doble bold", "Override opcional del offset de dobles bold; Auto conserva la fórmula actual.", *RANGES["bold_double_bond_offset_ratio"], optional=True),
    NumericControlSpec("triple_bond_offset_ratio", "Separación triple", "Override opcional de separación de las líneas de enlaces triples; Auto sigue el offset doble actual.", *RANGES["triple_bond_offset_ratio"], optional=True),
    NumericControlSpec("solid_wedge_width_ratio", "Anchura de cuña sólida", "Anchura máxima de la cuña sólida en el extremo ancho, relativa al enlace de referencia.", *RANGES["solid_wedge_width_ratio"]),
    NumericControlSpec("solid_wedge_tip_width_ratio", "Anchura del extremo estrecho", "Override opcional del ancho total de la punta en una unión con vecinos; Auto conserva su geometría contextual.", *RANGES["solid_wedge_tip_width_ratio"], optional=True),
    NumericControlSpec("hashed_wedge_count", "Número de rayas", "Cantidad real de trazos transversales de la cuña discontinua.", 3, 12, integer=True),
    NumericControlSpec("hashed_wedge_min_width_ratio", "Ancho mínimo de raya", "Ancho de la raya más estrecha al inicio de la cuña discontinua.", *RANGES["hashed_wedge_min_width_ratio"]),
    NumericControlSpec("hashed_wedge_max_width_ratio", "Ancho máximo de raya", "Ancho de la raya más ancha al final de la cuña discontinua.", *RANGES["hashed_wedge_max_width_ratio"]),
    NumericControlSpec("hashed_wedge_stroke_ratio", "Grosor de rayas", "Grosor Qt de las rayas hash; es distinto de su ancho longitudinal.", *RANGES["hashed_wedge_stroke_ratio"]),
    NumericControlSpec("hashed_wedge_terminal_shift_ratio", "Desplazamiento terminal hash", "Override opcional de la última raya junto a una unión sigma/pi; Auto conserva el desplazamiento actual.", *RANGES["hashed_wedge_terminal_shift_ratio"], optional=True),
    NumericControlSpec("label_clearance_ratio", "Separación enlace–etiqueta", "Override común y relativo de clearance medido desde el contorno visible del texto; Auto conserva la política por tipo actual.", *RANGES["label_clearance_ratio"], optional=True),
)

OPTIONAL_FALLBACKS = {
    "bold_double_bond_offset_ratio": lambda params: max(
        params.double_bond_offset_ratio,
        max(params.bond_stroke_ratio * 2.2, params.bond_stroke_ratio + 0.025) * 0.9,
    ),
    "triple_bond_offset_ratio": lambda params: params.double_bond_offset_ratio,
    "solid_wedge_tip_width_ratio": lambda _params: 0.24,
    "hashed_wedge_terminal_shift_ratio": lambda _params: 0.25,
    "label_clearance_ratio": lambda _params: 0.04,
}


class NumericControl(QWidget):
    """Accessible slider/spinbox pair with optional explicit override mode."""

    def __init__(self, spec: NumericControlSpec, parent=None) -> None:
        super().__init__(parent)
        self.spec = spec
        self._callback = None
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        self.label = QLabel(spec.label)
        self.label.setToolTip(spec.tooltip)
        row.addWidget(self.label, 2)
        self.override = None
        if spec.optional:
            self.override = QCheckBox("Override")
            self.override.setToolTip(spec.tooltip)
            row.addWidget(self.override)
        if spec.integer:
            self.slider = QSlider(Qt.Orientation.Horizontal)
            self.slider.setRange(int(spec.minimum), int(spec.maximum))
            self.spin = QSpinBox()
            self.spin.setRange(int(spec.minimum), int(spec.maximum))
            self.spin.setFixedWidth(78)
        else:
            self.slider = QSlider(Qt.Orientation.Horizontal)
            self.slider.setRange(0, 1000)
            self.spin = QDoubleSpinBox()
            self.spin.setRange(spec.minimum, spec.maximum)
            self.spin.setDecimals(3)
            self.spin.setSingleStep(0.005)
            self.spin.setFixedWidth(90)
        self.slider.setToolTip(spec.tooltip)
        self.spin.setToolTip(spec.tooltip)
        row.addWidget(self.slider, 4)
        row.addWidget(self.spin)
        self.reset_button = QToolButton()
        self.reset_button.setText("↺")
        self.reset_button.setToolTip(f"Restablecer: {spec.label}")
        self.reset_button.setAccessibleName(f"Restablecer {spec.label}")
        row.addWidget(self.reset_button)
        self.slider.valueChanged.connect(self._slider_changed)
        self.spin.valueChanged.connect(self._spin_changed)
        if self.override is not None:
            self.override.toggled.connect(self._override_changed)
        self._set_enabled(True)

    def bind(self, callback, reset_callback) -> None:
        self._callback = callback
        self.reset_button.clicked.connect(reset_callback)

    def _from_slider(self, raw: int) -> float | int:
        if self.spec.integer:
            return int(raw)
        fraction = raw / 1000.0
        return round(self.spec.minimum + fraction * (self.spec.maximum - self.spec.minimum), 3)

    def _to_slider(self, value: float | int) -> int:
        if self.spec.integer:
            return int(value)
        span = self.spec.maximum - self.spec.minimum
        if span <= 0:
            return 0
        return max(0, min(1000, round((float(value) - self.spec.minimum) / span * 1000)))

    def set_value(
        self,
        value: float | int | None,
        *,
        parameters: DepictionStyleParameters,
        optional_enabled: bool = True,
    ) -> None:
        is_override = value is not None
        shown = value
        if shown is None:
            shown = OPTIONAL_FALLBACKS[self.spec.key](parameters)
        if self.override is not None:
            blocker = QSignalBlocker(self.override)
            self.override.setChecked(bool(optional_enabled and is_override))
            del blocker
        if isinstance(self.spin, QDoubleSpinBox):
            spin_blocker = QSignalBlocker(self.spin)
            self.spin.setValue(float(shown))
        else:
            spin_blocker = QSignalBlocker(self.spin)
            self.spin.setValue(int(shown))
        slider_blocker = QSignalBlocker(self.slider)
        self.slider.setValue(self._to_slider(shown))
        del spin_blocker, slider_blocker
        self._set_enabled(self.override is None or self.override.isChecked())

    def _set_enabled(self, enabled: bool) -> None:
        self.slider.setEnabled(enabled)
        self.spin.setEnabled(enabled)

    def _emit(self, value: float | int) -> None:
        if self._callback is not None:
            self._callback(self.spec.key, value)

    def _slider_changed(self, value: int) -> None:
        converted = self._from_slider(value)
        blocker = QSignalBlocker(self.spin)
        self.spin.setValue(converted)
        del blocker
        if self.override is None or self.override.isChecked():
            self._emit(converted)

    def _spin_changed(self, value: float | int) -> None:
        blocker = QSignalBlocker(self.slider)
        self.slider.setValue(self._to_slider(value))
        del blocker
        if self.override is None or self.override.isChecked():
            self._emit(value)

    def _override_changed(self, enabled: bool) -> None:
        self._set_enabled(enabled)
        if self._callback is not None:
            self._callback(self.spec.key, self.spin.value() if enabled else None)


class ChoiceControl(QWidget):
    def __init__(self, label: str, options: dict[str, object], parent=None) -> None:
        super().__init__(parent)
        self.options = options
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        self.label = QLabel(label)
        row.addWidget(self.label, 2)
        self.combo = QComboBox()
        for key in options:
            self.combo.addItem(key.capitalize(), key)
        self.combo.setToolTip(f"Estilo de terminación {label.lower()}, aplicado por QPen en el renderer real.")
        row.addWidget(self.combo, 4)
        self.reset_button = QToolButton()
        self.reset_button.setText("↺")
        self.reset_button.setToolTip(f"Restablecer: {label}")
        row.addWidget(self.reset_button)

    def set_value(self, value: str) -> None:
        index = self.combo.findData(value)
        if index >= 0:
            blocker = QSignalBlocker(self.combo)
            self.combo.setCurrentIndex(index)
            del blocker


class RenderPanel(QWidget):
    def __init__(self, heading: str, parent=None) -> None:
        super().__init__(parent)
        self.heading = QLabel(heading)
        self.heading.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.heading.setStyleSheet("font-weight: 700; padding: 4px;")
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(3, 3, 3, 3)
        self.layout.addWidget(self.heading)
        self.adapter: ProductionCanvasAdapter | None = None
        self.diagnostics: GeometryDiagnostics | None = None
        self._diagnostics_enabled = False

    def set_case(
        self,
        graph,
        parameters: DepictionStyleParameters,
        theme: str,
    ) -> None:
        old_canvas = self.adapter.canvas if self.adapter is not None else None
        if old_canvas is not None:
            self.layout.removeWidget(old_canvas)
            old_canvas.close()
            old_canvas.deleteLater()
        self.adapter = ProductionCanvasAdapter.create(
            graph, parameters=parameters, theme=theme, parent=self
        )
        self.adapter.canvas.setMinimumSize(380, 300)
        self.layout.addWidget(self.adapter.canvas, 1)
        self.diagnostics = GeometryDiagnostics(self.adapter)
        if self._diagnostics_enabled:
            self.diagnostics.show()

    def set_diagnostics(self, enabled: bool) -> None:
        self._diagnostics_enabled = bool(enabled)
        if self.diagnostics is not None:
            self.diagnostics.set_enabled(enabled)


class DepictionLabWindow(QMainWindow):
    """Two-panel interactive parameter calibration interface."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("ChemUSON Molecular Depiction Visual Lab")
        self.resize(1480, 980)
        self.current_parameters = DepictionStyleParameters.current()
        self.working_parameters = DepictionStyleParameters.current()
        self.working_name = "Working"
        self.working_based_on = "current-renderer"
        self.case_id = CASES[0].case_id
        self.theme = "light"
        self.zoom = 1.0
        self.angle = 0.0
        self.current_panel = RenderPanel("CURRENT · renderer actual")
        self.working_panel = RenderPanel("WORKING · perfil experimental")

        root = QWidget()
        main = QVBoxLayout(root)
        main.setContentsMargins(10, 8, 10, 8)
        main.addWidget(self._build_toolbar())
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self.current_panel)
        splitter.addWidget(self.working_panel)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)
        main.addWidget(splitter, 5)
        self.controls = self._build_controls()
        main.addWidget(self.controls, 3)
        main.addWidget(self._build_actions())
        self.setCentralWidget(root)
        self._install_initial_case()
        self._connect_scrollbars()
        self._refresh_parameter_controls()
        apply_theme(self, self.theme)

    def _build_toolbar(self) -> QWidget:
        bar = QWidget()
        row = QHBoxLayout(bar)
        row.setContentsMargins(0, 0, 0, 4)
        row.addWidget(QLabel("Caso:"))
        self.case_combo = QComboBox()
        for case in CASES:
            case_label = f"{case.case_id} · {case.title}"
            if case.case_id == "B4":
                case_label += " [GEOMETRIC FIXTURE · NON-EDITORIAL · EXCLUDED FROM CALIBRATION]"
            self.case_combo.addItem(case_label, case.case_id)
        self.case_combo.currentIndexChanged.connect(self._case_selected)
        row.addWidget(self.case_combo, 3)
        row.addWidget(QLabel("Tema:"))
        self.theme_combo = QComboBox()
        self.theme_combo.addItem("Claro", "light")
        self.theme_combo.addItem("Oscuro", "dark")
        self.theme_combo.currentIndexChanged.connect(self._theme_selected)
        row.addWidget(self.theme_combo)
        row.addWidget(QLabel("Zoom:"))
        self.zoom_combo = QComboBox()
        for percent in (50, 75, 100, 150, 200, 300):
            self.zoom_combo.addItem(f"{percent}%", percent / 100.0)
        self.zoom_combo.setCurrentIndex(2)
        self.zoom_combo.currentIndexChanged.connect(self._zoom_selected)
        row.addWidget(self.zoom_combo)
        row.addWidget(QLabel("Rotación visual:"))
        self.angle_combo = QComboBox()
        for angle in (0, 30, 45, 60, 90, 120, 150, 180):
            self.angle_combo.addItem(f"{angle}°", float(angle))
        self.angle_combo.currentIndexChanged.connect(self._angle_selected)
        row.addWidget(self.angle_combo)
        self.diagnostics_toggle = QCheckBox("Diagnóstico")
        self.diagnostics_toggle.setToolTip("Muestra centros, endpoints, contornos, cajas de etiquetas y cruces geométricos; se excluye de exportaciones normales.")
        self.diagnostics_toggle.toggled.connect(self._diagnostics_toggled)
        row.addWidget(self.diagnostics_toggle)
        self.hide_controls_button = QToolButton()
        self.hide_controls_button.setText("Ocultar controles")
        self.hide_controls_button.setCheckable(True)
        self.hide_controls_button.toggled.connect(self._controls_visibility_changed)
        row.addWidget(self.hide_controls_button)
        return bar

    def _build_controls(self) -> QWidget:
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        body = QWidget()
        layout = QVBoxLayout(body)
        self.parameter_controls: dict[str, NumericControl] = {}
        self.choice_controls: dict[str, ChoiceControl] = {}
        groups = (
            ("A · Trazos y enlaces múltiples", NUMERIC_CONTROLS[:5]),
            ("B · Cuñas sólidas y discontinuas", NUMERIC_CONTROLS[5:12]),
            ("C · Etiquetas químicas", NUMERIC_CONTROLS[12:]),
        )
        for title, specs in groups:
            group = QGroupBox(title)
            form = QFormLayout(group)
            for spec in specs:
                control = NumericControl(spec)
                control.bind(
                    self._parameter_changed,
                    lambda checked=False, key=spec.key: self._reset_parameter(key),
                )
                self.parameter_controls[spec.key] = control
                form.addRow(control)
            layout.addWidget(group)
        line_group = QGroupBox("D · Terminaciones y uniones")
        line_form = QFormLayout(line_group)
        for key, label, options in (
            ("cap_style", "Terminación de línea", CAP_STYLES),
            ("join_style", "Unión de línea", JOIN_STYLES),
        ):
            control = ChoiceControl(label, options)
            control.combo.currentIndexChanged.connect(
                lambda _index, name=key, combo=control.combo: self._choice_changed(name, combo.currentData())
            )
            control.reset_button.clicked.connect(lambda _checked=False, name=key: self._reset_parameter(name))
            self.choice_controls[key] = control
            line_form.addRow(control)
        layout.addWidget(line_group)
        unavailable = QGroupBox("No disponibles todavía · controles intencionalmente omitidos")
        unavailable_note = QLabel(
            "Punta estrecha contextual de cuña; factores internos de miter/fork y límite de engrosamiento en bifurcaciones; "
            "offsets ópticos por elemento/dirección y side bearings tipográficos; trims de anillos dependientes de topología; "
            "orientación manual sigma/pi; antialiasing por panel. Estas propiedades no se aíslan de forma segura en esta fase, "
            "por lo que no se presentan como sliders. La caja de diagnóstico marca únicamente intersecciones geométricas Qt, "
            "no una puntuación estética."
        )
        unavailable_note.setWordWrap(True)
        unavailable.setToolTip(unavailable_note.text())
        unavailable_layout = QVBoxLayout(unavailable)
        unavailable_layout.addWidget(unavailable_note)
        layout.addWidget(unavailable)
        layout.addStretch(1)
        scroll.setWidget(body)
        return scroll

    def _build_actions(self) -> QWidget:
        bar = QWidget()
        layout = QHBoxLayout(bar)
        layout.setContentsMargins(0, 2, 0, 0)
        self.reset_all_button = QPushButton("Restablecer todo")
        self.reset_all_button.clicked.connect(self._reset_all)
        layout.addWidget(self.reset_all_button)
        self.editorial_preset_button = QPushButton("Cargar ChemUSON Editorial AB3")
        self.editorial_preset_button.setToolTip(
            "Carga el perfil experimental aprobado para A1–B3; B4 queda excluido y no cambia los defaults de producción."
        )
        self.editorial_preset_button.clicked.connect(self._load_editorial_ab3_preset)
        layout.addWidget(self.editorial_preset_button)
        layout.addWidget(QLabel("Perfil Working:"))
        self.profile_name = QLineEdit("Working")
        self.profile_name.setMaximumWidth(220)
        self.profile_name.setToolTip("Nombre del perfil experimental guardado sólo en la carpeta del laboratorio.")
        self.profile_name.textChanged.connect(self._profile_name_changed)
        layout.addWidget(self.profile_name)
        save_button = QPushButton("Guardar perfil")
        save_button.clicked.connect(self._save_profile)
        layout.addWidget(save_button)
        load_button = QPushButton("Abrir JSON…")
        load_button.clicked.connect(self._open_profile)
        layout.addWidget(load_button)
        self.compare_button = QPushButton("Comparar perfiles")
        self.compare_button.clicked.connect(self._compare_profiles)
        layout.addWidget(self.compare_button)
        export_button = QPushButton("Exportar comparación…")
        export_button.clicked.connect(self._export_case)
        layout.addWidget(export_button)
        gallery_button = QPushButton("Exportar galería…")
        gallery_button.clicked.connect(self._export_gallery)
        layout.addWidget(gallery_button)
        self.status = QLabel("Current es de solo lectura. Working modifica únicamente el panel experimental.")
        self.status.setWordWrap(True)
        layout.addWidget(self.status, 2)
        return bar

    def _install_initial_case(self) -> None:
        graph = case_graph(self.case_id)
        self.current_panel.set_case(graph, self.current_parameters, self.theme)
        self.working_panel.set_case(graph, self.working_parameters, self.theme)
        self._apply_views()
        self._assert_panels_match()

    def _assert_panels_match(self) -> None:
        if self.current_panel.adapter is None or self.working_panel.adapter is None:
            return
        left = graph_snapshot(self.current_panel.adapter.canvas.model)
        right = graph_snapshot(self.working_panel.adapter.canvas.model)
        if left != right:
            raise AssertionError("Current and Working panels do not contain identical molecular data")
        self.current_panel.adapter.assert_unchanged()
        self.working_panel.adapter.assert_unchanged()

    def _connect_scrollbars(self) -> None:
        for orientation in (Qt.Orientation.Horizontal, Qt.Orientation.Vertical):
            left = self.current_panel.adapter.canvas.horizontalScrollBar() if orientation == Qt.Orientation.Horizontal else self.current_panel.adapter.canvas.verticalScrollBar()
            right = self.working_panel.adapter.canvas.horizontalScrollBar() if orientation == Qt.Orientation.Horizontal else self.working_panel.adapter.canvas.verticalScrollBar()
            left.valueChanged.connect(lambda value, target=right: self._copy_scroll_value(target, value))
            right.valueChanged.connect(lambda value, target=left: self._copy_scroll_value(target, value))

    @staticmethod
    def _copy_scroll_value(target, value: int) -> None:
        if target.value() != value:
            target.setValue(value)

    def _apply_views(self) -> None:
        for panel in (self.current_panel, self.working_panel):
            if panel.adapter is not None:
                panel.adapter.set_view(self.zoom, self.angle)

    def _case_selected(self, index: int) -> None:
        case_id = self.case_combo.itemData(index)
        if not case_id or case_id == self.case_id:
            return
        self.case_id = str(case_id)
        graph = case_graph(self.case_id)
        self.current_panel.set_case(graph, self.current_parameters, self.theme)
        self.working_panel.set_case(graph, self.working_parameters, self.theme)
        self.current_panel.set_diagnostics(self.diagnostics_toggle.isChecked())
        self.working_panel.set_diagnostics(self.diagnostics_toggle.isChecked())
        self._apply_views()
        self._connect_scrollbars()
        self._assert_panels_match()
        case = next(item for item in CASES if item.case_id == self.case_id)
        if case.case_id == "B4":
            calibration_note = (
                "geometric fixture; non-editorial; excluded from preset selection / style calibration"
            )
        elif case.editorial_calibration_eligible:
            calibration_note = "editorial calibration sample"
        else:
            calibration_note = case.editorial_exclusion_reason or "diagnostic only"
        self.status.setText(
            f"Caso {self.case_id} · {len(graph.atoms)} átomos · {len(graph.bonds)} enlaces · "
            f"grafo intacto · {calibration_note}."
        )

    def _theme_selected(self, index: int) -> None:
        self.theme = str(self.theme_combo.itemData(index))
        for panel in (self.current_panel, self.working_panel):
            if panel.adapter is not None:
                panel.adapter.set_theme(self.theme)
        apply_theme(self, self.theme)
        self.status.setText(f"Tema {self.theme}; la escala y el grafo permanecen sincronizados.")

    def _zoom_selected(self, index: int) -> None:
        self.zoom = float(self.zoom_combo.itemData(index))
        self._apply_views()

    def _angle_selected(self, index: int) -> None:
        self.angle = float(self.angle_combo.itemData(index))
        self._apply_views()
        self._assert_panels_match()

    def _diagnostics_toggled(self, enabled: bool) -> None:
        self.current_panel.set_diagnostics(enabled)
        self.working_panel.set_diagnostics(enabled)

    def _controls_visibility_changed(self, hidden: bool) -> None:
        self.controls.setVisible(not hidden)
        self.hide_controls_button.setText("Mostrar controles" if hidden else "Ocultar controles")

    def _parameter_changed(self, key: str, value) -> None:
        try:
            candidate = replace(self.working_parameters, **{key: value})
            self.working_parameters = DepictionStyleParameters.from_mapping(candidate.to_mapping())
        except (TypeError, ParameterValidationError) as exc:
            self.status.setText(f"Parámetro rechazado: {exc}")
            self._refresh_parameter_controls()
            return
        if self.working_panel.adapter is not None:
            self.working_panel.adapter.apply_parameters(self.working_parameters)
        self._assert_panels_match()
        self._update_change_summary()

    def _choice_changed(self, key: str, value: str) -> None:
        self._parameter_changed(key, value)

    def _load_editorial_ab3_preset(self) -> None:
        preset = get_named_preset(EDITORIAL_AB3.name)
        self.working_name = preset.name
        self.working_based_on = preset.based_on
        self.profile_name.setText(preset.name)
        self.working_parameters = preset.parameters
        if self.working_panel.adapter is not None:
            self.working_panel.adapter.apply_parameters(preset.parameters)
        self._refresh_parameter_controls()
        self._assert_panels_match()
        self.status.setText(
            "Preset experimental ChemUSON Editorial AB3 cargado en Working; "
            "B4 está excluido de la calibración y los defaults globales no cambian."
        )

    def _reset_parameter(self, key: str) -> None:
        defaults = self.current_parameters.to_mapping()
        self._parameter_changed(key, defaults[key])
        self._refresh_parameter_controls()

    def _reset_all(self) -> None:
        self.working_name = "Working"
        self.working_based_on = "current-renderer"
        self.profile_name.setText(self.working_name)
        self.working_parameters = DepictionStyleParameters.current()
        if self.working_panel.adapter is not None:
            self.working_panel.adapter.apply_parameters(self.working_parameters)
        self._refresh_parameter_controls()
        self._update_change_summary()

    def _refresh_parameter_controls(self) -> None:
        values = self.working_parameters.to_mapping()
        for key, control in self.parameter_controls.items():
            control.set_value(values[key], parameters=self.working_parameters)
        for key, control in self.choice_controls.items():
            control.set_value(values[key])

    def _update_change_summary(self) -> None:
        changes = compare_profiles(self.current_parameters, self.working_parameters)
        if changes:
            names = ", ".join(changes)
            self.status.setText(f"{len(changes)} parámetro(s) distinto(s) de Current: {names}")
            self.status.setToolTip("\n".join(f"{key}: {value['current']} → {value['working']}" for key, value in changes.items()))
        else:
            self.status.setText("Working coincide con Current; no hay cambios experimentales.")
            self.status.setToolTip("")

    def _profile_name_changed(self, value: str) -> None:
        self.working_name = value.strip() or "Working"

    def _save_profile(self) -> None:
        try:
            path = profile_directory() / safe_profile_filename(self.working_name)
            save_profile(
                path,
                self.working_parameters,
                name=self.working_name,
                based_on=self.working_based_on,
            )
        except (OSError, ProfileError) as exc:
            self.status.setText(f"No se pudo guardar el perfil: {exc}")
            return
        self.status.setText(f"Perfil guardado en {path}")

    def _open_profile(self) -> None:
        path, _selected_filter = QFileDialog.getOpenFileName(
            self,
            "Abrir perfil del laboratorio",
            str(profile_directory()),
            "Perfiles de estilo (*.json)",
        )
        if not path:
            return
        try:
            name, parameters = load_profile(Path(path))
        except (OSError, ProfileError) as exc:
            # This is an explicit file operation, outside the no-dialog calibration loop.
            QMessageBox.warning(self, "Perfil inválido", str(exc))
            return
        self.working_name = name
        self.working_based_on = "imported lab profile"
        self.profile_name.setText(name)
        self.working_parameters = parameters
        if self.working_panel.adapter is not None:
            self.working_panel.adapter.apply_parameters(parameters)
        self._refresh_parameter_controls()
        self._update_change_summary()

    def _compare_profiles(self) -> None:
        changes = compare_profiles(self.current_parameters, self.working_parameters)
        if not changes:
            self.status.setText("Current y Working son idénticos.")
            return
        details = "; ".join(
            f"{key}: {value['current']} → {value['working']}"
            for key, value in changes.items()
        )
        self.status.setText(f"Comparación: {details}")
        self.status.setToolTip(details)

    def _export_case(self) -> None:
        directory = QFileDialog.getExistingDirectory(
            self, "Exportar evidencia de depiction", str(Path.home())
        )
        if not directory:
            return
        try:
            result = export_case_bundle(
                Path(directory),
                self.case_id,
                self.current_parameters,
                self.working_parameters,
                theme=self.theme,
                profile_name=self.working_name,
                profile_based_on=self.working_based_on,
            )
        except (OSError, RuntimeError, ValueError) as exc:
            QMessageBox.warning(self, "Error de exportación", str(exc))
            return
        self.status.setText(f"Evidencia exportada: {result['report']}")

    def _export_gallery(self) -> None:
        directory = QFileDialog.getExistingDirectory(
            self, "Exportar galería completa", str(Path.home())
        )
        if not directory:
            return
        try:
            result = export_gallery(
                Path(directory),
                self.current_parameters,
                self.working_parameters,
                theme=self.theme,
                profile_name=self.working_name,
                profile_based_on=self.working_based_on,
            )
        except (OSError, RuntimeError, ValueError) as exc:
            QMessageBox.warning(self, "Error de exportación", str(exc))
            return
        self.status.setText(f"Galería de {result['case_count']} casos exportada: {result['report']}")

    def closeEvent(self, event) -> None:
        for panel in (self.current_panel, self.working_panel):
            if panel.adapter is not None:
                if panel.diagnostics is not None:
                    panel.diagnostics.clear()
                panel.adapter.close()
                panel.adapter = None
        super().closeEvent(event)
