"""Interactive elemental analysis dialog."""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QApplication,
    QAbstractItemView,
    QDialog,
    QDoubleSpinBox,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from chemuson.core.elemental_analysis import (
    DEFAULT_TOLERANCE,
    SOLVENT_LIBRARY,
    FormulaError,
    SolvateCandidate,
    SolventEntry,
    compare_found_calculated,
    elemental_percentages,
    find_solvate_candidates,
    format_formula,
    format_publication_line,
    molecular_weight,
    parse_formula,
)


COMMON_COMPARISON_ELEMENTS = ("C", "H", "N", "O", "S")


class ElementalAnalysisDialog(QDialog):
    """Formula, found/calculated comparison, solvate search, and export tool."""

    def __init__(self, initial_formula: str = "", parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Elemental Analysis")
        self.resize(980, 700)
        self.setMinimumSize(760, 540)

        self._composition: dict[str, float] = {}
        self._percentages: dict[str, float] = {}
        self._found_values: dict[str, float | None] = {}
        self._selected_elements: set[str] = set()
        self._comparison_user_touched = False
        self._solvent_entries: list[SolventEntry] = list(SOLVENT_LIBRARY.values())
        self._custom_solvents: list[SolventEntry] = []
        self._solvate_candidates: list[SolvateCandidate] = []
        self._updating = False

        self.tabs = QTabWidget(self)
        self.formula_tab = self._build_formula_tab()
        self.comparison_tab = self._build_comparison_tab()
        self.solvate_tab = self._build_solvate_tab()
        self.export_tab = self._build_export_tab()
        self.tabs.addTab(self.formula_tab, "Formula")
        self.tabs.addTab(self.comparison_tab, "Experimental comparison")
        self.tabs.addTab(self.solvate_tab, "Hydrate/Solvate finder")
        self.tabs.addTab(self.export_tab, "Export")

        layout = QVBoxLayout(self)
        layout.addWidget(self.tabs)

        self._connect_signals()
        self._rebuild_solvent_table()
        if initial_formula:
            self.formula_edit.setText(initial_formula)
        else:
            self._recalculate_formula()

    def _build_formula_tab(self) -> QWidget:
        tab = QWidget(self)
        layout = QVBoxLayout(tab)

        formula_row = QHBoxLayout()
        self.formula_edit = QLineEdit()
        self.formula_edit.setPlaceholderText("C20H18N4O2\u00b70.5H2O")
        self.calculate_button = QPushButton("Calculate")
        formula_row.addWidget(self.formula_edit, 1)
        formula_row.addWidget(self.calculate_button)

        self.percent_decimals_spin = QSpinBox()
        self.percent_decimals_spin.setRange(0, 6)
        self.percent_decimals_spin.setValue(2)
        self.mw_decimals_spin = QSpinBox()
        self.mw_decimals_spin.setRange(0, 6)
        self.mw_decimals_spin.setValue(2)

        form = QFormLayout()
        form.addRow("Formula", formula_row)
        rounding_row = QHBoxLayout()
        rounding_row.addWidget(QLabel("Percent decimals"))
        rounding_row.addWidget(self.percent_decimals_spin)
        rounding_row.addSpacing(16)
        rounding_row.addWidget(QLabel("MW decimals"))
        rounding_row.addWidget(self.mw_decimals_spin)
        rounding_row.addStretch(1)
        form.addRow("Rounding", rounding_row)
        layout.addLayout(form)

        self.validation_label = QLabel()
        self.validation_label.setWordWrap(True)
        layout.addWidget(self.validation_label)

        summary = QGridLayout()
        self.normalized_formula_label = QLabel("N/D")
        self.molecular_weight_label = QLabel("N/D")
        summary.addWidget(QLabel("Normalized formula"), 0, 0)
        summary.addWidget(self.normalized_formula_label, 0, 1)
        summary.addWidget(QLabel("Molecular weight"), 1, 0)
        summary.addWidget(self.molecular_weight_label, 1, 1)
        summary.setColumnStretch(1, 1)
        layout.addLayout(summary)

        self.percent_table = QTableWidget(0, 2, self)
        self.percent_table.setHorizontalHeaderLabels(["Element", "Calculated %"])
        self.percent_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.percent_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.percent_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.percent_table, 1)
        return tab

    def _build_comparison_tab(self) -> QWidget:
        tab = QWidget(self)
        layout = QVBoxLayout(tab)

        controls = QHBoxLayout()
        controls.addWidget(QLabel("Tolerance (+/- percentage points)"))
        self.tolerance_spin = QDoubleSpinBox()
        self.tolerance_spin.setRange(0.01, 20.0)
        self.tolerance_spin.setDecimals(2)
        self.tolerance_spin.setSingleStep(0.05)
        self.tolerance_spin.setValue(DEFAULT_TOLERANCE)
        controls.addWidget(self.tolerance_spin)
        controls.addStretch(1)
        layout.addLayout(controls)

        self.comparison_table = QTableWidget(0, 7, self)
        self.comparison_table.setHorizontalHeaderLabels(
            ["Use", "Element", "Calculated %", "Found %", "Tolerance", "Delta", "Pass"]
        )
        self.comparison_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.comparison_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.comparison_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        layout.addWidget(self.comparison_table, 1)
        return tab

    def _build_solvate_tab(self) -> QWidget:
        tab = QWidget(self)
        layout = QVBoxLayout(tab)

        solvent_group = QGroupBox("Candidate solvents/adducts")
        solvent_layout = QVBoxLayout(solvent_group)
        self.solvent_table = QTableWidget(0, 3, self)
        self.solvent_table.setHorizontalHeaderLabels(["Use", "Name", "Formula"])
        self.solvent_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.solvent_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.solvent_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self.solvent_table.setMaximumHeight(210)
        solvent_layout.addWidget(self.solvent_table)

        custom_row = QHBoxLayout()
        self.custom_solvent_name_edit = QLineEdit()
        self.custom_solvent_name_edit.setPlaceholderText("Name")
        self.custom_solvent_formula_edit = QLineEdit()
        self.custom_solvent_formula_edit.setPlaceholderText("Custom formula, e.g. C5H12")
        self.add_custom_solvent_button = QPushButton("Add")
        custom_row.addWidget(self.custom_solvent_name_edit, 1)
        custom_row.addWidget(self.custom_solvent_formula_edit, 2)
        custom_row.addWidget(self.add_custom_solvent_button)
        solvent_layout.addLayout(custom_row)
        layout.addWidget(solvent_group)

        controls = QHBoxLayout()
        self.eq_min_spin = QDoubleSpinBox()
        self.eq_min_spin.setRange(0.0, 20.0)
        self.eq_min_spin.setDecimals(2)
        self.eq_min_spin.setSingleStep(0.25)
        self.eq_min_spin.setValue(0.25)
        self.eq_max_spin = QDoubleSpinBox()
        self.eq_max_spin.setRange(0.0, 20.0)
        self.eq_max_spin.setDecimals(2)
        self.eq_max_spin.setSingleStep(0.25)
        self.eq_max_spin.setValue(2.0)
        self.eq_step_spin = QDoubleSpinBox()
        self.eq_step_spin.setRange(0.01, 5.0)
        self.eq_step_spin.setDecimals(2)
        self.eq_step_spin.setSingleStep(0.25)
        self.eq_step_spin.setValue(0.25)
        self.max_components_spin = QSpinBox()
        self.max_components_spin.setRange(1, 4)
        self.max_components_spin.setValue(1)
        self.find_solvates_button = QPushButton("Find candidates")

        for label, widget in (
            ("Min eq", self.eq_min_spin),
            ("Max eq", self.eq_max_spin),
            ("Step", self.eq_step_spin),
            ("Max components", self.max_components_spin),
        ):
            controls.addWidget(QLabel(label))
            controls.addWidget(widget)
        controls.addStretch(1)
        controls.addWidget(self.find_solvates_button)
        layout.addLayout(controls)

        self.solvate_warning_label = QLabel("Candidate formulas only; verify chemically before assignment.")
        self.solvate_warning_label.setWordWrap(True)
        layout.addWidget(self.solvate_warning_label)

        self.solvate_results_table = QTableWidget(0, 7, self)
        self.solvate_results_table.setHorizontalHeaderLabels(
            ["Candidate formula", "MW", "Score", "Max |Delta|", "Pass", "Deltas", "Warning"]
        )
        self.solvate_results_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.solvate_results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.solvate_results_table, 1)
        return tab

    def _build_export_tab(self) -> QWidget:
        tab = QWidget(self)
        layout = QVBoxLayout(tab)

        buttons = QHBoxLayout()
        self.copy_calculated_button = QPushButton("Copy calculated table")
        self.copy_comparison_button = QPushButton("Copy comparison table")
        self.copy_publication_button = QPushButton("Copy manuscript line")
        buttons.addWidget(self.copy_calculated_button)
        buttons.addWidget(self.copy_comparison_button)
        buttons.addWidget(self.copy_publication_button)
        buttons.addStretch(1)
        layout.addLayout(buttons)

        self.publication_preview = QPlainTextEdit()
        self.publication_preview.setReadOnly(True)
        layout.addWidget(QLabel("Manuscript-ready line"))
        layout.addWidget(self.publication_preview, 1)
        return tab

    def _connect_signals(self) -> None:
        self.formula_edit.textChanged.connect(self._recalculate_formula)
        self.calculate_button.clicked.connect(self._recalculate_formula)
        self.percent_decimals_spin.valueChanged.connect(self._refresh_all_outputs)
        self.mw_decimals_spin.valueChanged.connect(self._refresh_all_outputs)
        self.tolerance_spin.valueChanged.connect(self._refresh_comparison_outputs)
        self.comparison_table.itemChanged.connect(self._on_comparison_item_changed)
        self.add_custom_solvent_button.clicked.connect(self._add_custom_solvent)
        self.find_solvates_button.clicked.connect(self._find_solvate_candidates)
        self.solvate_results_table.cellClicked.connect(self._load_solvate_candidate)
        self.copy_calculated_button.clicked.connect(self._copy_calculated_table)
        self.copy_comparison_button.clicked.connect(self._copy_comparison_table)
        self.copy_publication_button.clicked.connect(self._copy_publication_line)

    def _recalculate_formula(self, *_args) -> None:
        self._read_comparison_state()
        formula = self.formula_edit.text().strip()
        try:
            composition = parse_formula(formula)
            percentages = elemental_percentages(composition)
            mw = molecular_weight(composition)
        except FormulaError as exc:
            self._composition = {}
            self._percentages = {}
            self.validation_label.setText(f"Invalid formula: {exc}")
            self.validation_label.setStyleSheet("color: #b3261e;")
            self.normalized_formula_label.setText("N/D")
            self.molecular_weight_label.setText("N/D")
            self._rebuild_percent_table()
            self._rebuild_comparison_table()
            self._refresh_export_preview()
            return

        self._composition = composition
        self._percentages = percentages
        self.validation_label.setText("Formula parsed successfully.")
        self.validation_label.setStyleSheet("color: #146c2e;")
        self.normalized_formula_label.setText(format_formula(composition))
        self.molecular_weight_label.setText(f"{mw:.{self.mw_decimals_spin.value()}f}")
        self._rebuild_percent_table()
        self._rebuild_comparison_table()
        self._refresh_export_preview()

    def _refresh_all_outputs(self, *_args) -> None:
        if self._composition:
            self.normalized_formula_label.setText(format_formula(self._composition))
            self.molecular_weight_label.setText(
                f"{molecular_weight(self._composition):.{self.mw_decimals_spin.value()}f}"
            )
        self._rebuild_percent_table()
        self._rebuild_comparison_table()
        self._refresh_export_preview()
        self._rebuild_solvate_results_table()

    def _refresh_comparison_outputs(self, *_args) -> None:
        self._read_comparison_state()
        self._rebuild_comparison_table()
        self._refresh_export_preview()

    def _on_comparison_item_changed(self, _item: QTableWidgetItem) -> None:
        if self._updating:
            return
        self._comparison_user_touched = True
        self._read_comparison_state()
        self._rebuild_comparison_table()
        self._refresh_export_preview()

    def _read_comparison_state(self) -> None:
        if not hasattr(self, "comparison_table") or self._updating:
            return
        selected: set[str] = set()
        found: dict[str, float | None] = {}
        for row in range(self.comparison_table.rowCount()):
            element_item = self.comparison_table.item(row, 1)
            if element_item is None:
                continue
            element = element_item.text().strip()
            use_item = self.comparison_table.item(row, 0)
            if use_item is not None and use_item.checkState() == Qt.CheckState.Checked:
                selected.add(element)
            found_item = self.comparison_table.item(row, 3)
            found[element] = _item_float(found_item)
        if self._comparison_user_touched:
            self._selected_elements = selected
        self._found_values.update(found)

    def _rebuild_percent_table(self) -> None:
        self._updating = True
        try:
            rows = list(_element_order(self._percentages.keys()))
            self.percent_table.setRowCount(len(rows))
            decimals = self.percent_decimals_spin.value()
            for row, element in enumerate(rows):
                self.percent_table.setItem(row, 0, _readonly_item(element))
                self.percent_table.setItem(row, 1, _readonly_item(f"{self._percentages[element]:.{decimals}f}"))
        finally:
            self._updating = False

    def _rebuild_comparison_table(self) -> None:
        self._updating = True
        try:
            elements = _comparison_elements(self._percentages.keys(), self._found_values.keys())
            if (
                not self._comparison_user_touched
                and not self._selected_elements
                and self._percentages
            ):
                self._selected_elements = {
                    element for element in ("C", "H", "N") if element in self._percentages
                } or set(elements)
            rows_by_element = {
                row.element: row
                for row in compare_found_calculated(
                    self._percentages,
                    {element: self._found_values.get(element) for element in elements},
                    self.tolerance_spin.value(),
                )
            }
            self.comparison_table.setRowCount(len(elements))
            decimals = self.percent_decimals_spin.value()
            for row, element in enumerate(elements):
                comparison = rows_by_element[element]
                use_item = QTableWidgetItem("")
                use_item.setFlags(
                    Qt.ItemFlag.ItemIsUserCheckable
                    | Qt.ItemFlag.ItemIsEnabled
                    | Qt.ItemFlag.ItemIsSelectable
                )
                use_item.setCheckState(
                    Qt.CheckState.Checked
                    if element in self._selected_elements
                    else Qt.CheckState.Unchecked
                )
                self.comparison_table.setItem(row, 0, use_item)
                self.comparison_table.setItem(row, 1, _readonly_item(element))
                self.comparison_table.setItem(row, 2, _readonly_item(f"{comparison.calculated:.{decimals}f}"))

                found_item = QTableWidgetItem(
                    "" if comparison.found is None else f"{comparison.found:.{decimals}f}"
                )
                self.comparison_table.setItem(row, 3, found_item)
                self.comparison_table.setItem(row, 4, _readonly_item(f"+/- {comparison.tolerance:.2f}"))
                delta_text = "" if comparison.delta is None else f"{comparison.delta:+.{decimals}f}"
                self.comparison_table.setItem(row, 5, _readonly_item(delta_text))
                if comparison.passed is None:
                    passed_text = ""
                else:
                    passed_text = "Pass" if comparison.passed else "Fail"
                self.comparison_table.setItem(row, 6, _readonly_item(passed_text))
        finally:
            self._updating = False

    def _selected_comparison_elements_with_found(self) -> tuple[list[str], dict[str, float]]:
        self._read_comparison_state()
        elements = [element for element in _element_order(self._selected_elements) if element in self._percentages]
        found = {
            element: value
            for element in elements
            if (value := self._found_values.get(element)) is not None
        }
        return elements, found

    def _rebuild_solvent_table(self, extra_selected: set[tuple[str, str]] | None = None) -> None:
        entries = self._solvent_entries + self._custom_solvents
        selected_keys = self._selected_solvent_keys()
        if not selected_keys:
            selected_keys = {("H2O", "H2O")}
        if extra_selected:
            selected_keys.update(extra_selected)
        self.solvent_table.setRowCount(len(entries))
        for row, entry in enumerate(entries):
            use_item = QTableWidgetItem("")
            use_item.setFlags(
                Qt.ItemFlag.ItemIsUserCheckable
                | Qt.ItemFlag.ItemIsEnabled
                | Qt.ItemFlag.ItemIsSelectable
            )
            use_item.setCheckState(
                Qt.CheckState.Checked
                if (entry.name, entry.formula) in selected_keys
                else Qt.CheckState.Unchecked
            )
            self.solvent_table.setItem(row, 0, use_item)
            self.solvent_table.setItem(row, 1, _readonly_item(entry.name))
            self.solvent_table.setItem(row, 2, _readonly_item(entry.formula))

    def _selected_solvents(self) -> list[SolventEntry]:
        entries = self._solvent_entries + self._custom_solvents
        selected: list[SolventEntry] = []
        for row, entry in enumerate(entries):
            item = self.solvent_table.item(row, 0)
            if item is not None and item.checkState() == Qt.CheckState.Checked:
                selected.append(entry)
        return selected

    def _selected_solvent_keys(self) -> set[tuple[str, str]]:
        if not hasattr(self, "solvent_table"):
            return set()
        entries = self._solvent_entries + self._custom_solvents
        selected: set[tuple[str, str]] = set()
        for row, entry in enumerate(entries):
            item = self.solvent_table.item(row, 0)
            if item is not None and item.checkState() == Qt.CheckState.Checked:
                selected.add((entry.name, entry.formula))
        return selected

    def _add_custom_solvent(self) -> None:
        name = self.custom_solvent_name_edit.text().strip()
        formula = self.custom_solvent_formula_edit.text().strip()
        try:
            parse_formula(formula)
        except FormulaError as exc:
            self.solvate_warning_label.setText(f"Invalid custom solvent formula: {exc}")
            return
        entry = SolventEntry(name or formula, formula)
        self._custom_solvents.append(entry)
        self.custom_solvent_name_edit.clear()
        self.custom_solvent_formula_edit.clear()
        self._rebuild_solvent_table(extra_selected={(entry.name, entry.formula)})
        self.solvate_warning_label.setText("Custom solvent added.")

    def _find_solvate_candidates(self) -> None:
        if not self._composition:
            self.solvate_warning_label.setText("Enter a valid base formula before searching.")
            return
        elements, found = self._selected_comparison_elements_with_found()
        if not elements or not found:
            self.solvate_warning_label.setText(
                "Select comparison elements and enter found percentages before searching."
            )
            return
        if len(found) != len(elements):
            self.solvate_warning_label.setText(
                "Every selected comparison element needs a found percentage for ranking."
            )
            return
        solvents = self._selected_solvents()
        if not solvents:
            self.solvate_warning_label.setText("Select at least one candidate solvent/adduct.")
            return
        try:
            candidates = find_solvate_candidates(
                self.formula_edit.text().strip(),
                found,
                solvents,
                equivalent_range=(self.eq_min_spin.value(), self.eq_max_spin.value()),
                step_size=self.eq_step_spin.value(),
                selected_elements=elements,
                tolerances=self.tolerance_spin.value(),
                max_components=self.max_components_spin.value(),
                max_results=50,
            )
        except FormulaError as exc:
            self.solvate_warning_label.setText(f"Solvate search failed: {exc}")
            return
        self._solvate_candidates = candidates
        warning = "Candidate formulas only; verify chemically before assignment."
        if self.max_components_spin.value() > 1:
            warning += " Multiple solvent components can overfit elemental percentages."
        self.solvate_warning_label.setText(warning)
        self._rebuild_solvate_results_table()

    def _rebuild_solvate_results_table(self) -> None:
        self._updating = True
        try:
            decimals = self.percent_decimals_spin.value()
            mw_decimals = self.mw_decimals_spin.value()
            self.solvate_results_table.setRowCount(len(self._solvate_candidates))
            for row, candidate in enumerate(self._solvate_candidates):
                values = [
                    candidate.formula,
                    f"{candidate.molecular_weight:.{mw_decimals}f}",
                    f"{candidate.score:.3f}",
                    f"{candidate.max_absolute_delta:.{decimals}f}",
                    "Pass" if candidate.passed else "Fail",
                    "; ".join(
                        f"{element} {delta:+.{decimals}f}"
                        for element, delta in _ordered_delta_items(candidate.deltas)
                    ),
                    candidate.warning,
                ]
                for column, value in enumerate(values):
                    item = _readonly_item(value)
                    if column == 0:
                        item.setData(Qt.ItemDataRole.UserRole, row)
                    self.solvate_results_table.setItem(row, column, item)
        finally:
            self._updating = False

    def _load_solvate_candidate(self, row: int, _column: int) -> None:
        if row < 0 or row >= len(self._solvate_candidates):
            return
        candidate = self._solvate_candidates[row]
        self.formula_edit.setText(candidate.formula)
        self.tabs.setCurrentWidget(self.formula_tab)

    def _refresh_export_preview(self) -> None:
        if not hasattr(self, "publication_preview"):
            return
        if not self._percentages:
            self.publication_preview.setPlainText("")
            return
        self._read_comparison_state()
        elements = [element for element in _element_order(self._selected_elements) if element in self._percentages]
        if not elements:
            elements = _element_order(self._percentages.keys())
        formula = self.formula_edit.text().strip() or format_formula(self._composition)
        line = format_publication_line(
            formula,
            self._percentages,
            self._found_values,
            elements,
            decimals=self.percent_decimals_spin.value(),
        )
        self.publication_preview.setPlainText(line)

    def _copy_calculated_table(self) -> None:
        QApplication.clipboard().setText(self._calculated_table_text())

    def _copy_comparison_table(self) -> None:
        QApplication.clipboard().setText(self._comparison_table_text())

    def _copy_publication_line(self) -> None:
        QApplication.clipboard().setText(self.publication_preview.toPlainText())

    def _calculated_table_text(self) -> str:
        decimals = self.percent_decimals_spin.value()
        lines = ["Element\tCalculated %"]
        for element in _element_order(self._percentages.keys()):
            lines.append(f"{element}\t{self._percentages[element]:.{decimals}f}")
        return "\n".join(lines)

    def _comparison_table_text(self) -> str:
        self._read_comparison_state()
        decimals = self.percent_decimals_spin.value()
        selected = [element for element in _element_order(self._selected_elements) if element in self._percentages]
        rows = compare_found_calculated(
            self._percentages,
            {element: self._found_values.get(element) for element in selected},
            self.tolerance_spin.value(),
        )
        lines = ["Element\tCalculated %\tFound %\tDelta\tTolerance\tPass"]
        for row in rows:
            found = "" if row.found is None else f"{row.found:.{decimals}f}"
            delta = "" if row.delta is None else f"{row.delta:+.{decimals}f}"
            passed = "" if row.passed is None else ("Pass" if row.passed else "Fail")
            lines.append(
                f"{row.element}\t{row.calculated:.{decimals}f}\t{found}\t"
                f"{delta}\t{row.tolerance:.2f}\t{passed}"
            )
        return "\n".join(lines)


def _readonly_item(text: str) -> QTableWidgetItem:
    item = QTableWidgetItem(text)
    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
    return item


def _item_float(item: QTableWidgetItem | None) -> float | None:
    if item is None:
        return None
    text = item.text().strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _element_order(elements) -> list[str]:
    unique = dict.fromkeys(elements)
    return sorted(unique, key=lambda element: (0 if element == "C" else 1 if element == "H" else 2, element))


def _comparison_elements(percent_elements, found_elements) -> list[str]:
    combined = list(COMMON_COMPARISON_ELEMENTS)
    combined.extend(percent_elements)
    combined.extend(found_elements)
    return _element_order(combined)


def _ordered_delta_items(deltas: dict[str, float]) -> list[tuple[str, float]]:
    return [(element, deltas[element]) for element in _element_order(deltas.keys())]
