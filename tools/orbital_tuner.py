"""Herramienta manual para ajustar presets orbitales y ver el resultado en vivo."""

from __future__ import annotations

import argparse
import copy
import os
from pathlib import Path
import sys
from typing import Any

os.environ.setdefault("QT_QPA_PLATFORM", os.environ.get("QT_QPA_PLATFORM", ""))

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from PyQt6.QtCore import QFileSystemWatcher, Qt, QTimer
from PyQt6.QtGui import QAction, QKeySequence, QPixmap, QShortcut
from PyQt6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from chemuson.gui.orbitals import (
    ORBITAL_FAMILY_ORDER,
    ORBITAL_PALETTE_MODEL,
    ORBITAL_PRESET_CONFIG_PATH,
    build_orbital_renderer,
    default_orbital_presets_payload,
    load_orbital_presets_payload,
    reload_orbital_presets,
    save_orbital_presets_payload,
)


_CHOICE_FIELDS: dict[str, tuple[str, ...]] = {
    "orientation": ("up", "down", "left", "right"),
    "phase": ("positive", "negative", "neutral"),
    "gradient_mode": ("linear", "radial", "elliptical", "negative"),
    "torus_gradient_mode": ("linear", "radial", "elliptical", "negative"),
    "light_dir": (
        "north",
        "south",
        "east",
        "west",
        "northeast",
        "northwest",
        "southeast",
        "southwest",
        "center",
    ),
}

_FAMILY_TITLES = {
    "s": "Orbital s",
    "p": "Orbital p",
    "d": "Orbital d (clover)",
    "dz2": "Orbital dz2",
    "f": "Orbital f (experimental)",
    "fz3": "Orbital fz3 (experimental)",
    "sp3": "Orbital sp3",
    "sp_lobe": "Lobulo sp",
    "sigma_bonding": "Sigma bonding",
    "pi_bonding": "Pi bonding",
    "torus": "Torus",
}


def _copy_payload(payload: dict[str, Any]) -> dict[str, Any]:
    return copy.deepcopy(payload)


def _clear_layout(layout: QVBoxLayout | QFormLayout) -> None:
    while layout.count():
        item = layout.takeAt(0)
        widget = item.widget()
        child_layout = item.layout()
        if widget is not None:
            widget.deleteLater()
        elif child_layout is not None:
            _clear_layout(child_layout)  # type: ignore[arg-type]


def _friendly_name(key: str) -> str:
    return key.replace("_", " ")


class OrbitalTunerWindow(QMainWindow):
    def __init__(self, config_path: Path) -> None:
        super().__init__()
        self._config_path = config_path
        self._defaults = default_orbital_presets_payload(include_docs=True)
        self._payload = load_orbital_presets_payload(config_path, include_docs=True)
        self._building_form = False
        self._pending_external_reload = False

        self.setWindowTitle("ChemUSON Orbital Tuner")
        self.resize(1440, 920)

        self._watcher = QFileSystemWatcher(self)
        self._watcher.fileChanged.connect(self._on_file_changed)

        self._path_label = QLabel()
        self._path_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        self._status_label = QLabel()
        self._status_label.setWordWrap(True)

        self._family_list = QListWidget()
        self._family_list.currentItemChanged.connect(self._on_family_changed)

        self._triptych_label = QLabel()
        self._triptych_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._triptych_label.setMinimumHeight(300)

        self._sheet_label = QLabel()
        self._sheet_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._sheet_label.setMinimumHeight(220)

        self._editor_host = QWidget()
        self._editor_layout = QVBoxLayout(self._editor_host)
        self._editor_layout.setContentsMargins(8, 8, 8, 8)
        self._editor_layout.setSpacing(10)
        self._editor_layout.addStretch(1)

        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._scroll.setWidget(self._editor_host)

        self._build_ui()
        self._populate_family_list()
        self._ensure_watch_path()

        shortcut = QShortcut(QKeySequence("Ctrl+R"), self)
        shortcut.activated.connect(self.reload_from_disk)

        save_action = QAction(self)
        save_action.setShortcut(QKeySequence.StandardKey.Save)
        save_action.triggered.connect(self.save_to_disk)
        self.addAction(save_action)

        if self._family_list.count():
            self._family_list.setCurrentRow(0)
        self._refresh_previews()

    def _build_ui(self) -> None:
        central = QWidget()
        root = QVBoxLayout(central)
        root.setContentsMargins(10, 10, 10, 10)
        root.setSpacing(8)

        root.addWidget(self._path_label)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        root.addWidget(splitter, 1)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(8)
        left_layout.addWidget(QLabel("Familias orbitales"))
        left_layout.addWidget(self._family_list, 1)

        button_row = QHBoxLayout()
        save_button = QPushButton("Guardar JSON")
        save_button.clicked.connect(self.save_to_disk)
        reload_button = QPushButton("Load last saved presets")
        reload_button.clicked.connect(self.reload_from_disk)
        restore_button = QPushButton("Reset family to canonical defaults")
        restore_button.clicked.connect(self.restore_selected_family)
        export_button = QPushButton("Export current family preview")
        export_button.clicked.connect(self.export_triptych)
        button_row.addWidget(save_button)
        button_row.addWidget(reload_button)
        button_row.addWidget(restore_button)
        button_row.addWidget(export_button)
        left_layout.addLayout(button_row)
        left_layout.addWidget(self._status_label)

        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(10)
        right_layout.addWidget(QLabel("Preview outline / shaded / solid"))
        right_layout.addWidget(self._triptych_label)
        right_layout.addWidget(QLabel("Hoja completa actual"))
        right_layout.addWidget(self._sheet_label)
        right_layout.addWidget(QLabel("Parametros editables"))
        right_layout.addWidget(self._scroll, 1)

        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([280, 1100])

        self.setCentralWidget(central)

    def _populate_family_list(self) -> None:
        self._family_list.clear()
        for family in ORBITAL_FAMILY_ORDER:
            if family not in self._payload:
                continue
            label = _FAMILY_TITLES.get(family, family)
            item = QListWidgetItem(label)
            item.setData(Qt.ItemDataRole.UserRole, family)
            self._family_list.addItem(item)
        self._update_path_label()

    def _update_path_label(self) -> None:
        state = "existente" if self._config_path.exists() else "en memoria (aun no guardado)"
        self._path_label.setText(f"Preset file: {self._config_path} [{state}]")

    def _selected_family(self) -> str | None:
        item = self._family_list.currentItem()
        return None if item is None else str(item.data(Qt.ItemDataRole.UserRole))

    def _ensure_watch_path(self) -> None:
        tracked = set(self._watcher.files())
        if self._config_path.exists() and str(self._config_path) not in tracked:
            self._watcher.addPath(str(self._config_path))
        self._update_path_label()

    def _on_file_changed(self, _path: str) -> None:
        if self._pending_external_reload:
            return
        self._pending_external_reload = True
        QTimer.singleShot(150, self._reload_after_external_change)

    def _reload_after_external_change(self) -> None:
        self._pending_external_reload = False
        if not self._config_path.exists():
            return
        self.reload_from_disk()

    def _set_status(self, text: str) -> None:
        self._status_label.setText(text)

    def _on_family_changed(self, _current: QListWidgetItem | None, _previous: QListWidgetItem | None) -> None:
        self._rebuild_form()
        self._refresh_previews()

    def _nested_value(self, path: tuple[Any, ...]) -> Any:
        current: Any = self._payload
        for key in path:
            current = current[key]
        return current

    def _set_nested_value(self, path: tuple[Any, ...], value: Any) -> None:
        current: Any = self._payload
        for key in path[:-1]:
            current = current[key]
        current[path[-1]] = value

    def _rebuild_form(self) -> None:
        family = self._selected_family()
        self._building_form = True
        try:
            _clear_layout(self._editor_layout)
            if family is None:
                self._editor_layout.addStretch(1)
                return

            family_payload = self._payload[family]
            docs = family_payload.get("_docs", {})

            title = QLabel(f"{_FAMILY_TITLES.get(family, family)}  |  builder: {family_payload.get('builder', '')}")
            title.setStyleSheet("font-weight: 600;")
            self._editor_layout.addWidget(title)

            if family in {"f", "fz3"}:
                hint = QLabel("Para anadir o quitar lobulos, edita el JSON y usa Recargar archivo.")
                hint.setWordWrap(True)
                self._editor_layout.addWidget(hint)

            params_group = QGroupBox("params")
            params_layout = QVBoxLayout(params_group)
            params_layout.setContentsMargins(10, 10, 10, 10)
            self._build_value_widgets(
                parent_layout=params_layout,
                key="params",
                value=family_payload["params"],
                docs=docs,
                path=(family, "params"),
                group_title="params",
            )
            self._editor_layout.addWidget(params_group)
            self._editor_layout.addStretch(1)
        finally:
            self._building_form = False

    def _build_value_widgets(
        self,
        *,
        parent_layout: QVBoxLayout,
        key: str,
        value: Any,
        docs: Any,
        path: tuple[Any, ...],
        group_title: str,
    ) -> None:
        if isinstance(value, dict):
            for child_key, child_value in value.items():
                if str(child_key).startswith("_"):
                    continue
                child_docs = docs.get(child_key, "") if isinstance(docs, dict) else ""
                child_path = path + (child_key,)
                self._append_editor(parent_layout, child_key, child_value, child_docs, child_path)
            return
        self._append_editor(parent_layout, key, value, docs, path)

    def _append_editor(
        self,
        parent_layout: QVBoxLayout,
        key: str,
        value: Any,
        docs: Any,
        path: tuple[Any, ...],
    ) -> None:
        if isinstance(value, dict):
            box = QGroupBox(_friendly_name(key))
            box_layout = QVBoxLayout(box)
            box_layout.setContentsMargins(10, 10, 10, 10)
            if isinstance(docs, str) and docs:
                box.setToolTip(docs)
            for child_key, child_value in value.items():
                if str(child_key).startswith("_"):
                    continue
                child_docs = docs.get(child_key, "") if isinstance(docs, dict) else ""
                self._append_editor(box_layout, child_key, child_value, child_docs, path + (child_key,))
            parent_layout.addWidget(box)
            return

        if isinstance(value, list):
            box = QGroupBox(_friendly_name(key))
            box_layout = QVBoxLayout(box)
            box_layout.setContentsMargins(10, 10, 10, 10)
            if isinstance(docs, str) and docs:
                box.setToolTip(docs)
            for index, item in enumerate(value):
                child_docs = docs if isinstance(docs, dict) else ""
                item_box = QGroupBox(f"{key}[{index}]")
                item_layout = QVBoxLayout(item_box)
                item_layout.setContentsMargins(8, 8, 8, 8)
                if isinstance(item, dict):
                    for child_key, child_value in item.items():
                        if str(child_key).startswith("_"):
                            continue
                        nested_docs = child_docs.get(child_key, "") if isinstance(child_docs, dict) else ""
                        self._append_editor(item_layout, child_key, child_value, nested_docs, path + (index, child_key))
                else:
                    self._append_editor(item_layout, f"{key}[{index}]", item, child_docs, path + (index,))
                box_layout.addWidget(item_box)
            parent_layout.addWidget(box)
            return

        row = QWidget()
        row_layout = QFormLayout(row)
        row_layout.setContentsMargins(0, 0, 0, 0)
        row_layout.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        label = QLabel(_friendly_name(key))
        if isinstance(docs, str) and docs:
            label.setToolTip(docs)
        editor = self._editor_for_scalar(key, value, docs, path)
        row_layout.addRow(label, editor)
        parent_layout.addWidget(row)

    def _editor_for_scalar(self, key: str, value: Any, docs: Any, path: tuple[Any, ...]) -> QWidget:
        tooltip = docs if isinstance(docs, str) else ""

        if isinstance(value, bool):
            editor = QCheckBox()
            editor.setChecked(value)
            editor.setToolTip(tooltip)
            editor.toggled.connect(lambda checked, p=path: self._scalar_changed(p, bool(checked)))
            return editor

        if isinstance(value, (int, float)) and not isinstance(value, bool):
            editor = QDoubleSpinBox()
            editor.setRange(-9999.0, 9999.0)
            editor.setDecimals(4)
            editor.setSingleStep(0.1)
            editor.setValue(float(value))
            editor.setToolTip(tooltip)
            editor.valueChanged.connect(lambda number, p=path: self._scalar_changed(p, float(number)))
            return editor

        if key in _CHOICE_FIELDS:
            editor = QComboBox()
            editor.setEditable(True)
            editor.addItems(_CHOICE_FIELDS[key])
            editor.setCurrentText(str(value))
            editor.setToolTip(tooltip)
            editor.currentTextChanged.connect(lambda text, p=path: self._scalar_changed(p, str(text)))
            return editor

        editor = QLineEdit(str(value))
        editor.setToolTip(tooltip)
        editor.editingFinished.connect(lambda p=path, w=editor: self._scalar_changed(p, w.text()))
        return editor

    def _scalar_changed(self, path: tuple[Any, ...], value: Any) -> None:
        if self._building_form:
            return
        self._set_nested_value(path, value)
        self._refresh_previews()

    def _current_renderer(self):
        return build_orbital_renderer(self._payload)

    def _refresh_previews(self) -> None:
        family = self._selected_family()
        if family is None:
            return
        try:
            renderer = self._current_renderer()
            triptych = renderer.render_style_triptych(family, panel_size=220, gap=18)
            sheet = renderer.render_palette_image(ORBITAL_PALETTE_MODEL)
            self._triptych_label.setPixmap(QPixmap.fromImage(triptych))
            self._sheet_label.setPixmap(QPixmap.fromImage(sheet))
            self._set_status("Preview actualizada.")
        except Exception as exc:  # pragma: no cover - feedback de UI
            self._set_status(f"No se pudo regenerar el preview: {exc}")

    def save_to_disk(self) -> None:
        try:
            save_orbital_presets_payload(self._payload, self._config_path)
            reload_orbital_presets(self._config_path)
            self._ensure_watch_path()
            self._set_status(f"Presets guardados en {self._config_path}")
        except Exception as exc:  # pragma: no cover - feedback de UI
            QMessageBox.critical(self, "Error al guardar", str(exc))

    def reload_from_disk(self) -> None:
        family = self._selected_family()
        self._payload = load_orbital_presets_payload(self._config_path, include_docs=True)
        self._ensure_watch_path()
        self._rebuild_form()
        self._refresh_previews()
        self._set_status(f"Recargado desde {self._config_path}")
        if family is not None:
            for row in range(self._family_list.count()):
                item = self._family_list.item(row)
                if str(item.data(Qt.ItemDataRole.UserRole)) == family:
                    self._family_list.setCurrentRow(row)
                    break

    def restore_selected_family(self) -> None:
        family = self._selected_family()
        if family is None:
            return
        self._payload[family] = _copy_payload(self._defaults[family])
        self._rebuild_form()
        self._refresh_previews()
        self._set_status(f"Se restauraron los defaults de {family}")

    def export_triptych(self) -> None:
        family = self._selected_family()
        if family is None:
            return
        renderer = self._current_renderer()
        image = renderer.render_style_triptych(family, panel_size=220, gap=18)
        default_path = self._config_path.parent / f"{family}_preview.png"
        target, _ = QFileDialog.getSaveFileName(
            self,
            "Exportar preview orbital",
            str(default_path),
            "PNG (*.png)",
        )
        if not target:
            return
        if not image.save(target):
            QMessageBox.critical(self, "Error al exportar", f"No se pudo guardar {target}")
            return
        self._set_status(f"Preview exportada a {target}")


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=ORBITAL_PRESET_CONFIG_PATH,
        help=f"Ruta del archivo editable de presets (default: {ORBITAL_PRESET_CONFIG_PATH})",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    app = QApplication.instance() or QApplication(sys.argv)
    window = OrbitalTunerWindow(args.config)
    window.show()
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
