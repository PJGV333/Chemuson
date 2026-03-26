"""
Ventana principal de Chemuson.

Compone menús, barras de herramientas, docks y el lienzo central.
"""
from PyQt6.QtWidgets import (
    QApplication,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QFontDialog,
    QInputDialog,
    QMainWindow,
    QMenu,
    QMenuBar,
    QFileDialog,
    QMessageBox,
    QToolBar,
    QWidget,
    QFormLayout,
    QDoubleSpinBox,
    QComboBox,
    QAbstractItemView,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QPushButton,
    QVBoxLayout,
    QLabel,
    QHeaderView,
    QSizePolicy,
)
from PyQt6.QtCore import Qt, QSize, QSettings, QEvent, QTimer
from PyQt6.QtGui import QAction, QActionGroup, QKeySequence, QIcon, QPainter, QPixmap, QPen, QColor
from PyQt6.QtPrintSupport import QPrinter
from typing import Callable, Optional
from dataclasses import replace
from datetime import datetime
import json
import math
import os
import sys

from chemuson.gui.canvas import (
    BRANCH_ROTATION_STEP_DEG,
    FRAGMENT_ROTATION_STEP_DEG,
    ChemusonCanvas,
)
from chemuson.gui.periodic_table import PeriodicTableDialog
from chemuson.gui.orbitals import ORBITAL_MENU_ORDER, orbital_display_name, orbital_tool_id
from chemuson.gui.toolbar import ChemusonToolbar, SymbolPaletteToolbar
from chemuson.gui.styles import MAIN_STYLESHEET, TOOL_PALETTE_STYLESHEET
from chemuson.gui.icons import draw_generic_icon
from chemuson.gui.docks import PlantillasDock, InspectorDock, AppearanceDock
from chemuson.gui.dialogs import PreferencesDialog, QuickStartDialog, StyleDialog
from chemuson.gui.text_toolbar import TextFormatToolbar
from chemuson.gui.commands import ChangeAtomCommand
from chemuson.gui.template_library import TemplateLibrary, DEFAULT_CATEGORY_USER
from chemuson.gui.templates import (
    build_linear_chain_template,
)
from chemuson.chemio.rdkit_io import molfile_to_molgraph, molgraph_to_molfile
from chemuson.chemio.persistence import PersistenceManager
from chemuson.core.model import MolGraph
from chemuson.utils.autosave import AutosaveManager
from chemuson.utils import crash_reporter
from chemuson.version import get_app_version
from chemuson.update import (
    AutoUpdateCore,
    GitHubReleasesProvider,
    PortableUpdateContext,
    UpdateChannel,
    UpdateMode,
    UpdateSettings,
    UpdateTelemetryLogger,
    choose_windows_asset_flavor,
    detect_portable_update_context,
    detect_windows_install_context,
    is_portable_target_writable,
    is_windows_installer_asset,
    launch_portable_update_script,
    launch_inno_installer,
    mark_checked,
    should_check_now,
    write_portable_update_script,
)


FLATPAK_APP_ID = "io.github.PJGV333.Chemuson"


def _channel_display_name(channel: str) -> str:
    """Devuelve un nombre legible para el canal de updates."""
    value = str(channel or "").strip().lower()
    if value == "stable":
        return "estable"
    if value == "beta":
        return "beta"
    return value or "desconocido"


def _update_source_display_name(source: str) -> str:
    """Devuelve un nombre legible para el origen del feed de updates."""
    value = str(source or "").strip().lower()
    if value == "cache":
        return "caché local"
    if value == "remote":
        return "GitHub"
    if value == "error":
        return "error"
    return ""


def format_no_update_message(
    current_version: str,
    channel: str,
    reason: str = "",
    latest_version: str = "",
    source: str = "",
) -> str:
    """Construye mensaje explicativo cuando no hay update elegible."""
    lines = [
        "No hay una version publicada mas nueva para tu canal actual.",
        "",
        f"Version instalada: {str(current_version or '').strip() or 'desconocida'}",
        f"Canal: {_channel_display_name(channel)}",
    ]
    latest = str(latest_version or "").strip()
    if latest:
        lines.append(f"Ultima version consultada: {latest}")
    source_name = _update_source_display_name(source)
    if source_name:
        lines.append(f"Origen de datos: {source_name}")
    lines.extend(
        [
            "",
            "Este verificador compara releases publicadas; no distribuye commits sueltos.",
        ]
    )
    if str(source or "").strip().lower() == "cache":
        lines.extend(
            [
                "",
                "Aviso: el resultado proviene de la caché local. Si acabas de publicar una release, "
                "vuelve a intentarlo cuando GitHub esté accesible.",
            ]
        )
    detail = str(reason or "").strip()
    if detail:
        lines.extend(["", f"Detalle: {detail}"])
    return "\n".join(lines)


def is_running_in_flatpak() -> bool:
    """Detecta si Chemuson corre dentro de un sandbox Flatpak."""
    flatpak_id = str(os.getenv("FLATPAK_ID", "") or "").strip()
    if flatpak_id:
        return True
    return os.path.exists("/.flatpak-info")


def format_update_disabled_message(flatpak: bool = False, app_id: str = FLATPAK_APP_ID) -> str:
    """Construye mensaje cuando el chequeo interno de updates está deshabilitado."""
    if flatpak:
        return (
            "Esta edicion Flatpak no usa auto-update real dentro de Chemuson.\n"
            "Se actualiza con Flatpak, no desde la propia app.\n\n"
            f"Usa:\nflatpak update {app_id}\n\n"
            "Si instalaste desde un bundle local sin un remote configurado, "
            "instala el bundle mas reciente manualmente."
        )
    return "El chequeo de actualizaciones está deshabilitado por entorno."


class ChemusonWindow(QMainWindow):
    """
    Ventana principal del editor molecular Chemuson.
    Coordina acciones de menú, toolbars y el lienzo de dibujo.
    """
    def __init__(self) -> None:
        """Inicializa la instancia.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        super().__init__()
        self._app_version = get_app_version()
        self.setWindowTitle(f"Chemuson {self._app_version} - Editor Molecular Libre")
        self.resize(1200, 900)
        
        # Apply main stylesheet
        self.setStyleSheet(MAIN_STYLESHEET)

        # === CORE COMPONENTS ===
        self._create_actions()
        self._current_file_path: Optional[str] = None
        self._canvas_file_paths: dict[ChemusonCanvas, Optional[str]] = {}
        self._canvas_tab_titles: dict[ChemusonCanvas, str] = {}
        self._canvas_autosave_managers: dict[ChemusonCanvas, AutosaveManager] = {}
        self._untitled_counter = 1
        self._active_canvas_connected: Optional[ChemusonCanvas] = None
        self._numbering_default_mode = "atoms"
        self._numbering_default_include_export = True
        self._current_tool_id = "tool_select"
        self._settings = QSettings("Chemuson", "Chemuson")
        self._update_settings = self._load_update_preferences()
        self._name_advanced_default, self._name_rdkit_isolated_default = self._load_naming_preferences()
        self._windows_install_context = detect_windows_install_context()
        self._portable_update_context = detect_portable_update_context(
            windows_context=self._windows_install_context
        )
        self._pending_windows_installer_path = ""
        self._pending_windows_installer_version = ""
        self._pending_windows_download = None
        self._pending_portable_target_path = ""
        self._pending_portable_version = ""
        self._pending_portable_download = None
        self._pending_portable_relaunch = False
        self._pending_portable_context = PortableUpdateContext(
            is_portable=False,
            is_windows=False,
            is_appimage=False,
            target_path="",
            executable_path="",
            display_name="",
        )
        self._update_telemetry = UpdateTelemetryLogger()
        self._recent_files = self._load_recent_files()
        self.template_library = TemplateLibrary()
        self._template_icon_cache: dict[str, QIcon] = {}
        
        # === CENTRAL TABS/CANVAS ===
        self.tabs = QTabWidget(self)
        self.tabs.setDocumentMode(True)
        self.tabs.setTabsClosable(True)
        self.tabs.setMovable(True)
        self.setCentralWidget(self.tabs)
        self.canvas = self._create_document_tab(make_current=True)
        self.tabs.currentChanged.connect(self._on_tab_changed)
        self.tabs.tabCloseRequested.connect(self._on_tab_close_requested)
        self.action_aromatic_circles.setChecked(self.canvas.state.use_aromatic_circles)
        self.action_rules.setChecked(self.canvas.show_rulers)
        self.action_crosshair.setChecked(self.canvas.show_grid)
        self._load_numbering_preferences()
        self._sync_numbering_actions()
        
        # === DOCK WIDGETS ===
        self.templates_dock = PlantillasDock(self)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.templates_dock)
        self.templates_dock.hide()
        self.templates_dock.template_selected.connect(self._on_template_selected_from_gallery)
        self.templates_dock.new_category_requested.connect(self._on_new_template_category)
        self.templates_dock.import_requested.connect(self._on_import_template_library)
        self.templates_dock.export_requested.connect(self._on_export_template_library)
        self.templates_dock.rename_category_requested.connect(self._on_rename_template_category)
        self.templates_dock.delete_category_requested.connect(self._on_delete_template_category)
        self.templates_dock.rename_template_requested.connect(self._on_rename_template)
        self.templates_dock.delete_template_requested.connect(self._on_delete_template)
        
        self.inspector_dock = InspectorDock(self)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.inspector_dock)
        self.inspector_dock.hide()

        self.appearance_dock = AppearanceDock(self.canvas.drawing_style, self)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.appearance_dock)
        self.appearance_dock.hide()
        self.appearance_dock.appearance_changed.connect(self._apply_appearance_settings)

        # === LEFT TOOLBAR (Drawing tools) ===
        self.toolbar = ChemusonToolbar()
        self.toolbar.setStyleSheet(TOOL_PALETTE_STYLESHEET)
        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, self.toolbar)

        # === RIGHT SYMBOLS TOOLBAR ===
        self.symbols_toolbar = SymbolPaletteToolbar(action_group=self.toolbar.action_group)
        self.symbols_toolbar.setStyleSheet(TOOL_PALETTE_STYLESHEET)
        self.addToolBar(Qt.ToolBarArea.RightToolBarArea, self.symbols_toolbar)
        self.symbols_toolbar.set_text_menu(
            [
                self.action_label_font,
                self.action_label_size_set,
                None,
                self.action_label_bold,
                self.action_label_italic,
                self.action_label_underline,
                None,
                self.action_label_subscript,
                self.action_label_superscript,
                None,
                self.action_label_size_up,
                self.action_label_size_down,
            ],
            [self.action_label_color_element, self.action_label_color_black],
        )

        # === MENU AND TOOLBARS ===
        self._create_menu_bar()
        self._create_main_toolbar()
        self._sync_label_menu_state()
        self._migrate_legacy_templates()
        self._refresh_template_views()
        
        # === TEXT FORMAT TOOLBAR ===
        self.text_toolbar = TextFormatToolbar()
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.text_toolbar)
        # Initially hide it, or show it only when text tool is active?
        # User requested it to be available. We can leave it visible or toggle it.
        # For now, visible is fine.

        # Text Toolbar Connections
        self.text_toolbar.format_changed.connect(self._on_text_format_changed)
        self.text_toolbar.color_changed.connect(self._on_text_color_changed)
        self.text_toolbar.alignment_changed.connect(self._on_text_alignment_changed)
        
        # === SIGNAL CONNECTIONS ===
        self._connect_undo_redo()
        
        # Connect toolbar to canvas
        self.toolbar.tool_changed.connect(self._on_tool_changed)
        self.toolbar.bond_palette_changed.connect(self._handle_bond_palette)
        self.toolbar.ring_palette_changed.connect(self._handle_ring_palette)
        self.toolbar.element_palette_changed.connect(self._handle_element_palette)
        self.toolbar.periodic_table_requested.connect(self._show_periodic_table)

        # Connect symbols dock to canvas
        self.symbols_toolbar.tool_changed.connect(self._on_tool_changed)
        self.symbols_toolbar.tool_changed.connect(self._update_status)

        # Sync defaults selected during toolbar init
        self._apply_toolbar_defaults_to_canvas(self.canvas)
        
        # === STATUS BAR ===
        self._total_charge_label = QLabel()
        self._iupac_name_label = QLabel("Nombre IUPAC: N/D")
        self._iupac_name_label.setToolTip("Nombre IUPAC-lite del documento activo")
        self.statusBar().addPermanentWidget(self._iupac_name_label, 1)
        self.statusBar().addPermanentWidget(self._total_charge_label)
        self._update_total_charge_indicator()
        self._update_iupac_name_indicator()
        self._update_status(self.canvas.state.active_tool)
        
        # Update status bar when tool changes
        self.toolbar.tool_changed.connect(self._update_status)
        self._set_active_canvas(self.canvas)
        # Chequeo diferido para no bloquear arranque ni UI.
        QTimer.singleShot(1200, self._maybe_check_updates_startup)

    def _create_actions(self) -> None:
        """Initialize all QActions for menus and toolbars."""
        # --- File Actions ---
        self.action_new = QAction("Nuevo", self)
        self.action_new.setShortcut(QKeySequence.StandardKey.New)
        self.action_new.triggered.connect(self._on_file_new)
        
        self.action_open = QAction("Abrir...", self)
        self.action_open.setShortcut(QKeySequence.StandardKey.Open)
        self.action_open.triggered.connect(self._on_file_open)
        
        self.action_save = QAction("Guardar", self)
        self.action_save.setShortcuts(
            [QKeySequence.StandardKey.Save, QKeySequence("Ctrl+G")]
        )
        self.action_save.triggered.connect(self._on_file_save)

        self.action_recovery_center = QAction("Centro de recuperación...", self)
        self.action_recovery_center.triggered.connect(self._on_open_recovery_center)
        
        self.action_quit = QAction("Salir", self)
        self.action_quit.setShortcut(QKeySequence.StandardKey.Quit)
        self.action_quit.triggered.connect(self.close)
        
        # --- Edit Actions ---
        self.action_undo = QAction("Deshacer", self)
        self.action_undo.setShortcut(QKeySequence.StandardKey.Undo)
        
        self.action_redo = QAction("Rehacer", self)
        self.action_redo.setShortcut(QKeySequence.StandardKey.Redo)
        
        self.action_copy = QAction("Copiar", self)
        self.action_copy.setShortcut(QKeySequence.StandardKey.Copy)
        
        self.action_paste = QAction("Pegar", self)
        self.action_paste.setShortcut(QKeySequence.StandardKey.Paste)
        
        # --- View Actions ---
        self.action_show_carbons = QAction("Mostrar carbonos", self)
        self.action_show_carbons.setCheckable(True)
        self.action_show_carbons.setChecked(False)
        self.action_show_carbons.triggered.connect(self._on_toggle_carbons)
        
        self.action_show_hydrogens = QAction("Mostrar hidrógenos", self)
        self.action_show_hydrogens.setCheckable(True)
        self.action_show_hydrogens.setChecked(False)
        self.action_show_hydrogens.triggered.connect(self._on_toggle_hydrogens)
        
        self.action_aromatic_circles = QAction("Aromáticos como círculos", self)
        self.action_aromatic_circles.setCheckable(True)
        self.action_aromatic_circles.setChecked(False)
        self.action_aromatic_circles.triggered.connect(self._on_toggle_aromatic_circles)
        
        self.action_zoom_in = QAction("Zoom +", self)
        self.action_zoom_in.setShortcut(QKeySequence.StandardKey.ZoomIn)
        self.action_zoom_in.triggered.connect(self._on_zoom_in)
        
        self.action_zoom_out = QAction("Zoom -", self)
        self.action_zoom_out.setShortcut(QKeySequence.StandardKey.ZoomOut)
        self.action_zoom_out.triggered.connect(self._on_zoom_out)
        
        self.action_zoom_reset = QAction("Zoom 100%", self)
        self.action_zoom_reset.setShortcut("Ctrl+0")
        self.action_zoom_reset.triggered.connect(self._on_zoom_reset)
        
        self.action_show_main_toolbar_aux = QAction(
            "Mostrar copiar/pegar/zoom en barra superior",
            self,
        )
        self.action_show_main_toolbar_aux.setCheckable(True)
        self.action_show_main_toolbar_aux.setChecked(False)
        self.action_show_main_toolbar_aux.toggled.connect(self._on_toggle_main_toolbar_aux)

        self.action_rules = QAction("Reglas", self)
        self.action_rules.setCheckable(True)
        self.action_rules.setChecked(False)
        self.action_rules.triggered.connect(self._on_toggle_rules)

        self.action_crosshair = QAction("Cuadrícula", self)
        self.action_crosshair.setCheckable(True)
        self.action_crosshair.setChecked(False)
        self.action_crosshair.triggered.connect(self._on_toggle_crosshair)

        self.action_numbering_enabled = QAction("Mostrar numeración", self)
        self.action_numbering_enabled.setCheckable(True)
        self.action_numbering_enabled.setChecked(False)
        self.action_numbering_enabled.triggered.connect(self._on_toggle_numbering)

        self.action_numbering_mode_atoms = QAction("Numerar átomos", self)
        self.action_numbering_mode_atoms.setCheckable(True)
        self.action_numbering_mode_atoms.triggered.connect(
            lambda checked=False: self._on_set_numbering_mode("atoms")
        )
        self.action_numbering_mode_structures = QAction("Numerar estructuras", self)
        self.action_numbering_mode_structures.setCheckable(True)
        self.action_numbering_mode_structures.triggered.connect(
            lambda checked=False: self._on_set_numbering_mode("structures")
        )
        self.action_numbering_mode_both = QAction("Numerar ambos", self)
        self.action_numbering_mode_both.setCheckable(True)
        self.action_numbering_mode_both.triggered.connect(
            lambda checked=False: self._on_set_numbering_mode("both")
        )
        self._numbering_mode_group = QActionGroup(self)
        self._numbering_mode_group.setExclusive(True)
        self._numbering_mode_group.addAction(self.action_numbering_mode_atoms)
        self._numbering_mode_group.addAction(self.action_numbering_mode_structures)
        self._numbering_mode_group.addAction(self.action_numbering_mode_both)
        self.action_numbering_mode_atoms.setChecked(True)

        self.action_numbering_recalculate = QAction("Recalcular numeración", self)
        self.action_numbering_recalculate.triggered.connect(self._on_recalculate_numbering)

        self.action_numbering_export = QAction("Incluir numeración en exportación", self)
        self.action_numbering_export.setCheckable(True)
        self.action_numbering_export.setChecked(True)
        self.action_numbering_export.triggered.connect(self._on_toggle_numbering_export)

        self.action_clean_2d = QAction("Limpiar 2D", self)
        self.action_clean_2d.triggered.connect(self._on_clean_2d)
        self.action_clean_2d_full = QAction("Limpiar 2D (1 paso)", self)
        self.action_clean_2d_full.setShortcut(QKeySequence("Ctrl+K"))
        self.action_clean_2d_full.triggered.connect(self._on_clean_2d_full)

        self.action_template_linear_chain = QAction("Cadena lineal", self)
        self.action_template_linear_chain.triggered.connect(self._on_insert_linear_chain)
        self.action_template_new_category = QAction("Nueva categoría...", self)
        self.action_template_new_category.triggered.connect(self._on_new_template_category)
        self.action_template_import_library = QAction("Importar biblioteca...", self)
        self.action_template_import_library.triggered.connect(self._on_import_template_library)
        self.action_template_export_library = QAction("Exportar biblioteca...", self)
        self.action_template_export_library.triggered.connect(self._on_export_template_library)

        # --- Canvas Size Actions ---
        self.action_canvas_size_letter_portrait = QAction("Carta (vertical)", self)
        self.action_canvas_size_letter_portrait.triggered.connect(
            lambda: self._set_canvas_size(816, 1056)
        )
        self.action_canvas_size_letter_landscape = QAction("Carta (horizontal)", self)
        self.action_canvas_size_letter_landscape.triggered.connect(
            lambda: self._set_canvas_size(1056, 816)
        )

        self.action_canvas_size_a4_portrait = QAction("A4 (vertical)", self)
        self.action_canvas_size_a4_portrait.triggered.connect(
            lambda: self._set_canvas_size(794, 1123)
        )
        self.action_canvas_size_a4_landscape = QAction("A4 (horizontal)", self)
        self.action_canvas_size_a4_landscape.triggered.connect(
            lambda: self._set_canvas_size(1123, 794)
        )

        self.action_canvas_size_a3_portrait = QAction("A3 (vertical)", self)
        self.action_canvas_size_a3_portrait.triggered.connect(
            lambda: self._set_canvas_size(1123, 1587)
        )
        self.action_canvas_size_a3_landscape = QAction("A3 (horizontal)", self)
        self.action_canvas_size_a3_landscape.triggered.connect(
            lambda: self._set_canvas_size(1587, 1123)
        )

        self.action_canvas_size_custom = QAction("Personalizado...", self)
        self.action_canvas_size_custom.triggered.connect(self._on_canvas_custom_size)



        self.action_style = QAction("Dimensiones del dibujo...", self)
        self.action_style.triggered.connect(self._on_style_dialog)

        self.action_import_smiles = QAction("Importar SMILES...", self)
        self.action_import_smiles.triggered.connect(self._on_import_smiles)

        self.action_export_smiles = QAction("Exportar SMILES...", self)
        self.action_export_smiles.triggered.connect(self._on_export_smiles)

        self.action_draw_smiles = QAction("Dibujar desde SMILES...", self)
        self.action_draw_smiles.triggered.connect(self._on_import_smiles)

        # --- Analysis Actions ---
        self.action_analysis_name = QAction("Nombre (SMILES)", self)
        self.action_analysis_name.triggered.connect(lambda: self.canvas.run_analysis("name"))
        self.action_analysis_formula = QAction("Fórmula química", self)
        self.action_analysis_formula.triggered.connect(lambda: self.canvas.run_analysis("formula"))
        self.action_analysis_exact = QAction("Masa exacta", self)
        self.action_analysis_exact.triggered.connect(lambda: self.canvas.run_analysis("exact"))
        self.action_analysis_weight = QAction("Peso molecular", self)
        self.action_analysis_weight.triggered.connect(lambda: self.canvas.run_analysis("weight"))
        self.action_analysis_mz = QAction("m/z", self)
        self.action_analysis_mz.triggered.connect(lambda: self.canvas.run_analysis("mz"))
        self.action_analysis_elemental = QAction("Análisis elemental", self)
        self.action_analysis_elemental.triggered.connect(lambda: self.canvas.run_analysis("elemental"))
        self.action_analysis_all = QAction("Todo", self)
        self.action_analysis_all.triggered.connect(lambda: self.canvas.run_analysis("all"))

        # --- Rotation Actions ---
        self.action_rotate_left = QAction("Girar 90° a la izquierda", self)
        self.action_rotate_left.triggered.connect(lambda: self._on_rotate_selection(-90.0))

        self.action_rotate_right = QAction("Girar 90° a la derecha", self)
        self.action_rotate_right.triggered.connect(lambda: self._on_rotate_selection(90.0))

        self.action_flip_horizontal = QAction("Giro 180° horizontal", self)
        self.action_flip_horizontal.triggered.connect(self._on_flip_horizontal)

        self.action_flip_vertical = QAction("Giro 180° vertical", self)
        self.action_flip_vertical.triggered.connect(self._on_flip_vertical)

        self.action_branch_rotate_minus = QAction(
            f"Girar rama -{int(BRANCH_ROTATION_STEP_DEG)}°",
            self,
        )
        self.action_branch_rotate_minus.setShortcut(QKeySequence("Ctrl+Alt+Left"))
        self.action_branch_rotate_minus.triggered.connect(
            lambda: self._on_rotate_branch(-BRANCH_ROTATION_STEP_DEG)
        )

        self.action_branch_rotate_plus = QAction(
            f"Girar rama +{int(BRANCH_ROTATION_STEP_DEG)}°",
            self,
        )
        self.action_branch_rotate_plus.setShortcut(QKeySequence("Ctrl+Alt+Right"))
        self.action_branch_rotate_plus.triggered.connect(
            lambda: self._on_rotate_branch(BRANCH_ROTATION_STEP_DEG)
        )

        self.action_branch_invert = QAction("Invertir rama (180°)", self)
        self.action_branch_invert.setShortcut(QKeySequence("Ctrl+Alt+I"))
        self.action_branch_invert.triggered.connect(self._on_invert_branch)

        self.action_branch_auto_arrange = QAction("Autoacomodar rama", self)
        self.action_branch_auto_arrange.setShortcut(QKeySequence("Ctrl+Alt+A"))
        self.action_branch_auto_arrange.triggered.connect(self._on_auto_arrange_branch)

        self.action_fragment_pivot_set = QAction("Definir átomo pivote desde selección", self)
        self.action_fragment_pivot_set.triggered.connect(self._on_set_fragment_pivot)

        self.action_fragment_pivot_clear = QAction("Limpiar átomo pivote", self)
        self.action_fragment_pivot_clear.triggered.connect(self._on_clear_fragment_pivot)

        self.action_fragment_rotate_minus = QAction(
            f"Girar fragmento -{int(FRAGMENT_ROTATION_STEP_DEG)}°",
            self,
        )
        self.action_fragment_rotate_minus.triggered.connect(
            lambda: self._on_rotate_fragment(-FRAGMENT_ROTATION_STEP_DEG)
        )

        self.action_fragment_rotate_plus = QAction(
            f"Girar fragmento +{int(FRAGMENT_ROTATION_STEP_DEG)}°",
            self,
        )
        self.action_fragment_rotate_plus.triggered.connect(
            lambda: self._on_rotate_fragment(FRAGMENT_ROTATION_STEP_DEG)
        )

        self.action_fragment_invert = QAction("Invertir fragmento (180°)", self)
        self.action_fragment_invert.triggered.connect(self._on_invert_fragment)

        self.action_scale_selection = QAction("Redimensionar selección...", self)
        self.action_scale_selection.setShortcut(QKeySequence("Ctrl+Alt+S"))
        self.action_scale_selection.triggered.connect(
            lambda checked=False: self.canvas.open_selection_scale_dialog()
        )

        # --- Bond Thickness Actions ---
        self.action_bond_thickness_up = QAction("Aumentar grosor de enlace/flecha/corchete", self)
        self.action_bond_thickness_up.setShortcut(QKeySequence("Ctrl+Shift+Up"))
        self.action_bond_thickness_up.triggered.connect(self._on_bond_thickness_up)

        self.action_bond_thickness_down = QAction("Reducir grosor de enlace/flecha/corchete", self)
        self.action_bond_thickness_down.setShortcut(QKeySequence("Ctrl+Shift+Down"))
        self.action_bond_thickness_down.triggered.connect(self._on_bond_thickness_down)

        self.action_bond_thickness_reset = QAction("Restablecer grosor de enlace/flecha/corchete", self)
        self.action_bond_thickness_reset.setShortcut(QKeySequence("Ctrl+Shift+0"))
        self.action_bond_thickness_reset.triggered.connect(self._on_bond_thickness_reset)

        # --- Text Actions ---
        self.action_label_font = QAction("Fuente de etiquetas...", self)
        self.action_label_font.triggered.connect(self._on_label_font)

        self.action_label_size_set = QAction("Tamaño...", self)
        self.action_label_size_set.triggered.connect(self._on_label_font_size_dialog)

        self.action_label_bold = QAction("Negrita", self)
        self.action_label_bold.setCheckable(True)
        self.action_label_bold.triggered.connect(self._on_label_bold)

        self.action_label_italic = QAction("Cursiva", self)
        self.action_label_italic.setCheckable(True)
        self.action_label_italic.triggered.connect(self._on_label_italic)

        self.action_label_underline = QAction("Subrayado", self)
        self.action_label_underline.setCheckable(True)
        self.action_label_underline.triggered.connect(self._on_label_underline)

        self.action_label_subscript = QAction("Subíndice...", self)
        self.action_label_subscript.triggered.connect(self._on_label_subscript)

        self.action_label_superscript = QAction("Superíndice...", self)
        self.action_label_superscript.triggered.connect(self._on_label_superscript)

        self.action_label_size_up = QAction("Aumentar tamaño", self)
        self.action_label_size_up.triggered.connect(lambda: self._on_label_font_size(1.0))

        self.action_label_size_down = QAction("Reducir tamaño", self)
        self.action_label_size_down.triggered.connect(lambda: self._on_label_font_size(-1.0))

        self.action_text_align_left = QAction("Alinear a la izquierda", self)
        self.action_text_align_left.triggered.connect(
            lambda: self.canvas.update_text_alignment(Qt.AlignmentFlag.AlignLeft)
        )
        self.action_text_align_center = QAction("Centrar", self)
        self.action_text_align_center.triggered.connect(
            lambda: self.canvas.update_text_alignment(Qt.AlignmentFlag.AlignHCenter)
        )
        self.action_text_align_justify = QAction("Justificar", self)
        self.action_text_align_justify.triggered.connect(
            lambda: self.canvas.update_text_alignment(Qt.AlignmentFlag.AlignJustify)
        )

        self.action_label_color_element = QAction("Por elemento", self)
        self.action_label_color_element.setCheckable(True)
        self.action_label_color_element.triggered.connect(
            lambda checked=False: self._on_label_color_mode(True)
        )

        self.action_label_color_black = QAction("Negro", self)
        self.action_label_color_black.setCheckable(True)
        self.action_label_color_black.triggered.connect(
            lambda checked=False: self._on_label_color_mode(False)
        )

        self._label_color_group = QActionGroup(self)
        self._label_color_group.setExclusive(True)
        self._label_color_group.addAction(self.action_label_color_element)
        self._label_color_group.addAction(self.action_label_color_black)

        # --- Help / Updates ---
        self.action_check_updates_now = QAction("Buscar actualizaciones...", self)
        self.action_check_updates_now.triggered.connect(self._on_check_updates_now)

    def _create_menu_bar(self) -> None:
        """Create the main menu bar with all menus."""
        menubar = self.menuBar()
        
        # === Archivo (File) ===
        file_menu = menubar.addMenu("Archivo")
        file_menu.addAction(self.action_new)
        file_menu.addAction(self.action_open)
        file_menu.addAction(self.action_save)
        self.recent_menu = file_menu.addMenu("Archivos recientes")
        self._update_recent_menu()
        file_menu.addAction(self.action_recovery_center)
        file_menu.addSeparator()
        
        export_menu = file_menu.addMenu("Exportar como")
        self.action_export_png = QAction("PNG...", self)
        self.action_export_png.triggered.connect(lambda: self._on_export("png"))
        export_menu.addAction(self.action_export_png)
        
        self.action_export_svg = QAction("SVG...", self)
        self.action_export_svg.triggered.connect(lambda: self._on_export("svg"))
        export_menu.addAction(self.action_export_svg)
        
        self.action_export_pdf = QAction("PDF...", self)
        self.action_export_pdf.triggered.connect(lambda: self._on_export("pdf"))
        export_menu.addAction(self.action_export_pdf)
        
        file_menu.addSeparator()
        file_menu.addAction(self.action_quit)
        
        # === Editar (Edit) ===
        edit_menu = menubar.addMenu("Editar")
        edit_menu.addAction(self.action_undo)
        edit_menu.addAction(self.action_redo)
        edit_menu.addSeparator()
        edit_menu.addAction(self.action_copy)
        
        copy_as_menu = edit_menu.addMenu("Copiar como")
        self.action_copy_smiles = QAction("SMILES", self)
        self.action_copy_smiles.triggered.connect(lambda: self._on_copy_as("smiles"))
        copy_as_menu.addAction(self.action_copy_smiles)
        
        self.action_copy_molfile = QAction("Molfile", self)
        self.action_copy_molfile.triggered.connect(lambda: self._on_copy_as("molfile"))
        copy_as_menu.addAction(self.action_copy_molfile)
        
        self.action_copy_inchi = QAction("InChI", self)
        self.action_copy_inchi.triggered.connect(lambda: self._on_copy_as("inchi"))
        copy_as_menu.addAction(self.action_copy_inchi)
        
        edit_menu.addAction(self.action_paste)
        edit_menu.addSeparator()

        rotate_menu = edit_menu.addMenu("Rotar")
        rotate_menu.addAction(self.action_rotate_left)
        rotate_menu.addAction(self.action_rotate_right)
        rotate_menu.addSeparator()
        rotate_menu.addAction(self.action_flip_horizontal)
        rotate_menu.addAction(self.action_flip_vertical)
        rotate_menu.addSeparator()
        rotate_menu.addAction(self.action_branch_rotate_minus)
        rotate_menu.addAction(self.action_branch_rotate_plus)
        rotate_menu.addAction(self.action_branch_invert)
        rotate_menu.addAction(self.action_branch_auto_arrange)
        rotate_menu.addSeparator()
        self.fragment_rotate_menu = rotate_menu.addMenu("Fragmento con pivote")
        self.fragment_rotate_menu.addAction(self.action_fragment_pivot_set)
        self.fragment_rotate_menu.addAction(self.action_fragment_pivot_clear)
        self.fragment_rotate_menu.addSeparator()
        self.fragment_rotate_menu.addAction(self.action_fragment_rotate_minus)
        self.fragment_rotate_menu.addAction(self.action_fragment_rotate_plus)
        self.fragment_rotate_menu.addAction(self.action_fragment_invert)

        edit_menu.addAction(self.action_scale_selection)

        edit_menu.addSeparator()
        bond_thickness_menu = edit_menu.addMenu("Grosor de enlace/flecha/corchete")
        bond_thickness_menu.addAction(self.action_bond_thickness_up)
        bond_thickness_menu.addAction(self.action_bond_thickness_down)
        bond_thickness_menu.addAction(self.action_bond_thickness_reset)

        edit_menu.addSeparator()

        self.action_preferences = QAction("Preferencias...", self)
        self.action_preferences.triggered.connect(self._on_preferences)
        edit_menu.addAction(self.action_preferences)

        # === Ver (View) ===
        view_menu = menubar.addMenu("Ver")
        view_menu.addAction(self.action_show_carbons)
        view_menu.addAction(self.action_show_hydrogens)
        view_menu.addSeparator()
        view_menu.addAction(self.action_aromatic_circles)
        view_menu.addSeparator()
        view_menu.addAction(self.action_style)
        view_menu.addSeparator()
        view_menu.addAction(self.action_zoom_in)
        view_menu.addAction(self.action_zoom_out)
        view_menu.addAction(self.action_zoom_reset)
        view_menu.addSeparator()
        view_menu.addAction(self.action_show_main_toolbar_aux)
        view_menu.addSeparator()
        view_menu.addAction(self.action_rules)
        view_menu.addAction(self.action_crosshair)
        view_menu.addSeparator()
        numbering_menu = view_menu.addMenu("Numeración")
        numbering_menu.addAction(self.action_numbering_enabled)
        numbering_menu.addSeparator()
        numbering_menu.addAction(self.action_numbering_mode_atoms)
        numbering_menu.addAction(self.action_numbering_mode_structures)
        numbering_menu.addAction(self.action_numbering_mode_both)
        numbering_menu.addSeparator()
        numbering_menu.addAction(self.action_numbering_recalculate)
        numbering_menu.addAction(self.action_numbering_export)
        view_menu.addSeparator()
        canvas_size_menu = view_menu.addMenu("Tamaño de lienzo")
        canvas_size_menu.addAction(self.action_canvas_size_letter_portrait)
        canvas_size_menu.addAction(self.action_canvas_size_letter_landscape)
        canvas_size_menu.addSeparator()
        canvas_size_menu.addAction(self.action_canvas_size_a4_portrait)
        canvas_size_menu.addAction(self.action_canvas_size_a4_landscape)
        canvas_size_menu.addSeparator()
        canvas_size_menu.addAction(self.action_canvas_size_a3_portrait)
        canvas_size_menu.addAction(self.action_canvas_size_a3_landscape)
        canvas_size_menu.addSeparator()
        canvas_size_menu.addAction(self.action_canvas_size_custom)
        view_menu.addSeparator()
        
        # Docks visibility
        view_menu.addAction(self.symbols_toolbar.toggleViewAction())
        view_menu.addAction(self.templates_dock.toggleViewAction())
        view_menu.addAction(self.inspector_dock.toggleViewAction())
        view_menu.addAction(self.appearance_dock.toggleViewAction())
        
        # === Estructura (Structure) ===
        structure_menu = menubar.addMenu("Estructura")
        structure_menu.addAction(self.action_clean_2d)
        structure_menu.addAction(self.action_clean_2d_full)
        structure_menu.addSeparator()
        structure_menu.addSeparator()
        
        # Dynamic Templates Menu
        self.templates_menu = structure_menu.addMenu("Plantillas")
        
        # Save Template Action
        self.action_save_template = QAction("Guardar selección como plantilla...", self)
        self.action_save_template.triggered.connect(self._on_save_template)
        self._refresh_templates_menu()
        structure_menu.addAction(self.action_save_template)
        structure_menu.addAction(self.action_template_import_library)
        structure_menu.addAction(self.action_template_export_library)
        
        structure_menu.addSeparator()
        structure_menu.addAction(self.action_import_smiles)
        structure_menu.addAction(self.action_export_smiles)
        structure_menu.addSeparator()
        analysis_menu = structure_menu.addMenu("Análisis")
        analysis_menu.addAction(self.action_analysis_name)
        analysis_menu.addAction(self.action_analysis_formula)
        analysis_menu.addAction(self.action_analysis_exact)
        analysis_menu.addAction(self.action_analysis_weight)
        analysis_menu.addAction(self.action_analysis_mz)
        analysis_menu.addAction(self.action_analysis_elemental)
        analysis_menu.addSeparator()
        analysis_menu.addAction(self.action_analysis_all)

        # === Reacción (Reaction) ===
        reaction_menu = menubar.addMenu("Reacción")
        placeholder_reaction = QAction("Próximamente", self)
        placeholder_reaction.setEnabled(False)
        reaction_menu.addAction(placeholder_reaction)
        
        # === Ayuda (Help) ===
        help_menu = menubar.addMenu("Ayuda")
        
        self.action_quick_start = QAction("Guía rápida...", self)
        self.action_quick_start.triggered.connect(self._on_quick_start)
        help_menu.addAction(self.action_quick_start)
        help_menu.addAction(self.action_check_updates_now)
        
        help_menu.addSeparator()
        
        self.action_about = QAction("Acerca de Chemuson...", self)
        self.action_about.triggered.connect(self._on_about)
        help_menu.addAction(self.action_about)
    
    # -------------------------------------------------------------------------
    # Main Toolbar
    # -------------------------------------------------------------------------
    def _create_main_toolbar(self) -> None:
        """Create the main horizontal toolbar with common actions."""
        self.main_toolbar = QToolBar("Principal")
        self.main_toolbar.setMovable(False)
        self.main_toolbar.setFloatable(False)
        self.main_toolbar.setIconSize(QSize(24, 24))
        self.main_toolbar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.main_toolbar)
        
        # Set icons for actions with fallbacks where possible
        self.action_new.setIcon(QIcon.fromTheme("document-new", QIcon()))
        self.action_open.setIcon(QIcon.fromTheme("document-open", QIcon()))
        self.action_save.setIcon(QIcon.fromTheme("document-save", QIcon()))
        self.action_undo.setIcon(QIcon.fromTheme("edit-undo", QIcon()))
        self.action_redo.setIcon(QIcon.fromTheme("edit-redo", QIcon()))
        self.action_copy.setIcon(QIcon.fromTheme("edit-copy", QIcon()))
        self.action_paste.setIcon(QIcon.fromTheme("edit-paste", QIcon()))
        
        from chemuson.gui.icons import draw_generic_icon
        self.action_zoom_in.setIcon(draw_generic_icon("zoom_in"))
        self.action_zoom_out.setIcon(draw_generic_icon("zoom_out"))
        self.action_rotate_left.setIcon(draw_generic_icon("rotate_left"))
        self.action_rotate_right.setIcon(draw_generic_icon("rotate_right"))
        self.action_flip_horizontal.setIcon(draw_generic_icon("flip_horizontal"))
        self.action_flip_vertical.setIcon(draw_generic_icon("flip_vertical"))
        self.action_branch_rotate_minus.setIcon(draw_generic_icon("rotate_left"))
        self.action_branch_rotate_plus.setIcon(draw_generic_icon("rotate_right"))
        self.action_branch_invert.setIcon(draw_generic_icon("flip_horizontal"))
        self.action_branch_auto_arrange.setIcon(QIcon.fromTheme("edit-clear", QIcon()))
        self.action_clean_2d.setIcon(QIcon.fromTheme("edit-clear", QIcon()))
        from chemuson.gui.icons import draw_atom_icon
        self.action_draw_smiles.setIcon(draw_atom_icon("SMI"))
        
        # File actions
        self.main_toolbar.addAction(self.action_new)
        self.main_toolbar.addAction(self.action_open)
        self.main_toolbar.addAction(self.action_save)
        self.main_toolbar.addSeparator()
        
        # Edit actions
        self.main_toolbar.addAction(self.action_undo)
        self.main_toolbar.addAction(self.action_redo)
        self.main_toolbar.addSeparator()

        # Rotate actions
        self.main_toolbar.addAction(self.action_rotate_left)
        self.main_toolbar.addAction(self.action_rotate_right)
        self.main_toolbar.addAction(self.action_flip_horizontal)
        self.main_toolbar.addAction(self.action_flip_vertical)
        self.main_toolbar.addSeparator()
        
        # Structure actions
        self.main_toolbar.addAction(self.action_clean_2d)
        self.main_toolbar.addSeparator()
        self.main_toolbar.addAction(self.action_draw_smiles)

        # Optional top-toolbar actions (copy/paste/zoom +/-).
        self._main_toolbar_aux_separator_1 = QAction(self)
        self._main_toolbar_aux_separator_1.setSeparator(True)
        self._main_toolbar_aux_separator_2 = QAction(self)
        self._main_toolbar_aux_separator_2.setSeparator(True)
        self._set_main_toolbar_aux_visible(self.action_show_main_toolbar_aux.isChecked())

    def _set_main_toolbar_aux_visible(self, visible: bool) -> None:
        """Muestra/oculta copiar, pegar y zoom +/- en la barra superior."""
        if not hasattr(self, "main_toolbar"):
            return
        aux_actions = (
            self.action_copy,
            self.action_paste,
            self._main_toolbar_aux_separator_1,
            self.action_zoom_in,
            self.action_zoom_out,
            self._main_toolbar_aux_separator_2,
        )
        current = self.main_toolbar.actions()
        if visible:
            anchor = self.action_rotate_left if self.action_rotate_left in current else None
            for action in aux_actions:
                if action in self.main_toolbar.actions():
                    continue
                if anchor is not None:
                    self.main_toolbar.insertAction(anchor, action)
                else:
                    self.main_toolbar.addAction(action)
        else:
            for action in aux_actions:
                if action in self.main_toolbar.actions():
                    self.main_toolbar.removeAction(action)

    def _on_toggle_main_toolbar_aux(self, checked: bool) -> None:
        """Actualiza visibilidad de atajos auxiliares en barra superior."""
        self._set_main_toolbar_aux_visible(bool(checked))

    def _create_document_tab(self, make_current: bool = False) -> ChemusonCanvas:
        """Crea una nueva pestaña con su propio canvas independiente."""
        canvas = ChemusonCanvas()
        tab_title = (
            "Sin título"
            if self._untitled_counter == 1
            else f"Sin título {self._untitled_counter}"
        )
        self._untitled_counter += 1
        self._canvas_file_paths[canvas] = None
        self._canvas_tab_titles[canvas] = tab_title
        self._apply_default_numbering_to_canvas(canvas)
        self._apply_default_naming_to_canvas(canvas)
        index = self.tabs.addTab(canvas, tab_title)
        self.tabs.setTabToolTip(index, "Documento sin guardar")
        autosave_manager = AutosaveManager(self, canvas, canvas.undo_stack)
        autosave_manager.start()
        self._canvas_autosave_managers[canvas] = autosave_manager
        canvas.undo_stack.cleanChanged.connect(
            lambda clean, c=canvas: self._on_canvas_clean_state_changed(c, bool(clean))
        )
        self._update_tab_title(canvas)
        if make_current:
            self.tabs.setCurrentIndex(index)
        return canvas

    def _apply_default_numbering_to_canvas(self, canvas: ChemusonCanvas) -> None:
        """Aplica preferencias globales de numeración a un nuevo documento."""
        canvas.set_numbering_mode(self._numbering_default_mode)
        canvas.set_numbering_enabled(False)
        canvas.set_numbering_include_in_export(self._numbering_default_include_export)

    def _apply_default_naming_to_canvas(self, canvas: ChemusonCanvas) -> None:
        """Aplica preferencias globales de nomenclatura al documento."""
        canvas.name_advanced_enabled = bool(self._name_advanced_default)
        canvas.name_rdkit_isolated = bool(self._name_rdkit_isolated_default)

    def _set_canvas_file_path(self, canvas: ChemusonCanvas, filepath: Optional[str]) -> None:
        """Asigna ruta de archivo a una pestaña y actualiza su título."""
        clean_path = os.path.abspath(filepath) if filepath else None
        self._canvas_file_paths[canvas] = clean_path
        autosave_manager = self._canvas_autosave_managers.get(canvas)
        if autosave_manager is not None:
            autosave_manager.set_original_path(clean_path)
        if clean_path:
            self._canvas_tab_titles[canvas] = os.path.basename(clean_path)
        if canvas is self.canvas:
            self._current_file_path = clean_path
        self._update_tab_title(canvas)

    def _on_canvas_clean_state_changed(
        self,
        canvas: ChemusonCanvas,
        clean: Optional[bool] = None,
    ) -> None:
        """Actualiza el título de pestaña cuando cambia estado clean/dirty."""
        if not self._tabs_alive():
            return
        autosave_manager = self._canvas_autosave_managers.get(canvas)
        if autosave_manager is not None:
            is_clean = canvas.undo_stack.isClean() if clean is None else bool(clean)
            if is_clean:
                autosave_manager.cancel_debounce()
            else:
                autosave_manager.restart_debounce()
        self._update_tab_title(canvas)

    def _tabs_alive(self) -> bool:
        """Indica si el widget de pestañas sigue disponible."""
        tabs = getattr(self, "tabs", None)
        if tabs is None:
            return False
        try:
            tabs.count()
        except RuntimeError:
            return False
        return True

    def _canvas_from_tab_index(self, index: int) -> Optional[ChemusonCanvas]:
        """Obtiene el canvas asociado al índice de pestaña."""
        if not self._tabs_alive():
            return None
        if index < 0:
            return None
        widget = self.tabs.widget(index)
        return widget if isinstance(widget, ChemusonCanvas) else None

    def _update_tab_title(self, canvas: ChemusonCanvas) -> None:
        """Actualiza título y tooltip de la pestaña asociada al canvas."""
        if not self._tabs_alive():
            return
        try:
            index = self.tabs.indexOf(canvas)
            if index < 0:
                return
            base_title = self._canvas_tab_titles.get(canvas, "Sin título")
            dirty_suffix = " *" if not canvas.undo_stack.isClean() else ""
            self.tabs.setTabText(index, f"{base_title}{dirty_suffix}")
            path = self._canvas_file_paths.get(canvas)
            self.tabs.setTabToolTip(index, path or "Documento sin guardar")
        except RuntimeError:
            # Durante cierre, Qt puede destruir tabs/canvas antes de drenar señales.
            return

    def _on_tab_changed(self, index: int) -> None:
        """Activa el canvas correspondiente al cambiar de pestaña."""
        canvas = self._canvas_from_tab_index(index)
        if canvas is None:
            return
        self._set_active_canvas(canvas, clear_tool_selection=True)

    def _on_tab_close_requested(self, index: int) -> None:
        """Maneja cierre por botón X de una pestaña."""
        canvas = self._canvas_from_tab_index(index)
        if canvas is None:
            return
        self._close_canvas_tab(canvas)

    def _close_canvas_tab(self, canvas: ChemusonCanvas) -> bool:
        """Cierra una pestaña si no hay cambios pendientes o el usuario confirma."""
        if not self._confirm_discard_changes(canvas):
            return False
        index = self.tabs.indexOf(canvas)
        if index < 0:
            return False
        was_active = canvas is self.canvas
        self.tabs.removeTab(index)
        autosave_manager = self._canvas_autosave_managers.pop(canvas, None)
        if autosave_manager is not None:
            autosave_manager.stop()
        self._canvas_file_paths.pop(canvas, None)
        self._canvas_tab_titles.pop(canvas, None)
        if self._active_canvas_connected is canvas:
            self._disconnect_canvas_signals(canvas)
            self._active_canvas_connected = None
        canvas.deleteLater()

        if self.tabs.count() == 0:
            created = self._create_document_tab(make_current=True)
            self._set_active_canvas(created)
            self.statusBar().showMessage("Se creó un documento nuevo.")
            return True

        if was_active:
            current = self._canvas_from_tab_index(self.tabs.currentIndex())
            if current is not None:
                self._set_active_canvas(current)
        return True

    def _connect_canvas_signals(self, canvas: ChemusonCanvas) -> None:
        """Conecta señales del canvas activo a la UI principal."""
        canvas.selection_changed.connect(self._on_selection_changed)
        canvas.undo_stack.canUndoChanged.connect(self.action_undo.setEnabled)
        canvas.undo_stack.canRedoChanged.connect(self.action_redo.setEnabled)
        canvas.undo_stack.indexChanged.connect(self._on_active_undo_index_changed)

    def _disconnect_canvas_signals(self, canvas: ChemusonCanvas) -> None:
        """Desconecta señales del canvas activo anterior."""
        for signal, slot in (
            (canvas.selection_changed, self._on_selection_changed),
            (canvas.undo_stack.canUndoChanged, self.action_undo.setEnabled),
            (canvas.undo_stack.canRedoChanged, self.action_redo.setEnabled),
            (canvas.undo_stack.indexChanged, self._on_active_undo_index_changed),
        ):
            try:
                signal.disconnect(slot)
            except (TypeError, RuntimeError):
                pass

    def _sync_view_actions_from_canvas(self) -> None:
        """Sincroniza checks de menú Ver con el estado del canvas activo."""
        values = (
            (self.action_show_carbons, bool(self.canvas.state.show_implicit_carbons)),
            (self.action_show_hydrogens, bool(self.canvas.state.show_implicit_hydrogens)),
            (self.action_aromatic_circles, bool(self.canvas.state.use_aromatic_circles)),
            (self.action_rules, bool(self.canvas.show_rulers)),
            (self.action_crosshair, bool(self.canvas.show_grid)),
        )
        for action, checked in values:
            was_blocked = action.blockSignals(True)
            action.setChecked(checked)
            action.blockSignals(was_blocked)

    def _set_active_canvas(
        self,
        canvas: ChemusonCanvas,
        clear_tool_selection: bool = False,
    ) -> None:
        """Activa un canvas de pestaña y actualiza conexiones/estado de UI."""
        if canvas is None:
            return
        if self._active_canvas_connected is not None and self._active_canvas_connected is not canvas:
            self._disconnect_canvas_signals(self._active_canvas_connected)
        self.canvas = canvas
        self._current_file_path = self._canvas_file_paths.get(canvas)
        if self._active_canvas_connected is not canvas:
            self._connect_canvas_signals(canvas)
            self._active_canvas_connected = canvas

        self.action_undo.setEnabled(self.canvas.undo_stack.canUndo())
        self.action_redo.setEnabled(self.canvas.undo_stack.canRedo())
        self._sync_view_actions_from_canvas()
        self._sync_numbering_actions()
        self._sync_fragment_pivot_actions()
        self._sync_label_menu_state()
        if hasattr(self, "appearance_dock"):
            cap_mode = (
                "round"
                if self.canvas.drawing_style.cap_style == Qt.PenCapStyle.RoundCap
                else "flat"
            )
            self.appearance_dock.set_bond_caps(cap_mode)
        self._update_total_charge_indicator()
        self._update_iupac_name_indicator()
        if clear_tool_selection:
            self._clear_active_tool_selection()
        else:
            self._update_status(self.canvas.state.active_tool)
        if hasattr(self, "inspector_dock"):
            self.inspector_dock.update_selection(0, 0, 0, {})

    def _clear_active_tool_selection(self) -> None:
        """Limpia herramienta activa al cambiar de pestaña."""
        self._current_tool_id = "tool_none"
        self.toolbar.clear_tool_selection()
        self.canvas.set_current_tool(self._current_tool_id)
        self._update_status(self._current_tool_id)

    def _on_active_undo_index_changed(self, _index: int) -> None:
        """Actualiza widgets dependientes del estado del undo activo."""
        self._update_total_charge_indicator()
        self._update_iupac_name_indicator()
        self._sync_fragment_pivot_actions()
        self._update_tab_title(self.canvas)

    def _update_iupac_name_indicator(self) -> None:
        """Refresca el indicador de nombre IUPAC en la barra de estado."""
        if not hasattr(self, "_iupac_name_label"):
            return
        try:
            name = self.canvas.current_iupac_name()
        except Exception:
            name = "N/D"
        self._iupac_name_label.setText(f"Nombre IUPAC: {name or 'N/D'}")

    def _on_undo(self) -> None:
        """Deshace en la pestaña activa."""
        self.canvas.undo_stack.undo()

    def _on_redo(self) -> None:
        """Rehace en la pestaña activa."""
        self.canvas.undo_stack.redo()

    def _on_copy(self) -> None:
        """Copia desde la pestaña activa."""
        self.canvas.copy_to_clipboard()

    def _on_paste(self) -> None:
        """Pega en la pestaña activa."""
        self.canvas.paste_from_clipboard()

    def _on_text_format_changed(
        self,
        family: str,
        size: int,
        bold: bool,
        italic: bool,
        underline: bool,
        sub: bool,
        sup: bool,
        property_name: str,
    ) -> None:
        """Propaga cambios de formato de texto al canvas activo."""
        self.canvas.update_text_format(
            family,
            size,
            bold,
            italic,
            underline,
            sub,
            sup,
            property_name,
        )

    def _on_text_color_changed(self, color: QColor) -> None:
        """Propaga cambio de color de texto al canvas activo."""
        self.canvas.update_text_color(color)

    def _on_text_alignment_changed(self, alignment: Qt.AlignmentFlag) -> None:
        """Propaga cambio de alineación al canvas activo."""
        self.canvas.update_text_alignment(alignment)

    def _on_tool_changed(self, tool_id: str) -> None:
        """Actualiza herramienta activa en el documento de la pestaña actual."""
        self._current_tool_id = tool_id
        self.canvas.set_current_tool(tool_id)

    def _apply_toolbar_defaults_to_canvas(self, canvas: ChemusonCanvas) -> None:
        """Aplica selección actual de paletas al canvas indicado."""
        bond_spec = self.toolbar.current_bond_spec()
        canvas.state.active_bond_order = bond_spec.get("order", 1)
        canvas.state.active_bond_style = bond_spec.get("style", canvas.state.active_bond_style)
        canvas.state.active_bond_stereo = bond_spec.get("stereo", canvas.state.active_bond_stereo)
        canvas.state.active_bond_mode = bond_spec.get("mode", "increment")
        canvas.state.active_bond_aromatic = bond_spec.get("aromatic", False)

        ring_spec = self.toolbar.current_ring_spec()
        canvas.state.active_ring_size = ring_spec.get("size", canvas.state.active_ring_size)
        canvas.state.active_ring_aromatic = ring_spec.get("aromatic", False)
        canvas.state.active_ring_template = ring_spec.get("template")
        canvas.state.active_ring_anomeric = ring_spec.get("anomeric")

        canvas.state.default_element = self.toolbar.current_element()
        canvas.state.active_orbital_kind = self.symbols_toolbar.current_orbital_kind()
        canvas.set_current_tool(self._current_tool_id)

    def _connect_undo_redo(self) -> None:
        """Conecta acciones globales a la pestaña activa."""
        self.action_undo.triggered.connect(self._on_undo)
        self.action_redo.triggered.connect(self._on_redo)
        self.action_copy.triggered.connect(self._on_copy)
        self.action_paste.triggered.connect(self._on_paste)
        self.action_undo.setEnabled(False)
        self.action_redo.setEnabled(False)

    def _load_recent_files(self) -> list[str]:
        """Método auxiliar para  load recent files.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        value = self._settings.value("recent_files", [])
        if isinstance(value, str):
            return [value]
        if isinstance(value, list):
            return [str(p) for p in value if p]
        return []

    def _save_recent_files(self) -> None:
        """Método auxiliar para  save recent files.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self._settings.setValue("recent_files", self._recent_files)

    def _load_update_preferences(self) -> UpdateSettings:
        """Carga preferencias de actualización persistidas en QSettings."""
        enabled = self._setting_bool(self._settings.value("update/enabled", False), False)
        raw_channel = str(self._settings.value("update/channel", "stable") or "stable").strip().lower()
        raw_mode = str(self._settings.value("update/mode", "notify") or "notify").strip().lower()
        try:
            interval_hours = int(self._settings.value("update/check_interval_hours", 24) or 24)
        except Exception:
            interval_hours = 24
        require_sha256 = self._setting_bool(self._settings.value("update/require_sha256", True), True)
        require_signature = self._setting_bool(
            self._settings.value("update/require_signature", False),
            False,
        )
        highest_seen_version = str(
            self._settings.value("update/highest_seen_version", "") or ""
        ).strip()
        channel = UpdateChannel.STABLE
        mode = UpdateMode.NOTIFY
        try:
            channel = UpdateChannel(raw_channel)
        except Exception:
            channel = UpdateChannel.STABLE
        try:
            mode = UpdateMode(raw_mode)
        except Exception:
            mode = UpdateMode.NOTIFY
        last_check_iso = str(self._settings.value("update/last_check_iso", "") or "").strip()
        return UpdateSettings(
            enabled=bool(enabled),
            channel=channel,
            mode=mode,
            check_interval_hours=max(1, interval_hours),
            last_check_iso=last_check_iso,
            highest_seen_version=highest_seen_version,
            require_sha256=bool(require_sha256),
            require_signature=bool(require_signature),
        )

    def _save_update_preferences(self) -> None:
        """Persiste preferencias de actualización en QSettings."""
        self._settings.setValue("update/enabled", bool(self._update_settings.enabled))
        self._settings.setValue("update/channel", str(self._update_settings.channel.value))
        self._settings.setValue("update/mode", str(self._update_settings.mode.value))
        self._settings.setValue(
            "update/check_interval_hours",
            int(max(1, self._update_settings.check_interval_hours)),
        )
        self._settings.setValue("update/last_check_iso", str(self._update_settings.last_check_iso or ""))
        self._settings.setValue(
            "update/highest_seen_version",
            str(self._update_settings.highest_seen_version or ""),
        )
        self._settings.setValue("update/require_sha256", bool(self._update_settings.require_sha256))
        self._settings.setValue(
            "update/require_signature",
            bool(self._update_settings.require_signature),
        )

    def _update_settings_payload(self) -> dict:
        """Devuelve preferencias de update para precargar en diálogo."""
        return {
            "enabled": bool(self._update_settings.enabled),
            "channel": str(self._update_settings.channel.value),
            "mode": str(self._update_settings.mode.value),
            "check_interval_hours": int(max(1, self._update_settings.check_interval_hours)),
        }

    @staticmethod
    def _setting_bool(value, default: bool) -> bool:
        """Normaliza valores de QSettings a booleano."""
        if value is None:
            return bool(default)
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            lowered = value.strip().lower()
            if lowered in {"1", "true", "yes", "on", "si", "sí"}:
                return True
            if lowered in {"0", "false", "no", "off"}:
                return False
        try:
            return bool(int(value))
        except Exception:
            return bool(value)

    def _load_numbering_preferences(self) -> None:
        """Carga preferencias globales de numeración del usuario."""
        mode = str(self._settings.value("numbering/mode", "atoms") or "atoms").strip().lower()
        if mode not in {"atoms", "structures", "both"}:
            mode = "atoms"
        include_export = self._setting_bool(
            self._settings.value("numbering/include_export", True),
            True,
        )
        self._numbering_default_mode = mode
        self._numbering_default_include_export = bool(include_export)
        self._apply_default_numbering_to_canvas(self.canvas)

    def _load_naming_preferences(self) -> tuple[bool, bool]:
        """Carga preferencias globales de nomenclatura avanzada."""
        advanced = self._setting_bool(
            self._settings.value("naming/advanced_enabled", True),
            True,
        )
        isolated = self._setting_bool(
            self._settings.value("naming/rdkit_isolated", True),
            True,
        )
        return bool(advanced), bool(isolated)

    def _save_naming_preferences(self) -> None:
        """Persiste preferencias globales de nomenclatura avanzada."""
        self._settings.setValue("naming/advanced_enabled", bool(self._name_advanced_default))
        self._settings.setValue("naming/rdkit_isolated", bool(self._name_rdkit_isolated_default))

    def _naming_settings_payload(self) -> dict:
        """Payload para precargar preferencias de nomenclatura en diálogo."""
        return {
            "advanced_enabled": bool(self.canvas.name_advanced_enabled),
            "rdkit_isolated": bool(self.canvas.name_rdkit_isolated),
        }

    def _save_numbering_preferences(self) -> None:
        """Guarda preferencias globales de numeración del usuario."""
        # Se guarda modo/exports, pero no el estado "mostrar numeración" para
        # evitar que inicie activa en futuras aperturas.
        self._settings.remove("numbering/enabled")
        self._numbering_default_mode = str(self.canvas.state.numbering_mode)
        self._numbering_default_include_export = bool(self.canvas.state.numbering_include_in_export)
        self._settings.setValue("numbering/mode", str(self._numbering_default_mode))
        self._settings.setValue(
            "numbering/include_export",
            bool(self._numbering_default_include_export),
        )

    def _sync_numbering_actions(self) -> None:
        """Sincroniza acciones de menú de numeración con el estado del canvas."""
        mode = str(self.canvas.state.numbering_mode or "").strip().lower()
        if mode not in {"atoms", "structures", "both"}:
            mode = "atoms"
        self.canvas.state.numbering_mode = mode
        state_to_action = {
            "atoms": self.action_numbering_mode_atoms,
            "structures": self.action_numbering_mode_structures,
            "both": self.action_numbering_mode_both,
        }

        actions = [
            self.action_numbering_enabled,
            self.action_numbering_mode_atoms,
            self.action_numbering_mode_structures,
            self.action_numbering_mode_both,
            self.action_numbering_export,
        ]
        previous_blocks = []
        for action in actions:
            previous_blocks.append((action, action.blockSignals(True)))
        try:
            self.action_numbering_enabled.setChecked(bool(self.canvas.state.numbering_enabled))
            for key, action in state_to_action.items():
                action.setChecked(key == mode)
            self.action_numbering_export.setChecked(
                bool(self.canvas.state.numbering_include_in_export)
            )
            self.action_numbering_recalculate.setEnabled(bool(self.canvas.state.numbering_enabled))
        finally:
            for action, was_blocked in previous_blocks:
                action.blockSignals(was_blocked)

    def _maybe_check_updates_startup(self) -> None:
        """Chequea updates en inicio respetando política y frecuencia."""
        self._check_for_updates(force=False, interactive=False)

    def _on_check_updates_now(self) -> None:
        """Lanza chequeo manual de actualizaciones desde el menú Ayuda."""
        self._check_for_updates(force=True, interactive=True)

    def _preferred_update_asset_flavor(self) -> str | None:
        """Devuelve preferencia de asset según contexto de ejecución."""
        self._windows_install_context = detect_windows_install_context()
        self._portable_update_context = detect_portable_update_context(
            windows_context=self._windows_install_context
        )
        flavor = choose_windows_asset_flavor(self._windows_install_context)
        return flavor or None

    def _current_portable_update_context(self) -> PortableUpdateContext:
        """Recalcula contexto de auto-update para AppImage/binario portable."""
        self._windows_install_context = detect_windows_install_context()
        self._portable_update_context = detect_portable_update_context(
            windows_context=self._windows_install_context
        )
        return self._portable_update_context

    def _log_update_event(self, event: str, **fields) -> None:
        """Registra telemetría local mínima de update (sin datos sensibles)."""
        try:
            self._update_telemetry.log_event(event, **fields)
        except Exception:
            return

    def _check_for_updates(self, force: bool, interactive: bool) -> None:
        """Ejecuta chequeo de updates aplicando política y canal."""
        self._log_update_event(
            "check_start",
            force=bool(force),
            interactive=bool(interactive),
            channel=self._update_settings.channel.value,
            mode=self._update_settings.mode.value,
        )
        if os.getenv("CHEMUSON_DISABLE_UPDATE_CHECK", "").strip().lower() in {"1", "true", "yes"}:
            flatpak_runtime = is_running_in_flatpak()
            self._log_update_event(
                "check_skipped",
                reason_code="disabled_flatpak" if flatpak_runtime else "disabled_env",
                channel=self._update_settings.channel.value,
                mode=self._update_settings.mode.value,
            )
            if interactive:
                QMessageBox.information(
                    self,
                    "Actualizaciones",
                    format_update_disabled_message(flatpak=flatpak_runtime),
                )
            return
        if not force and not should_check_now(self._update_settings):
            self._log_update_event(
                "check_skipped",
                reason_code="interval_not_elapsed",
                channel=self._update_settings.channel.value,
                mode=self._update_settings.mode.value,
            )
            return
        mark_checked(self._update_settings)
        self._save_update_preferences()
        provider = None
        try:
            manual_check = bool(force and interactive)
            provider = GitHubReleasesProvider(
                "PJGV333",
                "Chemuson",
                timeout=8.0 if manual_check else 4.0,
                allow_cached_fallback=not manual_check,
            )
            updater = AutoUpdateCore(provider, self._update_settings)
            result = updater.check_for_updates(
                get_app_version(),
                preferred_asset_flavor=self._preferred_update_asset_flavor(),
            )
        except Exception as exc:
            self._log_update_event(
                "check_error",
                error_type=exc.__class__.__name__,
                channel=self._update_settings.channel.value,
                mode=self._update_settings.mode.value,
            )
            if self._update_settings.mode == UpdateMode.NOTIFY or interactive:
                self.statusBar().showMessage(f"No se pudo comprobar actualización: {exc}", 10000)
                QMessageBox.warning(
                    self,
                    "Actualizaciones",
                    f"No se pudo comprobar actualización:\n{exc}",
                )
            return
        source = ""
        if provider is not None:
            source = str(getattr(provider, "last_fetch_source", "") or "")
        if source == "cache":
            self.statusBar().showMessage(
                "GitHub no disponible: se usó caché local de actualizaciones.",
                12000,
            )
        self._log_update_event(
            "check_result",
            available=bool(getattr(result, "available", False)),
            current_version=str(getattr(result, "current_version", "") or ""),
            latest_version=str(getattr(result, "latest_version", "") or ""),
            source=source or "unknown",
            channel=self._update_settings.channel.value,
            mode=self._update_settings.mode.value,
        )
        # Persiste metadatos de seguridad (p. ej. highest_seen_version anti-downgrade).
        self._save_update_preferences()
        self._handle_update_check_result(result, interactive=interactive, source=source)

    def _download_update_candidate(self, candidate):
        """Descarga y verifica un candidato de update; devuelve `DownloadedUpdate`."""
        version = str(getattr(candidate.release, "version", "latest") or "latest")
        self._log_update_event(
            "download_start",
            latest_version=version,
            channel=self._update_settings.channel.value,
            mode=self._update_settings.mode.value,
        )
        provider = GitHubReleasesProvider("PJGV333", "Chemuson", timeout=12.0)
        updater = AutoUpdateCore(provider, self._update_settings)
        download_dir = os.path.join(os.path.expanduser("~"), ".chemuson", "updates", version)
        try:
            downloaded = updater.download_candidate(candidate, download_dir)
            verification = updater.verify_download(downloaded)
            if not verification.ok:
                reason = verification.reason or "Falló verificación del paquete descargado."
                self._log_update_event(
                    "download_failed",
                    latest_version=version,
                    reason_code="verification_failed",
                )
                raise RuntimeError(reason)
            self._log_update_event("download_ok", latest_version=version)
            return downloaded
        except Exception as exc:
            self._log_update_event(
                "download_failed",
                latest_version=version,
                error_type=exc.__class__.__name__,
            )
            raise

    def _queue_windows_installer_update(self, candidate, show_errors: bool = True) -> bool:
        """Descarga instalador Windows y lo deja en cola para aplicar al salir."""
        candidate_version = str(getattr(candidate.release, "version", "") or "")
        if (
            self._pending_windows_installer_path
            and self._pending_windows_installer_version == candidate_version
            and os.path.exists(self._pending_windows_installer_path)
            and self._pending_windows_download is not None
        ):
            self._log_update_event(
                "queue_reused",
                latest_version=candidate_version,
                context="windows_installer",
            )
            return True
        try:
            downloaded = self._download_update_candidate(candidate)
        except Exception as exc:
            self._log_update_event(
                "queue_failed",
                latest_version=candidate_version,
                context="windows_installer",
                error_type=exc.__class__.__name__,
            )
            if show_errors:
                QMessageBox.warning(
                    self,
                    "Actualización",
                    f"No se pudo preparar el instalador de actualización:\n{exc}",
                )
            else:
                self.statusBar().showMessage(
                    f"No se pudo preparar la instalación silenciosa de actualización: {exc}",
                    12000,
                )
            return False
        self._pending_windows_download = downloaded
        self._pending_windows_installer_path = str(getattr(downloaded, "artifact_path", "") or "")
        self._pending_windows_installer_version = candidate_version
        self._log_update_event(
            "queue_ok",
            latest_version=candidate_version,
            context="windows_installer",
        )
        return True

    def _queue_portable_binary_update(
        self,
        candidate,
        context: PortableUpdateContext,
        show_errors: bool = True,
    ) -> bool:
        """Descarga un binario portable y lo deja listo para reemplazo al salir."""
        candidate_version = str(getattr(candidate.release, "version", "") or "")
        target_path = str(getattr(context, "target_path", "") or "").strip()
        if (
            self._pending_portable_target_path
            and self._pending_portable_target_path == target_path
            and self._pending_portable_version == candidate_version
            and self._pending_portable_download is not None
            and os.path.exists(str(getattr(self._pending_portable_download, "artifact_path", "") or ""))
        ):
            self._log_update_event(
                "queue_reused",
                latest_version=candidate_version,
                context="portable_binary",
            )
            return True
        try:
            if not context.can_self_update:
                raise RuntimeError("No se pudo determinar la ruta del ejecutable actual.")
            if not is_portable_target_writable(target_path):
                raise PermissionError(
                    "Chemuson no tiene permisos para reemplazar el ejecutable actual."
                )
            downloaded = self._download_update_candidate(candidate)
        except Exception as exc:
            self._log_update_event(
                "queue_failed",
                latest_version=candidate_version,
                context="portable_binary",
                error_type=exc.__class__.__name__,
            )
            if show_errors:
                QMessageBox.warning(
                    self,
                    "Actualización",
                    f"No se pudo preparar la actualización portable:\n{exc}",
                )
            else:
                self.statusBar().showMessage(
                    f"No se pudo preparar el auto-update real: {exc}",
                    12000,
                )
            return False
        self._pending_portable_download = downloaded
        self._pending_portable_target_path = target_path
        self._pending_portable_version = candidate_version
        self._pending_portable_relaunch = False
        self._pending_portable_context = context
        self._log_update_event(
            "queue_ok",
            latest_version=candidate_version,
            context="portable_binary",
            target_kind="appimage" if context.is_appimage else "portable",
        )
        return True

    def _offer_windows_installer_update(self, result, interactive: bool) -> None:
        """Gestiona oferta/cola de update para instalaciones Windows."""
        candidate = getattr(result, "candidate", None)
        if candidate is None:
            return
        version = str(getattr(result, "latest_version", "") or "")
        if self._update_settings.mode == UpdateMode.SILENT and not interactive:
            if self._queue_windows_installer_update(candidate, show_errors=False):
                self.statusBar().showMessage(
                    f"Instalación silenciosa {version} lista para aplicarse al cerrar Chemuson.",
                    20000,
                )
            return

        reply = QMessageBox.question(
            self,
            "Actualización disponible",
            (
                f"Hay una nueva versión disponible ({version}).\n\n"
                "Esta edición no usa auto-update real.\n"
                "Chemuson descargará el instalador oficial y lo ejecutará en silencio al cerrar.\n\n"
                "¿Quieres prepararlo ahora?"
            ),
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        if not self._queue_windows_installer_update(candidate):
            return
        self.statusBar().showMessage(
            f"Instalador {version} preparado. Se ejecutará en silencio al cerrar Chemuson.",
            20000,
        )
        if interactive:
            close_now = QMessageBox.question(
                self,
                "Aplicar actualización",
                "El instalador de actualización está listo.\n"
                "¿Deseas cerrar Chemuson ahora para ejecutar la instalación silenciosa?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if close_now == QMessageBox.StandardButton.Yes:
                self.close()

    def _offer_portable_binary_update(
        self,
        result,
        interactive: bool,
        context: PortableUpdateContext,
    ) -> None:
        """Gestiona oferta/cola de update para AppImage o binario portable."""
        candidate = getattr(result, "candidate", None)
        if candidate is None:
            return
        version = str(getattr(result, "latest_version", "") or "")
        target_label = "AppImage" if context.is_appimage else "ejecutable portable"
        if self._update_settings.mode == UpdateMode.SILENT and not interactive:
            if self._queue_portable_binary_update(candidate, context, show_errors=False):
                self.statusBar().showMessage(
                    f"Auto-update real {version} listo: Chemuson reemplazará el {target_label} al cerrar.",
                    20000,
                )
            return

        reply = QMessageBox.question(
            self,
            "Actualización disponible",
            (
                f"Hay una nueva versión disponible ({version}).\n\n"
                "Esto es auto-update real.\n"
                f"Chemuson descargará la actualización y reemplazará el {target_label} actual al cerrar.\n\n"
                "¿Quieres prepararla ahora?"
            ),
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        if not self._queue_portable_binary_update(candidate, context):
            return
        self.statusBar().showMessage(
            f"Auto-update real {version} preparado. Se aplicará al cerrar Chemuson.",
            20000,
        )
        if interactive:
            close_now = QMessageBox.question(
                self,
                "Aplicar actualización",
                "El auto-update real está listo.\n"
                "¿Deseas cerrar Chemuson ahora para completar el reemplazo?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if close_now == QMessageBox.StandardButton.Yes:
                self._pending_portable_relaunch = True
                if not self.close():
                    self._pending_portable_relaunch = False

    def _handle_update_check_result(self, result, interactive: bool = False, source: str = "") -> None:
        """Notifica updates y habilita flujo de instalador en Windows instalado."""
        if not getattr(result, "available", False):
            if interactive:
                QMessageBox.information(
                    self,
                    "Actualizaciones",
                    format_no_update_message(
                        str(getattr(result, "current_version", "") or self._app_version),
                        str(getattr(result, "channel", "") or self._update_settings.channel.value),
                        str(getattr(result, "reason", "") or ""),
                        str(getattr(result, "latest_version", "") or ""),
                        source,
                    ),
                )
            return
        candidate = getattr(result, "candidate", None)
        artifact_name = ""
        if candidate is not None and getattr(candidate, "artifact", None) is not None:
            artifact_name = str(candidate.artifact.name)

        # En instalaciones Windows se prioriza setup y se aplica al cierre.
        if (
            self._windows_install_context.is_windows
            and self._windows_install_context.installed
            and artifact_name
            and is_windows_installer_asset(artifact_name)
        ):
            self._offer_windows_installer_update(result, interactive=interactive)
            return

        portable_context = self._current_portable_update_context()
        if candidate is not None and portable_context.can_self_update:
            self._offer_portable_binary_update(
                result,
                interactive=interactive,
                context=portable_context,
            )
            return

        message = (
            f"Actualización disponible: {result.latest_version}"
            + (f" ({artifact_name})" if artifact_name else "")
        )
        self.statusBar().showMessage(message, 15000)
        if self._update_settings.mode == UpdateMode.NOTIFY or interactive:
            release_url = ""
            if candidate is not None and getattr(candidate, "release", None) is not None:
                release_url = str(candidate.release.html_url or "")
            source_line = f"\n\nOrigen de datos: {'caché local' if source == 'cache' else 'GitHub'}"
            extra = f"\n\nRelease: {release_url}" if release_url else ""
            QMessageBox.information(
                self,
                "Actualización disponible",
                f"Hay una nueva versión disponible ({result.latest_version}).{extra}{source_line}",
            )

    def _apply_pending_windows_update_on_exit(self) -> bool:
        """Ejecuta instalador pendiente antes de salir cuando corresponde."""
        installer_path = str(self._pending_windows_installer_path or "").strip()
        if not installer_path:
            return True
        clear_pending = False
        try:
            downloaded = self._pending_windows_download
            if downloaded is None:
                raise RuntimeError("No hay metadata de verificación del instalador pendiente.")
            provider = GitHubReleasesProvider("PJGV333", "Chemuson", timeout=6.0)
            updater = AutoUpdateCore(provider, self._update_settings)
            verify_result = updater.verify_download(downloaded)
            if not verify_result.ok:
                raise RuntimeError(
                    f"No se ejecutará artefacto no verificado: {verify_result.reason}"
                )
            launch_inno_installer(installer_path, silent=True)
            self._log_update_event(
                "apply_queued",
                latest_version=self._pending_windows_installer_version or "",
                context="windows_installer",
                result="launched",
            )
            clear_pending = True
            return True
        except Exception as exc:
            self._log_update_event(
                "apply_failed",
                latest_version=self._pending_windows_installer_version or "",
                context="windows_installer",
                error_type=exc.__class__.__name__,
            )
            reply = QMessageBox.question(
                self,
                "Actualización pendiente",
                (
                    "No se pudo lanzar el instalador de actualización pendiente.\n"
                    f"Error: {exc}\n\n"
                    "¿Deseas salir de todas formas sin aplicar la actualización?"
                ),
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if reply == QMessageBox.StandardButton.Yes:
                clear_pending = True
                return True
            return False
        finally:
            if clear_pending:
                self._pending_windows_download = None
                self._pending_windows_installer_path = ""
                self._pending_windows_installer_version = ""

    def _apply_pending_portable_update_on_exit(self) -> bool:
        """Lanza helper para reemplazar AppImage/binario portable al cerrar."""
        target_path = str(self._pending_portable_target_path or "").strip()
        if not target_path:
            return True
        clear_pending = False
        try:
            downloaded = self._pending_portable_download
            if downloaded is None:
                raise RuntimeError("No hay metadata de verificación del binario pendiente.")
            provider = GitHubReleasesProvider("PJGV333", "Chemuson", timeout=6.0)
            updater = AutoUpdateCore(provider, self._update_settings)
            verify_result = updater.verify_download(downloaded)
            if not verify_result.ok:
                raise RuntimeError(
                    f"No se reemplazará binario no verificado: {verify_result.reason}"
                )
            version = str(self._pending_portable_version or "latest").strip() or "latest"
            script_dir = os.path.dirname(str(getattr(downloaded, "artifact_path", "") or ""))
            rollback_dir = os.path.join(os.path.expanduser("~"), ".chemuson", "rollback")
            os.makedirs(rollback_dir, exist_ok=True)
            backup_name = (
                f"{os.path.basename(target_path)}."
                f"{datetime.utcnow().strftime('%Y%m%d_%H%M%S')}.bak"
            )
            backup_path = os.path.join(rollback_dir, backup_name)
            log_path = os.path.join(script_dir, "portable-update.log")
            script_name = (
                "apply-portable-update.cmd"
                if self._pending_portable_context.is_windows
                else "apply-portable-update.sh"
            )
            script_path = os.path.join(script_dir, script_name)
            write_portable_update_script(
                script_path,
                source_path=str(getattr(downloaded, "artifact_path", "") or ""),
                target_path=target_path,
                backup_path=backup_path,
                log_path=log_path,
                pid=os.getpid(),
                relaunch=bool(self._pending_portable_relaunch),
                platform_name="windows" if self._pending_portable_context.is_windows else "linux",
            )
            launch_portable_update_script(
                script_path,
                platform_name="windows" if self._pending_portable_context.is_windows else "linux",
            )
            self._log_update_event(
                "apply_queued",
                latest_version=version,
                context="portable_binary",
                result="helper_launched",
                target_kind="appimage" if self._pending_portable_context.is_appimage else "portable",
            )
            clear_pending = True
            return True
        except Exception as exc:
            self._log_update_event(
                "apply_failed",
                latest_version=self._pending_portable_version or "",
                context="portable_binary",
                error_type=exc.__class__.__name__,
            )
            reply = QMessageBox.question(
                self,
                "Actualización pendiente",
                (
                    "No se pudo preparar el reemplazo del ejecutable pendiente.\n"
                    f"Error: {exc}\n\n"
                    "¿Deseas salir de todas formas sin aplicar la actualización?"
                ),
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if reply == QMessageBox.StandardButton.Yes:
                clear_pending = True
                return True
            return False
        finally:
            if clear_pending:
                self._pending_portable_download = None
                self._pending_portable_target_path = ""
                self._pending_portable_version = ""
                self._pending_portable_relaunch = False
                self._pending_portable_context = PortableUpdateContext(
                    is_portable=False,
                    is_windows=False,
                    is_appimage=False,
                    target_path="",
                    executable_path="",
                    display_name="",
                )

    def _add_recent_file(self, filepath: str) -> None:
        """Método auxiliar para  add recent file.

        Args:
            filepath: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        if not filepath:
            return
        path = os.path.abspath(filepath)
        self._recent_files = [p for p in self._recent_files if p != path]
        self._recent_files.insert(0, path)
        self._recent_files = self._recent_files[:10]
        self._save_recent_files()
        self._update_recent_menu()

    def _update_recent_menu(self) -> None:
        """Método auxiliar para  update recent menu.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        if not hasattr(self, "recent_menu"):
            return
        self.recent_menu.clear()
        existing = [p for p in self._recent_files if os.path.exists(p)]
        self._recent_files = existing
        if not existing:
            empty = QAction("Sin recientes", self)
            empty.setEnabled(False)
            self.recent_menu.addAction(empty)
            return
        for path in existing:
            action = QAction(path, self)
            action.triggered.connect(lambda checked=False, p=path: self._open_recent_file(p))
            self.recent_menu.addAction(action)

    def _open_recent_file(self, filepath: str) -> None:
        """Método auxiliar para  open recent file.

        Args:
            filepath: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        if not filepath or not os.path.exists(filepath):
            QMessageBox.warning(self, "Archivo no encontrado", "El archivo no existe.")
            self._update_recent_menu()
            return
        self._open_file_path(filepath)
    
    # -------------------------------------------------------------------------
    # File Menu Handlers
    # -------------------------------------------------------------------------
    def _on_file_new(self) -> None:
        """Create a new empty canvas."""
        canvas = self._create_document_tab(make_current=True)
        self._apply_toolbar_defaults_to_canvas(canvas)
        self._set_active_canvas(canvas)
        self._update_total_charge_indicator()
        self.statusBar().showMessage("Nuevo documento creado")
    
    def _on_file_open(self) -> None:
        """Open a molecule file (.cmsn or .mol)."""
        filepaths, _selected_filter = QFileDialog.getOpenFileNames(
            self,
            "Abrir archivo(s)",
            "",
            "Archivos de Chemuson (*.cmsn);;Archivos MOL (*.mol *.sdf);;Todos los archivos (*.*)"
        )
        for filepath in filepaths:
            self._open_file_path(filepath)

    def _load_file_into_canvas(self, filepath: str, canvas: ChemusonCanvas) -> None:
        """Carga un archivo químico dentro de un canvas específico."""
        if filepath.lower().endswith(".cmsn"):
            canvas.clear_canvas()
            PersistenceManager.load_from_file(filepath, canvas)
            return

        with open(filepath, "r", encoding="utf-8") as f:
            molfile = f.read()
        graph = molfile_to_molgraph(molfile)
        canvas.clear_canvas()
        canvas._insert_molgraph(graph)

    @staticmethod
    def _read_autosave_metadata(filepath: str) -> Optional[dict]:
        """Lee metadatos mínimos de un autosave para la tabla de recuperación."""
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                payload = json.load(f)
        except Exception:
            return None
        if payload.get("application") != "Chemuson":
            return None
        metadata = payload.get("autosave_metadata")
        if not isinstance(metadata, dict):
            metadata = {}
        raw_path = metadata.get("original_path")
        raw_timestamp = metadata.get("timestamp")
        return {
            "autosave_path": filepath,
            "original_path": str(raw_path) if raw_path else None,
            "timestamp": str(raw_timestamp) if raw_timestamp else "Desconocida",
        }

    @staticmethod
    def _list_autosave_entries(directory: str) -> list[dict]:
        """Lista autosaves válidos de un directorio ordenados por fecha desc."""
        if not os.path.isdir(directory):
            return []
        entries: list[dict] = []
        for name in sorted(os.listdir(directory)):
            if not name.endswith(".json"):
                continue
            filepath = os.path.join(directory, name)
            if not os.path.isfile(filepath):
                continue
            metadata = ChemusonWindow._read_autosave_metadata(filepath)
            if metadata is None:
                continue
            metadata["filename"] = name
            entries.append(metadata)
        entries.sort(key=lambda entry: os.path.getmtime(entry["autosave_path"]), reverse=True)
        return entries

    @staticmethod
    def _archive_autosave(path: str, autosave_dir: str) -> str:
        """Mueve un autosave a la carpeta `old` y devuelve su nueva ruta."""
        old_dir = os.path.join(autosave_dir, "old")
        os.makedirs(old_dir, exist_ok=True)
        basename = os.path.basename(path)
        target = os.path.join(old_dir, basename)
        if os.path.exists(target):
            root, ext = os.path.splitext(basename)
            target = os.path.join(
                old_dir,
                f"{root}_{datetime.now().strftime('%Y%m%d_%H%M%S_%f')}{ext}",
            )
        os.replace(path, target)
        return target

    def _open_autosave_document(
        self,
        autosave_path: str,
        original_path: Optional[str] = None,
    ) -> bool:
        """Abre un autosave como pestaña nueva."""
        canvas = self._create_document_tab(make_current=True)
        self._apply_toolbar_defaults_to_canvas(canvas)
        try:
            PersistenceManager.load_from_file(autosave_path, canvas)
            canvas.undo_stack.setClean()
            self._set_canvas_file_path(canvas, original_path)
            self.tabs.setCurrentWidget(canvas)
            self._set_active_canvas(canvas)
            self._update_total_charge_indicator()
            return True
        except Exception as exc:
            index = self.tabs.indexOf(canvas)
            if index >= 0:
                self.tabs.removeTab(index)
            autosave_manager = self._canvas_autosave_managers.pop(canvas, None)
            if autosave_manager is not None:
                autosave_manager.stop()
            self._canvas_file_paths.pop(canvas, None)
            self._canvas_tab_titles.pop(canvas, None)
            canvas.deleteLater()
            if self.tabs.count() == 0:
                replacement = self._create_document_tab(make_current=True)
                self._apply_toolbar_defaults_to_canvas(replacement)
                self._set_active_canvas(replacement)
            QMessageBox.warning(
                self,
                "Error al abrir autosave",
                f"No se pudo abrir el autosave:\n{exc}",
            )
            return False

    def _on_open_recovery_center(self) -> None:
        """Abre manualmente el centro de recuperación y archivos recientes."""
        self._show_recovery_center(show_only_if_pending=False)

    def _show_recovery_center(self, show_only_if_pending: bool = False) -> None:
        """Muestra diálogo con pendientes, recientes y recuperados."""
        autosave_dir = AutosaveManager.default_autosave_dir()
        pending_entries = self._list_autosave_entries(autosave_dir)
        recovered_entries = self._list_autosave_entries(os.path.join(autosave_dir, "old"))
        recent_entries = [p for p in self._recent_files if os.path.exists(p)]
        if recent_entries != self._recent_files:
            self._recent_files = recent_entries
            self._save_recent_files()
            self._update_recent_menu()

        if show_only_if_pending and not pending_entries:
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Centro de recuperación")
        dialog.resize(980, 560)
        layout = QVBoxLayout(dialog)
        layout.addWidget(
            QLabel(
                "Gestiona autosaves pendientes, sesiones recuperadas "
                "y tus archivos recientes."
            )
        )

        tabs = QTabWidget(dialog)
        layout.addWidget(tabs)

        def _setup_table(table: QTableWidget, action_column: int, action_width: int) -> None:
            """Configura tabla con filas legibles y columna de acciones estable."""
            table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
            table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
            table.verticalHeader().setVisible(False)
            table.verticalHeader().setDefaultSectionSize(48)
            table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
            if table.columnCount() > 2:
                table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
            table.horizontalHeader().setSectionResizeMode(action_column, QHeaderView.ResizeMode.Fixed)
            table.setColumnWidth(action_column, action_width)

        def _set_actions_cell(
            table: QTableWidget,
            row: int,
            buttons: list[tuple[str, Callable]],
        ) -> None:
            """Inserta botones de acción en una celda respetando márgenes."""
            action_widget = QWidget(table)
            action_layout = QHBoxLayout(action_widget)
            action_layout.setContentsMargins(8, 6, 8, 6)
            action_layout.setSpacing(8)
            button_widgets: list[QPushButton] = []
            for text, callback in buttons:
                button = QPushButton(text, action_widget)
                # Garantiza legibilidad incluso con escalado alto del sistema.
                button.setMinimumHeight(38)
                button.setMinimumWidth(136)
                button.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
                button.clicked.connect(callback)
                action_layout.addWidget(button)
                button_widgets.append(button)
            table.setCellWidget(row, table.columnCount() - 1, action_widget)
            margins = action_layout.contentsMargins()
            required_width = margins.left() + margins.right()
            if button_widgets:
                required_width += action_layout.spacing() * (len(button_widgets) - 1)
            for button in button_widgets:
                required_width += max(button.minimumSizeHint().width(), button.sizeHint().width())
            action_column = table.columnCount() - 1
            if required_width > table.columnWidth(action_column):
                table.setColumnWidth(action_column, required_width)
            required_height = margins.top() + margins.bottom()
            required_height += max(button.minimumSizeHint().height() for button in button_widgets)
            table.setRowHeight(
                row,
                max(
                    table.rowHeight(row),
                    required_height + 2,
                ),
            )

        pending_tab = QWidget()
        pending_layout = QVBoxLayout(pending_tab)
        if not pending_entries:
            pending_layout.addWidget(QLabel("No hay autosaves pendientes."))
        else:
            pending_table = QTableWidget(len(pending_entries), 3, pending_tab)
            pending_table.setHorizontalHeaderLabels(["Archivo original", "Timestamp", "Acciones"])
            _setup_table(pending_table, action_column=2, action_width=360)

            def _mark_pending_done(row: int, status: str) -> None:
                pending_table.setCellWidget(row, 2, QLabel(status))
                for col in (0, 1):
                    item = pending_table.item(row, col)
                    if item is not None:
                        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEnabled)

            for row, entry in enumerate(pending_entries):
                original_path_value = entry.get("original_path")
                original_path = original_path_value or "Sin nombre"
                timestamp = entry.get("timestamp") or "Desconocida"
                autosave_path = entry["autosave_path"]

                pending_table.setItem(row, 0, QTableWidgetItem(original_path))
                pending_table.setItem(row, 1, QTableWidgetItem(timestamp))

                def _recover_handler(
                    _checked: bool = False,
                    r: int = row,
                    p: str = autosave_path,
                    original: Optional[str] = original_path_value,
                ) -> None:
                    if not self._open_autosave_document(p, original_path=original):
                        return
                    try:
                        self._archive_autosave(p, autosave_dir)
                        _mark_pending_done(r, "Recuperado")
                        self.statusBar().showMessage("Autosave recuperado")
                    except Exception as exc:
                        QMessageBox.warning(
                            self,
                            "Error al mover autosave",
                            f"No se pudo archivar el autosave:\n{exc}",
                        )

                def _discard_handler(
                    _checked: bool = False,
                    r: int = row,
                    p: str = autosave_path,
                ) -> None:
                    try:
                        self._archive_autosave(p, autosave_dir)
                        _mark_pending_done(r, "Descartado")
                    except Exception as exc:
                        QMessageBox.warning(
                            self,
                            "Error al descartar",
                            f"No se pudo mover el autosave:\n{exc}",
                        )

                _set_actions_cell(
                    pending_table,
                    row,
                    [
                        ("Recuperar", _recover_handler),
                        ("Descartar", _discard_handler),
                    ],
                )

            pending_layout.addWidget(pending_table)
        tabs.addTab(pending_tab, f"Pendientes ({len(pending_entries)})")

        recent_tab = QWidget()
        recent_layout = QVBoxLayout(recent_tab)
        if not recent_entries:
            recent_layout.addWidget(QLabel("No hay archivos recientes disponibles."))
        else:
            recent_table = QTableWidget(len(recent_entries), 2, recent_tab)
            recent_table.setHorizontalHeaderLabels(["Archivo reciente", "Acciones"])
            _setup_table(recent_table, action_column=1, action_width=170)
            for row, filepath in enumerate(recent_entries):
                recent_table.setItem(row, 0, QTableWidgetItem(filepath))

                def _open_recent_handler(
                    _checked: bool = False,
                    p: str = filepath,
                ) -> None:
                    if not os.path.exists(p):
                        QMessageBox.warning(self, "Archivo no encontrado", "El archivo no existe.")
                        self._update_recent_menu()
                        return
                    self._open_file_path(p)

                _set_actions_cell(recent_table, row, [("Abrir", _open_recent_handler)])
            recent_layout.addWidget(recent_table)
        tabs.addTab(recent_tab, f"Recientes ({len(recent_entries)})")

        recovered_tab = QWidget()
        recovered_layout = QVBoxLayout(recovered_tab)
        if not recovered_entries:
            recovered_layout.addWidget(QLabel("No hay autosaves recuperados/archivados."))
        else:
            recovered_table = QTableWidget(len(recovered_entries), 3, recovered_tab)
            recovered_table.setHorizontalHeaderLabels(["Archivo original", "Timestamp", "Acciones"])
            _setup_table(recovered_table, action_column=2, action_width=170)
            for row, entry in enumerate(recovered_entries):
                original_path_value = entry.get("original_path")
                original_path = original_path_value or "Sin nombre"
                timestamp = entry.get("timestamp") or "Desconocida"
                autosave_path = entry["autosave_path"]
                recovered_table.setItem(row, 0, QTableWidgetItem(original_path))
                recovered_table.setItem(row, 1, QTableWidgetItem(timestamp))

                def _open_recovered_handler(
                    _checked: bool = False,
                    p: str = autosave_path,
                    original: Optional[str] = original_path_value,
                ) -> None:
                    self._open_autosave_document(p, original_path=original)

                _set_actions_cell(recovered_table, row, [("Abrir", _open_recovered_handler)])
            recovered_layout.addWidget(recovered_table)
        tabs.addTab(recovered_tab, f"Recuperados ({len(recovered_entries)})")

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(dialog.reject)
        buttons.accepted.connect(dialog.accept)
        layout.addWidget(buttons)
        dialog.exec()

    @staticmethod
    def check_autosaves(window: "ChemusonWindow") -> None:
        """Busca autosaves pendientes y abre el centro de recuperación al iniciar."""
        window._show_recovery_center(show_only_if_pending=True)

    def _open_file_path(self, filepath: str) -> None:
        """Método auxiliar para  open file path.

        Args:
            filepath: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        if not filepath:
            return
        canvas = self._create_document_tab(make_current=True)
        self._apply_toolbar_defaults_to_canvas(canvas)
        try:
            self._load_file_into_canvas(filepath, canvas)
            canvas.undo_stack.setClean()
            self._set_canvas_file_path(canvas, filepath)
            self.tabs.setCurrentWidget(canvas)
            self._set_active_canvas(canvas)
            self._add_recent_file(filepath)
            self._update_total_charge_indicator()
            self.statusBar().showMessage(f"Abierto: {filepath}")
        except Exception as e:
            index = self.tabs.indexOf(canvas)
            if index >= 0:
                self.tabs.removeTab(index)
            autosave_manager = self._canvas_autosave_managers.pop(canvas, None)
            if autosave_manager is not None:
                autosave_manager.stop()
            self._canvas_file_paths.pop(canvas, None)
            self._canvas_tab_titles.pop(canvas, None)
            canvas.deleteLater()
            if self.tabs.count() == 0:
                replacement = self._create_document_tab(make_current=True)
                self._apply_toolbar_defaults_to_canvas(replacement)
                self._set_active_canvas(replacement)
            QMessageBox.critical(self, "Error", f"No se pudo abrir el archivo:\n{e}")
    
    def _on_file_save(self) -> None:
        """Save the current work in .cmsn format."""
        filepath = self._canvas_file_paths.get(self.canvas)
        selected_filter = ""
        if not filepath:
            filepath, selected_filter = QFileDialog.getSaveFileName(
                self,
                "Guardar archivo",
                "",
                "Archivo de Chemuson (*.cmsn);;Archivo MOL (*.mol);;Todos los archivos (*.*)"
            )
        if filepath:
            try:
                autosave_manager = self._canvas_autosave_managers.get(self.canvas)
                if filepath.lower().endswith(".mol") or filepath.lower().endswith(".sdf") or "MOL" in selected_filter:
                    # Export as .mol if explicitly requested
                    from chemuson.chemio.rdkit_io import molgraph_to_molfile
                    molfile = molgraph_to_molfile(self.canvas.graph)
                    with open(filepath, "w") as f:
                        f.write(molfile)
                else:
                    # Save as .cmsn (native)
                    if not filepath.lower().endswith(".cmsn"):
                        filepath += ".cmsn"
                    PersistenceManager.save_to_file(filepath, self.canvas)
                self._set_canvas_file_path(self.canvas, filepath)
                self.canvas.undo_stack.setClean()
                if autosave_manager is not None:
                    autosave_manager.cleanup_after_save()
                self._add_recent_file(filepath)
                self.statusBar().showMessage(f"Guardado: {filepath}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"No se pudo guardar:\n{e}")
    
    def _on_export(self, format: str) -> None:
        """Export the canvas in the specified format."""
        if format == "png":
            image = self.canvas._render_scene_image()
            if image:
                filepath, _ = QFileDialog.getSaveFileName(self, "Exportar PNG", "", "Imagen PNG (*.png)")
                if filepath:
                    image.save(filepath, "PNG")
                    self.statusBar().showMessage(f"Exportado: {filepath}")
        elif format == "svg":
            svg_data = self.canvas._render_scene_svg()
            filepath, _ = QFileDialog.getSaveFileName(self, "Exportar SVG", "", "Imagen SVG (*.svg)")
            if filepath:
                if svg_data:
                    with open(filepath, "w") as f:
                        f.write(svg_data.decode("utf-8", errors="replace"))
                else:
                    from chemuson.chemio.rdkit_io import molgraph_to_svg
                    svg_text = molgraph_to_svg(self.canvas.graph)
                    with open(filepath, "w") as f:
                        f.write(svg_text)
                    if self.canvas.state.numbering_enabled and self.canvas.state.numbering_include_in_export:
                        QMessageBox.information(
                            self,
                            "SVG sin numeración",
                            "Este entorno no permite render SVG desde la escena.\n"
                            "Se exportó SVG químico base sin overlay de numeración.",
                        )
                self.statusBar().showMessage(f"Exportado: {filepath}")
        elif format == "pdf":
            filepath, _ = QFileDialog.getSaveFileName(self, "Exportar PDF", "", "Documento PDF (*.pdf)")
            if filepath:
                try:
                    printer = QPrinter(QPrinter.PrinterMode.HighResolution)
                    printer.setOutputFormat(QPrinter.OutputFormat.PdfFormat)
                    printer.setOutputFileName(filepath)

                    bounds = self.canvas._render_scene_bounds(selected_only=False)
                    if bounds is None:
                        QMessageBox.warning(self, "Error", "No hay contenido para exportar.")
                        return

                    def render_pdf() -> bool:
                        painter = QPainter(printer)
                        if not painter.isActive():
                            return False
                        self.canvas.scene.render(
                            painter,
                            printer.pageRect(QPrinter.Unit.Point),
                            bounds,
                        )
                        painter.end()
                        return True

                    ok = self.canvas._with_hidden_render_items(render_pdf)
                    if not ok:
                        QMessageBox.warning(self, "Error", "No se pudo exportar PDF.")
                        return
                    self.statusBar().showMessage(f"Exportado: {filepath}")
                except Exception as e:
                    QMessageBox.warning(self, "Error", f"No se pudo exportar PDF:\n{e}")

    def _confirm_discard_changes(self, canvas: Optional[ChemusonCanvas] = None) -> bool:
        """Método auxiliar para  confirm discard changes.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        target = canvas or self.canvas
        if target.undo_stack.isClean():
            return True
        index = self.tabs.indexOf(target)
        title = self._canvas_tab_titles.get(target, "Documento")
        if index >= 0:
            title = self.tabs.tabText(index).replace(" *", "")
        reply = QMessageBox.question(
            self,
            "Cambios sin guardar",
            f"'{title}' tiene cambios sin guardar. ¿Deseas guardar antes de cerrar?",
            QMessageBox.StandardButton.Save
            | QMessageBox.StandardButton.Discard
            | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Save,
        )
        if reply == QMessageBox.StandardButton.Save:
            previous_index = self.tabs.currentIndex()
            if index >= 0 and previous_index != index:
                self.tabs.setCurrentIndex(index)
                self._set_active_canvas(target)
            self._on_file_save()
            saved = target.undo_stack.isClean()
            if (
                previous_index >= 0
                and previous_index < self.tabs.count()
                and self.tabs.currentIndex() != previous_index
            ):
                self.tabs.setCurrentIndex(previous_index)
            return saved
        if reply == QMessageBox.StandardButton.Discard:
            return True
        return False

    def closeEvent(self, event) -> None:
        """Método auxiliar para closeEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        canvases = [
            canvas
            for i in range(self.tabs.count())
            if (canvas := self._canvas_from_tab_index(i)) is not None
        ]
        for canvas in canvases:
            if not self._confirm_discard_changes(canvas):
                event.ignore()
                return
        if not self._apply_pending_portable_update_on_exit():
            event.ignore()
            return
        if not self._apply_pending_windows_update_on_exit():
            event.ignore()
            return
        event.accept()

    def changeEvent(self, event) -> None:
        """Método auxiliar para changeEvent.

        Args:
            event: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        if event.type() == QEvent.Type.ActivationChange and self.isActiveWindow():
            self.canvas.restore_text_edit_focus()
        super().changeEvent(event)
    
    # -------------------------------------------------------------------------
    # Edit Menu Handlers
    # -------------------------------------------------------------------------
    def _on_copy_as(self, format: str) -> None:
        """Copy molecule in specified format to clipboard."""
        try:
            from PyQt6.QtWidgets import QApplication
            clipboard = QApplication.clipboard()
            graph, _bbox = self.canvas._analysis_graph_and_bbox()
            target_graph = graph if graph is not None else self.canvas.graph
            
            if format == "smiles":
                from chemuson.chemio.rdkit_io import molgraph_to_smiles
                text = molgraph_to_smiles(target_graph) if target_graph.atoms else ""
            elif format == "molfile":
                from chemuson.chemio.rdkit_io import molgraph_to_molfile
                text = molgraph_to_molfile(target_graph) if target_graph.atoms else ""
            elif format == "inchi":
                # InChI requires additional RDKit import
                from chemuson.chemio.rdkit_io import molgraph_to_rdkit
                from rdkit.Chem.inchi import MolToInchi
                if target_graph.atoms:
                    mol = molgraph_to_rdkit(target_graph)
                    text = MolToInchi(mol)
                else:
                    text = ""
            else:
                text = ""
            
            clipboard.setText(text)
            self.statusBar().showMessage(f"Copiado como {format.upper()}")
        except Exception as e:
            self.statusBar().showMessage(f"Error: {e}")
    
    def _on_preferences(self) -> None:
        """Open preferences dialog."""
        dialog = PreferencesDialog(
            self.canvas.state,
            self.canvas.drawing_style,
            update_settings=self._update_settings_payload(),
            naming_settings=self._naming_settings_payload(),
            parent=self,
        )
        dialog.preferences_changed.connect(self._apply_preferences)
        dialog.exec()

    def _on_style_dialog(self) -> None:
        """Open drawing style dialog."""
        dialog = StyleDialog(
            self.canvas.drawing_style,
            self.canvas.state.bond_length,
            self.canvas.state.label_font_size,
            self.canvas.state.numbering_font_size,
            self,
        )
        if dialog.exec() == QDialog.DialogCode.Accepted:
            result = dialog.selected_dimensions()
            self.canvas.apply_document_dimensions(
                style=result["style"],
                label_font_size=result["label_font_size"],
                numbering_font_size=result["numbering_font_size"],
                scale_existing=result["scale_existing"],
                scale_factor=result["scale_factor"],
            )

    def _apply_preferences(self, prefs: dict) -> None:
        """Aplica preferences.

        Args:
            prefs: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.state.show_implicit_carbons = prefs.get("show_carbons", False)
        self.canvas.state.show_implicit_hydrogens = prefs.get("show_hydrogens", False)
        self.canvas.state.use_aromatic_circles = prefs.get("aromatic_circles", False)
        bond_caps = prefs.get("bond_caps")
        if bond_caps:
            self._apply_bond_caps(bond_caps)

        self.action_show_carbons.setChecked(self.canvas.state.show_implicit_carbons)
        self.action_show_hydrogens.setChecked(self.canvas.state.show_implicit_hydrogens)
        self.action_aromatic_circles.setChecked(self.canvas.state.use_aromatic_circles)

        self.canvas.refresh_atom_visibility()
        self.canvas.refresh_aromatic_circles()
        self._sync_label_menu_state()

        if bond_caps:
            self.appearance_dock.set_bond_caps(bond_caps)

        update_enabled = self._setting_bool(prefs.get("update_enabled", self._update_settings.enabled), self._update_settings.enabled)
        update_channel_raw = str(prefs.get("update_channel", self._update_settings.channel.value) or self._update_settings.channel.value).strip().lower()
        update_mode_raw = str(prefs.get("update_mode", self._update_settings.mode.value) or self._update_settings.mode.value).strip().lower()
        try:
            update_interval = int(prefs.get("update_check_interval_hours", self._update_settings.check_interval_hours))
        except Exception:
            update_interval = self._update_settings.check_interval_hours
        try:
            channel = UpdateChannel(update_channel_raw)
        except Exception:
            channel = UpdateChannel.STABLE
        try:
            mode = UpdateMode(update_mode_raw)
        except Exception:
            mode = UpdateMode.NOTIFY
        self._update_settings.enabled = bool(update_enabled)
        self._update_settings.channel = channel
        self._update_settings.mode = mode
        self._update_settings.check_interval_hours = max(1, int(update_interval))
        self._save_update_preferences()

        advanced_enabled = self._setting_bool(
            prefs.get("name_advanced_enabled", self.canvas.name_advanced_enabled),
            self.canvas.name_advanced_enabled,
        )
        rdkit_isolated = self._setting_bool(
            prefs.get("name_rdkit_isolated", self.canvas.name_rdkit_isolated),
            self.canvas.name_rdkit_isolated,
        )
        self.canvas.name_advanced_enabled = bool(advanced_enabled)
        self.canvas.name_rdkit_isolated = bool(rdkit_isolated)
        self._name_advanced_default = bool(advanced_enabled)
        self._name_rdkit_isolated_default = bool(rdkit_isolated)
        self._save_naming_preferences()
        self._update_iupac_name_indicator()

    def _apply_appearance_settings(self, prefs: dict) -> None:
        """Aplica appearance settings.

        Args:
            prefs: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        bond_caps = prefs.get("bond_caps")
        if bond_caps:
            self._apply_bond_caps(bond_caps)

    def _apply_bond_caps(self, bond_caps: str) -> None:
        """Aplica bond caps.

        Args:
            bond_caps: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        if bond_caps == "round":
            cap_style = Qt.PenCapStyle.RoundCap
            join_style = Qt.PenJoinStyle.RoundJoin
        else:
            cap_style = Qt.PenCapStyle.FlatCap
            join_style = Qt.PenJoinStyle.MiterJoin
        if (
            self.canvas.drawing_style.cap_style != cap_style
            or self.canvas.drawing_style.join_style != join_style
        ):
            style = replace(
                self.canvas.drawing_style,
                cap_style=cap_style,
                join_style=join_style,
            )
            self.canvas.apply_drawing_style(style)
    
    # -------------------------------------------------------------------------
    # View Menu Handlers
    # -------------------------------------------------------------------------
    def _on_toggle_numbering(self, checked: bool) -> None:
        """Activa/desactiva numeración de átomos/estructuras."""
        self.canvas.set_numbering_enabled(checked)
        self._sync_numbering_actions()
        self._save_numbering_preferences()
        self.statusBar().showMessage(
            "Numeración: visible" if checked else "Numeración: oculta"
        )

    def _on_set_numbering_mode(self, mode: str) -> None:
        """Cambia el modo de numeración mostrado en el lienzo."""
        self.canvas.set_numbering_mode(mode)
        self._sync_numbering_actions()
        self._save_numbering_preferences()
        labels = {
            "atoms": "átomos",
            "structures": "estructuras",
            "both": "átomos y estructuras",
        }
        active_mode = self.canvas.state.numbering_mode
        self.statusBar().showMessage(f"Numeración: {labels.get(active_mode, active_mode)}")

    def _on_recalculate_numbering(self) -> None:
        """Recalcula explícitamente la numeración visual."""
        self.canvas.recompute_numbering(force_reset=True)
        self.statusBar().showMessage("Numeración recalculada")

    def _on_toggle_numbering_export(self, checked: bool) -> None:
        """Define si la numeración se incluye al exportar."""
        self.canvas.set_numbering_include_in_export(checked)
        self._sync_numbering_actions()
        self._save_numbering_preferences()
        self.statusBar().showMessage(
            "Exportación: numeración incluida"
            if checked
            else "Exportación: numeración excluida"
        )

    def _on_toggle_carbons(self, checked: bool) -> None:
        """Toggle visibility of implicit carbons."""
        self.canvas.state.show_implicit_carbons = checked
        self.canvas.refresh_atom_visibility()
        self.statusBar().showMessage(
            "Carbonos: visibles" if checked else "Carbonos: ocultos"
        )
    
    def _on_toggle_hydrogens(self, checked: bool) -> None:
        """Toggle visibility of implicit hydrogens."""
        self.canvas.state.show_implicit_hydrogens = checked
        self.canvas.refresh_atom_visibility()
        self.statusBar().showMessage(
            "Hidrógenos: visibles" if checked else "Hidrógenos: ocultos"
        )
    
    def _on_toggle_aromatic_circles(self, checked: bool) -> None:
        """Toggle aromatic circle display mode."""
        self.canvas.state.use_aromatic_circles = checked
        self.canvas.refresh_aromatic_circles()
        self.statusBar().showMessage(
            "Aromáticos: círculos" if checked else "Aromáticos: Kekulé"
        )

    def _on_toggle_rules(self, checked: bool) -> None:
        """Toggle rulers on the canvas."""
        self.canvas.set_show_rulers(checked)
        self.statusBar().showMessage(
            "Reglas: visibles" if checked else "Reglas: ocultas"
        )

    def _on_toggle_crosshair(self, checked: bool) -> None:
        """Toggle grid display on the canvas."""
        self.canvas.set_show_grid(checked)
        self.statusBar().showMessage(
            "Cuadrícula: visible" if checked else "Cuadrícula: oculta"
        )
    
    def _on_zoom_in(self) -> None:
        """Zoom in the canvas."""
        self.canvas.zoom_in()
    
    def _on_zoom_out(self) -> None:
        """Zoom out the canvas."""
        self.canvas.zoom_out()
    
    def _on_zoom_reset(self) -> None:
        """Reset zoom to 100% and center view on paper."""
        self.canvas.resetTransform()
        self.canvas.center_on_paper()
        self.statusBar().showMessage("Zoom: 100%")

    def _on_rotate_selection(self, angle_deg: float) -> None:
        """Maneja rotate selection.

        Args:
            angle_deg: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.rotate_selection_degrees(angle_deg)

    def _on_flip_horizontal(self) -> None:
        """Maneja flip horizontal.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.flip_selection_horizontal()

    def _on_flip_vertical(self) -> None:
        """Maneja flip vertical.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.flip_selection_vertical()

    def _on_rotate_branch(self, angle_deg: float) -> None:
        """Gira la rama menor del enlace seleccionado con snap químico."""
        self.canvas.rotate_selected_branch_degrees(angle_deg)

    def _on_invert_branch(self) -> None:
        """Invierte 180° la rama menor del enlace seleccionado."""
        self.canvas.invert_selected_branch()

    def _on_auto_arrange_branch(self) -> None:
        """Autoacomoda la rama menor del enlace seleccionado."""
        self.canvas.auto_arrange_selected_branch()

    def _selected_single_atom_id(self) -> Optional[int]:
        """Devuelve el único átomo seleccionado si existe y es válido."""
        atom_ids = [
            atom_id
            for atom_id in self.canvas.state.selected_atoms
            if atom_id in self.canvas.model.atoms
        ]
        if self.canvas.state.selected_bonds or len(atom_ids) != 1:
            return None
        return int(atom_ids[0])

    def _sync_fragment_pivot_actions(self) -> None:
        """Sincroniza acciones del menú de rotación de fragmentos."""
        pivot_atom_id = self.canvas.fragment_pivot_atom_id()
        selected_atom_id = self._selected_single_atom_id()
        has_transform_selection = bool(
            self.canvas.state.selected_atoms or self.canvas.state.selected_bonds
        )

        if selected_atom_id is None:
            self.action_fragment_pivot_set.setText("Definir átomo pivote desde selección")
        else:
            atom = self.canvas.model.get_atom(selected_atom_id)
            self.action_fragment_pivot_set.setText(
                f"Usar átomo {selected_atom_id} ({atom.element}) como pivote"
            )

        if pivot_atom_id is None:
            self.action_fragment_pivot_clear.setText("Limpiar átomo pivote")
            if hasattr(self, "fragment_rotate_menu"):
                self.fragment_rotate_menu.setTitle("Fragmento con pivote")
        else:
            atom = self.canvas.model.get_atom(pivot_atom_id)
            self.action_fragment_pivot_clear.setText(
                f"Limpiar pivote (átomo {pivot_atom_id}, {atom.element})"
            )
            if hasattr(self, "fragment_rotate_menu"):
                self.fragment_rotate_menu.setTitle(
                    f"Fragmento con pivote (átomo {pivot_atom_id})"
                )

        self.action_fragment_pivot_set.setEnabled(selected_atom_id is not None)
        self.action_fragment_pivot_clear.setEnabled(pivot_atom_id is not None)
        can_rotate = pivot_atom_id is not None and has_transform_selection
        self.action_fragment_rotate_minus.setEnabled(can_rotate)
        self.action_fragment_rotate_plus.setEnabled(can_rotate)
        self.action_fragment_invert.setEnabled(can_rotate)

    def _on_set_fragment_pivot(self) -> None:
        """Usa el átomo seleccionado como pivote para rotar fragmentos."""
        atom_id = self._selected_single_atom_id()
        if atom_id is None:
            self.canvas._show_status_message(
                "Selecciona un solo átomo para definir el pivote del fragmento."
            )
            self._sync_fragment_pivot_actions()
            return
        self.canvas.set_fragment_pivot_atom(atom_id)
        self._sync_fragment_pivot_actions()

    def _on_clear_fragment_pivot(self) -> None:
        """Limpia el pivote activo de rotación de fragmentos."""
        self.canvas.set_fragment_pivot_atom(None)
        self._sync_fragment_pivot_actions()

    def _on_rotate_fragment(self, angle_deg: float) -> None:
        """Gira el fragmento seleccionado alrededor del pivote activo."""
        self.canvas.rotate_selected_fragment_around_pivot(angle_deg)
        self._sync_fragment_pivot_actions()

    def _on_invert_fragment(self) -> None:
        """Invierte 180° el fragmento seleccionado alrededor del pivote activo."""
        self._on_rotate_fragment(180.0)

    def _on_bond_thickness_up(self) -> None:
        """Maneja bond thickness up.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.increase_selected_bond_thickness()

    def _on_bond_thickness_down(self) -> None:
        """Maneja bond thickness down.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.decrease_selected_bond_thickness()

    def _on_bond_thickness_reset(self) -> None:
        """Maneja bond thickness reset.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.reset_selected_bond_thickness()

    # -------------------------------------------------------------------------
    # Text Menu Handlers
    # -------------------------------------------------------------------------
    def _sync_label_menu_state(self) -> None:
        """Sincroniza label menu state.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.action_label_bold.setChecked(self.canvas.state.label_font_bold)
        self.action_label_italic.setChecked(self.canvas.state.label_font_italic)
        self.action_label_underline.setChecked(self.canvas.state.label_font_underline)
        self.action_label_color_element.setChecked(self.canvas.state.use_element_colors)
        self.action_label_color_black.setChecked(not self.canvas.state.use_element_colors)

    def _on_label_font(self) -> None:
        """Maneja label font.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        font, ok = QFontDialog.getFont(
            self.canvas.label_font(),
            self,
            "Fuente de etiquetas",
        )
        if ok:
            self.canvas.apply_label_font(font)
            self._sync_label_menu_state()

    def _on_label_font_size_dialog(self) -> None:
        """Maneja label font size dialog.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        size, ok = QInputDialog.getDouble(
            self,
            "Tamaño de etiquetas",
            "Tamaño (pt):",
            float(self.canvas.current_label_size_value()),
            6.0,
            72.0,
            1,
        )
        if not ok:
            return
        self.canvas.apply_label_font_size(float(size))
        self._sync_label_menu_state()

    def _on_label_bold(self, checked: bool) -> None:
        """Maneja label bold.

        Args:
            checked: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        font = self.canvas.label_font()
        font.setBold(checked)
        self.canvas.apply_label_font(font)
        self._sync_label_menu_state()

    def _on_label_italic(self, checked: bool) -> None:
        """Maneja label italic.

        Args:
            checked: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        font = self.canvas.label_font()
        font.setItalic(checked)
        self.canvas.apply_label_font(font)
        self._sync_label_menu_state()

    def _on_label_underline(self, checked: bool) -> None:
        """Maneja label underline.

        Args:
            checked: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        font = self.canvas.label_font()
        font.setUnderline(checked)
        self.canvas.apply_label_font(font)
        self._sync_label_menu_state()

    def _on_label_font_size(self, delta: float) -> None:
        """Maneja label font size.

        Args:
            delta: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.adjust_label_font_size(delta)

    def _set_canvas_size(self, width: int, height: int) -> None:
        """Método auxiliar para  set canvas size.

        Args:
            width: Descripción del parámetro.
            height: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.set_paper_size(width, height)
        self.statusBar().showMessage(f"Lienzo: {width} x {height} px")

    def _on_canvas_custom_size(self) -> None:
        """Maneja canvas custom size.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        dialog = QDialog(self)
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
            """Método auxiliar para apply unit settings.

            Args:
                unit: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la interfaz.
            """
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

        apply_unit_settings("cm")
        width_spin.setValue(self.canvas.paper_width / px_per_cm)
        height_spin.setValue(self.canvas.paper_height / px_per_cm)

        def on_unit_changed(text: str) -> None:
            """Método auxiliar para on unit changed.

            Args:
                text: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la interfaz.
            """
            old_unit = "cm"
            if unit_combo.property("last_unit"):
                old_unit = unit_combo.property("last_unit")
            unit_combo.setProperty("last_unit", text)

            def to_px(value: float, unit: str) -> float:
                """Método auxiliar para to px.

                Args:
                    value: Descripción del parámetro.
                    unit: Descripción del parámetro.

                Returns:
                    Resultado de la operación o None.

                Side Effects:
                    Puede modificar el estado interno o la interfaz.
                """
                if unit == "px":
                    return value
                if unit == "in":
                    return value * px_per_in
                return value * px_per_cm

            def from_px(value: float, unit: str) -> float:
                """Método auxiliar para from px.

                Args:
                    value: Descripción del parámetro.
                    unit: Descripción del parámetro.

                Returns:
                    Resultado de la operación o None.

                Side Effects:
                    Puede modificar el estado interno o la interfaz.
                """
                if unit == "px":
                    return value
                if unit == "in":
                    return value / px_per_in
                return value / px_per_cm

            current_px_w = to_px(width_spin.value(), old_unit)
            current_px_h = to_px(height_spin.value(), old_unit)
            apply_unit_settings(text)
            width_spin.setValue(from_px(current_px_w, text))
            height_spin.setValue(from_px(current_px_h, text))

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

        self._set_canvas_size(width_px, height_px)

    def _apply_label_script(self, marker: str, title: str) -> None:
        """Aplica label script.

        Args:
            marker: Descripción del parámetro.
            title: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        if self.canvas.state.selected_bonds or len(self.canvas.state.selected_atoms) != 1:
            self.statusBar().showMessage("Selecciona un átomo para aplicar el formato.")
            return
        value, ok = QInputDialog.getText(self, title, "Texto:")
        if not ok:
            return
        cleaned = value.strip()
        if not cleaned:
            return
        atom_id = next(iter(self.canvas.state.selected_atoms))
        atom = self.canvas.model.get_atom(atom_id)
        label = atom.element
        charge = ""
        if label and label[-1] in "+-":
            charge = label[-1]
            label = label[:-1]
        new_label = f"{label}{marker}{cleaned}{charge}"
        cmd = ChangeAtomCommand(self.canvas.model, self.canvas, atom_id, new_label)
        self.canvas.undo_stack.push(cmd)

    def _on_label_subscript(self) -> None:
        """Maneja label subscript.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self._apply_label_script("_", "Subíndice")

    def _on_label_superscript(self) -> None:
        """Maneja label superscript.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self._apply_label_script("^", "Superíndice")

    def _on_label_color_mode(self, use_element_colors: bool) -> None:
        """Maneja label color mode.

        Args:
            use_element_colors: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.set_use_element_colors(use_element_colors)
        self._sync_label_menu_state()
        self.statusBar().showMessage(
            "Etiquetas: por elemento" if use_element_colors else "Etiquetas: negro"
        )
    
    # -------------------------------------------------------------------------
    # Structure Menu Handlers
    # -------------------------------------------------------------------------
    @staticmethod
    def _coords_center(coords: dict[int, tuple[float, float]]) -> tuple[float, float]:
        """Calcula el centro geométrico de un conjunto de coordenadas."""
        xs = [x for x, _ in coords.values()]
        ys = [y for _, y in coords.values()]
        if not xs:
            return 0.0, 0.0
        return sum(xs) / len(xs), sum(ys) / len(ys)

    @staticmethod
    def _average_bond_length(coords: dict[int, tuple[float, float]], bonds) -> float:
        """Promedio de longitudes de enlaces disponibles en `coords`."""
        lengths: list[float] = []
        for bond in bonds:
            p1 = coords.get(bond.a1_id)
            p2 = coords.get(bond.a2_id)
            if p1 is None or p2 is None:
                continue
            dx = p2[0] - p1[0]
            dy = p2[1] - p1[1]
            dist = math.hypot(dx, dy)
            if dist > 1e-6:
                lengths.append(dist)
        if not lengths:
            return 0.0
        return sum(lengths) / len(lengths)

    @classmethod
    def _rescale_coords_to_bond_length(
        cls,
        coords: dict[int, tuple[float, float]],
        bonds,
        target_length: float,
    ) -> dict[int, tuple[float, float]]:
        """Reescala coordenadas alrededor de su centro para ajustar longitud media."""
        if not coords:
            return {}
        target = float(target_length)
        if target <= 1e-6:
            return dict(coords)
        current = cls._average_bond_length(coords, bonds)
        if current <= 1e-6:
            return dict(coords)
        scale = target / current
        if not math.isfinite(scale):
            return dict(coords)
        # Evita explosiones numéricas por geometrías degeneradas.
        scale = max(0.05, min(200.0, scale))
        cx, cy = cls._coords_center(coords)
        return {
            atom_id: (cx + (x - cx) * scale, cy + (y - cy) * scale)
            for atom_id, (x, y) in coords.items()
        }

    @classmethod
    def _align_coords_to_reference(
        cls,
        reference: dict[int, tuple[float, float]],
        coords: dict[int, tuple[float, float]],
    ) -> dict[int, tuple[float, float]]:
        """Alinea `coords` a la pose actual usando rotación rígida o reflexión."""
        if not reference or not coords:
            return dict(coords)

        common_ids = [atom_id for atom_id in reference if atom_id in coords]
        if not common_ids:
            return dict(coords)

        ref_common = {atom_id: reference[atom_id] for atom_id in common_ids}
        coord_common = {atom_id: coords[atom_id] for atom_id in common_ids}
        ref_cx, ref_cy = cls._coords_center(ref_common)
        src_cx, src_cy = cls._coords_center(coord_common)

        if len(common_ids) == 1:
            dx = ref_cx - src_cx
            dy = ref_cy - src_cy
            return {
                atom_id: (x + dx, y + dy)
                for atom_id, (x, y) in coords.items()
            }

        def _candidate(mirror_x: bool) -> tuple[dict[int, tuple[float, float]], float]:
            sum_cos = 0.0
            sum_sin = 0.0
            for atom_id in common_ids:
                qx, qy = coords[atom_id]
                px, py = reference[atom_id]
                qx -= src_cx
                qy -= src_cy
                px -= ref_cx
                py -= ref_cy
                if mirror_x:
                    qx = -qx
                sum_cos += qx * px + qy * py
                sum_sin += qx * py - qy * px

            theta = math.atan2(sum_sin, sum_cos)
            cos_t = math.cos(theta)
            sin_t = math.sin(theta)

            transformed: dict[int, tuple[float, float]] = {}
            error = 0.0
            for atom_id, (x, y) in coords.items():
                qx = x - src_cx
                qy = y - src_cy
                if mirror_x:
                    qx = -qx
                rx = qx * cos_t - qy * sin_t + ref_cx
                ry = qx * sin_t + qy * cos_t + ref_cy
                transformed[atom_id] = (rx, ry)
                if atom_id in reference:
                    px, py = reference[atom_id]
                    error += (rx - px) ** 2 + (ry - py) ** 2
            return transformed, error

        direct, direct_error = _candidate(False)
        mirrored, mirrored_error = _candidate(True)
        return mirrored if mirrored_error < direct_error else direct

    @classmethod
    def _project_missing_hydrogen_coords(
        cls,
        before: dict[int, tuple[float, float]],
        after: dict[int, tuple[float, float]],
        bonds,
        atom_elements: dict[int, str],
    ) -> dict[int, tuple[float, float]]:
        """Projeta H faltantes trasladando su vector local desde átomos ancla."""
        projected = dict(after)
        if not before:
            return projected

        adjacency: dict[int, list[int]] = {}
        for bond in bonds:
            adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
            adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)

        for atom_id, element in atom_elements.items():
            if element != "H" or atom_id in projected:
                continue
            before_pos = before.get(atom_id)
            if before_pos is None:
                continue
            candidate_positions: list[tuple[float, float]] = []
            for anchor_id in adjacency.get(atom_id, []):
                anchor_before = before.get(anchor_id)
                anchor_after = projected.get(anchor_id)
                if anchor_before is None or anchor_after is None:
                    continue
                dx = before_pos[0] - anchor_before[0]
                dy = before_pos[1] - anchor_before[1]
                candidate_positions.append((anchor_after[0] + dx, anchor_after[1] + dy))
            if candidate_positions:
                avg_x = sum(pos[0] for pos in candidate_positions) / len(candidate_positions)
                avg_y = sum(pos[1] for pos in candidate_positions) / len(candidate_positions)
                projected[atom_id] = (avg_x, avg_y)

        for atom_id, coord in before.items():
            projected.setdefault(atom_id, coord)
        return projected

    @staticmethod
    def _is_acyclic_structure(atom_ids: set[int], bonds) -> bool:
        """Indica si el subgrafo seleccionado no contiene ciclos."""
        if not atom_ids:
            return True
        adjacency: dict[int, list[int]] = {}
        for bond in bonds:
            if bond.a1_id not in atom_ids or bond.a2_id not in atom_ids:
                continue
            adjacency.setdefault(bond.a1_id, []).append(bond.a2_id)
            adjacency.setdefault(bond.a2_id, []).append(bond.a1_id)

        visited: set[int] = set()
        for start in atom_ids:
            if start in visited:
                continue
            stack: list[tuple[int, int]] = [(start, -1)]
            while stack:
                node, parent = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                for neigh in adjacency.get(node, []):
                    if neigh == parent:
                        continue
                    if neigh in visited:
                        return False
                    stack.append((neigh, node))
        return True

    def _on_clean_2d(self) -> None:
        """Maneja clean 2d.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self._run_clean_2d(step_ratio=0.35, fallback_iterations=30, status_suffix="(paso)")

    def _on_clean_2d_full(self) -> None:
        """Maneja clean 2d full.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self._run_clean_2d(step_ratio=1.0, fallback_iterations=200, status_suffix="(1 paso)")

    def _run_clean_2d(self, step_ratio: float, fallback_iterations: int, status_suffix: str) -> None:
        """Clean 2D coordinates using RDKit or fallback."""
        atom_ids, bonds = self.canvas._selected_structure_ids()
        target_ids = atom_ids if atom_ids else set(self.canvas.model.atoms.keys())
        if not target_ids:
            return
        scale_bonds = bonds if atom_ids else list(self.canvas.model.bonds.values())

        if self._is_acyclic_structure(target_ids, scale_bonds):
            # Para cadenas acíclicas, el layout geométrico interno mantiene mejor
            # ángulos ideales (sp3/sp2/sp) que la depicción RDKit con H explícitos.
            self.canvas.clean_2d_fallback(target_ids, iterations=max(40, fallback_iterations))
            msg = "Selección 2D limpiada" if atom_ids else "Estructura 2D limpiada"
            self.statusBar().showMessage(f"{msg} {status_suffix} (acíclico)")
            return

        try:
            from chemuson.chemio.rdkit_io import molgraph_to_rdkit_with_map
            from rdkit import Chem
            from rdkit.Chem import AllChem

            graph = self.canvas._build_selection_graph(atom_ids, bonds) if atom_ids else self.canvas.graph
            mol, id_map = molgraph_to_rdkit_with_map(graph)
            for aid, rd_idx in id_map.items():
                try:
                    mol.GetAtomWithIdx(rd_idx).SetIntProp("_chemuson_id", int(aid))
                except Exception:
                    continue

            before = {
                aid: (self.canvas.model.get_atom(aid).x, self.canvas.model.get_atom(aid).y)
                for aid in target_ids
            }
            before_cx, before_cy = self._coords_center(before)
            before_avg_len = self._average_bond_length(before, scale_bonds)

            raw_after = {}
            used_no_h_layout = False
            try:
                mol_no_h = Chem.RemoveHs(Chem.Mol(mol), sanitize=True)
                if 0 < mol_no_h.GetNumAtoms() < mol.GetNumAtoms():
                    AllChem.Compute2DCoords(mol_no_h)
                    conf_no_h = mol_no_h.GetConformer()
                    for atom in mol_no_h.GetAtoms():
                        if not atom.HasProp("_chemuson_id"):
                            continue
                        aid = int(atom.GetIntProp("_chemuson_id"))
                        if aid not in target_ids:
                            continue
                        pos = conf_no_h.GetAtomPosition(atom.GetIdx())
                        raw_after[aid] = (pos.x, pos.y)
                    if raw_after:
                        atom_elements = {
                            aid: self.canvas.model.get_atom(aid).element for aid in target_ids
                        }
                        raw_after = self._project_missing_hydrogen_coords(
                            before=before,
                            after=raw_after,
                            bonds=scale_bonds,
                            atom_elements=atom_elements,
                        )
                        used_no_h_layout = True
            except Exception:
                used_no_h_layout = False

            if not used_no_h_layout:
                AllChem.Compute2DCoords(mol)
                conf = mol.GetConformer()
                for aid, rd_idx in id_map.items():
                    if aid not in target_ids:
                        continue
                    pos = conf.GetAtomPosition(rd_idx)
                    raw_after[aid] = (pos.x, pos.y)

            target_bond_len = before_avg_len if before_avg_len > 1e-6 else float(self.canvas.state.bond_length)
            raw_after = self._rescale_coords_to_bond_length(raw_after, scale_bonds, target_bond_len)
            raw_after = self._align_coords_to_reference(before, raw_after)
            after_cx, after_cy = self._coords_center(raw_after)
            after = {}
            for aid, (x, y) in raw_after.items():
                x = x - after_cx + before_cx
                y = y - after_cy + before_cy
                if step_ratio < 1.0:
                    bx, by = before[aid]
                    x = bx + (x - bx) * step_ratio
                    y = by + (y - by) * step_ratio
                after[aid] = (x, y)
            # La interpolación lineal entre dos orientaciones del mismo esqueleto
            # puede reducir temporalmente el tamaño percibido (especialmente anillos).
            # Re-normalizamos la longitud media para mantener escala visual.
            after = self._rescale_coords_to_bond_length(after, scale_bonds, target_bond_len)
            final_cx, final_cy = self._coords_center(after)
            if abs(final_cx - before_cx) > 1e-9 or abs(final_cy - before_cy) > 1e-9:
                after = {
                    aid: (x - final_cx + before_cx, y - final_cy + before_cy)
                    for aid, (x, y) in after.items()
                }

            from chemuson.gui.commands import MoveAtomsCommand
            cmd = MoveAtomsCommand(self.canvas.model, self.canvas, before, after)
            self.canvas.undo_stack.push(cmd)
            self.canvas._update_selection_overlay()
            msg = "Selección 2D limpiada" if atom_ids else "Estructura 2D limpiada"
            self.statusBar().showMessage(f"{msg} {status_suffix}".strip())
            return
        except Exception as e:
            message = str(e)
            if "No module named" in message and "rdkit" in message:
                self.canvas.clean_2d_fallback(target_ids, iterations=fallback_iterations)
                msg = "Selección 2D limpiada" if atom_ids else "Estructura 2D limpiada"
                self.statusBar().showMessage(f"{msg} {status_suffix} (básico)")
                return
            self.statusBar().showMessage(f"Error: {e}")

    def _insert_template(self, label: str, graph) -> None:
        """Método auxiliar para  insert template.

        Args:
            label: Descripción del parámetro.
            graph: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.begin_template_insert_mode(graph, label)
        self.statusBar().showMessage(
            f"Plantilla '{label}' lista. Haz clic para colocar; clic sobre átomo para conectar."
        )

    def _get_templates_dir(self) -> str:
        """Ruta legada para plantillas en disco (`src/templates`)."""
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        templates_dir = os.path.join(base_dir, "templates")
        if not os.path.exists(templates_dir):
            os.makedirs(templates_dir)
        return templates_dir

    def _migrate_legacy_templates(self) -> None:
        """Importa plantillas `.mol` legadas a la nueva biblioteca JSON."""
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
                    str(tpl.get("name", "")).strip().lower(),
                    str(tpl.get("molblock", "")).strip(),
                )
                for tpl in self.template_library.list_templates()
            }
            changed = False
            for filename in files:
                filepath = os.path.join(templates_dir, filename)
                try:
                    with open(filepath, "r", encoding="utf-8") as fh:
                        molblock = fh.read().strip()
                except Exception:
                    continue
                if not molblock:
                    continue
                name = os.path.splitext(filename)[0].replace("_", " ").strip() or "Plantilla"
                signature = (name.lower(), molblock)
                if signature in existing:
                    continue
                self.template_library.add_template(
                    name=name,
                    category=DEFAULT_CATEGORY_USER,
                    molblock=molblock,
                )
                existing.add(signature)
                changed = True
            if changed:
                self.statusBar().showMessage("Plantillas legadas importadas a la biblioteca.")
        except Exception:
            # Migración best-effort: no bloquear inicio por fallos en archivos legados.
            return

    def _template_preview_icon(self, template_id: str) -> QIcon:
        """Genera miniatura de una plantilla para la galería lateral."""
        cache = self._template_icon_cache.get(template_id)
        if cache is not None:
            return cache
        icon = QIcon()
        try:
            graph = self.template_library.graph_from_template(template_id)
            if not graph.atoms:
                self._template_icon_cache[template_id] = icon
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
        self._template_icon_cache[template_id] = icon
        return icon

    def _refresh_template_views(self) -> None:
        """Sincroniza menú y dock con la biblioteca de plantillas."""
        self._template_icon_cache.clear()
        grouped = self.template_library.grouped_templates()
        grouped_with_icons: list[dict] = []
        for group in grouped:
            templates = []
            for template in group.get("templates", []):
                entry = dict(template)
                template_id = str(entry.get("id", "")).strip()
                if template_id:
                    entry["icon"] = self._template_preview_icon(template_id)
                templates.append(entry)
            grouped_with_icons.append({"name": group.get("name", ""), "templates": templates})
        self.templates_dock.set_templates(grouped_with_icons)
        self._refresh_templates_menu()

    def _refresh_templates_menu(self) -> None:
        """Construye menú dinámico de plantillas por categoría."""
        self.templates_menu.clear()
        grouped = self.template_library.grouped_templates()
        total_templates = 0
        for group in grouped:
            category = str(group.get("name", "")).strip()
            if not category:
                continue
            submenu = self.templates_menu.addMenu(category)
            templates = list(group.get("templates", []))
            if not templates:
                empty = QAction("(Vacío)", self)
                empty.setEnabled(False)
                submenu.addAction(empty)
                continue
            for template in templates:
                template_id = str(template.get("id", "")).strip()
                label = str(template.get("name", "")).strip() or "Plantilla"
                if not template_id:
                    continue
                action = QAction(label, self)
                action.triggered.connect(
                    lambda checked=False, tid=template_id: self._start_template_insert_by_id(tid)
                )
                submenu.addAction(action)
                total_templates += 1
        if total_templates == 0:
            empty = QAction("(Sin plantillas)", self)
            empty.setEnabled(False)
            self.templates_menu.addAction(empty)
        self.templates_menu.addSeparator()
        self.templates_menu.addAction(self.action_save_template)
        self.templates_menu.addAction(self.action_template_new_category)
        self.templates_menu.addAction(self.action_template_import_library)
        self.templates_menu.addAction(self.action_template_export_library)

    def _start_template_insert_by_id(self, template_id: str, *, place_now: bool = False) -> None:
        """Carga plantilla desde biblioteca e inicia inserción.

        Args:
            template_id: ID de plantilla en biblioteca.
            place_now: Si es `True`, inserta inmediatamente en el lienzo.
        """
        try:
            template = self.template_library.get_template(template_id)
            graph = self.template_library.graph_from_template(template_id)
            label = str(template.get("name", "Plantilla")).strip() or "Plantilla"
            if place_now:
                target = self.canvas._last_scene_pos
                if target is None:
                    target = self.canvas.mapToScene(self.canvas.viewport().rect().center())
                self.canvas._insert_molgraph_at(graph, target)
                self.canvas.cancel_template_insert_mode()
                self.statusBar().showMessage(f"Plantilla '{label}' insertada")
            else:
                self._insert_template(label, graph)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"No se pudo cargar la plantilla:\n{e}")

    def _on_template_selected_from_gallery(self, payload: dict) -> None:
        """Maneja selección de plantilla desde el dock lateral."""
        template_id = str(payload.get("id", "")).strip()
        if not template_id:
            return
        self._start_template_insert_by_id(template_id)

    def _prompt_template_destination(self) -> Optional[tuple[str, str]]:
        """Solicita nombre y categoría para guardar plantilla."""
        name, ok = QInputDialog.getText(self, "Guardar plantilla", "Nombre de la plantilla:")
        if not ok:
            return None
        clean_name = " ".join(name.strip().split())
        if not clean_name:
            return None
        categories = self.template_library.categories()
        if not categories:
            categories = [DEFAULT_CATEGORY_USER]
        if DEFAULT_CATEGORY_USER not in categories:
            categories.append(DEFAULT_CATEGORY_USER)
        choices = categories + ["+ Nueva categoría..."]
        default_idx = choices.index(DEFAULT_CATEGORY_USER) if DEFAULT_CATEGORY_USER in choices else 0
        category_choice, ok = QInputDialog.getItem(
            self,
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
            new_category, ok = QInputDialog.getText(self, "Nueva categoría", "Nombre de la categoría:")
            if not ok:
                return None
            category = self.template_library.add_category(new_category)
        else:
            category = self.template_library.add_category(category_choice)
        return clean_name, category

    def _on_save_template(self) -> None:
        """Guarda selección (o molécula completa) como plantilla reusable."""
        try:
            atom_ids, bonds = self.canvas._selected_structure_ids()
            graph_to_save = (
                self.canvas._build_selection_graph(atom_ids, bonds)
                if atom_ids
                else self.canvas.graph
            )
            if not graph_to_save.atoms:
                QMessageBox.warning(self, "Aviso", "No hay nada para guardar.")
                return
            target = self._prompt_template_destination()
            if target is None:
                return
            name, category = target
            self.template_library.add_template_from_graph(name, category, graph_to_save)
            self._refresh_template_views()
            self.statusBar().showMessage(f"Plantilla guardada: {name} ({category})")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error al guardar plantilla:\n{e}")

    def _on_new_template_category(self) -> None:
        """Crea categoría personalizada de plantillas."""
        name, ok = QInputDialog.getText(self, "Nueva categoría", "Nombre de la categoría:")
        if not ok:
            return
        category = " ".join(name.strip().split())
        if not category:
            return
        self.template_library.add_category(category)
        self._refresh_template_views()
        self.statusBar().showMessage(f"Categoría creada: {category}")

    def _on_rename_template_category(self, old_name: str) -> None:
        """Renombra una categoría existente."""
        if not old_name:
            return
        new_name, ok = QInputDialog.getText(
            self,
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
            self.template_library.rename_category(old_name, clean_new)
            self._refresh_template_views()
            self.statusBar().showMessage(f"Categoría renombrada: {old_name} -> {clean_new}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"No se pudo renombrar la categoría:\n{e}")

    def _on_delete_template_category(self, name: str) -> None:
        """Elimina una categoría y conserva sus plantillas en categoría de respaldo."""
        if not name:
            return
        reply = QMessageBox.question(
            self,
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
            self.template_library.delete_category(name, fallback_category=DEFAULT_CATEGORY_USER)
            self._refresh_template_views()
            self.statusBar().showMessage(f"Categoría eliminada: {name}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"No se pudo eliminar la categoría:\n{e}")

    def _on_rename_template(self, template_id: str) -> None:
        """Renombra plantilla individual."""
        if not template_id:
            return
        try:
            template = self.template_library.get_template(template_id)
        except Exception:
            return
        current_name = str(template.get("name", "Plantilla"))
        new_name, ok = QInputDialog.getText(
            self,
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
            self.template_library.rename_template(template_id, clean_new)
            self._refresh_template_views()
            self.statusBar().showMessage(f"Plantilla renombrada: {clean_new}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"No se pudo renombrar la plantilla:\n{e}")

    def _on_delete_template(self, template_id: str) -> None:
        """Elimina plantilla de la biblioteca."""
        if not template_id:
            return
        try:
            template = self.template_library.get_template(template_id)
        except Exception:
            return
        name = str(template.get("name", "Plantilla"))
        reply = QMessageBox.question(
            self,
            "Eliminar plantilla",
            f"¿Eliminar la plantilla '{name}'?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if reply != QMessageBox.StandardButton.Yes:
            return
        try:
            self.template_library.delete_template(template_id)
            self._refresh_template_views()
            self.statusBar().showMessage(f"Plantilla eliminada: {name}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"No se pudo eliminar la plantilla:\n{e}")

    def _on_import_template_library(self) -> None:
        """Importa bibliotecas de plantillas desde JSON."""
        filepath, _ = QFileDialog.getOpenFileName(
            self,
            "Importar biblioteca de plantillas",
            "",
            "Template Library (*.json)",
        )
        if not filepath:
            return
        choice = QMessageBox.question(
            self,
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
            added = self.template_library.import_from_file(filepath, merge=merge)
            self._refresh_template_views()
            if merge:
                self.statusBar().showMessage(f"Biblioteca importada: {added} plantilla(s) agregadas.")
            else:
                self.statusBar().showMessage("Biblioteca importada (reemplazo completo).")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"No se pudo importar la biblioteca:\n{e}")

    def _on_export_template_library(self) -> None:
        """Exporta biblioteca de plantillas a JSON."""
        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Exportar biblioteca de plantillas",
            "chemuson_templates.json",
            "Template Library (*.json)",
        )
        if not filepath:
            return
        try:
            self.template_library.export_to_file(filepath)
            self.statusBar().showMessage("Biblioteca de plantillas exportada.")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"No se pudo exportar la biblioteca:\n{e}")

    def _insert_template_from_file(self, filepath: str) -> None:
        """Método auxiliar para  insert template from file.

        Args:
            filepath: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                molblock = f.read()
            graph = molfile_to_molgraph(molblock)
            name = os.path.splitext(os.path.basename(filepath))[0]
            self._insert_template(name, graph)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error al cargar plantilla:\n{e}")

    def _on_insert_linear_chain(self) -> None:
        """Maneja insert linear chain.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        graph = build_linear_chain_template(self.canvas.state.bond_length)
        self._insert_template("Cadena lineal", graph)

    def _on_import_smiles(self) -> None:
        """Import a molecule from a SMILES string."""
        smiles, ok = QInputDialog.getText(self, "Importar SMILES", "SMILES:")
        if not ok or not smiles.strip():
            return
        try:
            from chemuson.chemio.rdkit_io import smiles_to_molgraph
            graph = smiles_to_molgraph(smiles.strip())
            self.canvas._insert_molgraph(graph)
            self.statusBar().showMessage("SMILES insertado")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"No se pudo importar SMILES:\n{e}")

    def _on_export_smiles(self) -> None:
        """Export the current molecule as SMILES."""
        try:
            from chemuson.chemio.rdkit_io import molgraph_to_smiles
            smiles = molgraph_to_smiles(self.canvas.graph)
            QMessageBox.information(self, "SMILES", smiles)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"No se pudo exportar SMILES:\n{e}")
    
    # -------------------------------------------------------------------------
    # Help Menu Handlers
    # -------------------------------------------------------------------------
    def _on_quick_start(self) -> None:
        """Show quick start dialog."""
        dialog = QuickStartDialog(self)
        dialog.exec()
    
    def _on_about(self) -> None:
        """Show about dialog."""
        version = get_app_version()
        QMessageBox.about(
            self,
            "Acerca de Chemuson",
            "<h2>Chemuson</h2>"
            "<p>Editor Molecular Libre</p>"
            f"<p>Versión {version}</p>"
            "<p>Un editor de estructuras químicas de código abierto "
            "inspirado en ChemDoodle.</p>"
        )
    
    # -------------------------------------------------------------------------
    # Toolbar Handlers
    # -------------------------------------------------------------------------
    def _handle_bond_palette(self, bond_spec: dict) -> None:
        """Método auxiliar para  handle bond palette.

        Args:
            bond_spec: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.set_active_bond(bond_spec)
        self._update_status("tool_bond")

    def _handle_ring_palette(self, ring_spec: dict) -> None:
        """Método auxiliar para  handle ring palette.

        Args:
            ring_spec: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.set_active_ring(ring_spec)
        self._update_status("tool_ring")

    def _handle_element_palette(self, element: str) -> None:
        """Método auxiliar para  handle element palette.

        Args:
            element: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self.canvas.set_active_element(element)
        self._update_status("tool_atom")

    def _show_periodic_table(self) -> None:
        """Método auxiliar para  show periodic table.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        dialog = PeriodicTableDialog(self)
        dialog.element_selected.connect(self.toolbar.select_element)
        dialog.exec()
    
    def _update_status(self, tool_id: str) -> None:
        """Update status bar with current tool."""
        ring_label = f"Anillo {self.canvas.state.active_ring_size}"
        orbital_label = orbital_display_name(self.canvas.state.active_orbital_kind)
        if self.canvas.state.active_ring_template:
            template_name = {
                "haworth": "Haworth",
                "chair": "Silla",
            }.get(self.canvas.state.active_ring_template, "Anillo")
            anomeric = self.canvas.state.active_ring_anomeric
            suffix = ""
            if anomeric == "alpha":
                suffix = " α"
            elif anomeric == "beta":
                suffix = " β"
            ring_label = f"{template_name}{suffix}".strip()
        tool_names = {
            "tool_none": "Ninguna",
            "tool_select": "Seleccionar",
            "tool_select_lasso": "Seleccion (lazo)",
            "tool_bond": "Enlace",
            "tool_coordination_bond": "Enlace coordinativo",
            "tool_rotate_3d_precise": "Rotación 3D precisa",
            "tool_ring": ring_label,
            "tool_atom": f"Elemento {self.canvas.state.default_element}",
            "tool_orbital": orbital_label,
            "tool_coordination_center": "Centro de coordinación (esfera)",
            "tool_chain": "Cadena",
            "tool_arrow_forward": "Flecha directa",
            "tool_arrow_forward_open": "Flecha directa abierta",
            "tool_arrow_forward_dashed": "Flecha directa discontinua",
            "tool_arrow_retro": "Flecha retro",
            "tool_arrow_retro_open": "Flecha retro abierta",
            "tool_arrow_retro_dashed": "Flecha retro discontinua",
            "tool_arrow_both": "Flecha doble",
            "tool_arrow_both_open": "Flecha doble abierta",
            "tool_arrow_both_dashed": "Flecha doble discontinua",
            "tool_arrow_equilibrium": "Equilibrio",
            "tool_arrow_equilibrium_dashed": "Equilibrio discontinuo",
            "tool_arrow_retrosynthetic": "Flecha retrosintesis",
            "tool_arrow_curved": "Flecha curva",
            "tool_arrow_curved_fishhook": "Flecha curva (1 e-)",
            "tool_brackets_round": "Parentesis",
            "tool_brackets_square": "Corchetes",
            "tool_brackets_square_left": "Corchete izquierdo",
            "tool_brackets_square_right": "Corchete derecho",
            "tool_brackets_corner": "Esquinas",
            "tool_brackets_curly": "Llaves",
            "tool_brackets_curly_left": "Llave izquierda",
            "tool_brackets_curly_right": "Llave derecha",
            "tool_brackets_frame": "Marco",
            "tool_brackets_frame_rounded": "Marco redondeado",
            "tool_charge_plus": "Carga positiva",
            "tool_charge_minus": "Carga negativa",
            "tool_charge": "Carga alterna",
            "tool_symbol_plus": "Signo más",
            "tool_symbol_minus": "Signo menos",
            "tool_symbol_radical": "Electrón desapareado",
            "tool_symbol_lone_pair": "Par solitario",
            "tool_symbol_wavy_anchor": "Ancla ondulada",
            "tool_symbol_radical_cation": "Radical catión",
            "tool_symbol_radical_anion": "Radical anión",
            "tool_symbol_partial_plus": "Carga parcial (+)",
            "tool_symbol_partial_minus": "Carga parcial (-)",
        }
        tool_names.update(
            {
                orbital_tool_id(kind): orbital_display_name(kind)
                for kind in ORBITAL_MENU_ORDER
            }
        )
        name = tool_names.get(tool_id, tool_id)
        self.statusBar().showMessage(f"Herramienta: {name}")

    def _on_selection_changed(self, num_atoms: int, num_bonds: int, num_text: int, details: dict):
        """Handle selection change to update UI components."""
        self.inspector_dock.update_selection(num_atoms, num_bonds, num_text, details)
        self._update_total_charge_indicator()
        self._sync_fragment_pivot_actions()
        
        # Sync Text Toolbar if a single text item is selected
        if num_text == 1 and details.get("type") == "text":
            font = details.get("font")
            settings = {
                "color": details.get("color"),
                "sub": details.get("sub"),
                "sup": details.get("sup")
            }
            self.text_toolbar.set_state(font, settings)

    def _update_total_charge_indicator(self) -> None:
        """Actualiza el indicador de carga total en la barra de estado."""
        charge = int(self.canvas.model.total_formal_charge())
        if charge > 0:
            charge_text = f"+{charge}"
        else:
            charge_text = str(charge)
        self._total_charge_label.setText(f"Carga total: {charge_text}")


def run_app() -> None:
    """Punto de entrada de la GUI de Chemuson."""
    crash_reporter.install()
    try:
        app = QApplication(sys.argv)
        app.setApplicationName("Chemuson")
        app.setApplicationVersion(get_app_version())
        window = ChemusonWindow()
        ChemusonWindow.check_autosaves(window)
        window.show()
        exit_code = app.exec()
    except Exception as exc:
        log_path = crash_reporter.write_crash_log(exc)
        if QApplication.instance() is not None:
            QMessageBox.critical(
                None,
                "Chemuson - Error crítico",
                "No se pudo iniciar la aplicación.\n"
                f"Se guardó un reporte en:\n{log_path}",
            )
        else:
            sys.stderr.write(
                "No se pudo iniciar la aplicación.\n"
                f"Se guardó un reporte en: {log_path}\n"
            )
        return
    sys.exit(exit_code)
