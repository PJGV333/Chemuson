"""Composición de alto nivel de :class:`ChemusonWindow`."""

from __future__ import annotations

from typing import Optional

from PyQt6.QtCore import QSettings, QTimer, Qt
from PyQt6.QtWidgets import QLabel, QTabWidget, QTextEdit

from chemuson.gui.actions import (
    create_edit_actions,
    create_file_actions,
    create_project_actions,
    create_update_actions,
)
from chemuson.gui.canvas import ChemusonCanvas
from chemuson.gui.controllers import (
    Clean2DController,
    CompChem3DController,
    DocumentController,
    ExportController,
    FileController,
    RecoveryController,
    TemplateController,
    TextFormatController,
    UpdateController,
    ValidationController,
    ViewController,
)
from chemuson.gui.docks import (
    AppearanceDock,
    ChemicalPropertiesDock,
    CompChemDock,
    InspectorDock,
    PlantillasDock,
    SpectroscopyDock,
    ValidationDock,
)
from chemuson.gui.main_window_ui_builder import MainWindowUiBuilder
from chemuson.gui.semantic_diagram_workflow import SemanticDiagramWorkflow
from chemuson.gui.tab_manager import CanvasTabManager
from chemuson.gui.template_browser_service import TemplateBrowserService
from chemuson.gui.template_library import TemplateLibrary
from chemuson.gui.text_toolbar import TextFormatToolbar
from chemuson.gui.toolbar import ChemusonToolbar, SymbolPaletteToolbar
from chemuson.version import get_app_version


def assemble_application_shell(self) -> None:
    """Ensambla regiones, colaboradores y conexiones de la ventana principal."""
    self._ui_builder = MainWindowUiBuilder()
    self._app_version = get_app_version()
    self.current_theme = "light"
    self.setWindowTitle(f"Chemuson {self._app_version} - Editor Molecular Libre")
    self.resize(1200, 900)

    # === CORE COMPONENTS ===
    create_file_actions(self)
    create_edit_actions(self)
    create_project_actions(self)
    self._ui_builder.create_local_actions(self)
    create_update_actions(self)
    self._current_file_path: Optional[str] = None
    self._active_canvas_connected: Optional[ChemusonCanvas] = None
    self._numbering_default_mode = "atoms"
    self._numbering_default_include_export = True
    self._current_tool_id = "tool_select"
    self._settings = QSettings("Chemuson", "Chemuson")
    self._document_controller = DocumentController()
    self._file_controller = FileController()
    self._update_controller = UpdateController(self._settings, self._app_version)
    self._update_settings = self._update_controller.settings
    self._name_advanced_default, self._name_rdkit_isolated_default = (
        self._load_naming_preferences()
    )
    self._recent_files = self._load_recent_files()
    self.template_library = TemplateLibrary()
    self._template_icon_cache: dict[str, object] = {}
    self._template_browser = TemplateBrowserService()
    self._semantic_diagram_workflow = SemanticDiagramWorkflow()
    self._view_controller = ViewController()
    self._export_controller = ExportController()
    self._clean2d_controller = Clean2DController()
    self._validation_controller = ValidationController()
    self._compchem3d_controller = CompChem3DController(self)
    self._text_format_controller = TextFormatController()
    self._recovery_controller = RecoveryController()
    self._template_controller = TemplateController()
    self._compchem_coordset = None
    self._latest_compchem_job_id = 0
    self._compchem_job_backends: dict[int, str] = {}
    self._compchem_job_operations: dict[int, str] = {}
    self._compchem3d_controller.frame_ready.connect(self._on_compchem_frame_ready)
    self._compchem3d_controller.job_finished.connect(self._on_compchem_job_finished)

    # === CENTRAL TABS/CANVAS ===
    self.tabs = QTabWidget(self)
    self.tabs.setDocumentMode(True)
    self.tabs.setTabsClosable(True)
    self.tabs.setMovable(True)
    self.setCentralWidget(self.tabs)
    self._tab_manager = CanvasTabManager(self.tabs, autosave_parent=self)
    self._canvas_file_paths = self._tab_manager.file_paths
    self._canvas_tab_titles = self._tab_manager.tab_titles
    self._canvas_autosave_managers = self._tab_manager.autosave_managers
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
    self.templates_dock.template_selected.connect(
        self._on_template_selected_from_gallery
    )
    self.templates_dock.new_category_requested.connect(self._on_new_template_category)
    self.templates_dock.import_requested.connect(self._on_import_template_library)
    self.templates_dock.export_requested.connect(self._on_export_template_library)
    self.templates_dock.rename_category_requested.connect(
        self._on_rename_template_category
    )
    self.templates_dock.delete_category_requested.connect(
        self._on_delete_template_category
    )
    self.templates_dock.rename_template_requested.connect(self._on_rename_template)
    self.templates_dock.delete_template_requested.connect(self._on_delete_template)

    self.inspector_dock = InspectorDock(self)
    self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.inspector_dock)
    self.inspector_dock.hide()

    self.validation_dock = ValidationDock(self)
    self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.validation_dock)
    self.validation_dock.hide()
    self.validation_dock.issue_selected.connect(self._select_validation_issue_from_dock)
    self.validation_dock.correction_requested.connect(
        self._on_validation_correction_requested
    )
    self.validation_dock.refresh_requested.connect(self._on_validate_structure)
    self.validation_dock.next_requested.connect(
        lambda: self._on_navigate_validation_issue(1)
    )
    self.validation_dock.previous_requested.connect(
        lambda: self._on_navigate_validation_issue(-1)
    )

    self.chemical_properties_dock = ChemicalPropertiesDock(self)
    self.addDockWidget(
        Qt.DockWidgetArea.RightDockWidgetArea, self.chemical_properties_dock
    )
    self.chemical_properties_dock.hide()

    self.spectroscopy_dock = SpectroscopyDock(self)
    self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.spectroscopy_dock)
    self.spectroscopy_dock.hide()
    self.spectroscopy_dock.peak_atom_selected.connect(self._select_atom_from_spectrum)

    self.compchem_dock = CompChemDock(self)
    self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.compchem_dock)
    self.compchem_dock.hide()
    self.compchem_dock.generate_requested.connect(self._on_compchem_generate)
    self.compchem_dock.optimize_requested.connect(self._on_compchem_optimize)
    self.compchem_dock.project_requested.connect(self._on_compchem_project_to_2d)
    self.compchem_dock.reset_requested.connect(self._on_compchem_reset)
    self.compchem_dock.export_xyz_requested.connect(self._on_compchem_export_xyz)
    self.compchem_dock.export_input_requested.connect(self._on_compchem_export_input)

    self._properties_update_timer = QTimer(self)
    self._properties_update_timer.setSingleShot(True)
    self._properties_update_timer.setInterval(120)
    self._properties_update_timer.timeout.connect(
        self._refresh_chemical_properties_dock
    )
    self._descriptor_jobs = {}
    self._next_descriptor_job_id = 1
    self._latest_descriptor_job_id = 0
    self._name2structure_jobs = {}
    self._cancelled_name2structure_jobs: set[int] = set()
    self._next_name2structure_job_id = 1

    self.appearance_dock = AppearanceDock(self.canvas.drawing_style, self)
    self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.appearance_dock)
    self.appearance_dock.hide()
    self.appearance_dock.appearance_changed.connect(self._apply_appearance_settings)

    # === LEFT TOOLBAR (Drawing tools) ===
    self.toolbar = ChemusonToolbar()
    self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, self.toolbar)

    # === RIGHT SYMBOLS TOOLBAR ===
    self.symbols_toolbar = SymbolPaletteToolbar(action_group=self.toolbar.action_group)
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
    self._template_browser.migrate_legacy_templates(self._template_browser_context())
    self._template_browser.refresh_template_views(self._template_browser_context())

    # === TEXT FORMAT TOOLBAR ===
    self.text_toolbar = TextFormatToolbar()
    self.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.text_toolbar)
    self._external_text_editor: QTextEdit | None = None
    self._external_text_cursor_state: tuple[int, int] | None = None
    self._external_text_selected_range: tuple[int, int] | None = None
    self.text_toolbar.format_changed.connect(self._on_text_format_changed)
    self.text_toolbar.color_changed.connect(self._on_text_color_changed)
    self.text_toolbar.alignment_changed.connect(self._on_text_alignment_changed)
    self.text_toolbar.opacity_changed.connect(self._on_opacity_changed)

    # === SIGNAL CONNECTIONS ===
    self._connect_undo_redo()
    self.toolbar.tool_changed.connect(self._on_tool_changed)
    self.toolbar.bond_palette_changed.connect(self._handle_bond_palette)
    self.toolbar.ring_palette_changed.connect(self._handle_ring_palette)
    self.toolbar.element_palette_changed.connect(self._handle_element_palette)
    self.toolbar.periodic_table_requested.connect(self._show_periodic_table)
    self.symbols_toolbar.tool_changed.connect(self._on_tool_changed)
    self.symbols_toolbar.tool_changed.connect(self._update_status)
    self.symbols_toolbar.atomic_diagram_requested.connect(
        self._open_atomic_diagram_dialog
    )
    self.symbols_toolbar.diatomic_mo_diagram_requested.connect(
        self._open_diatomic_mo_diagram_dialog
    )
    self.symbols_toolbar.ligand_field_diagram_requested.connect(
        self._open_ligand_field_diagram_dialog
    )
    self.symbols_toolbar.electronic_diagram_preset_requested.connect(
        self._insert_semantic_preset
    )
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
    self.toolbar.tool_changed.connect(self._update_status)
    self._set_active_canvas(self.canvas)
    self._apply_theme()
    QTimer.singleShot(1200, self._maybe_check_updates_startup)
