"""
Ventana principal de Chemuson.

Compone menús, barras de herramientas, docks y el lienzo central.
"""
from PyQt6.QtWidgets import (
    QApplication,
    QDialog,
    QMainWindow,
    QFileDialog,
    QMessageBox,
    QTabWidget,
    QLabel,
    QTextEdit,
)
from PyQt6.QtCore import Qt, QSettings, QEvent, QTimer, QPointF
from PyQt6.QtGui import QAction, QColor, QTextCursor
from typing import Optional
import math
import os
import sys

from chemuson.gui.canvas import (
    ChemusonCanvas,
)
from chemuson.gui.elemental_analysis_dialog import ElementalAnalysisDialog
from chemuson.gui.periodic_table import PeriodicTableDialog
from chemuson.gui.toolbar import ChemusonToolbar, SymbolPaletteToolbar
from chemuson.gui.styles import get_main_stylesheet, get_tool_palette_stylesheet
from chemuson.gui.icons import set_icon_theme
from chemuson.gui.docks import (
    AppearanceDock,
    ChemicalPropertiesDock,
    InspectorDock,
    PlantillasDock,
)
from chemuson.gui.dialogs import PreferencesDialog, QuickStartDialog, StyleDialog
from chemuson.gui.text_toolbar import TextFormatToolbar
from chemuson.gui.template_library import TemplateLibrary, DEFAULT_CATEGORY_USER
from chemuson.gui.template_browser_service import (
    TemplateBrowserContext,
    TemplateBrowserService,
)
from chemuson.gui.main_window_ui_builder import MainWindowUiBuilder
from chemuson.gui.semantic_diagram_workflow import (
    SemanticDiagramWorkflow,
    SemanticDiagramWorkflowContext,
)
from chemuson.gui.actions import (
    create_edit_actions,
    create_file_actions,
    create_project_actions,
    create_update_actions,
)
from chemuson.gui.controllers import (
    Clean2DController,
    DocumentController,
    DocumentDiscardContext,
    DocumentTabsContext,
    ExportController,
    FileController,
    FileWorkflowContext,
    RecentFilesContext,
    RecoveryController,
    TemplateController,
    TemplateControllerContext,
    TextFormatController,
    UpdateController,
    UpdateControllerContext,
    ViewController,
)
from chemuson.gui.controllers.update_controller import (
    FLATPAK_APP_ID as FLATPAK_APP_ID,
    format_no_update_message as format_no_update_message,
    format_update_disabled_message as format_update_disabled_message,
)
from chemuson.gui.rich_text_dialog_service import (
    open_rich_text_value_dialog,
    rich_text_editor_value,
)
from chemuson.gui.tab_manager import CanvasTabManager
from chemuson.update import UpdateSettings
from chemuson.utils import crash_reporter
from chemuson.version import get_app_version

__all__ = [
    "ChemusonWindow",
    "FLATPAK_APP_ID",
    "format_no_update_message",
    "format_update_disabled_message",
]


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
        self._ui_builder = MainWindowUiBuilder()
        self._app_version = get_app_version()
        self.current_theme = "light"
        self.setWindowTitle(f"Chemuson {self._app_version} - Editor Molecular Libre")
        self.resize(1200, 900)

        # === CORE COMPONENTS ===
        self._create_actions()
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
        self._name_advanced_default, self._name_rdkit_isolated_default = self._load_naming_preferences()
        self._recent_files = self._load_recent_files()
        self.template_library = TemplateLibrary()
        self._template_icon_cache: dict[str, object] = {}
        self._template_browser = TemplateBrowserService()
        self._semantic_diagram_workflow = SemanticDiagramWorkflow()
        self._view_controller = ViewController()
        self._export_controller = ExportController()
        self._clean2d_controller = Clean2DController()
        self._text_format_controller = TextFormatController()
        self._recovery_controller = RecoveryController()
        self._template_controller = TemplateController()
        
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

        self.chemical_properties_dock = ChemicalPropertiesDock(self)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.chemical_properties_dock)
        self.chemical_properties_dock.hide()
        self._properties_update_timer = QTimer(self)
        self._properties_update_timer.setSingleShot(True)
        self._properties_update_timer.setInterval(120)
        self._properties_update_timer.timeout.connect(self._refresh_chemical_properties_dock)

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
        # Initially hide it, or show it only when text tool is active?
        # User requested it to be available. We can leave it visible or toggle it.
        # For now, visible is fine.

        # Text Toolbar Connections
        self.text_toolbar.format_changed.connect(self._on_text_format_changed)
        self.text_toolbar.color_changed.connect(self._on_text_color_changed)
        self.text_toolbar.alignment_changed.connect(self._on_text_alignment_changed)
        self.text_toolbar.opacity_changed.connect(self._on_opacity_changed)
        
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
        self.symbols_toolbar.atomic_diagram_requested.connect(self._open_atomic_diagram_dialog)
        self.symbols_toolbar.diatomic_mo_diagram_requested.connect(self._open_diatomic_mo_diagram_dialog)
        self.symbols_toolbar.ligand_field_diagram_requested.connect(self._open_ligand_field_diagram_dialog)
        self.symbols_toolbar.electronic_diagram_preset_requested.connect(
            self._insert_semantic_preset
        )

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
        self._apply_theme()
        # Chequeo diferido para no bloquear arranque ni UI.
        QTimer.singleShot(1200, self._maybe_check_updates_startup)

    def _create_actions(self) -> None:
        """Inicializa acciones y grupos de UI de alto nivel."""
        create_file_actions(self)
        create_edit_actions(self)
        create_project_actions(self)
        self._ui_builder.create_local_actions(self)
        create_update_actions(self)

    def _apply_theme(self) -> None:
        """Aplica el tema actual a la ventana y a las barras con estilo propio."""
        set_icon_theme(self.current_theme)
        self.setStyleSheet(get_main_stylesheet(self.current_theme))
        toolbar_stylesheet = get_tool_palette_stylesheet(self.current_theme)
        for toolbar in (getattr(self, "toolbar", None), getattr(self, "symbols_toolbar", None)):
            if toolbar is not None:
                toolbar.setStyleSheet(toolbar_stylesheet)
        if hasattr(self, "_ui_builder"):
            self._ui_builder.refresh_main_toolbar_icons(self)
        for widget_name in ("toolbar", "symbols_toolbar", "text_toolbar"):
            widget = getattr(self, widget_name, None)
            if widget is not None and hasattr(widget, "refresh_icons"):
                widget.refresh_icons()

    def toggle_theme(self, checked: bool | None = None) -> None:
        """Cambia entre modo claro y oscuro."""
        if checked is None:
            self.current_theme = "dark" if self.current_theme == "light" else "light"
        else:
            self.current_theme = "dark" if checked else "light"
        self._apply_theme()

    def _create_menu_bar(self) -> None:
        """Construye menús principales delegando el wiring repetitivo."""
        self._ui_builder.build_menu_bar(self)
    
    # -------------------------------------------------------------------------
    # Main Toolbar
    # -------------------------------------------------------------------------
    def _create_main_toolbar(self) -> None:
        """Construye la barra superior principal."""
        self._ui_builder.build_main_toolbar(self)

    def _set_main_toolbar_aux_visible(self, visible: bool) -> None:
        """Muestra/oculta copiar, pegar y zoom +/- en la barra superior."""
        self._ui_builder.set_main_toolbar_aux_visible(self, bool(visible))

    def _on_toggle_main_toolbar_aux(self, checked: bool) -> None:
        """Actualiza visibilidad de atajos auxiliares en barra superior."""
        self._set_main_toolbar_aux_visible(bool(checked))

    def _create_document_tab(self, make_current: bool = False) -> ChemusonCanvas:
        """Crea una nueva pestaña con su propio canvas independiente."""
        return self._tab_manager.create_document_tab(
            make_current=make_current,
            prepare_canvas=lambda canvas: (
                self._apply_default_numbering_to_canvas(canvas),
                self._apply_default_naming_to_canvas(canvas),
            ),
            clean_state_changed=self._on_canvas_clean_state_changed,
        )

    def _apply_default_numbering_to_canvas(self, canvas: ChemusonCanvas) -> None:
        """Aplica preferencias globales de numeración a un nuevo documento."""
        self._view_controller.apply_default_numbering_to_canvas(self, canvas)

    def _apply_default_naming_to_canvas(self, canvas: ChemusonCanvas) -> None:
        """Aplica preferencias globales de nomenclatura al documento."""
        canvas.name_advanced_enabled = bool(self._name_advanced_default)
        canvas.name_rdkit_isolated = bool(self._name_rdkit_isolated_default)

    def _set_canvas_file_path(self, canvas: ChemusonCanvas, filepath: Optional[str]) -> None:
        """Asigna ruta de archivo a una pestaña y actualiza su título."""
        self._document_controller.set_canvas_file_path(
            self._document_tabs_context(),
            canvas,
            filepath,
        )

    def _on_canvas_clean_state_changed(
        self,
        canvas: ChemusonCanvas,
        clean: Optional[bool] = None,
    ) -> None:
        """Actualiza el título de pestaña cuando cambia estado clean/dirty."""
        self._tab_manager.on_canvas_clean_state_changed(canvas, clean)

    def _tabs_alive(self) -> bool:
        """Indica si el widget de pestañas sigue disponible."""
        return self._tab_manager.tabs_alive()

    def _canvas_from_tab_index(self, index: int) -> Optional[ChemusonCanvas]:
        """Obtiene el canvas asociado al índice de pestaña."""
        return self._tab_manager.canvas_from_tab_index(index)

    def _update_tab_title(self, canvas: ChemusonCanvas) -> None:
        """Actualiza título y tooltip de la pestaña asociada al canvas."""
        self._document_controller.update_tab_title(self._document_tabs_context(), canvas)

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
        had_single_tab = self.tabs.count() == 1
        closed = self._tab_manager.close_canvas_tab(
            canvas,
            confirm_discard=self._confirm_discard_changes,
            before_discard=self._before_canvas_discard,
            create_replacement=lambda: self._create_document_tab(make_current=True),
            activate_canvas=self._set_active_canvas,
        )
        if closed and had_single_tab and self.tabs.count() == 1:
            self.statusBar().showMessage("Se creó un documento nuevo.")
        return closed

    def _before_canvas_discard(self, canvas: ChemusonCanvas) -> None:
        """Desconecta el canvas activo antes de retirarlo del tab widget."""
        if self._active_canvas_connected is canvas:
            self._disconnect_canvas_signals(canvas)
            self._active_canvas_connected = None

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
        self._current_file_path = self._tab_manager.file_path_for(canvas)
        if self._active_canvas_connected is not canvas:
            self._connect_canvas_signals(canvas)
            self._active_canvas_connected = canvas

        self.action_undo.setEnabled(self.canvas.undo_stack.canUndo())
        self.action_redo.setEnabled(self.canvas.undo_stack.canRedo())
        self._sync_clipboard_actions()
        self._sync_semantic_diagram_actions()
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
        self._schedule_chemical_properties_update()
        self.text_toolbar.set_opacity_percent(self.canvas.current_opacity_percent())
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
        self._schedule_chemical_properties_update()
        self._sync_fragment_pivot_actions()
        self._update_tab_title(self.canvas)
        self.text_toolbar.set_opacity_percent(self.canvas.current_opacity_percent())

    def _update_iupac_name_indicator(self) -> None:
        """Refresca el indicador de nombre IUPAC en la barra de estado."""
        if not hasattr(self, "_iupac_name_label"):
            return
        try:
            name = self.canvas.current_iupac_name()
        except Exception:
            name = "N/D"
        self._iupac_name_label.setText(f"Nombre IUPAC: {name or 'N/D'}")

    def _schedule_chemical_properties_update(self) -> None:
        """Programa recálculo ligero del dock de propiedades químicas."""
        timer = getattr(self, "_properties_update_timer", None)
        if timer is None:
            return
        timer.start()

    def _refresh_chemical_properties_dock(self) -> None:
        """Refresca propiedades calculadas del documento activo."""
        dock = getattr(self, "chemical_properties_dock", None)
        if dock is None:
            return
        dock.update_properties(self._chemical_properties_rows())

    def _chemical_properties_rows(self) -> list[tuple[str, str]]:
        """Calcula propiedades químicas rápidas sin bloquear con RDKit."""
        canvas = getattr(self, "canvas", None)
        if canvas is None or not getattr(canvas.model, "atoms", None):
            return []
        graph = canvas.model
        try:
            counts = canvas._analysis_atom_counts(graph)
            formula = canvas._analysis_formula(counts)
            exact_mass = canvas._analysis_exact_mass(counts)
            molecular_weight = canvas._analysis_molecular_weight(counts)
            issues = graph.validate_detailed()
        except Exception:
            return [("Estado", "N/D")]

        rows = [
            ("Fórmula", formula or "N/D"),
            (
                "Masa exacta",
                f"{exact_mass:.4f}" if exact_mass is not None else "N/D",
            ),
            (
                "Peso molecular",
                f"{molecular_weight:.4f}" if molecular_weight is not None else "N/D",
            ),
            ("Carga total", str(graph.total_formal_charge())),
            ("Átomos", str(len(graph.atoms))),
            ("Enlaces", str(len(graph.bonds))),
            ("Errores de valencia", str(len(issues))),
        ]
        elemental_line = canvas._analysis_elemental_line(counts, molecular_weight)
        if elemental_line:
            rows.append(("Análisis elemental", elemental_line.replace("Elemental Analysis: ", "")))
        return rows

    def _on_undo(self) -> None:
        """Deshace en la pestaña activa."""
        self.canvas.undo_stack.undo()

    def _on_redo(self) -> None:
        """Rehace en la pestaña activa."""
        self.canvas.undo_stack.redo()

    def _on_copy(self) -> None:
        """Copia desde la pestaña activa."""
        self.canvas.copy_to_clipboard()

    def _on_cut(self) -> None:
        """Corta desde la pestaña activa."""
        self.canvas.cut_to_clipboard()

    def _on_paste(self) -> None:
        """Pega en la pestaña activa."""
        self.canvas.paste_from_clipboard()

    def _on_duplicate(self) -> None:
        """Duplica la seleccion en la pestaña activa."""
        self.canvas.duplicate_selection()

    def _on_delete(self) -> None:
        """Elimina la seleccion en la pestaña activa."""
        self.canvas.delete_selection()

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
        self._text_format_controller.on_text_format_changed(
            self,
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
        self._text_format_controller.on_text_color_changed(self, color)

    def _on_text_alignment_changed(self, alignment: Qt.AlignmentFlag) -> None:
        """Propaga cambio de alineación al canvas activo."""
        self._text_format_controller.on_text_alignment_changed(self, alignment)

    def _on_opacity_changed(self, value: int) -> None:
        """Propaga cambio de opacidad al canvas activo."""
        if getattr(self, "_external_text_editor", None) is not None:
            return
        self.canvas.apply_opacity_percent(float(value))

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
        canvas.state.active_energy_diagram_kind = self.symbols_toolbar.current_energy_diagram_kind()
        canvas.state.active_orbital_kind = self.symbols_toolbar.current_orbital_kind()
        canvas.set_current_tool(self._current_tool_id)

    def _visible_canvas_center_scene_pos(self) -> QPointF:
        """Devuelve el centro visible actual del canvas activo."""
        return self.canvas.mapToScene(self.canvas.viewport().rect().center())

    def _set_external_text_editor(self, editor: QTextEdit | None) -> None:
        """Registra un editor de texto enriquecido temporal para el toolbar superior."""
        self._text_format_controller.set_external_text_editor(self, editor)

    def _on_external_text_copy_available(self, _available: bool) -> None:
        self._text_format_controller.on_external_text_copy_available(self, _available)

    def _remember_external_text_cursor(self) -> None:
        """Guarda anchor/posición del editor temporal para no perder selección al cambiar foco."""
        self._text_format_controller.remember_external_text_cursor(self)

    def _external_text_cursor_for_formatting(self, editor: QTextEdit) -> QTextCursor:
        """Obtiene el cursor adecuado para formatear, priorizando selección."""
        return self._text_format_controller.external_text_cursor_for_formatting(self, editor)

    def _sync_text_toolbar_from_external_editor(self) -> None:
        """Sincroniza el toolbar desde un QTextEdit temporal."""
        self._text_format_controller.sync_text_toolbar_from_external_editor(self)

    def _apply_text_format_to_external_editor(
        self,
        family: str,
        size: int,
        bold: bool,
        italic: bool,
        underline: bool,
        sub: bool,
        sup: bool,
        property_name: str,
    ) -> bool:
        """Aplica cambios del toolbar a un editor temporal si existe."""
        return self._text_format_controller.apply_text_format_to_external_editor(
            self, family, size, bold, italic, underline, sub, sup, property_name
        )

    def _apply_text_color_to_external_editor(self, color: QColor) -> bool:
        return self._text_format_controller.apply_text_color_to_external_editor(self, color)

    def _apply_text_alignment_to_external_editor(self, alignment: Qt.AlignmentFlag) -> bool:
        return self._text_format_controller.apply_text_alignment_to_external_editor(self, alignment)

    def _rich_text_editor_value(self, editor: QTextEdit) -> str:
        """Compatibilidad para el flujo de diálogo enriquecido."""
        return rich_text_editor_value(editor)

    def _open_rich_text_value_dialog(
        self,
        *,
        title: str,
        label: str,
        initial_text: str = "",
    ) -> tuple[str, bool]:
        """Abre un editor enriquecido no modal conectado al toolbar de texto."""
        return open_rich_text_value_dialog(
            self,
            title=title,
            label=label,
            initial_text=initial_text,
            set_external_text_editor=self._set_external_text_editor,
            sync_text_toolbar=self._sync_text_toolbar_from_external_editor,
        )

    def _insert_semantic_preset(self, preset_name: str) -> None:
        self._semantic_diagram_workflow.insert_preset(
            self._semantic_diagram_workflow_context(),
            preset_name,
        )

    def _open_atomic_diagram_dialog(
        self,
        existing_item=None,
    ) -> None:
        self._semantic_diagram_workflow.open_atomic_diagram_dialog(
            self._semantic_diagram_workflow_context(),
            existing_item=existing_item,
        )

    def _open_atomic_species_diagram_dialog(
        self,
        existing_item=None,
    ) -> None:
        self._semantic_diagram_workflow.open_atomic_species_diagram_dialog(
            self._semantic_diagram_workflow_context(),
            existing_item=existing_item,
        )

    def _open_diatomic_mo_diagram_dialog(
        self,
        existing_item=None,
    ) -> None:
        self._semantic_diagram_workflow.open_diatomic_mo_diagram_dialog(
            self._semantic_diagram_workflow_context(),
            existing_item=existing_item,
        )

    def _open_ligand_field_diagram_dialog(
        self,
        existing_item=None,
    ) -> None:
        self._semantic_diagram_workflow.open_ligand_field_diagram_dialog(
            self._semantic_diagram_workflow_context(),
            existing_item=existing_item,
        )

    def _open_elemental_analysis_dialog(self) -> None:
        """Open the interactive elemental analysis tool."""
        existing = getattr(self, "_elemental_analysis_dialog", None)
        if existing is not None and existing.isVisible():
            existing.raise_()
            existing.activateWindow()
            return
        initial_formula = ""
        try:
            counts = self.canvas._analysis_atom_counts(self.canvas.model)
            initial_formula = self.canvas._analysis_formula(counts)
        except Exception:
            initial_formula = ""
        dialog = ElementalAnalysisDialog(initial_formula=initial_formula, parent=self)
        dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self._elemental_analysis_dialog = dialog
        dialog.destroyed.connect(lambda *_args: setattr(self, "_elemental_analysis_dialog", None))
        dialog.show()

    def _edit_selected_semantic_diagram(self) -> None:
        self._semantic_diagram_workflow.edit_selected_semantic_diagram(
            self._semantic_diagram_workflow_context()
        )

    def _connect_undo_redo(self) -> None:
        """Conecta acciones globales a la pestaña activa."""
        self.action_undo.triggered.connect(self._on_undo)
        self.action_redo.triggered.connect(self._on_redo)
        self.action_copy.triggered.connect(self._on_copy)
        self.action_cut.triggered.connect(self._on_cut)
        self.action_paste.triggered.connect(self._on_paste)
        self.action_duplicate.triggered.connect(self._on_duplicate)
        self.action_delete.triggered.connect(self._on_delete)
        self.action_edit_electronic_diagram.triggered.connect(self._edit_selected_semantic_diagram)
        self.action_undo.setEnabled(False)
        self.action_redo.setEnabled(False)
        self.action_copy.setEnabled(False)
        self.action_cut.setEnabled(False)
        self.action_paste.setEnabled(
            bool(getattr(self.canvas, "can_paste_from_clipboard", lambda: False)())
        )
        self.action_duplicate.setEnabled(False)
        self.action_delete.setEnabled(False)
        self.action_edit_electronic_diagram.setEnabled(False)
        try:
            QApplication.clipboard().dataChanged.connect(self._sync_clipboard_actions)
        except Exception:
            pass

    def _load_recent_files(self) -> list[str]:
        """Método auxiliar para  load recent files."""
        return self._document_controller.load_recent_files(self._recent_files_context([]))

    def _save_recent_files(self) -> None:
        """Método auxiliar para  save recent files."""
        self._document_controller.save_recent_files(self._recent_files_context())

    def _load_update_preferences(self) -> UpdateSettings:
        """Compatibilidad con el estado centralizado en UpdateController."""
        return self._update_controller.settings

    def _save_update_preferences(self) -> None:
        """Persiste preferencias de actualización en QSettings."""
        self._update_controller.save_preferences()

    def _update_settings_payload(self) -> dict:
        """Devuelve preferencias de update para precargar en diálogo."""
        return self._update_controller.settings_payload()

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
        self._view_controller.load_numbering_preferences(self)

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
        self._view_controller.save_numbering_preferences(self)

    def _sync_numbering_actions(self) -> None:
        """Sincroniza acciones de menú de numeración con el estado del canvas."""
        self._view_controller.sync_numbering_actions(self)

    def _maybe_check_updates_startup(self) -> None:
        """Chequea updates en inicio respetando política y frecuencia."""
        self._check_for_updates(force=False, interactive=False)

    def _on_check_updates_now(self) -> None:
        """Lanza chequeo manual de actualizaciones desde el menú Ayuda."""
        self._check_for_updates(force=True, interactive=True)

    def _update_controller_context(self) -> UpdateControllerContext:
        """Dependencias mínimas que UpdateController necesita desde la ventana."""
        return UpdateControllerContext(
            parent=self,
            show_status=lambda message, timeout=0: self.statusBar().showMessage(message, timeout),
            close_window=self.close,
        )

    def _check_for_updates(self, force: bool, interactive: bool) -> None:
        """Ejecuta chequeo de updates delegando en UpdateController."""
        self._update_controller.check_for_updates(
            self._update_controller_context(),
            force=force,
            interactive=interactive,
        )







    def _apply_pending_windows_update_on_exit(self) -> bool:
        """Aplica instalador pendiente usando UpdateController."""
        return self._update_controller.apply_pending_windows_update_on_exit(
            self._update_controller_context()
        )

    def _apply_pending_portable_update_on_exit(self) -> bool:
        """Aplica reemplazo portable pendiente usando UpdateController."""
        return self._update_controller.apply_pending_portable_update_on_exit(
            self._update_controller_context()
        )

    def _add_recent_file(self, filepath: str) -> None:
        """Método auxiliar para  add recent file."""
        self._document_controller.add_recent_file(self._recent_files_context(), filepath)

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
        """Abre una entrada reciente usando el workflow desacoplado."""
        self._file_controller.open_recent_file(self._file_workflow_context(), filepath)
    
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
        self._file_controller.load_file_into_canvas(filepath, canvas)

    @staticmethod
    def _read_autosave_metadata(filepath: str) -> Optional[dict]:
        """Lee metadatos mínimos de un autosave para la tabla de recuperación."""
        return RecoveryController.read_autosave_metadata(filepath)

    @staticmethod
    def _list_autosave_entries(directory: str) -> list[dict]:
        """Lista autosaves válidos de un directorio ordenados por fecha desc."""
        return RecoveryController.list_autosave_entries(directory)

    @staticmethod
    def _archive_autosave(path: str, autosave_dir: str) -> str:
        """Mueve un autosave a la carpeta `old` y devuelve su nueva ruta."""
        return RecoveryController.archive_autosave(path, autosave_dir)

    def _open_autosave_document(
        self,
        autosave_path: str,
        original_path: Optional[str] = None,
    ) -> bool:
        """Abre un autosave como pestaña nueva."""
        return self._recovery_controller.open_autosave_document(self, autosave_path, original_path)

    def _on_open_recovery_center(self) -> None:
        """Abre manualmente el centro de recuperación y archivos recientes."""
        self._recovery_controller.on_open_recovery_center(self)

    def _show_recovery_center(self, show_only_if_pending: bool = False) -> None:
        """Muestra diálogo con pendientes, recientes y recuperados."""
        self._recovery_controller.show_recovery_center(self, show_only_if_pending=show_only_if_pending)

    @staticmethod
    def check_autosaves(window: "ChemusonWindow") -> None:
        """Busca autosaves pendientes y abre el centro de recuperación al iniciar."""
        window._recovery_controller.check_autosaves(window)

    def _open_file_path(self, filepath: str) -> None:
        """Abre un archivo físico usando el workflow desacoplado."""
        self._file_controller.open_file_path(self._file_workflow_context(), filepath)

    def _on_file_save(self) -> None:
        """Guarda el documento activo usando el workflow desacoplado."""
        self._file_controller.save_file(self._file_workflow_context())
    
    def _on_export(self, format: str) -> None:
        """Export the canvas in the specified format."""
        self._export_controller.export(self, format)

    def _confirm_discard_changes(self, canvas: ChemusonCanvas) -> bool:
        """Confirma descarte del canvas especificado."""
        return self._document_controller.confirm_discard_changes(
            self._document_discard_context(),
            canvas,
        )

    def _file_workflow_context(self) -> FileWorkflowContext:
        """Construye el contrato mínimo requerido por FileController."""
        return FileWorkflowContext(
            parent=self,
            canvas=self.canvas,
            tabs=self.tabs,
            tab_manager=self._tab_manager,
            create_document_tab=self._create_document_tab,
            apply_toolbar_defaults=self._apply_toolbar_defaults_to_canvas,
            set_active_canvas=self._set_active_canvas,
            before_canvas_discard=self._before_canvas_discard,
            add_recent_file=self._add_recent_file,
            refresh_recent_menu=self._update_recent_menu,
            update_total_charge_indicator=self._update_total_charge_indicator,
            show_status=lambda message: self.statusBar().showMessage(message),
        )

    def _recent_files_context(
        self,
        recent_files: Optional[list[str]] = None,
    ) -> RecentFilesContext:
        """Construye el contexto mínimo requerido por DocumentController para recientes."""
        return RecentFilesContext(
            settings=self._settings,
            recent_files=self._recent_files if recent_files is None else recent_files,
            persist_recent_files=self._save_recent_files,
            refresh_recent_menu=self._update_recent_menu,
        )

    def _document_tabs_context(self) -> DocumentTabsContext:
        """Contrato explícito para operaciones de pestañas/rutas."""
        return DocumentTabsContext(
            tab_manager=self._tab_manager,
            current_canvas=self.canvas,
            current_file_path_setter=self._set_current_file_path,
        )

    def _set_current_file_path(self, filepath: Optional[str]) -> None:
        self._current_file_path = filepath

    def _document_discard_context(self) -> DocumentDiscardContext:
        """Construye el contexto mínimo requerido para confirmar descarte."""
        return DocumentDiscardContext(
            parent=self,
            canvas=self.canvas,
            tab_manager=self._tab_manager,
            save_canvas=self._on_file_save,
            activate_canvas=self._set_active_canvas,
        )

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
        """Copia la estructura en el formato solicitado usando FileController."""
        self._file_controller.copy_as(self._file_workflow_context(), format)
    
    def _on_preferences(self) -> None:
        """Open preferences dialog."""
        dialog = PreferencesDialog(
            self.canvas.state,
            self.canvas.drawing_style,
            current_theme=self.current_theme,
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
        theme = str(prefs.get("theme", self.current_theme))
        if bond_caps:
            self._apply_bond_caps(bond_caps)
        if theme != self.current_theme:
            self.current_theme = "dark" if theme == "dark" else "light"
            self._apply_theme()

        self.action_show_carbons.setChecked(self.canvas.state.show_implicit_carbons)
        self.action_show_hydrogens.setChecked(self.canvas.state.show_implicit_hydrogens)
        self.action_aromatic_circles.setChecked(self.canvas.state.use_aromatic_circles)

        self.canvas.refresh_atom_visibility()
        self.canvas.refresh_aromatic_circles()
        self._sync_label_menu_state()

        if bond_caps:
            self.appearance_dock.set_bond_caps(bond_caps)

        self._update_controller.apply_preferences(prefs)
        self._update_settings = self._update_controller.settings

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
        """Aplica preferencias visuales del documento."""
        self._view_controller.apply_appearance_settings(self, prefs)

    def _apply_bond_caps(self, bond_caps: str) -> None:
        """Compatibilidad para aplicar estilo de terminación de enlace."""
        self._view_controller.apply_bond_caps(self, bond_caps)
    
    # -------------------------------------------------------------------------
    # View Menu Handlers
    # -------------------------------------------------------------------------
    def _on_toggle_numbering(self, checked: bool) -> None:
        self._view_controller.on_toggle_numbering(self, checked)

    def _on_set_numbering_mode(self, mode: str) -> None:
        self._view_controller.on_set_numbering_mode(self, mode)

    def _on_recalculate_numbering(self) -> None:
        self._view_controller.on_recalculate_numbering(self)

    def _on_toggle_numbering_export(self, checked: bool) -> None:
        self._view_controller.on_toggle_numbering_export(self, checked)

    def _on_toggle_carbons(self, checked: bool) -> None:
        self._view_controller.on_toggle_carbons(self, checked)
    
    def _on_toggle_hydrogens(self, checked: bool) -> None:
        self._view_controller.on_toggle_hydrogens(self, checked)
    
    def _on_toggle_aromatic_circles(self, checked: bool) -> None:
        self._view_controller.on_toggle_aromatic_circles(self, checked)

    def _on_toggle_rules(self, checked: bool) -> None:
        self._view_controller.on_toggle_rules(self, checked)

    def _on_toggle_crosshair(self, checked: bool) -> None:
        self._view_controller.on_toggle_crosshair(self, checked)
    
    def _on_zoom_in(self) -> None:
        self._view_controller.on_zoom_in(self)
    
    def _on_zoom_out(self) -> None:
        self._view_controller.on_zoom_out(self)
    
    def _on_zoom_reset(self) -> None:
        self._view_controller.on_zoom_reset(self)

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
        self._text_format_controller.sync_label_menu_state(self)

    def _on_label_font(self) -> None:
        self._text_format_controller.on_label_font(self)

    def _on_label_font_size_dialog(self) -> None:
        self._text_format_controller.on_label_font_size_dialog(self)

    def _on_label_bold(self, checked: bool) -> None:
        self._text_format_controller.on_label_bold(self, checked)

    def _on_label_italic(self, checked: bool) -> None:
        self._text_format_controller.on_label_italic(self, checked)

    def _on_label_underline(self, checked: bool) -> None:
        self._text_format_controller.on_label_underline(self, checked)

    def _on_label_font_size(self, delta: float) -> None:
        self._text_format_controller.on_label_font_size(self, delta)

    def _set_canvas_size(self, width: int, height: int) -> None:
        self._view_controller.set_canvas_size(self, width, height)

    def _on_canvas_custom_size(self) -> None:
        self._view_controller.open_canvas_custom_size_dialog(self)

    def _on_label_subscript(self) -> None:
        self._text_format_controller.on_label_subscript(self)

    def _on_label_superscript(self) -> None:
        self._text_format_controller.on_label_superscript(self)

    def _on_label_color_mode(self, use_element_colors: bool) -> None:
        self._text_format_controller.on_label_color_mode(self, use_element_colors)
    
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
        """Maneja clean 2d."""
        self._run_clean_2d(step_ratio=0.35, fallback_iterations=30, status_suffix="(paso)")

    def _on_clean_2d_full(self) -> None:
        """Maneja clean 2d full.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        self._run_clean_2d(step_ratio=1.0, fallback_iterations=200, status_suffix="(1 paso)")

    def _on_validate_structure(self) -> None:
        """Ejecuta validación química detallada sobre el documento activo."""
        errors = self.canvas.validate_structure()
        if errors:
            self.statusBar().showMessage(
                f"Validación: {len(errors)} error(es) de valencia. Use F8 para navegar.",
                9000,
            )
        else:
            self.statusBar().showMessage("Validación: sin errores de valencia.", 7000)

    def _on_navigate_validation_issue(self, step: int) -> None:
        """Navega al siguiente/anterior diagnóstico de valencia."""
        self.canvas.navigate_validation_issue(step)

    def _run_clean_2d(self, step_ratio: float, fallback_iterations: int, status_suffix: str) -> None:
        """Clean 2D coordinates using RDKit or fallback."""
        self._clean2d_controller.run_clean_2d(self, step_ratio, fallback_iterations, status_suffix)

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

    def _refresh_template_views(self) -> None:
        """Sincroniza menú y dock con la biblioteca de plantillas."""
        self._template_browser.refresh_template_views(self._template_browser_context())

    def _refresh_templates_menu(self) -> None:
        """Construye menú dinámico de plantillas por categoría."""
        self._template_browser.refresh_templates_menu(self._template_browser_context())

    def _start_template_insert_by_id(self, template_id: str, *, place_now: bool = False) -> None:
        """Carga plantilla desde biblioteca e inicia inserción."""
        self._template_controller.start_template_insert_by_id(
            self._template_controller_context(),
            template_id,
            place_now=place_now,
        )

    def _on_template_selected_from_gallery(self, payload: dict) -> None:
        """Maneja selección de plantilla desde el dock lateral."""
        template_id = str(payload.get("id", "")).strip()
        if not template_id:
            return
        self._start_template_insert_by_id(template_id)

    def _prompt_template_destination(self) -> Optional[tuple[str, str]]:
        """Solicita nombre y categoría para guardar plantilla."""
        return self._template_controller.prompt_template_destination(
            self._template_controller_context()
        )

    def _on_save_template(self) -> None:
        """Guarda selección (o molécula completa) como plantilla reusable."""
        self._template_controller.on_save_template(self._template_controller_context())

    def _on_new_template_category(self) -> None:
        """Crea categoría personalizada de plantillas."""
        self._template_controller.on_new_template_category(self._template_controller_context())

    def _on_rename_template_category(self, old_name: str) -> None:
        """Renombra una categoría existente."""
        self._template_controller.on_rename_template_category(
            self._template_controller_context(),
            old_name,
        )

    def _on_delete_template_category(self, name: str) -> None:
        """Elimina una categoría y conserva sus plantillas en categoría de respaldo."""
        self._template_controller.on_delete_template_category(
            self._template_controller_context(),
            name,
        )

    def _on_rename_template(self, template_id: str) -> None:
        """Renombra plantilla individual."""
        self._template_controller.on_rename_template(
            self._template_controller_context(),
            template_id,
        )

    def _on_delete_template(self, template_id: str) -> None:
        """Elimina plantilla de la biblioteca."""
        self._template_controller.on_delete_template(
            self._template_controller_context(),
            template_id,
        )

    def _on_import_template_library(self) -> None:
        """Importa bibliotecas de plantillas desde JSON."""
        self._template_controller.on_import_template_library(
            self._template_controller_context()
        )

    def _on_export_template_library(self) -> None:
        """Exporta biblioteca de plantillas a JSON."""
        self._template_controller.on_export_template_library(
            self._template_controller_context()
        )

    def _insert_template_from_file(self, filepath: str) -> None:
        """Carga e inserta una plantilla legada en formato MOL."""
        self._template_browser.insert_template_from_file(
            self._template_browser_context(),
            filepath,
        )

    def _on_insert_linear_chain(self) -> None:
        """Maneja insert linear chain."""
        self._template_controller.on_insert_linear_chain(
            self._template_controller_context()
        )

    def _on_import_smiles(self) -> None:
        """Import a molecule from a SMILES string."""
        self._template_controller.on_import_smiles(self._template_controller_context())

    def _on_export_smiles(self) -> None:
        """Export the current molecule as SMILES."""
        self._template_controller.on_export_smiles(self._template_controller_context())

    def _template_controller_context(self) -> TemplateControllerContext:
        """Construye el contexto mínimo requerido por TemplateController."""
        return TemplateControllerContext(
            parent=self,
            canvas=self.canvas,
            template_library=self.template_library,
            show_status=self.statusBar().showMessage,
            refresh_template_views=self._refresh_template_views,
            insert_template=self._insert_template,
        )

    def _template_browser_context(self) -> TemplateBrowserContext:
        """Dependencias de presentación/migración de plantillas."""
        return TemplateBrowserContext(
            parent=self,
            template_library=self.template_library,
            templates_menu=self.templates_menu,
            templates_dock=self.templates_dock,
            action_save_template=self.action_save_template,
            action_template_new_category=self.action_template_new_category,
            action_template_import_library=self.action_template_import_library,
            action_template_export_library=self.action_template_export_library,
            preview_cache=self._template_icon_cache,
            show_status=self.statusBar().showMessage,
            start_template_insert_by_id=self._start_template_insert_by_id,
            insert_template=self._insert_template,
            default_category_user=DEFAULT_CATEGORY_USER,
        )

    def _semantic_diagram_workflow_context(self) -> SemanticDiagramWorkflowContext:
        """Dependencias mínimas para editar/insertar diagramas semánticos."""
        return SemanticDiagramWorkflowContext(
            parent=self,
            canvas=self.canvas,
            visible_center=self._visible_canvas_center_scene_pos,
        )

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
        """Actualiza la barra de estado con la herramienta activa."""
        name = self._ui_builder.tool_status_label(self, tool_id)
        self.statusBar().showMessage(f"Herramienta: {name}")

    def _on_selection_changed(self, num_atoms: int, num_bonds: int, num_text: int, details: dict):
        """Handle selection change to update UI components."""
        self.inspector_dock.update_selection(num_atoms, num_bonds, num_text, details)
        self._update_total_charge_indicator()
        self._sync_fragment_pivot_actions()
        self._sync_clipboard_actions()
        self._sync_semantic_diagram_actions()
        if getattr(self, "_external_text_editor", None) is not None:
            return
        self.text_toolbar.set_opacity_percent(self.canvas.current_opacity_percent())
        
        # Sync Text Toolbar if a single text item is selected
        if num_text == 1 and details.get("type") == "text":
            font = details.get("font")
            settings = {
                "color": details.get("color"),
                "sub": details.get("sub"),
                "sup": details.get("sup")
            }
            self.text_toolbar.set_state(font, settings)

    def _sync_clipboard_actions(self) -> None:
        """Sincroniza Copiar/Pegar con la selección y el portapapeles activos."""
        if not hasattr(self, "action_copy") or not hasattr(self, "action_paste"):
            return
        canvas = getattr(self, "canvas", None)
        if canvas is None:
            self.action_copy.setEnabled(False)
            self.action_paste.setEnabled(False)
            return
        self.action_copy.setEnabled(bool(canvas.has_copyable_selection()))
        self.action_cut.setEnabled(bool(canvas.has_copyable_selection()))
        self.action_paste.setEnabled(bool(canvas.can_paste_from_clipboard()))
        self.action_duplicate.setEnabled(bool(canvas.has_copyable_selection()))
        self.action_delete.setEnabled(bool(canvas.has_copyable_selection()))

    def _sync_semantic_diagram_actions(self) -> None:
        """Sincroniza la edición completa de diagramas semánticos con la selección."""
        if not hasattr(self, "action_edit_electronic_diagram"):
            return
        canvas = getattr(self, "canvas", None)
        if canvas is None:
            self.action_edit_electronic_diagram.setEnabled(False)
            return
        self.action_edit_electronic_diagram.setEnabled(
            canvas.selected_semantic_diagram_item() is not None
        )

    def _update_total_charge_indicator(self) -> None:
        self._view_controller.update_total_charge_indicator(self)


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
