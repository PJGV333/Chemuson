"""
Ventana principal de Chemuson.

Compone menús, barras de herramientas, docks y el lienzo central.
"""
from PyQt6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFontDialog,
    QInputDialog,
    QMainWindow,
    QMenu,
    QMenuBar,
    QFileDialog,
    QMessageBox,
    QToolBar,
    QFormLayout,
    QDoubleSpinBox,
    QComboBox,
    QVBoxLayout,
    QLabel,
)
from PyQt6.QtCore import Qt, QSize, QSettings, QEvent
from PyQt6.QtGui import QAction, QActionGroup, QKeySequence, QIcon, QPainter, QPixmap, QPen, QColor
from PyQt6.QtPrintSupport import QPrinter
from typing import Optional
from dataclasses import replace

from gui.canvas import ChemusonCanvas
from gui.periodic_table import PeriodicTableDialog
from gui.toolbar import ChemusonToolbar, SymbolPaletteToolbar
from gui.styles import MAIN_STYLESHEET, TOOL_PALETTE_STYLESHEET
from gui.icons import draw_generic_icon
from gui.docks import PlantillasDock, InspectorDock, AppearanceDock
from gui.dialogs import PreferencesDialog, QuickStartDialog, StyleDialog
from gui.text_toolbar import TextFormatToolbar
from gui.commands import ChangeAtomCommand
from gui.template_library import TemplateLibrary, DEFAULT_CATEGORY_USER
from gui.templates import (
    build_linear_chain_template,
)
import os
from chemio.rdkit_io import molfile_to_molgraph, molgraph_to_molfile
from chemio.persistence import PersistenceManager
from core.model import MolGraph


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
        self.setWindowTitle("Chemuson - Editor Molecular Libre")
        self.resize(1200, 900)
        
        # Apply main stylesheet
        self.setStyleSheet(MAIN_STYLESHEET)

        # === CORE COMPONENTS ===
        self._create_actions()
        self._current_file_path: Optional[str] = None
        self._settings = QSettings("Chemuson", "Chemuson")
        self._recent_files = self._load_recent_files()
        self.template_library = TemplateLibrary()
        self._template_icon_cache: dict[str, QIcon] = {}
        
        # === CENTRAL CANVAS ===
        self.canvas = ChemusonCanvas()
        self.setCentralWidget(self.canvas)
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
        self.toolbar.set_text_menu(
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

        # === RIGHT SYMBOLS TOOLBAR ===
        self.symbols_toolbar = SymbolPaletteToolbar(action_group=self.toolbar.action_group)
        self.symbols_toolbar.setStyleSheet(TOOL_PALETTE_STYLESHEET)
        self.addToolBar(Qt.ToolBarArea.RightToolBarArea, self.symbols_toolbar)

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
        self.text_toolbar.format_changed.connect(self.canvas.update_text_format)
        self.text_toolbar.color_changed.connect(self.canvas.update_text_color)
        self.text_toolbar.alignment_changed.connect(self.canvas.update_text_alignment)
        
        # === SIGNAL CONNECTIONS ===
        self._connect_undo_redo()
        
        # Connect toolbar to canvas
        self.toolbar.tool_changed.connect(self.canvas.set_current_tool)
        self.toolbar.bond_palette_changed.connect(self._handle_bond_palette)
        self.toolbar.ring_palette_changed.connect(self._handle_ring_palette)
        self.toolbar.element_palette_changed.connect(self._handle_element_palette)
        self.toolbar.periodic_table_requested.connect(self._show_periodic_table)
        
        # Connect selection updates
        self.canvas.selection_changed.connect(self._on_selection_changed)

        # Connect symbols dock to canvas
        self.symbols_toolbar.tool_changed.connect(self.canvas.set_current_tool)
        self.symbols_toolbar.tool_changed.connect(self._update_status)

        # Sync defaults selected during toolbar init
        bond_spec = self.toolbar.current_bond_spec()
        self.canvas.state.active_bond_order = bond_spec.get("order", 1)
        self.canvas.state.active_bond_style = bond_spec.get("style", self.canvas.state.active_bond_style)
        self.canvas.state.active_bond_stereo = bond_spec.get("stereo", self.canvas.state.active_bond_stereo)
        self.canvas.state.active_bond_mode = bond_spec.get("mode", "increment")
        self.canvas.state.active_bond_aromatic = bond_spec.get("aromatic", False)

        ring_spec = self.toolbar.current_ring_spec()
        self.canvas.state.active_ring_size = ring_spec.get("size", self.canvas.state.active_ring_size)
        self.canvas.state.active_ring_aromatic = ring_spec.get("aromatic", False)
        self.canvas.state.active_ring_template = ring_spec.get("template")
        self.canvas.state.active_ring_anomeric = ring_spec.get("anomeric")

        self.canvas.state.default_element = self.toolbar.current_element()
        
        # === STATUS BAR ===
        self._total_charge_label = QLabel()
        self.statusBar().addPermanentWidget(self._total_charge_label)
        self._update_total_charge_indicator()
        self._update_status(self.canvas.state.active_tool)
        
        # Update status bar when tool changes
        self.toolbar.tool_changed.connect(self._update_status)
        self.canvas.undo_stack.indexChanged.connect(
            lambda _index: self._update_total_charge_indicator()
        )

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



        self.action_style = QAction("Estilo de dibujo...", self)
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

        # --- Bond Thickness Actions ---
        self.action_bond_thickness_up = QAction("Aumentar grosor de enlace", self)
        self.action_bond_thickness_up.setShortcut(QKeySequence("Ctrl+Shift+Up"))
        self.action_bond_thickness_up.triggered.connect(self._on_bond_thickness_up)

        self.action_bond_thickness_down = QAction("Reducir grosor de enlace", self)
        self.action_bond_thickness_down.setShortcut(QKeySequence("Ctrl+Shift+Down"))
        self.action_bond_thickness_down.triggered.connect(self._on_bond_thickness_down)

        self.action_bond_thickness_reset = QAction("Restablecer grosor de enlace", self)
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

        edit_menu.addSeparator()
        bond_thickness_menu = edit_menu.addMenu("Grosor de enlace")
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
        
        from gui.icons import draw_generic_icon
        self.action_zoom_in.setIcon(draw_generic_icon("zoom_in"))
        self.action_zoom_out.setIcon(draw_generic_icon("zoom_out"))
        self.action_rotate_left.setIcon(draw_generic_icon("rotate_left"))
        self.action_rotate_right.setIcon(draw_generic_icon("rotate_right"))
        self.action_flip_horizontal.setIcon(draw_generic_icon("flip_horizontal"))
        self.action_flip_vertical.setIcon(draw_generic_icon("flip_vertical"))
        self.action_clean_2d.setIcon(QIcon.fromTheme("edit-clear", QIcon()))
        from gui.icons import draw_atom_icon
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
    
    def _connect_undo_redo(self) -> None:
        """Connect undo/redo actions to the canvas undo stack."""
        self.action_undo.triggered.connect(self.canvas.undo_stack.undo)
        self.action_redo.triggered.connect(self.canvas.undo_stack.redo)
        
        self.canvas.undo_stack.canUndoChanged.connect(self.action_undo.setEnabled)
        self.canvas.undo_stack.canRedoChanged.connect(self.action_redo.setEnabled)
        
        # Initial state
        self.action_undo.setEnabled(self.canvas.undo_stack.canUndo())
        self.action_redo.setEnabled(self.canvas.undo_stack.canRedo())
        
        # Connect copy/paste
        self.action_copy.triggered.connect(self.canvas.copy_to_clipboard)
        self.action_paste.triggered.connect(self.canvas.paste_from_clipboard)

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
        mode = str(self._settings.value("numbering/mode", "atoms") or "atoms")
        include_export = self._setting_bool(
            self._settings.value("numbering/include_export", True),
            True,
        )
        self.canvas.set_numbering_mode(mode)
        # Arranque por defecto: numeración desactivada (no persistir "encendido").
        self.canvas.set_numbering_enabled(False)
        self.canvas.set_numbering_include_in_export(include_export)

    def _save_numbering_preferences(self) -> None:
        """Guarda preferencias globales de numeración del usuario."""
        # Se guarda modo/exports, pero no el estado "mostrar numeración" para
        # evitar que inicie activa en futuras aperturas.
        self._settings.remove("numbering/enabled")
        self._settings.setValue("numbering/mode", str(self.canvas.state.numbering_mode))
        self._settings.setValue(
            "numbering/include_export",
            bool(self.canvas.state.numbering_include_in_export),
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
        self.canvas.clear_canvas()
        self._current_file_path = None
        self.canvas.undo_stack.setClean()
        self._update_total_charge_indicator()
        self.statusBar().showMessage("Nuevo documento creado")
    
    def _on_file_open(self) -> None:
        """Open a molecule file (.cmsn or .mol)."""
        filepath, selected_filter = QFileDialog.getOpenFileName(
            self,
            "Abrir archivo",
            "",
            "Archivos de Chemuson (*.cmsn);;Archivos MOL (*.mol *.sdf);;Todos los archivos (*.*)"
        )
        if filepath:
            self._open_file_path(filepath)

    def _open_file_path(self, filepath: str) -> None:
        """Método auxiliar para  open file path.

        Args:
            filepath: Descripción del parámetro.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        try:
            if filepath.lower().endswith(".cmsn"):
                self.canvas.clear_canvas() # Ensure clean slate
                PersistenceManager.load_from_file(filepath, self.canvas)
            else:
                # Fallback to RDKit for .mol/.sdf
                with open(filepath, "r") as f:
                    molfile = f.read()
                graph = molfile_to_molgraph(molfile)
                self.canvas.clear_canvas()
                self.canvas._insert_molgraph(graph)
            self._current_file_path = filepath
            self.canvas.undo_stack.setClean()
            self._add_recent_file(filepath)
            self._update_total_charge_indicator()
            self.statusBar().showMessage(f"Abierto: {filepath}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"No se pudo abrir el archivo:\n{e}")
    
    def _on_file_save(self) -> None:
        """Save the current work in .cmsn format."""
        filepath = self._current_file_path
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
                if filepath.lower().endswith(".mol") or filepath.lower().endswith(".sdf") or "MOL" in selected_filter:
                    # Export as .mol if explicitly requested
                    from chemio.rdkit_io import molgraph_to_molfile
                    molfile = molgraph_to_molfile(self.canvas.graph)
                    with open(filepath, "w") as f:
                        f.write(molfile)
                else:
                    # Save as .cmsn (native)
                    if not filepath.lower().endswith(".cmsn"):
                        filepath += ".cmsn"
                    PersistenceManager.save_to_file(filepath, self.canvas)
                self._current_file_path = filepath
                self.canvas.undo_stack.setClean()
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
                    from chemio.rdkit_io import molgraph_to_svg
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

    def _confirm_discard_changes(self) -> bool:
        """Método auxiliar para  confirm discard changes.

        Returns:
            Resultado de la operación o None.

        Side Effects:
            Puede modificar el estado interno o la interfaz.
        """
        if self.canvas.undo_stack.isClean():
            return True
        reply = QMessageBox.question(
            self,
            "Cambios sin guardar",
            "Hay cambios sin guardar. ¿Deseas guardar antes de salir?",
            QMessageBox.StandardButton.Save
            | QMessageBox.StandardButton.Discard
            | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Save,
        )
        if reply == QMessageBox.StandardButton.Save:
            self._on_file_save()
            return self.canvas.undo_stack.isClean()
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
        if self._confirm_discard_changes():
            event.accept()
        else:
            event.ignore()

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
            
            if format == "smiles":
                from chemio.rdkit_io import molgraph_to_smiles
                text = molgraph_to_smiles(self.canvas.graph)
            elif format == "molfile":
                from chemio.rdkit_io import molgraph_to_molfile
                text = molgraph_to_molfile(self.canvas.graph)
            elif format == "inchi":
                # InChI requires additional RDKit import
                from chemio.rdkit_io import molgraph_to_rdkit
                from rdkit.Chem.inchi import MolToInchi
                mol = molgraph_to_rdkit(self.canvas.graph)
                text = MolToInchi(mol)
            else:
                text = ""
            
            clipboard.setText(text)
            self.statusBar().showMessage(f"Copiado como {format.upper()}")
        except Exception as e:
            self.statusBar().showMessage(f"Error: {e}")
    
    def _on_preferences(self) -> None:
        """Open preferences dialog."""
        dialog = PreferencesDialog(self.canvas.state, self.canvas.drawing_style, self)
        dialog.preferences_changed.connect(self._apply_preferences)
        dialog.exec()

    def _on_style_dialog(self) -> None:
        """Open drawing style dialog."""
        dialog = StyleDialog(self.canvas.drawing_style, self.canvas.state.bond_length, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            style, _bond_length = dialog.selected_style()
            self.canvas.apply_drawing_style(style)

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
            float(self.canvas.state.label_font_size),
            6.0,
            72.0,
            1,
        )
        if not ok:
            return
        font = self.canvas.label_font()
        font.setPointSizeF(float(size))
        self.canvas.apply_label_font(font)
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
        font = self.canvas.label_font()
        size = font.pointSizeF()
        if size <= 0:
            size = font.pointSize()
        if size <= 0:
            size = 10.0
        size = max(6.0, size + delta)
        font.setPointSizeF(size)
        self.canvas.apply_label_font(font)

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

        def center(coords: dict[int, tuple[float, float]]) -> tuple[float, float]:
            """Método auxiliar para center.

            Args:
                coords: Descripción del parámetro.

            Returns:
                Resultado de la operación o None.

            Side Effects:
                Puede modificar el estado interno o la interfaz.
            """
            xs = [x for x, _ in coords.values()]
            ys = [y for _, y in coords.values()]
            if not xs:
                return 0.0, 0.0
            return sum(xs) / len(xs), sum(ys) / len(ys)

        try:
            from chemio.rdkit_io import molgraph_to_rdkit_with_map
            from rdkit.Chem import AllChem

            graph = self.canvas._build_selection_graph(atom_ids, bonds) if atom_ids else self.canvas.graph
            mol, id_map = molgraph_to_rdkit_with_map(graph)
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()

            before = {
                aid: (self.canvas.model.get_atom(aid).x, self.canvas.model.get_atom(aid).y)
                for aid in target_ids
            }
            before_cx, before_cy = center(before)

            raw_after = {}
            for aid, rd_idx in id_map.items():
                if aid not in target_ids:
                    continue
                pos = conf.GetAtomPosition(rd_idx)
                raw_after[aid] = (pos.x, pos.y)

            after_cx, after_cy = center(raw_after)
            after = {}
            for aid, (x, y) in raw_after.items():
                x = x - after_cx + before_cx
                y = y - after_cy + before_cy
                if step_ratio < 1.0:
                    bx, by = before[aid]
                    x = bx + (x - bx) * step_ratio
                    y = by + (y - by) * step_ratio
                after[aid] = (x, y)

            from gui.commands import MoveAtomsCommand
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

    def _start_template_insert_by_id(self, template_id: str) -> None:
        """Activa inserción por clic para una plantilla de biblioteca."""
        try:
            template = self.template_library.get_template(template_id)
            graph = self.template_library.graph_from_template(template_id)
            label = str(template.get("name", "Plantilla")).strip() or "Plantilla"
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
            from chemio.rdkit_io import smiles_to_molgraph
            graph = smiles_to_molgraph(smiles.strip())
            self.canvas.clear_canvas()
            self.canvas._insert_molgraph(graph)
            self.statusBar().showMessage("SMILES importado")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"No se pudo importar SMILES:\n{e}")

    def _on_export_smiles(self) -> None:
        """Export the current molecule as SMILES."""
        try:
            from chemio.rdkit_io import molgraph_to_smiles
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
        QMessageBox.about(
            self,
            "Acerca de Chemuson",
            "<h2>Chemuson</h2>"
            "<p>Editor Molecular Libre</p>"
            "<p>Versión 0.2.0</p>"
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
            "tool_select": "Seleccionar",
            "tool_select_lasso": "Seleccion (lazo)",
            "tool_bond": "Enlace",
            "tool_coordination_bond": "Enlace coordinativo",
            "tool_rotate_3d_precise": "Rotación 3D precisa",
            "tool_ring": ring_label,
            "tool_atom": f"Elemento {self.canvas.state.default_element}",
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
            "tool_brackets_curly": "Llaves",
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
        name = tool_names.get(tool_id, tool_id)
        self.statusBar().showMessage(f"Herramienta: {name}")

    def _on_selection_changed(self, num_atoms: int, num_bonds: int, num_text: int, details: dict):
        """Handle selection change to update UI components."""
        self.inspector_dock.update_selection(num_atoms, num_bonds, num_text, details)
        self._update_total_charge_indicator()
        
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
