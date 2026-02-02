"""
Chemuson Main Window
Page-based molecular editor with menu bar, toolbars, and paper canvas.
"""
from PyQt6.QtWidgets import (
    QDialog,
    QFontDialog,
    QInputDialog,
    QMainWindow,
    QMenu,
    QMenuBar,
    QFileDialog,
    QMessageBox,
    QToolBar,
)
from PyQt6.QtCore import Qt, QSize
from PyQt6.QtGui import QAction, QActionGroup, QKeySequence, QIcon, QPainter
from PyQt6.QtPrintSupport import QPrinter

from gui.canvas import ChemusonCanvas
from gui.periodic_table import PeriodicTableDialog
from gui.toolbar import ChemusonToolbar
from gui.styles import MAIN_STYLESHEET, TOOL_PALETTE_STYLESHEET
from gui.icons import draw_generic_icon
from gui.docks import PlantillasDock, InspectorDock
from gui.dialogs import PreferencesDialog, QuickStartDialog, StyleDialog
from gui.text_toolbar import TextFormatToolbar
from gui.commands import ChangeAtomCommand
from gui.templates import (
    build_linear_chain_template,
)
import os
import glob
from chemio.rdkit_io import molfile_to_molgraph, molgraph_to_molfile
from core.model import MolGraph


class ChemusonWindow(QMainWindow):
    """
    Main window for the Chemuson molecular editor.
    Features a page-based editing interface similar to ChemDoodle.
    """
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Chemuson - Editor Molecular Libre")
        self.resize(1200, 900)
        
        # Apply main stylesheet
        self.setStyleSheet(MAIN_STYLESHEET)
        
        # === CORE COMPONENTS ===
        self._create_actions()
        
        # === CENTRAL CANVAS ===
        self.canvas = ChemusonCanvas()
        self.setCentralWidget(self.canvas)
        self.action_aromatic_circles.setChecked(self.canvas.state.use_aromatic_circles)
        self.action_rules.setChecked(self.canvas.show_rulers)
        self.action_crosshair.setChecked(self.canvas.show_grid)
        
        # === DOCK WIDGETS ===
        self.templates_dock = PlantillasDock(self)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.templates_dock)
        self.templates_dock.hide()
        
        self.inspector_dock = InspectorDock(self)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.inspector_dock)
        self.inspector_dock.hide()
        
        # === MENU AND TOOLBARS ===
        self._create_menu_bar()
        self._create_main_toolbar()
        self._sync_label_menu_state()
        
        # === LEFT TOOLBAR (Drawing tools) ===
        self.toolbar = ChemusonToolbar()
        self.toolbar.setStyleSheet(TOOL_PALETTE_STYLESHEET)
        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, self.toolbar)
        
        # === TEXT FORMAT TOOLBAR ===
        self.text_toolbar = TextFormatToolbar()
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.text_toolbar)
        # Initially hide it, or show it only when text tool is active?
        # User requested it to be available. We can leave it visible or toggle it.
        # For now, visible is fine.

        # Text Toolbar Connections
        self.text_toolbar.format_changed.connect(self.canvas.update_text_format)
        self.text_toolbar.color_changed.connect(self.canvas.update_text_color)
        
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

        # Sync defaults selected during toolbar init
        self._handle_bond_palette(self.toolbar.current_bond_spec())
        self._handle_ring_palette(self.toolbar.current_ring_spec())
        self._handle_element_palette(self.toolbar.current_element())
        
        # === STATUS BAR ===
        self.statusBar().showMessage("Herramienta: Carbono (C)")
        
        # Update status bar when tool changes
        self.toolbar.tool_changed.connect(self._update_status)

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
        self.action_save.setShortcut(QKeySequence.StandardKey.Save)
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

        self.action_rules = QAction("Reglas", self)
        self.action_rules.setCheckable(True)
        self.action_rules.setChecked(False)
        self.action_rules.triggered.connect(self._on_toggle_rules)

        self.action_crosshair = QAction("Cuadrícula", self)
        self.action_crosshair.setCheckable(True)
        self.action_crosshair.setChecked(False)
        self.action_crosshair.triggered.connect(self._on_toggle_crosshair)

        self.action_clean_2d = QAction("Limpiar 2D", self)
        self.action_clean_2d.triggered.connect(self._on_clean_2d)

        self.action_template_linear_chain = QAction("Cadena lineal", self)
        self.action_template_linear_chain.triggered.connect(self._on_insert_linear_chain)

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

        # --- Rotation Actions ---
        self.action_rotate_left = QAction("Girar 90° a la izquierda", self)
        self.action_rotate_left.triggered.connect(lambda: self._on_rotate_selection(-90.0))

        self.action_rotate_right = QAction("Girar 90° a la derecha", self)
        self.action_rotate_right.triggered.connect(lambda: self._on_rotate_selection(90.0))

        self.action_flip_horizontal = QAction("Giro 180° horizontal", self)
        self.action_flip_horizontal.triggered.connect(self._on_flip_horizontal)

        self.action_flip_vertical = QAction("Giro 180° vertical", self)
        self.action_flip_vertical.triggered.connect(self._on_flip_vertical)

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

        self.action_preferences = QAction("Preferencias...", self)
        self.action_preferences.triggered.connect(self._on_preferences)
        edit_menu.addAction(self.action_preferences)

        # === Texto (Text) ===
        text_menu = menubar.addMenu("Texto")
        text_menu.addAction(self.action_label_font)
        text_menu.addAction(self.action_label_size_set)
        text_menu.addSeparator()
        text_menu.addAction(self.action_label_bold)
        text_menu.addAction(self.action_label_italic)
        text_menu.addAction(self.action_label_underline)
        text_menu.addSeparator()
        text_menu.addAction(self.action_label_subscript)
        text_menu.addAction(self.action_label_superscript)
        text_menu.addSeparator()
        text_menu.addAction(self.action_label_size_up)
        text_menu.addAction(self.action_label_size_down)
        text_menu.addSeparator()
        label_color_menu = text_menu.addMenu("Color de etiquetas")
        label_color_menu.addAction(self.action_label_color_element)
        label_color_menu.addAction(self.action_label_color_black)
        
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
        view_menu.addAction(self.action_rules)
        view_menu.addAction(self.action_crosshair)
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
        view_menu.addAction(self.templates_dock.toggleViewAction())
        view_menu.addAction(self.inspector_dock.toggleViewAction())
        
        # === Estructura (Structure) ===
        structure_menu = menubar.addMenu("Estructura")
        structure_menu.addAction(self.action_clean_2d)
        structure_menu.addSeparator()
        structure_menu.addSeparator()
        
        # Dynamic Templates Menu
        self.templates_menu = structure_menu.addMenu("Plantillas")
        self._refresh_templates_menu()
        
        # Save Template Action
        self.action_save_template = QAction("Guardar selección como plantilla...", self)
        self.action_save_template.triggered.connect(self._on_save_template)
        structure_menu.addAction(self.action_save_template)
        
        structure_menu.addSeparator()
        structure_menu.addAction(self.action_import_smiles)
        structure_menu.addAction(self.action_export_smiles)

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
        
        # Clipboard actions
        self.main_toolbar.addAction(self.action_copy)
        self.main_toolbar.addAction(self.action_paste)
        self.main_toolbar.addSeparator()
        
        # View actions
        self.main_toolbar.addAction(self.action_zoom_in)
        self.main_toolbar.addAction(self.action_zoom_out)
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
    
    # -------------------------------------------------------------------------
    # File Menu Handlers
    # -------------------------------------------------------------------------
    def _on_file_new(self) -> None:
        """Create a new empty canvas."""
        self.canvas.clear_canvas()
        self.statusBar().showMessage("Nuevo documento creado")
    
    def _on_file_open(self) -> None:
        """Open a molecule file."""
        filepath, _ = QFileDialog.getOpenFileName(
            self,
            "Abrir archivo",
            "",
            "Archivos MOL (*.mol *.sdf);;Todos los archivos (*.*)"
        )
        if filepath:
            try:
                with open(filepath, "r") as f:
                    molfile = f.read()
                from chemio.rdkit_io import molfile_to_molgraph
                graph = molfile_to_molgraph(molfile)
                self.canvas.clear_canvas()
                self.canvas._insert_molgraph(graph)
                self.statusBar().showMessage(f"Abierto: {filepath}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"No se pudo abrir el archivo:\n{e}")
    
    def _on_file_save(self) -> None:
        """Save the current molecule."""
        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Guardar archivo",
            "",
            "Archivo MOL (*.mol);;Todos los archivos (*.*)"
        )
        if filepath:
            try:
                from chemio.rdkit_io import molgraph_to_molfile
                molfile = molgraph_to_molfile(self.canvas.graph)
                with open(filepath, "w") as f:
                    f.write(molfile)
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
            from chemio.rdkit_io import molgraph_to_svg
            svg_text = molgraph_to_svg(self.canvas.graph)
            filepath, _ = QFileDialog.getSaveFileName(self, "Exportar SVG", "", "Imagen SVG (*.svg)")
            if filepath:
                with open(filepath, "w") as f:
                    f.write(svg_text)
                self.statusBar().showMessage(f"Exportado: {filepath}")
        elif format == "pdf":
            filepath, _ = QFileDialog.getSaveFileName(self, "Exportar PDF", "", "Documento PDF (*.pdf)")
            if filepath:
                try:
                    printer = QPrinter(QPrinter.PrinterMode.HighResolution)
                    printer.setOutputFormat(QPrinter.OutputFormat.PdfFormat)
                    printer.setOutputFileName(filepath)
                    
                    # Target rectangle for the paper on the PDF page
                    painter = QPainter(printer)
                    paper_rect = self.canvas.paper.sceneBoundingRect()
                    self.canvas.scene.render(painter, printer.pageRect(QPrinter.Unit.Point), paper_rect)
                    painter.end()
                    self.statusBar().showMessage(f"Exportado: {filepath}")
                except Exception as e:
                    QMessageBox.warning(self, "Error", f"No se pudo exportar PDF:\n{e}")
    
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
        dialog = PreferencesDialog(self.canvas.state, self)
        dialog.preferences_changed.connect(self._apply_preferences)
        dialog.exec()

    def _on_style_dialog(self) -> None:
        """Open drawing style dialog."""
        dialog = StyleDialog(self.canvas.drawing_style, self.canvas.state.bond_length, self)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            style, _bond_length = dialog.selected_style()
            self.canvas.apply_drawing_style(style)

    def _apply_preferences(self, prefs: dict) -> None:
        self.canvas.state.show_implicit_carbons = prefs.get("show_carbons", False)
        self.canvas.state.show_implicit_hydrogens = prefs.get("show_hydrogens", False)
        self.canvas.state.use_aromatic_circles = prefs.get("aromatic_circles", False)

        self.action_show_carbons.setChecked(self.canvas.state.show_implicit_carbons)
        self.action_show_hydrogens.setChecked(self.canvas.state.show_implicit_hydrogens)
        self.action_aromatic_circles.setChecked(self.canvas.state.use_aromatic_circles)

        self.canvas.refresh_atom_visibility()
        self.canvas.refresh_aromatic_circles()
        self._sync_label_menu_state()
    
    # -------------------------------------------------------------------------
    # View Menu Handlers
    # -------------------------------------------------------------------------
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
        self.canvas.rotate_selection_degrees(angle_deg)

    def _on_flip_horizontal(self) -> None:
        self.canvas.flip_selection_horizontal()

    def _on_flip_vertical(self) -> None:
        self.canvas.flip_selection_vertical()

    # -------------------------------------------------------------------------
    # Text Menu Handlers
    # -------------------------------------------------------------------------
    def _sync_label_menu_state(self) -> None:
        self.action_label_bold.setChecked(self.canvas.state.label_font_bold)
        self.action_label_italic.setChecked(self.canvas.state.label_font_italic)
        self.action_label_underline.setChecked(self.canvas.state.label_font_underline)
        self.action_label_color_element.setChecked(self.canvas.state.use_element_colors)
        self.action_label_color_black.setChecked(not self.canvas.state.use_element_colors)

    def _on_label_font(self) -> None:
        font, ok = QFontDialog.getFont(
            self.canvas.label_font(),
            self,
            "Fuente de etiquetas",
        )
        if ok:
            self.canvas.apply_label_font(font)
            self._sync_label_menu_state()

    def _on_label_font_size_dialog(self) -> None:
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
        font = self.canvas.label_font()
        font.setBold(checked)
        self.canvas.apply_label_font(font)
        self._sync_label_menu_state()

    def _on_label_italic(self, checked: bool) -> None:
        font = self.canvas.label_font()
        font.setItalic(checked)
        self.canvas.apply_label_font(font)
        self._sync_label_menu_state()

    def _on_label_underline(self, checked: bool) -> None:
        font = self.canvas.label_font()
        font.setUnderline(checked)
        self.canvas.apply_label_font(font)
        self._sync_label_menu_state()

    def _on_label_font_size(self, delta: float) -> None:
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
        self.canvas.set_paper_size(width, height)
        self.statusBar().showMessage(f"Lienzo: {width} x {height} px")

    def _on_canvas_custom_size(self) -> None:
        width, ok = QInputDialog.getInt(
            self,
            "Tamaño de lienzo",
            "Ancho (px):",
            int(self.canvas.paper_width),
            200,
            20000,
            1,
        )
        if not ok:
            return
        height, ok = QInputDialog.getInt(
            self,
            "Tamaño de lienzo",
            "Alto (px):",
            int(self.canvas.paper_height),
            200,
            20000,
            1,
        )
        if not ok:
            return
        self._set_canvas_size(width, height)

    def _apply_label_script(self, marker: str, title: str) -> None:
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
        self._apply_label_script("_", "Subíndice")

    def _on_label_superscript(self) -> None:
        self._apply_label_script("^", "Superíndice")

    def _on_label_color_mode(self, use_element_colors: bool) -> None:
        self.canvas.set_use_element_colors(use_element_colors)
        self._sync_label_menu_state()
        self.statusBar().showMessage(
            "Etiquetas: por elemento" if use_element_colors else "Etiquetas: negro"
        )
    
    # -------------------------------------------------------------------------
    # Structure Menu Handlers
    # -------------------------------------------------------------------------
    def _on_clean_2d(self) -> None:
        """Clean 2D coordinates using RDKit."""
        try:
            from chemio.rdkit_io import molgraph_to_rdkit, rdkit_to_molgraph
            from rdkit.Chem import AllChem
            
            mol = molgraph_to_rdkit(self.canvas.graph)
            AllChem.Compute2DCoords(mol)
            cleaned = rdkit_to_molgraph(mol)
            
            # Preserve current center position
            self.canvas.clear_canvas()
            self.canvas._insert_molgraph(cleaned)
            self.statusBar().showMessage("Estructura 2D limpiada")
        except Exception as e:
            self.statusBar().showMessage(f"Error: {e}")

    def _insert_template(self, label: str, graph) -> None:
        self.canvas._insert_molgraph(graph)
        self.statusBar().showMessage(f"Plantilla: {label}")

    def _get_templates_dir(self) -> str:
        """Get the directory where user templates are stored."""
        # Use a local templates directory relative to the source
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        templates_dir = os.path.join(base_dir, "templates")
        if not os.path.exists(templates_dir):
            os.makedirs(templates_dir)
        return templates_dir

    def _refresh_templates_menu(self) -> None:
        """Scan templates directory and populate the menu."""
        self.templates_menu.clear()
        
        # Add standard linear chain (generated)
        self.action_linear = QAction("Cadena lineal (Fischer)", self)
        self.action_linear.triggered.connect(self._on_insert_linear_chain)
        self.templates_menu.addAction(self.action_linear)
        self.templates_menu.addSeparator()

        templates_dir = self._get_templates_dir()
        files = glob.glob(os.path.join(templates_dir, "*.mol"))
        
        if not files:
            no_templates = QAction("(Sin plantillas guardadas)", self)
            no_templates.setEnabled(False)
            self.templates_menu.addAction(no_templates)
            return

        for filepath in sorted(files):
            filename = os.path.basename(filepath)
            name = os.path.splitext(filename)[0].replace("_", " ").title()
            action = QAction(name, self)
            # Use default argument capture for lambda loop
            action.triggered.connect(lambda checked=False, f=filepath: self._insert_template_from_file(f))
            self.templates_menu.addAction(action)

    def _on_save_template(self) -> None:
        """Save the current selection (or whole canvas) as a new template."""
        # 1. Determine what to save
        if self.canvas.state.selected_atoms:
            # Create a subgraph from selection
            # This requires extracting specific atoms and bonds
            # For simplicity in this task, let's warn if multiple disjoint, 
            # but currently we just allow saving the subgraph.
            
            # Since MolGraph doesn't have a 'subgraph' method easily exposed here,
            # we can clone the graph and remove unselected.
            # Or simplified: Serialize the Whole Canvas if selection is empty?
            pass
        
        # Strategy: Save the WHOLE canvas if nothing is selected? 
        # Or just the selection. let's support Selection Only if atoms selected.
        
        try:
            graph_to_save = MolGraph()
            
            if self.canvas.state.selected_atoms:
                # Copy selected atoms
                visible_atoms = self.canvas.state.selected_atoms
                # We need to copy them with new IDs or keep IDs? 
                # molgraph_to_molfile handles arbitrary IDs usually.
                
                # Copy atoms
                old_to_new = {}
                for atom_id in visible_atoms:
                    old_atom = self.canvas.model.get_atom(atom_id)
                    new_atom = graph_to_save.add_atom(
                        old_atom.element, 
                        old_atom.x, 
                        old_atom.y, 
                        charge=old_atom.charge,
                        explicit=old_atom.is_explicit
                    )
                    old_to_new[atom_id] = new_atom.id
                
                # Copy bonds if both ends are selected
                for bond in self.canvas.model.bonds.values():
                    if bond.a1_id in visible_atoms and bond.a2_id in visible_atoms:
                        graph_to_save.add_bond(
                            old_to_new[bond.a1_id],
                            old_to_new[bond.a2_id],
                            order=bond.order,
                            style=bond.style,
                            stereo=bond.stereo
                        )
            else:
                # Save entire canvas
                # We can just use deepcopy or serialize/deserialize
                # Simplest: use rdkit io to roundtrip or just copy logic
                # Let's use the canvas graph directly for serialization
                graph_to_save = self.canvas.graph

            if not graph_to_save.atoms:
                QMessageBox.warning(self, "Aviso", "No hay nada para guardar.")
                return

            name, ok = QInputDialog.getText(self, "Nueva Plantilla", "Nombre de la plantilla:")
            if not ok or not name.strip():
                return
            
            safe_name = "".join(c for c in name if c.isalnum() or c in (' ', '_', '-')).strip()
            filename = f"{safe_name}.mol"
            path = os.path.join(self._get_templates_dir(), filename)
            
            molblock = molgraph_to_molfile(graph_to_save)
            with open(path, "w") as f:
                f.write(molblock)
            
            self.statusBar().showMessage(f"Plantilla guardada: {name}")
            self._refresh_templates_menu()
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error al guardar plantilla:\n{e}")

    def _insert_template_from_file(self, filepath: str) -> None:
        try:
            with open(filepath, "r") as f:
                molblock = f.read()
            graph = molfile_to_molgraph(molblock)
            name = os.path.splitext(os.path.basename(filepath))[0]
            self._insert_template(name, graph)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Error al cargar plantilla:\n{e}")

    def _on_insert_linear_chain(self) -> None:
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
        self.canvas.set_active_bond(bond_spec)
        self._update_status("tool_bond")

    def _handle_ring_palette(self, ring_spec: dict) -> None:
        self.canvas.set_active_ring(ring_spec)
        self._update_status("tool_ring")

    def _handle_element_palette(self, element: str) -> None:
        self.canvas.set_active_element(element)
        self._update_status("tool_atom")

    def _show_periodic_table(self) -> None:
        dialog = PeriodicTableDialog(self)
        dialog.element_selected.connect(self.toolbar.select_element)
        dialog.exec()
    
    def _update_status(self, tool_id: str) -> None:
        """Update status bar with current tool."""
        tool_names = {
            "tool_select": "Seleccionar",
            "tool_select_lasso": "Seleccion (lazo)",
            "tool_erase": "Borrar",
            "tool_bond": "Enlace",
            "tool_ring": f"Anillo {self.canvas.state.active_ring_size}",
            "tool_atom": f"Elemento {self.canvas.state.default_element}",
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
        }
        name = tool_names.get(tool_id, tool_id)
        self.statusBar().showMessage(f"Herramienta: {name}")

    def _on_selection_changed(self, num_atoms: int, num_bonds: int, num_text: int, details: dict):
        """Handle selection change to update UI components."""
        self.inspector_dock.update_selection(num_atoms, num_bonds, num_text, details)
        
        # Sync Text Toolbar if a single text item is selected
        if num_text == 1 and details.get("type") == "text":
            font = details.get("font")
            settings = {
                "color": details.get("color"),
                "sub": details.get("sub"),
                "sup": details.get("sup")
            }
            self.text_toolbar.set_state(font, settings)
