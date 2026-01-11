"""
Chemuson Main Window
Page-based molecular editor with vertical toolbar and paper canvas.
"""
from PyQt6.QtWidgets import QMainWindow, QToolBar
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction, QKeySequence

from gui.canvas import ChemusonCanvas
from gui.periodic_table import PeriodicTableDialog
from gui.toolbar import ChemusonToolbar


class ChemusonWindow(QMainWindow):
    """
    Main window for the Chemuson molecular editor.
    Features a page-based editing interface similar to ChemDoodle.
    """
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Chemuson - Editor Molecular Libre")
        self.resize(1200, 900)
        
        # Set gray "desktop" background
        self.setStyleSheet("""
            QMainWindow {
                background-color: #E0E0E0;
            }
        """)
        
        # === TOP TOOLBAR (File actions - placeholder) ===
        self.file_toolbar = QToolBar("Archivo")
        self.file_toolbar.setMovable(False)
        self.file_toolbar.setFloatable(False)
        self.file_toolbar.setStyleSheet("""
            QToolBar {
                background-color: #F5F5F5;
                border-bottom: 1px solid #CCCCCC;
                spacing: 4px;
                padding: 4px;
            }
        """)
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, self.file_toolbar)
        
        # Placeholder file actions
        self._add_file_action("📄 Nuevo", "file_new")
        self._add_file_action("📂 Abrir", "file_open")
        self._add_file_action("💾 Guardar", "file_save")
        self.file_toolbar.addSeparator()
        self._add_file_action("🔍+ Zoom +", "zoom_in")
        self._add_file_action("🔍- Zoom -", "zoom_out")
        
        # === LEFT TOOLBAR (Drawing tools) ===
        self.toolbar = ChemusonToolbar()
        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, self.toolbar)
        
        # === CENTRAL CANVAS ===
        self.canvas = ChemusonCanvas()
        self.setCentralWidget(self.canvas)

        # Undo/redo actions from the canvas stack
        self.file_toolbar.addSeparator()
        undo_action = self.canvas.undo_stack.createUndoAction(self, "Deshacer")
        undo_action.setShortcut(QKeySequence.StandardKey.Undo)
        redo_action = self.canvas.undo_stack.createRedoAction(self, "Rehacer")
        redo_action.setShortcut(QKeySequence.StandardKey.Redo)
        self.file_toolbar.addAction(undo_action)
        self.file_toolbar.addAction(redo_action)

        # Clipboard actions
        self.file_toolbar.addSeparator()
        copy_action = QAction("Copiar", self)
        copy_action.setShortcut(QKeySequence.StandardKey.Copy)
        copy_action.triggered.connect(self.canvas.copy_to_clipboard)
        paste_action = QAction("Pegar", self)
        paste_action.setShortcut(QKeySequence.StandardKey.Paste)
        paste_action.triggered.connect(self.canvas.paste_from_clipboard)
        self.file_toolbar.addAction(copy_action)
        self.file_toolbar.addAction(paste_action)
        
        # Connect toolbar to canvas
        self.toolbar.tool_changed.connect(self.canvas.set_current_tool)
        self.toolbar.bond_palette_changed.connect(self._handle_bond_palette)
        self.toolbar.ring_palette_changed.connect(self._handle_ring_palette)
        self.toolbar.element_palette_changed.connect(self._handle_element_palette)
        self.toolbar.periodic_table_requested.connect(self._show_periodic_table)
        # Sync defaults selected during toolbar init
        self._handle_bond_palette(self.toolbar.current_bond_spec())
        self._handle_ring_palette(self.toolbar.current_ring_spec())
        self._handle_element_palette(self.toolbar.current_element())
        
        # === STATUS BAR ===
        self.statusBar().showMessage("Herramienta: Carbono (C)")
        self.statusBar().setStyleSheet("""
            QStatusBar {
                background-color: #F5F5F5;
                border-top: 1px solid #CCCCCC;
                color: #666666;
            }
        """)
        
        # Update status bar when tool changes
        self.toolbar.tool_changed.connect(self._update_status)
    
    def _add_file_action(self, text: str, action_id: str) -> QAction:
        """Add a file action to the top toolbar."""
        action = QAction(text, self)
        action.setObjectName(action_id)
        self.file_toolbar.addAction(action)
        # Connect actions (placeholder - will implement later)
        action.triggered.connect(lambda: self._handle_file_action(action_id))
        return action
    
    def _handle_file_action(self, action_id: str) -> None:
        """Handle file toolbar actions."""
        if action_id == "zoom_in":
            self.canvas.zoom_in()
        elif action_id == "zoom_out":
            self.canvas.zoom_out()
        elif action_id == "file_new":
            self.canvas.clear_canvas()
        else:
            self.statusBar().showMessage(f"Acción '{action_id}' no implementada aún")

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
            'tool_select': 'Seleccionar',
            'tool_erase': 'Borrar',
            'tool_bond': 'Enlace',
            'tool_ring': f'Anillo {self.canvas.state.active_ring_size}',
            'tool_atom': f'Elemento {self.canvas.state.default_element}',
        }
        name = tool_names.get(tool_id, tool_id)
        self.statusBar().showMessage(f"Herramienta: {name}")
