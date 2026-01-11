"""
Chemuson Main Window
Page-based molecular editor with vertical toolbar and paper canvas.
"""
from PyQt6.QtWidgets import QMainWindow, QToolBar, QLabel
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction, QKeySequence

from gui.canvas import ChemusonCanvas
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
        
        # Connect toolbar to canvas
        self.toolbar.tool_changed.connect(self.canvas.set_current_tool)
        
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
    
    def _update_status(self, tool_id: str) -> None:
        """Update status bar with current tool."""
        tool_names = {
            'tool_select': 'Seleccionar',
            'tool_erase': 'Borrar',
            'bond_single': 'Enlace Sencillo',
            'bond_double': 'Enlace Doble',
            'ring_benzene': 'Anillo de Benceno',
            'atom_C': 'Carbono (C)',
            'atom_N': 'Nitrógeno (N)',
            'atom_O': 'Oxígeno (O)',
            'atom_S': 'Azufre (S)',
            'atom_P': 'Fósforo (P)',
            'atom_F': 'Flúor (F)',
            'atom_Cl': 'Cloro (Cl)',
            'atom_Br': 'Bromo (Br)',
        }
        name = tool_names.get(tool_id, tool_id)
        self.statusBar().showMessage(f"Herramienta: {name}")
