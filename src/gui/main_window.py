from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QWidget, QLabel
from PyQt6.QtCore import Qt
from gui.canvas import ChemusonCanvas
from gui.toolbar import ChemusonToolbar

class ChemusonWindow(QMainWindow):
    """
    Ventana principal del editor molecular Chemuson.
    """
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Chemuson - Editor Molecular Libre")
        self.resize(1024, 768)
        
        # Barra de herramientas
        self.toolbar = ChemusonToolbar()
        self.addToolBar(self.toolbar)
        
        # Widget central
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        
        self.layout = QVBoxLayout(self.central_widget)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        # Canvas interactivo
        self.canvas = ChemusonCanvas()
        self.layout.addWidget(self.canvas)

        # Conectar Toolbar con Canvas
        self.toolbar.tool_changed.connect(self.canvas.set_current_tool)

        # Label de estado (opcional, para feedback de herramientas)
        self.info_label = QLabel("Herramienta: Carbono (C)")
        self.info_label.setAlignment(Qt.AlignmentFlag.AlignRight)
        self.info_label.setStyleSheet("background-color: #f8f8f8; padding: 2px; color: #888; border-top: 1px solid #ddd;")
        self.layout.addWidget(self.info_label)
        
        # Actualizar label al cambiar herramienta
        self.toolbar.tool_changed.connect(lambda t: self.info_label.setText(f"Herramienta: {t}"))
        
        # Aquí luego añadiremos la barra de herramientas y el canvas
