import sys
import os

# Aseguramos que Python encuentre los módulos dentro de 'src'
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from PyQt6.QtWidgets import QApplication
from gui.main_window import ChemusonWindow

def main():
    """
    Punto de entrada principal de la aplicación Chemuson.
    """
    app = QApplication(sys.argv)
    app.setApplicationName("Chemuson")
    
    window = ChemusonWindow()
    window.show()
    
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
