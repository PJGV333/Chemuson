"""Punto de entrada de la aplicación de dibujo molecular Chemuson.

Este módulo inicializa PyQt6, carga la ventana principal y arranca el bucle
de eventos. Es el archivo que se ejecuta al iniciar la aplicación.
"""

import sys
import os

# Aseguramos que Python encuentre los módulos dentro de `src` al ejecutar
# el archivo directamente.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from PyQt6.QtWidgets import QApplication
from gui.main_window import ChemusonWindow

def main():
    """
    Arranca la aplicación Qt y muestra la ventana principal.

    Returns:
        Código de salida del proceso de Qt.

    Side Effects:
        Crea la instancia de `QApplication`, muestra la ventana y entra en el
        bucle de eventos de Qt.
    """
    app = QApplication(sys.argv)
    app.setApplicationName("Chemuson")
    
    window = ChemusonWindow()
    window.show()
    
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
