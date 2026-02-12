"""
Diálogo de tabla periódica para seleccionar elementos.
"""
from __future__ import annotations

from PyQt6.QtWidgets import QDialog, QGridLayout, QPushButton, QVBoxLayout, QLabel
from PyQt6.QtCore import pyqtSignal, Qt


# Representación de la tabla periódica como cuadrícula de símbolos.
ELEMENT_GRID = [
    ["H", None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, "He"],
    ["Li", "Be", None, None, None, None, None, None, None, None, None, None, "B", "C", "N", "O", "F", "Ne"],
    ["Na", "Mg", None, None, None, None, None, None, None, None, None, None, "Al", "Si", "P", "S", "Cl", "Ar"],
    ["K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"],
    ["Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe"],
    ["Cs", "Ba", "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"],
    ["Fr", "Ra", "Ac", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"],
    [None, None, None, "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"],
    [None, None, None, "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"],
]


class PeriodicTableDialog(QDialog):
    """Diálogo modal con una tabla periódica simplificada."""
    element_selected = pyqtSignal(str)

    def __init__(self, parent=None) -> None:
        """Inicializa el diálogo y construye la interfaz."""
        super().__init__(parent)
        self.setWindowTitle("Tabla Periódica")
        self.setModal(True)
        self._build_ui()

    def _build_ui(self) -> None:
        """Construye la cuadrícula de botones de elementos."""
        layout = QVBoxLayout(self)
        title = QLabel("Selecciona un elemento")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title)

        grid = QGridLayout()
        grid.setSpacing(2)

        for row_idx, row in enumerate(ELEMENT_GRID):
            for col_idx, symbol in enumerate(row):
                if symbol is None:
                    continue
                button = QPushButton(symbol)
                button.setFixedSize(32, 28)
                button.clicked.connect(lambda checked=False, s=symbol: self._select(s))
                grid.addWidget(button, row_idx, col_idx)

        layout.addLayout(grid)

    def _select(self, symbol: str) -> None:
        """Emite el elemento seleccionado y cierra el diálogo."""
        self.element_selected.emit(symbol)
        self.accept()
