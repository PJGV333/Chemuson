"""Onboarding nativo de la UI moderna (Fase 7 del plan de modernización).

``OnboardingOverlay`` es un hijo de la ventana principal que, en la primera
ejecución, guía al usuario por las tres zonas primarias (rail de
herramientas, lienzo y panel lateral) mediante:

- una máscara translúcida que oscurece la ventana, con un "agujero"
  (región transparente) que resalta la zona del paso actual; y
- una tarjeta con título, texto, botones ``Anterior``/``Siguiente``/
  ``Cerrar`` y una casilla ``No volver a mostrar``.

No toca la escena ni el canvas: el agujero se posiciona sobre geometrías
públicas (``tool_rail``, el widget central y ``side_panel``) mapeadas a las
coordenadas del overlay. La persistencia usa ``platform.settings`` (el mismo
``QSettings`` que el resto de la GUI) con la clave
``ui/onboarding/completed``.

El overlay emite :signal:`finished(bool)` al cerrarse: ``True`` si se
completaron los tres pasos o si se cerró anticipadamente con "No volver a
mostrar" marcado; ``False`` si se cerró anticipadamente sin marcarlo (el
onboarding debe volver a aparecer en el siguiente arranque). El ensamblado
decide la persistencia de ``ui/onboarding/completed``.
"""
from __future__ import annotations

from PyQt6.QtCore import QPoint, QRect, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QPainter
from PyQt6.QtWidgets import (
    QCheckBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

__all__ = ["OnboardingOverlay"]

#: Pasos fijos: (título, texto, clave de zona).
#: El texto es el del brief de la Fase 7.
_STEPS: tuple[tuple[str, str, str], ...] = (
    ("Rail de herramientas",
     "Elige aquí las herramientas de dibujo y anotación.", "rail"),
    ("Lienzo",
     "Dibuja, selecciona y edita tus estructuras en el lienzo.", "canvas"),
    ("Panel lateral",
     "Inspector, validación, propiedades, plantillas y apariencia están aquí.",
     "sidepanel"),
)

#: Color de la máscara (slate-900 ~55 %): oscurece sin tapar el contenido.
_MASK_COLOR = QColor(15, 23, 42, 140)
#: Radio del borde del "agujero" que resalta la zona.
_HOLE_RADIUS = 10
#: Margen entre el agujero y el borde del widget objetivo.
_HOLE_MARGIN = 6
#: Métricas de la tarjeta.
_CARD_WIDTH = 320
_CARD_BG = QColor("#F8FAFC")
_CARD_TITLE = QColor("#0F172A")
_CARD_BODY = QColor("#475569")


class _Card(QWidget):
    """Tarjeta del paso: título, texto, navegación y ``No volver a mostrar``."""

    previous = pyqtSignal()
    next = pyqtSignal()
    close = pyqtSignal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("onboardCard")
        self.setFixedWidth(_CARD_WIDTH)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 14, 18, 14)
        layout.setSpacing(8)

        self._title = QLabel(self)
        self._title.setStyleSheet(f"color: {_CARD_TITLE.name()}; font-weight: 600;")
        layout.addWidget(self._title)

        self._body = QLabel(self)
        self._body.setWordWrap(True)
        self._body.setStyleSheet(f"color: {_CARD_BODY.name()};")
        layout.addWidget(self._body)

        self._no_more = QCheckBox("No volver a mostrar", self)
        layout.addWidget(self._no_more)

        buttons = QHBoxLayout()
        self._prev = QPushButton("Anterior", self)
        self._next = QPushButton("Siguiente", self)
        self._close = QPushButton("Cerrar", self)
        self._prev.clicked.connect(self.previous.emit)
        self._next.clicked.connect(self.next.emit)
        self._close.clicked.connect(self.close.emit)
        buttons.addWidget(self._prev)
        buttons.addWidget(self._close)
        buttons.addStretch(1)
        buttons.addWidget(self._next)
        layout.addLayout(buttons)

    def set_step(self, index: int) -> None:
        title, body, _key = _STEPS[index]
        self._title.setText(f"{index + 1} de {len(_STEPS)} · {title}")
        self._body.setText(body)
        self._prev.setEnabled(index > 0)
        self._next.setText("Siguiente")

    @property
    def no_more(self) -> bool:
        return self._no_more.isChecked()

    @no_more.setter
    def no_more(self, value: bool) -> None:
        self._no_more.setChecked(value)

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setBrush(QColor(_CARD_BG))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawRoundedRect(self.rect().adjusted(1, 1, -1, -1), 10, 10)
        painter.end()


class OnboardingOverlay(QWidget):
    """Overlay de onboarding: máscara + agujero sobre la zona + tarjeta.

    Args:
        parent: Ventana principal que contiene el overlay.
        targets: Lista de ``QWidget`` objetivo por paso (en el orden de
            :data:`_STEPS``: rail, lienzo/central, panel lateral). ``None``
            en una posición desactiva el resalte de ese paso (tarjeta central).
    """

    finished = pyqtSignal(bool)

    def __init__(
        self,
        parent: QWidget | None,
        targets: list[QWidget | None] | None = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("onboardOverlay")
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground)
        self._targets = list(targets) if targets is not None else [None, None, None]
        self._step = 0
        self._card = _Card(self)
        self._card.previous.connect(self._on_previous)
        self._card.next.connect(self._on_next)
        self._card.close.connect(self._on_close)
        self._card.set_step(self._step)

    # ------------------------------------------------------------------
    # Ciclo de vida
    # ------------------------------------------------------------------
    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._relayout()

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._relayout()

    # ------------------------------------------------------------------
    # Navegación
    # ------------------------------------------------------------------
    def _on_previous(self) -> None:
        if self._step > 0:
            self._step -= 1
            self._card.set_step(self._step)
            self._relayout()

    def _on_next(self) -> None:
        if self._step < len(_STEPS) - 1:
            self._step += 1
            self._card.set_step(self._step)
            self._relayout()
        else:
            # Completar los tres pasos cuenta como onboarding completado.
            self._finish(completed=True)

    def _on_close(self) -> None:
        # Cerrar anticipadamente: solo cuenta como completado si el usuario
        # marcó "No volver a mostrar".
        self._finish(completed=self._card.no_more)

    # ------------------------------------------------------------------
    # API pública de navegación (también usada por tests)
    # ------------------------------------------------------------------
    def advance(self) -> None:
        """Avanza al paso siguiente; en el último paso completa (``finished(True)``)."""
        self._on_next()

    def go_back(self) -> None:
        """Vuelve al paso anterior (sin efecto en el primero)."""
        self._on_previous()

    def request_close(self) -> None:
        """Cierra el onboarding (equivalente a pulsar ``Cerrar``).

        Cuenta como completado solo si el usuario marcó "No volver a
        mostrar"; en caso contrario emite ``finished(False)`` y el onboarding
        se ofrece de nuevo en el siguiente arranque.
        """
        self._on_close()

    def set_no_more(self, value: bool) -> None:
        """Marca o desmarca "No volver a mostrar" (equivalente al usuario)."""
        self._card.no_more = value

    @property
    def card(self) -> _Card:
        """La tarjeta del paso (acceso para tests/presentación)."""
        return self._card

    def _finish(self, completed: bool) -> None:
        """Cierra el overlay y emite ``finished`` con el resultado.

        Args:
            completed: ``True`` si el onboarding debe considerarse completado
                (3 pasos completados o cierre con "No volver a mostrar");
                ``False`` si se cerró anticipadamente sin marcarlo (debe
                repetirse en el siguiente arranque).
        """
        self.finished.emit(completed)
        self.hide()

    # ------------------------------------------------------------------
    # Geometría (máscara + agujero + tarjeta)
    # ------------------------------------------------------------------
    def _target_rect(self, index: int) -> QRect:
        widget = self._targets[index] if index < len(self._targets) else None
        if widget is None or not widget.isVisible():
            return QRect()
        top_left = widget.mapTo(self, QPoint(0, 0))
        return QRect(top_left, widget.size())

    def _relayout(self) -> None:
        self.setGeometry(self.parentWidget().rect() if self.parentWidget() else self.rect())
        hole = self._target_rect(self._step)
        if not hole.isValid():
            hole = self.rect().adjusted(
                self.width() // 4, self.height() // 4,
                -self.width() // 4, -self.height() // 4,
            )
        hole = hole.adjusted(_HOLE_MARGIN, _HOLE_MARGIN, -_HOLE_MARGIN, -_HOLE_MARGIN)
        self._hole = hole
        self._position_card()
        self.update()

    def _position_card(self) -> None:
        card = self._card
        card.show()
        card.adjustSize()
        cw, ch = card.width(), card.height()
        x = max(16, min((self.width() - cw) // 2, self.width() - cw - 16))
        y = self._hole.bottom() + 16
        if y + ch > self.height() - 16:
            y = self._hole.top() - ch - 16
        y = max(16, min(y, self.height() - ch - 16))
        card.move(x, y)
        card.raise_()

    # ------------------------------------------------------------------
    # Pintado (máscara + agujero transparente)
    # ------------------------------------------------------------------
    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.fillRect(self.rect(), _MASK_COLOR)
        # El "agujero" se limpia para dejar ver la ventana bajo la máscara.
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_Clear)
        painter.drawRoundedRect(self._hole, _HOLE_RADIUS, _HOLE_RADIUS)
        painter.setCompositionMode(QPainter.CompositionMode.CompositionMode_SourceOver)
        painter.end()

    # ------------------------------------------------------------------
    # Acceso (para tests)
    # ------------------------------------------------------------------
    @property
    def step(self) -> int:
        return self._step

    @property
    def no_more(self) -> bool:
        return self._card.no_more
