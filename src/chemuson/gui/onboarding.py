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

**Render de la máscara (corrección post-push, gate manual Fase 7):** el
agujero se obtiene restando caminos con :class:`QPainterPath`
(``outer - inner``) y rellenando solo esa diferencia. Antes se usaba
``CompositionMode_Clear``, que en KDE/Wayland real dejaba franjas/bordes
negros alrededor de la zona resaltada (a veces hasta la barra inferior):
limpiar píxeles del backing store de un widget translúcido no es un contrato
estable en ese backend. Con la resta de caminos la máscara se rellena y el
agujero simplemente **no se pinta**, de modo que la ventana bajo el overlay
se ve sin artefactos en ToolRail, Canvas y SidePanel, en light y dark.

**Presentación de la tarjeta:** la tarjeta es theme-aware. Ya no usa colores
hardcodeados ni ``setStyleSheet`` por widget: se resuelve con el sistema QSS
de tokens (``#onboardCard`` y sus hijos en ``theme/qss.py``), de modo que
light y dark son coherentes y los hijos (``QCheckBox``, ``QPushButton``) no
heredan un tema distinto al de la tarjeta. La altura del cuerpo se fija al
máximo de los tres pasos para que la tarjeta no "salte" al cambiar de paso.

El overlay emite :signal:`finished(bool)` al cerrarse: ``True`` si se
completaron los tres pasos o si se cerró anticipadamente con "No volver a
mostrar" marcado; ``False`` si se cerró anticipadamente sin marcarlo (el
onboarding debe volver a aparecer en el siguiente arranque). El ensamblado
decide la persistencia de ``ui/onboarding/completed``.
"""
from __future__ import annotations

from PyQt6.QtCore import QEvent, QPoint, QRect, QRectF, QTimer, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QFontMetrics, QPainter, QPainterPath
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
#: Es intencionadamente oscuro en ambos temas (spotlight de onboarding).
_MASK_COLOR = QColor(15, 23, 42, 140)
#: Radio del borde del "agujero" que resalta la zona.
_HOLE_RADIUS = 10
#: Margen entre el agujero y el borde del widget objetivo.
_HOLE_MARGIN = 6
#: Métricas de la tarjeta. Los colores NO viven aquí: los resuelve el QSS de
#: tokens (``#onboardCard`` en ``theme/qss.py``) para que light y dark sean
#: coherentes.
_CARD_WIDTH = 320
_CARD_PADDING_X = 18
_CARD_PADDING_Y = 14


class _Card(QWidget):
    """Tarjeta del paso: título, texto, navegación y ``No volver a mostrar``.

    La presentación es theme-aware vía objectName + QSS de tokens
    (``#onboardCard`` y sus hijos en ``theme/qss.py``); la tarjeta no pinta
    colores propios ni hereda un tema distinto al de sus hijos.
    """

    previous = pyqtSignal()
    next = pyqtSignal()
    close = pyqtSignal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("onboardCard")
        # El fondo/borde de la tarjeta viene del QSS de tokens: un QWidget
        # necesita WA_StyledBackground para que QSS le pinte el fondo.
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.setFixedWidth(_CARD_WIDTH)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(_CARD_PADDING_X, _CARD_PADDING_Y,
                                  _CARD_PADDING_X, _CARD_PADDING_Y)
        layout.setSpacing(8)

        self._title = QLabel(self)
        self._title.setObjectName("onboardTitle")
        self._title.setWordWrap(True)
        layout.addWidget(self._title)

        self._body = QLabel(self)
        self._body.setObjectName("onboardBody")
        self._body.setWordWrap(True)
        # Altura reservada = la que necesita el texto más largo de los tres
        # pasos: así la tarjeta no cambia de tamaño al navegar (sin saltos).
        self._body.setMinimumHeight(self._max_body_height())
        layout.addWidget(self._body)

        self._no_more = QCheckBox("No volver a mostrar", self)
        self._no_more.setObjectName("onboardCheck")
        layout.addWidget(self._no_more)

        # Botones centrados (stretch a ambos lados): con el QSS de la tarjeta
        # (``min-width: 0``) los tres caben sin clipping.
        buttons = QHBoxLayout()
        buttons.setContentsMargins(0, 0, 0, 0)
        buttons.setSpacing(8)
        self._prev = QPushButton("Anterior", self)
        self._prev.setObjectName("onboardPrev")
        self._close = QPushButton("Cerrar", self)
        self._close.setObjectName("onboardClose")
        self._next = QPushButton("Siguiente", self)
        self._next.setObjectName("onboardNext")
        self._prev.clicked.connect(self.previous.emit)
        self._next.clicked.connect(self.next.emit)
        self._close.clicked.connect(self.close.emit)
        buttons.addStretch(1)
        buttons.addWidget(self._prev)
        buttons.addWidget(self._close)
        buttons.addWidget(self._next)
        buttons.addStretch(1)
        layout.addLayout(buttons)

    def set_step(self, index: int) -> None:
        title, body, _key = _STEPS[index]
        self._title.setText(f"{index + 1} de {len(_STEPS)} · {title}")
        self._body.setText(body)
        self._prev.setEnabled(index > 0)
        self._next.setText("Siguiente")

    def _max_body_height(self) -> int:
        """Altura (px) que necesita el texto más largo de los tres pasos."""
        metrics = QFontMetrics(self._body.font())
        available = _CARD_WIDTH - 2 * _CARD_PADDING_X
        return max(
            metrics.boundingRect(
                QRect(0, 0, available, 0),
                Qt.AlignmentFlag.AlignLeft | Qt.TextFlag.TextWordWrap,
                body,
            ).height()
            for _title, body, _key in _STEPS
        )

    @property
    def no_more(self) -> bool:
        return self._no_more.isChecked()

    @no_more.setter
    def no_more(self, value: bool) -> None:
        self._no_more.setChecked(value)

    @property
    def title_label(self) -> QLabel:
        return self._title

    @property
    def body_label(self) -> QLabel:
        return self._body

    @property
    def no_more_checkbox(self) -> QCheckBox:
        return self._no_more

    @property
    def buttons(self) -> tuple[QPushButton, QPushButton, QPushButton]:
        return (self._prev, self._close, self._next)


class OnboardingOverlay(QWidget):
    """Overlay de onboarding: máscara + agujero transparente + tarjeta.

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
        self._hole = QRect()
        self._relayout_timer = QTimer(self)
        self._relayout_timer.setSingleShot(True)
        self._relayout_timer.timeout.connect(self._relayout)
        self._geometry_event_types = {
            QEvent.Type.Show,
            QEvent.Type.Hide,
            QEvent.Type.Resize,
            QEvent.Type.Move,
            QEvent.Type.LayoutRequest,
        }
        for event_name in ("ScreenChangeInternal", "DevicePixelRatioChange"):
            event_type = getattr(QEvent.Type, event_name, None)
            if event_type is not None:
                self._geometry_event_types.add(event_type)
        parent_widget = self.parentWidget()
        if parent_widget is not None:
            parent_widget.installEventFilter(self)
        for target in self._targets:
            if target is not None:
                target.installEventFilter(self)
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

    def eventFilter(self, watched, event) -> bool:  # noqa: N802
        """Recalcula el spotlight cuando cambian ventana, layout o pantalla."""
        if (
            (watched is self.parentWidget() or watched in self._targets)
            and event.type() in self._geometry_event_types
            and not self._relayout_timer.isActive()
        ):
            # Espera a que Qt termine de propagar el resize/layout del padre y
            # sus hijos; así las coordenadas globales ya reflejan el layout final.
            self._relayout_timer.start(0)
        return super().eventFilter(watched, event)

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
        """Rect del widget objetivo en coordenadas del overlay.

        El overlay es hermano (no ancestro) de las zonas resaltadas, así que
        ``mapTo`` no es válido: se mapea globalmente y se trae al overlay, que
        cubre exactamente la rect de la ventana. Sin esto Qt avisa
        (``QWidget::mapTo(): parent must be in parent hierarchy``) y el agujero
        queda mal posicionado.
        """
        widget = self._targets[index] if index < len(self._targets) else None
        if widget is None or not widget.isVisible():
            return QRect()
        top_left = self.mapFromGlobal(widget.mapToGlobal(QPoint(0, 0)))
        return QRect(top_left, widget.size())

    def _relayout(self) -> None:
        parent = self.parentWidget()
        if parent is not None and self.geometry() != parent.rect():
            self.setGeometry(parent.rect())
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
    # Pintado (máscara por resta de caminos; el agujero NO se pinta)
    # ------------------------------------------------------------------
    def _mask_path(self) -> QPainterPath:
        """Camino rellenable de la máscara: ``rect completo - agujero``.

        El contrato es geométrico, no de "limpieza de píxeles": se rellena
        únicamente ``outer - inner`` y el agujero queda sin pintar (visible a
        través del overlay). Esto evita los artefactos negros que producía
        ``CompositionMode_Clear`` en KDE/Wayland.
        """
        outer = QPainterPath()
        outer.addRect(QRectF(self.rect()))
        inner = QPainterPath()
        inner.addRoundedRect(QRectF(self._hole), _HOLE_RADIUS, _HOLE_RADIUS)
        return outer.subtracted(inner)

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.fillPath(self._mask_path(), _MASK_COLOR)
        painter.end()

    # ------------------------------------------------------------------
    # Acceso (para tests)
    # ------------------------------------------------------------------
    @property
    def step(self) -> int:
        return self._step

    @property
    def hole(self) -> QRect:
        """Rectángulo del agujero actual (en coordenadas del overlay)."""
        return self._hole

    @property
    def no_more(self) -> bool:
        return self._card.no_more
