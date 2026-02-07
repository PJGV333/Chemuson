"""
Elementos de escena de Chemuson.

Define subclases de QGraphicsItem para átomos, enlaces y anotaciones con
renderizado profesional.
"""
from __future__ import annotations

import math
from typing import Optional
from PyQt6.QtWidgets import (
    QGraphicsEllipseItem,
    QGraphicsPathItem,
    QGraphicsTextItem,
    QGraphicsItem,
    QGraphicsRectItem,
    QStyle,
)
from PyQt6.QtGui import QColor, QFont, QFontMetrics, QPainterPath, QPen, QBrush
from PyQt6.QtCore import Qt, QRectF, QPointF


from core.model import Atom, Bond, BondStyle
from gui.style import DrawingStyle, CHEMDOODLE_LIKE
from gui.wedge_geometry import compute_wedge_points


# Colores de elementos (esquema CPK simplificado) para etiquetas.
ELEMENT_COLORS = {
    'C': '#333333',   # Carbon - dark gray
    'N': '#3050F8',   # Nitrogen - blue
    'O': '#FF0D0D',   # Oxygen - red
    'S': '#FFFF30',   # Sulfur - yellow
    'P': '#FF8000',   # Phosphorus - orange
    'F': '#90E050',   # Fluorine - light green
    'Cl': '#1FF01F',  # Chlorine - green
    'Br': '#A62929',  # Bromine - dark red
    'I': '#940094',   # Iodine - purple
    'H': '#000000',   # Hydrogen - black
}

# Colores específicos para etiquetas (H en negro).
LABEL_ELEMENT_COLORS = {**ELEMENT_COLORS, "H": "#000000"}
# Etiquetas abreviadas que se muestran como texto directo.
ABBREVIATION_LABELS = {"Me", "Et", "Pr", "iPr", "tBu", "Bu", "Ph", "R", "TBS", "Si"}

# Fine-tuning knobs for wedge integration at the wide terminal when there is
# exactly one outgoing simple bond. Keep these at module scope so they are
# easy to tweak without touching geometry logic.
WEDGE_SINGLE_END_OVERLAP_STROKE_MULT = 4.0
WEDGE_SINGLE_END_OVERLAP_WIDTH_MULT = 3.0
# Single fine-tuning knob for end_deg == 1 terminal join.
# >1.0 pushes the wide corner more into the simple bond; <1.0 retracts it.
WEDGE_SINGLE_END_JOIN_TUNE = 1
# Terminal cap for end_deg == 1: keeps a short flat base segment colinear
# with the outgoing simple bond (ChemDraw-like trapezoid ending).
WEDGE_SINGLE_END_TERMINAL_CLIP_STROKE_MULT = 10000
WEDGE_SINGLE_END_TERMINAL_CLIP_WIDTH_MULT = 50000
# Small forward overlay so wedge cap sits on top of the simple bond and
# hides round line caps at the junction.
WEDGE_SINGLE_END_CAP_OVERLAY_STROKE_MULT = 0.1
WEDGE_SINGLE_END_CAP_OVERLAY_WIDTH_MULT = 10000
WEDGE_SINGLE_END_CAP_OVERLAY_SPAN_MULT = 0.12
WEDGE_SINGLE_END_MIN_SEP_MULT = 0.22
WEDGE_MULTI_END_MIN_SEP_MULT = 1.6
WEDGE_MULTI_END_CORNER_JOIN_STROKE_MULT = 0.90
WEDGE_MULTI_END_CORNER_JOIN_WIDTH_MULT = 0.22
WEDGE_MULTI_END_CORNER_BLEND = 0.72
# Extra outward reach at the wide end so wedge corners touch adjacent
# simple bonds without changing fork/bifurcation geometry.
WEDGE_MULTI_END_CONTACT_STROKE_MULT = 0.6
WEDGE_MULTI_END_CONTACT_WIDTH_MULT = 0.14
WEDGE_MULTI_END_EXTENSION_STROKE_MULT = 0.87 #0.85
WEDGE_MULTI_END_EXTENSION_WIDTH_MULT = 0.22 #0.2
WEDGE_MULTI_END_FORK_DEPTH_STROKE_MULT = 1.8 #1.8
WEDGE_MULTI_END_FORK_DEPTH_WIDTH_MULT = 0.45 #0.45
WEDGE_MULTI_END_FORK_BLEND = 0.85


def _build_wavy_path(
    start: QPointF,
    end: QPointF,
    amplitude: float,
    wavelength: float,
    min_cycles: int = 2,
    segments_per_cycle: int = 12,
) -> QPainterPath:
    """Construye un trazo ondulado entre dos puntos.

    Args:
        start: Punto inicial en escena.
        end: Punto final en escena.
        amplitude: Amplitud de la onda.
        wavelength: Longitud de onda.
        min_cycles: Número mínimo de ciclos.
        segments_per_cycle: Segmentos por ciclo.

    Returns:
        Ruta QPainterPath con la onda calculada.

    Side Effects:
        No tiene efectos laterales.
    """
    path = QPainterPath()
    dx = end.x() - start.x()
    dy = end.y() - start.y()
    length = math.hypot(dx, dy)
    if length <= 1e-6:
        return path

    base_wavelength = max(wavelength, 1e-3)
    cycles = max(min_cycles, int(length / base_wavelength))
    cycles = max(1, cycles)
    actual_wavelength = length / cycles
    amplitude = min(amplitude, actual_wavelength * 0.3)

    nx = -dy / length
    ny = dx / length
    segments = max(24, cycles * segments_per_cycle)
    for i in range(segments + 1):
        t = i / segments
        offset = math.sin(t * 2.0 * math.pi * cycles) * amplitude
        px = start.x() + dx * t + nx * offset
        py = start.y() + dy * t + ny * offset
        if i == 0:
            path.moveTo(px, py)
        else:
            path.lineTo(px, py)
    return path


class AtomItem(QGraphicsEllipseItem):
    """Elemento gráfico que representa un átomo (con ocultación implícita)."""
    
    def __init__(
        self,
        atom: Atom,
        radius: float = 12.0,
        show_carbon: bool = True,
        show_hydrogen: bool = True,
        label_font: Optional[QFont] = None,
        style: DrawingStyle = CHEMDOODLE_LIKE,
        use_element_colors: bool = True,
    ) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Args:
            atom: Átomo del modelo asociado.
            radius: Radio visual en píxeles.
            show_carbon: Si se muestran etiquetas de carbono implícito.
            show_hydrogen: Si se muestran etiquetas de hidrógeno implícito.
            label_font: Fuente base para la etiqueta.
            style: Estilo de dibujo aplicado.
            use_element_colors: Si se colorea la etiqueta según el elemento.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__(-radius, -radius, radius * 2, radius * 2)
        self.atom_id = atom.id
        self.element = atom.element
        self._radius = radius
        self._show_carbon = show_carbon
        self._show_hydrogen = show_hydrogen
        self._is_explicit = atom.is_explicit
        self._is_selected = False
        self._is_hover = False
        self._valence_error = False
        self._use_element_colors = use_element_colors
        self._element_color = QColor(self._label_color_for_element(atom.element))
        self._charge = atom.charge
        self._style = style
        self._display_label = atom.element
        self._label_anchor_override: Optional[str] = None
        self._label_anchor_width: Optional[float] = None
        self._label_anchor_offset = 0.0
        self._label_offset = QPointF(0.0, 0.0)
        
        self.setPos(atom.x, atom.y)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        
        # Create text label
        self.label = QGraphicsTextItem("", self)
        base_font = label_font or QFont("Arial", 10, QFont.Weight.Bold)
        self.label.setFont(base_font)
        self._user_underline = base_font.underline()
        self.label.document().setDocumentMargin(0)
        self._set_label_text(self._display_label)

        # Charge label
        self.charge_label = QGraphicsTextItem("", self)
        charge_font = QFont(base_font)
        charge_font.setPointSizeF(max(base_font.pointSizeF() * 0.8, 6.0))
        self.charge_label.setFont(charge_font)
        self.charge_label.setDefaultTextColor(QColor("#333333"))

        self._update_charge_label()
        
        # Set color based on element
        self._set_element_color()
        
        # Apply visibility
        self._update_visibility()
    
    def _center_label(self) -> None:
        """Centra la etiqueta del átomo.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        rect = self.label.boundingRect()
        x = -rect.width() / 2
        if self._label_anchor_width is not None:
            x = -self._label_anchor_width / 2 - self._label_anchor_offset
        self.label.setPos(x + self._label_offset.x(), -rect.height() / 2 + self._label_offset.y())
        if hasattr(self, "charge_label"):
            self._position_charge_label()

    def _position_charge_label(self) -> None:
        """Posiciona la etiqueta de carga.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        rect = self.label.boundingRect()
        label_pos = self.label.pos()
        x = label_pos.x() + rect.width() / 2 + 2
        y = label_pos.y() - rect.height() / 2 - 2
        self.charge_label.setPos(x, y)
    
    def _base_label_color(self) -> str:
        """Calcula el color base de la etiqueta.

        Returns:
            Color base (hex) para etiquetas.

        Side Effects:
            No tiene efectos laterales.
        """
        return "#000000" if not self._use_element_colors else "#333333"

    def _label_color_for_element(self, element: str) -> str:
        """Calcula el color de etiqueta para un elemento.

        Args:
            element: Símbolo del elemento químico.

        Returns:
            Color (hex) según el elemento.

        Side Effects:
            No tiene efectos laterales.
        """
        if not self._use_element_colors:
            return "#000000"
        return LABEL_ELEMENT_COLORS.get(element, "#333333")

    def _set_element_color(self) -> None:
        """Actualiza el color del elemento.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._element_color = QColor(self._label_color_for_element(self.element))
        self._apply_label_style()

    def _apply_label_style(self) -> None:
        """Aplica el estilo visual de la etiqueta.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        font = self.label.font()
        font.setUnderline(self._user_underline or self._valence_error)
        self.label.setFont(font)
        if self._valence_error:
            color = "#C0392B"
        else:
            color = self._element_color.name()
        self.label.setDefaultTextColor(QColor(color))
        charge_color = "#C0392B" if self._valence_error else self._base_label_color()
        if hasattr(self, "charge_label"):
            self.charge_label.setDefaultTextColor(QColor(charge_color))
        self._set_label_text(self._display_label)
    
    def _should_hide_element(self) -> bool:
        """Determina si hide elemento.

        Returns:
            True si debe ocultarse la etiqueta/círculo.

        Side Effects:
            No tiene efectos laterales.
        """
        if self._is_explicit:
            return False
        if self.element == "C" and not self._show_carbon:
            return True
        if self.element == "H" and not self._show_hydrogen:
            return True
        return False

    def _should_draw_circle(self) -> bool:
        """Determina si draw círculo.

        Returns:
            True si debe dibujarse el círculo del átomo.

        Side Effects:
            No tiene efectos laterales.
        """
        if self.element == "H" or self.element == "C":
            return False
        return not (self._is_explicit and self.element not in {"C", "H"})

    def _update_visibility(self) -> None:
        """Actualiza la visibilidad del átomo.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        if self._should_hide_element():
            # Make circle invisible but keep selectable (minimal hit area)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self.setPen(QPen(Qt.PenStyle.NoPen))
            self.label.setVisible(False)
            return

        self.label.setVisible(True)
        self._apply_normal_style()
    
    def _apply_normal_style(self) -> None:
        """Aplica el estilo normal según el estado.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        if not self._should_draw_circle():
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self.setPen(QPen(Qt.PenStyle.NoPen))
            return
        if self._is_hover:
            self.setBrush(QBrush(QColor("#F0F0F0")))
            pen = QPen(QColor("#333333"), self._style.stroke_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
        else:
            self.setBrush(QBrush(QColor(self._style.atom_fill_color)))
            pen = QPen(QColor(self._style.atom_stroke_color), self._style.stroke_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)

    def paint(self, painter, option, widget=None) -> None:
        """Dibuja el elemento en la escena con el estilo actual.

        Args:
            painter: Pintor Qt usado para dibujar.
            option: Opciones de estilo del render.
            widget: Widget contenedor opcional.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        painter.setPen(self.pen())
        painter.setBrush(self.brush())
        painter.drawEllipse(self.rect())
    
    def set_visibility_flags(self, show_carbon: bool, show_hydrogen: bool) -> None:
        """Actualiza banderas de visibilidad.

        Args:
            show_carbon: Si se muestran etiquetas de carbono implícito.
            show_hydrogen: Si se muestran etiquetas de hidrógeno implícito.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._show_carbon = show_carbon
        self._show_hydrogen = show_hydrogen
        self._update_visibility()
    
    def set_selected(self, selected: bool) -> None:
        """Actualiza selección.

        Args:
            selected: Estado de selección (True/False).

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._is_selected = selected
        self._update_visibility()
    
    def set_hover(self, hover: bool) -> None:
        """Actualiza estado hover.

        Args:
            hover: Estado de hover (True/False).

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._is_hover = hover
        if not self._should_hide_element():
            self._apply_normal_style()
    
    def set_element(self, element: str, is_explicit: Optional[bool] = None) -> None:
        """Actualiza elemento.

        Args:
            element: Símbolo del elemento químico.
            is_explicit: Si el símbolo debe mostrarse explícito.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.element = element
        if is_explicit is not None:
            self._is_explicit = is_explicit
        self._display_label = element
        self._label_anchor_override = None
        self._set_element_color()
        self._update_visibility()

    def set_display_label(self, label: str, anchor: Optional[str] = None) -> None:
        """Actualiza etiqueta mostrada.

        Args:
            label: Texto de etiqueta a mostrar.
            anchor: Ancla opcional para centrar la etiqueta.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._display_label = label
        self._label_anchor_override = anchor
        self._set_label_text(label)

    def set_label_offset(self, offset: QPointF) -> None:
        """Actualiza desplazamiento de etiqueta.

        Args:
            offset: Desplazamiento en coordenadas de escena.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._label_offset = QPointF(offset)
        self._center_label()

    def set_label_font(self, font: QFont) -> None:
        """Actualiza fuente de etiqueta.

        Args:
            font: Fuente tipográfica a aplicar.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.label.setFont(font)
        self._user_underline = font.underline()
        charge_font = QFont(font)
        charge_font.setPointSizeF(max(font.pointSizeF() * 0.8, 6.0))
        self.charge_label.setFont(charge_font)
        self._apply_label_style()

    def set_valence_error(self, has_error: bool) -> None:
        """Actualiza error de valencia.

        Args:
            has_error: Indica si hay error de valencia.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._valence_error = has_error
        self._apply_label_style()

    def set_charge(self, charge: int) -> None:
        """Actualiza carga.

        Args:
            charge: Carga formal del átomo.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._charge = charge
        if charge == 0:
            self.charge_label.setVisible(False)
            return
        sign = "+" if charge > 0 else "-"
        magnitude = abs(charge)
        if magnitude == 1:
            label = "⊕" if charge > 0 else "⊖"
        else:
            label = f"{sign}{magnitude}"
        self.charge_label.setPlainText(label)
        self._position_charge_label()
        self.charge_label.setVisible(True)

    def _update_charge_label(self) -> None:
        """Actualiza la etiqueta de carga.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.set_charge(self._charge)

    def _set_label_text(self, text: str) -> None:
        """Establece el texto de la etiqueta.

        Args:
            text: Texto a mostrar.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        anchor = self._label_anchor_override or self._first_label_token(text)
        if anchor:
            metrics = QFontMetrics(self.label.font())
            self._label_anchor_width = metrics.horizontalAdvance(anchor)
            self._label_anchor_offset = self._anchor_prefix_width(text, anchor)
        else:
            self._label_anchor_width = None
            self._label_anchor_offset = 0.0
        html = self._format_label_html(text)
        if html:
            self.label.setHtml(html)
        else:
            self.label.setPlainText(text)
        self._center_label()

    def _anchor_prefix_width(self, text: str, anchor: str) -> float:
        """Calcula el ancho del prefijo antes del ancla.

        Args:
            text: Texto a mostrar.
            anchor: Ancla opcional para centrar la etiqueta.

        Returns:
            Ancho del prefijo antes del ancla en píxeles.

        Side Effects:
            No tiene efectos laterales.
        """
        if not text or not anchor:
            return 0.0
        index = text.rfind(anchor)
        if index <= 0:
            return 0.0
        prefix = text[:index]
        return self._measure_label_width(prefix)

    def _measure_label_width(self, text: str) -> float:
        """Calcula etiqueta width.

        Args:
            text: Texto a mostrar.

        Returns:
            Ancho calculado del texto en píxeles.

        Side Effects:
            No tiene efectos laterales.
        """
        if not text:
            return 0.0
        base_font = self.label.font()
        base_metrics = QFontMetrics(base_font)
        base_size = base_font.pointSizeF()
        if base_size <= 0:
            base_size = 10.0
        sub_size = max(base_size - 2.0, 6.0)
        sub_font = QFont(base_font)
        sub_font.setPointSizeF(sub_size)
        sub_metrics = QFontMetrics(sub_font)
        width = 0.0
        i = 0
        script_mode: Optional[str] = None
        while i < len(text):
            ch = text[i]
            if ch in "_^":
                script_mode = "sub" if ch == "_" else "sup"
                i += 1
                continue
            if script_mode is not None:
                j = i
                while j < len(text) and text[j].isalnum():
                    j += 1
                if j > i:
                    width += sub_metrics.horizontalAdvance(text[i:j])
                    i = j
                    script_mode = None
                    continue
                script_mode = None
            if ch.isdigit():
                j = i + 1
                while j < len(text) and text[j].isdigit():
                    j += 1
                width += sub_metrics.horizontalAdvance(text[i:j])
                i = j
                continue
            width += base_metrics.horizontalAdvance(ch)
            i += 1
        return width

    def _first_label_token(self, text: str) -> str | None:
        """Obtiene el primer token de la etiqueta.

        Args:
            text: Texto a mostrar.

        Returns:
            Primer token relevante o None si no existe.

        Side Effects:
            No tiene efectos laterales.
        """
        cleaned = text.strip()
        if not cleaned:
            return None
        if cleaned[-1] in "+-":
            cleaned = cleaned[:-1]
        cleaned = cleaned.lstrip("_^")
        if not cleaned:
            return None
        abbrs = sorted(ABBREVIATION_LABELS, key=len, reverse=True)
        for abbr in abbrs:
            if cleaned.startswith(abbr):
                return abbr
        ch = cleaned[0]
        if ch.isalpha():
            if ch.isupper():
                if len(cleaned) > 1 and cleaned[1].islower():
                    candidate = cleaned[:2]
                    if candidate in LABEL_ELEMENT_COLORS:
                        return candidate
                if ch in LABEL_ELEMENT_COLORS:
                    return ch
                j = 1
                while j < len(cleaned) and cleaned[j].islower():
                    j += 1
                return cleaned[:j]
            j = 1
            while j < len(cleaned) and cleaned[j].islower():
                j += 1
            return cleaned[:j]
        if ch.isdigit():
            return ch
        return ch

    def _format_label_html(self, text: str) -> str | None:
        """Formatea etiqueta html.

        Args:
            text: Texto a mostrar.

        Returns:
            HTML formateado o None si no aplica.

        Side Effects:
            No tiene efectos laterales.
        """
        cleaned = text.strip()
        if not cleaned:
            return None

        charge = None
        if cleaned[-1] in "+-":
            charge = cleaned[-1]
            cleaned = cleaned[:-1]

        tokens: list[tuple[str, str]] = []
        i = 0
        abbrs = sorted(ABBREVIATION_LABELS, key=len, reverse=True)
        script_mode: Optional[str] = None
        while i < len(cleaned):
            ch = cleaned[i]
            if ch in "_^":
                script_mode = "sub" if ch == "_" else "sup"
                i += 1
                continue
            if script_mode is not None:
                j = i
                while j < len(cleaned) and cleaned[j].isalnum():
                    j += 1
                if j > i:
                    tokens.append((script_mode, cleaned[i:j]))
                    i = j
                    script_mode = None
                    continue
                script_mode = None
            matched = None
            for abbr in abbrs:
                if cleaned.startswith(abbr, i):
                    matched = abbr
                    break
            if matched:
                tokens.append(("abbr", matched))
                i += len(matched)
                continue
            ch = cleaned[i]
            if ch.isdigit():
                j = i + 1
                while j < len(cleaned) and cleaned[j].isdigit():
                    j += 1
                tokens.append(("sub", cleaned[i:j]))
                i = j
                continue
            if ch.isalpha():
                if ch.isupper():
                    symbol = None
                    if i + 1 < len(cleaned) and cleaned[i + 1].islower():
                        candidate = cleaned[i : i + 2]
                        if candidate in LABEL_ELEMENT_COLORS:
                            symbol = candidate
                            i += 2
                    if symbol is None and ch in LABEL_ELEMENT_COLORS:
                        symbol = ch
                        i += 1
                    if symbol is not None:
                        tokens.append(("elem", symbol))
                        continue
                    j = i + 1
                    while j < len(cleaned) and cleaned[j].islower():
                        j += 1
                    tokens.append(("abbr", cleaned[i:j]))
                    i = j
                    continue
                j = i + 1
                while j < len(cleaned) and cleaned[j].islower():
                    j += 1
                tokens.append(("abbr", cleaned[i:j]))
                i = j
                continue
            if ch in "()":
                tokens.append(("text", ch))
                i += 1
                continue
            tokens.append(("text", ch))
            i += 1

        if not tokens:
            return None

        base_size = self.label.font().pointSizeF()
        if base_size <= 0:
            base_size = 10.0
        sub_size = max(base_size - 2.0, 6.0)
        parts: list[str] = []
        error_color = "#C0392B"
        use_error = self._valence_error
        base_color = self._base_label_color()
        last_color = error_color if use_error else base_color
        for kind, value in tokens:
            if kind == "elem":
                color = error_color if use_error else self._label_color_for_element(value)
                parts.append(f'<span style="color:{color};">{value}</span>')
                last_color = color
                continue
            if kind in {"abbr", "text"}:
                color = error_color if use_error else base_color
                parts.append(f'<span style="color:{color};">{value}</span>')
                last_color = color
                continue
            if kind == "sub":
                parts.append(
                    f'<sub><span style="font-size:{sub_size}pt; color:{last_color};">{value}</span></sub>'
                )
                continue
            if kind == "sup":
                parts.append(
                    f'<sup><span style="font-size:{sub_size}pt; color:{last_color};">{value}</span></sup>'
                )
                continue

        if charge:
            parts.append(
                f'<sup><span style="font-size:{sub_size}pt; color:{last_color};">{charge}</span></sup>'
            )

        return "".join(parts)

    def set_style(self, style: DrawingStyle) -> None:
        """Actualiza estilo.

        Args:
            style: Estilo de dibujo aplicado.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._style = style
        if not self._should_hide_element():
            self._apply_normal_style()

    def set_use_element_colors(self, use_element_colors: bool) -> None:
        """Actualiza uso de colores por elemento.

        Args:
            use_element_colors: Si se colorea la etiqueta según el elemento.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._use_element_colors = use_element_colors
        self._set_element_color()


class BondItem(QGraphicsPathItem):
    """Elemento gráfico que representa un enlace químico."""
    
    def __init__(
        self,
        bond: Bond,
        atom1: Atom,
        atom2: Atom,
        render_aromatic_as_circle: bool = True,
        style: DrawingStyle = CHEMDOODLE_LIKE,
    ) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Args:
            bond: Enlace del modelo asociado.
            atom1: Átomo inicial del enlace.
            atom2: Átomo final del enlace.
            render_aromatic_as_circle: Si los aromáticos se dibujan con círculo.
            style: Estilo de dibujo aplicado.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__()
        self.bond_id = bond.id
        self.a1_id = bond.a1_id
        self.a2_id = bond.a2_id
        self.order = bond.order
        self.style = bond.style
        self.stereo = bond.stereo
        self.is_aromatic = bond.is_aromatic
        self.display_order = bond.display_order
        self.ring_id = bond.ring_id
        self.length_px = bond.length_px
        self._stroke_px = bond.stroke_px
        self._color = bond.color
        self.render_aromatic_as_circle = render_aromatic_as_circle
        self._style = style
        self._label_shrink_start = 0.0
        self._label_shrink_end = 0.0
        self._endpoint_trim_start = 0.0
        self._endpoint_trim_end = 0.0
        self._endpoint_extend_start = 0.0
        self._endpoint_extend_end = 0.0
        self._offset_sign = 1
        self._ring_center: QPointF | None = None
        self._bond_in_ring = False
        self._prefer_full_length = False
        self._symmetric_double = False
        self._wedge_join_start: list[tuple[float, float, float]] = []
        self._wedge_join_end: list[tuple[float, float, float]] = []
        self.setZValue(-5)
        pen_color = QColor(self._style.bond_color)
        if self._color:
            candidate = QColor(self._color)
            if candidate.isValid():
                pen_color = candidate
        pen = QPen(pen_color, self._style.stroke_px)
        pen.setCapStyle(self._style.cap_style)
        pen.setJoinStyle(self._style.join_style)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.update_positions(atom1, atom2)

    def paint(self, painter, option, widget=None) -> None:
        """Dibuja el elemento en la escena con el estilo actual.

        Args:
            painter: Pintor Qt usado para dibujar.
            option: Opciones de estilo del render.
            widget: Widget contenedor opcional.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        painter.setPen(self.pen())
        painter.setBrush(self.brush())
        painter.drawPath(self.path())

    def set_bond(self, bond: Bond, atom1: Atom, atom2: Atom) -> None:
        """Actualiza enlace.

        Args:
            bond: Enlace del modelo asociado.
            atom1: Átomo inicial del enlace.
            atom2: Átomo final del enlace.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.order = bond.order
        self.style = bond.style
        self.stereo = bond.stereo
        self.is_aromatic = bond.is_aromatic
        self.display_order = bond.display_order
        self.ring_id = bond.ring_id
        self.length_px = bond.length_px
        self._stroke_px = bond.stroke_px
        self._color = bond.color
        self.update_positions(atom1, atom2)

    def set_render_aromatic_as_circle(self, enabled: bool) -> None:
        """Actualiza renderizado aromático con círculo.

        Args:
            enabled: Habilita o deshabilita el modo.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.render_aromatic_as_circle = enabled
        # Canvas refresh will call update_positions with current atom positions.

    def set_ring_context(self, ring_center: QPointF | None) -> None:
        """Actualiza contexto de anillo.

        Args:
            ring_center: Centro del anillo en coordenadas de escena.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._ring_center = ring_center

    def set_bond_in_ring(self, in_ring: bool) -> None:
        """Actualiza enlace en anillo.

        Args:
            in_ring: Indica si el enlace pertenece a un anillo.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._bond_in_ring = bool(in_ring)

    def set_offset_sign(self, sign: int) -> None:
        """Actualiza signo de desplazamiento.

        Args:
            sign: Signo de desplazamiento (+1/-1).

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._offset_sign = 1 if sign >= 0 else -1

    def set_multibond_rendering(self, prefer_full_length: bool, symmetric_double: bool) -> None:
        """Actualiza renderizado de multienlaces.

        Args:
            prefer_full_length: Si se prefiere longitud completa en multienlaces.
            symmetric_double: Si el doble enlace se dibuja simétrico.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._prefer_full_length = bool(prefer_full_length)
        self._symmetric_double = bool(symmetric_double)

    def set_label_shrink(self, start: float, end: float) -> None:
        """Actualiza recorte de etiquetas.

        Args:
            start: Punto inicial en escena.
            end: Punto final en escena.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._label_shrink_start = max(0.0, float(start))
        self._label_shrink_end = max(0.0, float(end))

    def set_endpoint_trim(self, start: float, end: float) -> None:
        """Actualiza recorte de extremos.

        Args:
            start: Punto inicial en escena.
            end: Punto final en escena.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._endpoint_trim_start = max(0.0, float(start))
        self._endpoint_trim_end = max(0.0, float(end))

    def set_endpoint_extend(self, start: float, end: float) -> None:
        """Actualiza extensión de extremos.

        Args:
            start: Punto inicial en escena.
            end: Punto final en escena.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._endpoint_extend_start = max(0.0, float(start))
        self._endpoint_extend_end = max(0.0, float(end))

    def set_wedge_join_neighbors(
        self,
        start_neighbors: list[tuple[float, float, float]],
        end_neighbors: list[tuple[float, float, float]],
    ) -> None:
        """Configura contexto de enlaces vecinos para adaptar la cuña.

        Cada vecino se representa como `(ux, uy, width_px)` en coordenadas de
        escena, donde `(ux, uy)` es el vector unitario desde el átomo del
        extremo de la cuña hacia el vecino.
        """
        cleaned_start: list[tuple[float, float, float]] = []
        for ux, uy, width in start_neighbors:
            length = math.hypot(ux, uy)
            if length <= 1e-6:
                continue
            cleaned_start.append((ux / length, uy / length, max(0.0, float(width))))
        cleaned_end: list[tuple[float, float, float]] = []
        for ux, uy, width in end_neighbors:
            length = math.hypot(ux, uy)
            if length <= 1e-6:
                continue
            cleaned_end.append((ux / length, uy / length, max(0.0, float(width))))
        self._wedge_join_start = cleaned_start
        self._wedge_join_end = cleaned_end

    def _adapt_wedge_corner_for_neighbors(
        self,
        base_cx: float,
        base_cy: float,
        ux: float,
        uy: float,
        nx: float,
        ny: float,
        side_sign: float,
        default_corner: tuple[float, float],
        neighbors: list[tuple[float, float, float]],
        stroke_px: float,
    ) -> tuple[float, float]:
        """Adapta una esquina de base al enlace vecino del mismo lado."""
        dx = default_corner[0] - base_cx
        dy = default_corner[1] - base_cy
        default_x = dx * ux + dy * uy
        default_y = side_sign * (dx * nx + dy * ny)
        if default_y <= 1e-6:
            return default_corner

        best_neighbor: tuple[float, float, float] | None = None
        best_neighbor_any: tuple[float, float, float] | None = None
        best_score = -1e9
        best_score_any = -1e9
        for nux, nuy, nwidth in neighbors:
            cross = ux * nuy - uy * nux
            dot = ux * nux + uy * nuy
            any_score = abs(cross) + max(0.0, dot) * 10 #0.28
            if any_score > best_score_any:
                best_score_any = any_score
                best_neighbor_any = (nux, nuy, nwidth)
            if side_sign * cross <= 0.03:
                continue
            score = abs(cross) + max(0.0, dot) * 15  # 0.3
            if score > best_score:
                best_score = score
                best_neighbor = (nux, nuy, nwidth)
        # If no neighbor falls on this side, keep default for single-neighbor
        # terminals; fallback to strongest only when multiple neighbors exist.
        if best_neighbor is None:
            if len(neighbors) <= 1 or best_neighbor_any is None:
                return default_corner
            best_neighbor = best_neighbor_any

        nux, nuy, nwidth = best_neighbor
        nnx = -nuy
        nny = nux
        half_neighbor = max(nwidth * 0.5, stroke_px * 0.45)
        e1x = base_cx + nnx * half_neighbor
        e1y = base_cy + nny * half_neighbor
        e2x = base_cx - nnx * half_neighbor
        e2y = base_cy - nny * half_neighbor
        e1_side = side_sign * ((e1x - base_cx) * nx + (e1y - base_cy) * ny)
        e2_side = side_sign * ((e2x - base_cx) * nx + (e2y - base_cy) * ny)
        edge_x, edge_y = (e1x, e1y) if e1_side >= e2_side else (e2x, e2y)

        join = max(stroke_px * 1.05, min(default_y * 1.35, nwidth * 1.10))
        target_x = edge_x + nux * join
        target_y = edge_y + nuy * join
        tx = (target_x - base_cx) * ux + (target_y - base_cy) * uy
        ty = side_sign * ((target_x - base_cx) * nx + (target_y - base_cy) * ny)

        # Constrain adaptation: allow elongated joins without collapsing width.
        max_forward = max(stroke_px * 1.70, default_y * 2.10)
        max_back = stroke_px * 0.16
        tx = max(default_x - max_back, min(default_x + max_forward, tx))
        min_y = max(stroke_px * 0.45, default_y * 0.42)
        max_y = max(default_y * 1.25, min_y + 0.1)
        ty = max(min_y, min(max_y, ty))

        blend = 0.95
        fx = default_x * (1.0 - blend) + tx * blend
        fy = default_y * (1.0 - blend) + ty * blend
        return (
            base_cx + ux * fx + nx * side_sign * fy,
            base_cy + uy * fx + ny * side_sign * fy,
        )

    def _adapt_wedge_tip_for_neighbors(
        self,
        tip_x: float,
        tip_y: float,
        ux: float,
        uy: float,
        nx: float,
        ny: float,
        width: float,
        stroke_px: float,
        neighbors: list[tuple[float, float, float]],
    ) -> tuple[tuple[float, float], tuple[float, float]]:
        """Adapta la punta: punto agudo o borde corto (trapezoidal)."""
        best: tuple[float, float, float] | None = None
        best_score = -1e9
        for nux, nuy, nwidth in neighbors:
            # Prefer neighbor opposite to wedge axis (chain continuation).
            oppose = -(ux * nux + uy * nuy)
            if oppose <= 0.05:
                continue
            cross = abs(ux * nuy - uy * nux)
            score = oppose + cross * 0.2
            if score > best_score:
                best_score = score
                best = (nux, nuy, nwidth)

        if best is None:
            return (tip_x, tip_y), (tip_x, tip_y)

        _, _, nwidth = best
        tip_half = max(stroke_px * 0.42, min(nwidth * 0.32, width * 0.14))
        tip_back = max(stroke_px * 0.10, min(stroke_px * 0.36, width * 0.06))
        cx = tip_x - ux * tip_back
        cy = tip_y - uy * tip_back
        return (cx + nx * tip_half, cy + ny * tip_half), (cx - nx * tip_half, cy - ny * tip_half)

    def _extend_line_endpoints(
        self,
        p1x: float,
        p1y: float,
        p2x: float,
        p2y: float,
        ux: float,
        uy: float,
        trim_start: float,
        trim_end: float,
        pen_width: float,
    ) -> tuple[float, float, float, float]:
        """Extiende extremos de línea para trazos planos.

        Args:
            p1x: Coordenada X inicial.
            p1y: Coordenada Y inicial.
            p2x: Coordenada X final.
            p2y: Coordenada Y final.
            ux: Componente X del vector unitario.
            uy: Componente Y del vector unitario.
            trim_start: Recorte desde el inicio.
            trim_end: Recorte desde el final.
            pen_width: Grosor del trazo actual.

        Returns:
            Resultado calculado.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        if self._style.cap_style != Qt.PenCapStyle.FlatCap:
            return p1x, p1y, p2x, p2y
        extend_start = self._endpoint_extend_start
        extend_end = self._endpoint_extend_end
        if extend_start <= 0.0 and extend_end <= 0.0:
            return p1x, p1y, p2x, p2y
        if trim_start <= 0.0:
            p1x -= ux * extend_start
            p1y -= uy * extend_start
        if trim_end <= 0.0:
            p2x += ux * extend_end
            p2y += uy * extend_end
        return p1x, p1y, p2x, p2y

    def update_positions(self, atom1: Atom, atom2: Atom) -> None:
        """Actualiza posiciones.

        Args:
            atom1: Átomo inicial del enlace.
            atom2: Átomo final del enlace.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        x1, y1 = atom1.x, atom1.y
        x2, y2 = atom2.x, atom2.y
        path = QPainterPath()
        color = QColor(self._style.bond_color)
        if self._color:
            candidate = QColor(self._color)
            if candidate.isValid():
                color = candidate
        dx = x2 - x1
        dy = y2 - y1
        length = math.hypot(dx, dy)
        if length <= 1e-6:
            self.setPath(path)
            return
        nx = -dy / length  # Normal vector perpendicular to bond
        ny = dx / length
        ux = dx / length
        uy = dy / length
        render_length = length
        p1x, p1y = x1, y1
        p2x, p2y = x2, y2
        trim_start = self._label_shrink_start + self._endpoint_trim_start
        trim_end = self._label_shrink_end + self._endpoint_trim_end
        if trim_start + trim_end > 0:
            min_length = max(1.0, self._style.stroke_px)
            max_trim = max(0.0, render_length - min_length)
            if trim_start + trim_end > max_trim and max_trim > 0:
                scale = max_trim / (trim_start + trim_end)
                trim_start *= scale
                trim_end *= scale
            p1x += ux * trim_start
            p1y += uy * trim_start
            p2x -= ux * trim_end
            p2y -= uy * trim_end
            render_length = math.hypot(p2x - p1x, p2y - p1y)

        offset_sign = self._offset_sign
        if self._ring_center is not None:
            midx = (x1 + x2) / 2
            midy = (y1 + y2) / 2
            vx = self._ring_center.x() - midx
            vy = self._ring_center.y() - midy
            offset_sign = 1 if (nx * vx + ny * vy) >= 0 else -1

        # Aromatic bonds: if circle mode, draw as single line
        if self.is_aromatic and self.render_aromatic_as_circle:
            stroke_px = self._stroke_px if self._stroke_px is not None else self._style.stroke_px
            e1x, e1y, e2x, e2y = self._extend_line_endpoints(
                p1x, p1y, p2x, p2y, ux, uy, trim_start, trim_end, stroke_px
            )
            path.moveTo(e1x, e1y)
            path.lineTo(e2x, e2y)
            pen = QPen(color, stroke_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self.setPath(path)
            return

        effective_order = self.order
        if self.is_aromatic and self.display_order is not None:
            effective_order = self.display_order

        stroke_px = self._stroke_px if self._stroke_px is not None else self._style.stroke_px
        stroke_scale = stroke_px / self._style.stroke_px if self._style.stroke_px > 1e-6 else 1.0

        if self.style == BondStyle.PLAIN:
            e1x, e1y, e2x, e2y = self._extend_line_endpoints(
                p1x, p1y, p2x, p2y, ux, uy, trim_start, trim_end, stroke_px
            )
            if effective_order == 1:
                path.moveTo(e1x, e1y)
                path.lineTo(e2x, e2y)
            else:
                offset = self._style.double_offset_px
                if effective_order == 2:
                    if self._symmetric_double:
                        half_offset = offset * 0.5
                        path.moveTo(e1x + nx * half_offset, e1y + ny * half_offset)
                        path.lineTo(e2x + nx * half_offset, e2y + ny * half_offset)
                        path.moveTo(e1x - nx * half_offset, e1y - ny * half_offset)
                        path.lineTo(e2x - nx * half_offset, e2y - ny * half_offset)
                    else:
                        path.moveTo(e1x, e1y)
                        path.lineTo(e2x, e2y)
                        use_inner_trim = (not self._symmetric_double) and (not self._prefer_full_length)
                        q1x = e1x + nx * offset * offset_sign
                        q1y = e1y + ny * offset * offset_sign
                        q2x = e2x + nx * offset * offset_sign
                        q2y = e2y + ny * offset * offset_sign
                        if use_inner_trim:
                            q1x += ux * self._style.inner_trim_px
                            q1y += uy * self._style.inner_trim_px
                            q2x -= ux * self._style.inner_trim_px
                            q2y -= uy * self._style.inner_trim_px
                        path.moveTo(q1x, q1y)
                        path.lineTo(q2x, q2y)
                else:  # Triple bond
                    path.moveTo(e1x, e1y)
                    path.lineTo(e2x, e2y)
                    path.moveTo(e1x + nx * offset, e1y + ny * offset)
                    path.lineTo(e2x + nx * offset, e2y + ny * offset)
                    path.moveTo(e1x - nx * offset, e1y - ny * offset)
                    path.lineTo(e2x - nx * offset, e2y - ny * offset)
            pen = QPen(color, stroke_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            
        elif self.style == BondStyle.BOLD:
            bold_px = max(stroke_px * 2.2, stroke_px + 1.0)
            e1x, e1y, e2x, e2y = self._extend_line_endpoints(
                p1x, p1y, p2x, p2y, ux, uy, trim_start, trim_end, bold_px
            )
            path.moveTo(e1x, e1y)
            path.lineTo(e2x, e2y)
            pen = QPen(color, bold_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))

        elif self.style == BondStyle.INTERACTION:
            e1x, e1y, e2x, e2y = self._extend_line_endpoints(
                p1x, p1y, p2x, p2y, ux, uy, trim_start, trim_end, stroke_px
            )
            path.moveTo(e1x, e1y)
            path.lineTo(e2x, e2y)
            pen = QPen(color, stroke_px, Qt.PenStyle.DotLine)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))

        elif self.style == BondStyle.WEDGE:
            # ChemDraw-like wedge proportions: moderate base width and
            # sub-linear growth when custom stroke is increased.
            width = self._style.wedge_width_px * (0.72 + 0.28 * math.sqrt(max(stroke_scale, 1e-6)))
            width = max(width, stroke_px * 2.3)
            width = min(width, render_length * 0.34)
            wedge_trim_start = 0.0 if self._bond_in_ring else trim_start
            wedge_trim_end = 0.0 if self._bond_in_ring else trim_end
            tip, base1, base2 = compute_wedge_points(
                (x1, y1),
                (x2, y2),
                width,
                trim_start=wedge_trim_start,
                trim_end=wedge_trim_end,
            )
            base_cx = (base1[0] + base2[0]) * 0.5
            base_cy = (base1[1] + base2[1]) * 0.5
            tip_cx, tip_cy = tip
            wedge_dx = base_cx - tip_cx
            wedge_dy = base_cy - tip_cy
            wedge_len = math.hypot(wedge_dx, wedge_dy)
            if wedge_len <= 1e-6:
                self.setPath(path)
                return
            w_ux = wedge_dx / wedge_len
            w_uy = wedge_dy / wedge_len
            w_nx = -w_uy
            w_ny = w_ux
            half_w = width * 0.5

            # Tip adapts to incoming bond(s): point, short edge (trapezoid),
            # or elongated arrow-like profile depending on connectivity.
            tip_pos, tip_neg = self._adapt_wedge_tip_for_neighbors(
                tip_cx,
                tip_cy,
                w_ux,
                w_uy,
                w_nx,
                w_ny,
                width,
                stroke_px,
                self._wedge_join_start,
            )
            if wedge_trim_start <= 1e-6 and self._style.cap_style == Qt.PenCapStyle.RoundCap:
                tip_overlap = min(stroke_px * 0.18, width * 0.035)
                if tip_overlap > 0.0:
                    tip_pos = (tip_pos[0] - w_ux * tip_overlap, tip_pos[1] - w_uy * tip_overlap)
                    tip_neg = (tip_neg[0] - w_ux * tip_overlap, tip_neg[1] - w_uy * tip_overlap)

            # Base adaptation: each side can bend toward a neighbor bond edge.
            base_pos_default = (base_cx + w_nx * half_w, base_cy + w_ny * half_w)
            base_neg_default = (base_cx - w_nx * half_w, base_cy - w_ny * half_w)
            base_pos = self._adapt_wedge_corner_for_neighbors(
                base_cx,
                base_cy,
                w_ux,
                w_uy,
                w_nx,
                w_ny,
                1.0,
                base_pos_default,
                self._wedge_join_end,
                stroke_px,
            )
            base_neg = self._adapt_wedge_corner_for_neighbors(
                base_cx,
                base_cy,
                w_ux,
                w_uy,
                w_nx,
                w_ny,
                -1.0,
                base_neg_default,
                self._wedge_join_end,
                stroke_px,
            )
            # Prevent corner inversion around the base axis.
            pos_proj = (base_pos[0] - base_cx) * w_nx + (base_pos[1] - base_cy) * w_ny
            neg_proj = (base_neg[0] - base_cx) * w_nx + (base_neg[1] - base_cy) * w_ny
            if pos_proj < neg_proj:
                base_pos, base_neg = base_neg, base_pos
            # Keep a minimum terminal width while allowing strong integration.
            end_deg = len(self._wedge_join_end)
            anchored_corner = 0  # 1 => base_pos locked, -1 => base_neg locked
            single_end_dir: tuple[float, float] | None = None
            if end_deg == 1:
                # Lock only the corner that points toward the outgoing bond
                # to the shared atom vertex (tiny overlap for a seamless join).
                nux, nuy, nwidth = self._wedge_join_end[0]
                single_end_dir = (nux, nuy)
                overlap = min(
                    stroke_px * WEDGE_SINGLE_END_OVERLAP_STROKE_MULT,
                    width * WEDGE_SINGLE_END_OVERLAP_WIDTH_MULT,
                )
                overlap *= max(0.0, WEDGE_SINGLE_END_JOIN_TUNE)
                # Snap to the nearest edge of the outgoing simple-bond stroke
                # (not to its centerline) to eliminate tiny terminal nubs.
                joint_cx = base_cx + nux * overlap
                joint_cy = base_cy + nuy * overlap
                nnx, nny = -nuy, nux
                half_neighbor = max(nwidth * 0.5, stroke_px * 0.50)
                edge_a = (joint_cx + nnx * half_neighbor, joint_cy + nny * half_neighbor)
                edge_b = (joint_cx - nnx * half_neighbor, joint_cy - nny * half_neighbor)
                pos_along = (base_pos[0] - base_cx) * nux + (base_pos[1] - base_cy) * nuy
                neg_along = (base_neg[0] - base_cx) * nux + (base_neg[1] - base_cy) * nuy
                if pos_along >= neg_along:
                    d_a = (base_pos[0] - edge_a[0]) ** 2 + (base_pos[1] - edge_a[1]) ** 2
                    d_b = (base_pos[0] - edge_b[0]) ** 2 + (base_pos[1] - edge_b[1]) ** 2
                    base_pos = edge_a if d_a <= d_b else edge_b
                    anchored_corner = 1
                else:
                    d_a = (base_neg[0] - edge_a[0]) ** 2 + (base_neg[1] - edge_a[1]) ** 2
                    d_b = (base_neg[0] - edge_b[0]) ** 2 + (base_neg[1] - edge_b[1]) ** 2
                    base_neg = edge_a if d_a <= d_b else edge_b
                    anchored_corner = -1
            sep = (base_pos[0] - base_neg[0]) * w_nx + (base_pos[1] - base_neg[1]) * w_ny
            if end_deg >= 2:
                min_sep = width * WEDGE_MULTI_END_MIN_SEP_MULT
            elif end_deg == 1:
                min_sep = width * WEDGE_SINGLE_END_MIN_SEP_MULT
            else:
                min_sep = width * 0.50
            if abs(sep) < min_sep:
                corr = (min_sep - abs(sep)) * 0.5
                sign = 1.0 if sep >= 0.0 else -1.0
                if end_deg == 1 and anchored_corner != 0:
                    # Preserve the bonded corner exactly on the atom/bond join;
                    # widen only the opposite corner to avoid tiny protrusions.
                    full_corr = min_sep - abs(sep)
                    if anchored_corner > 0:
                        base_neg = (
                            base_neg[0] - w_nx * full_corr * sign,
                            base_neg[1] - w_ny * full_corr * sign,
                        )
                    else:
                        base_pos = (
                            base_pos[0] + w_nx * full_corr * sign,
                            base_pos[1] + w_ny * full_corr * sign,
                        )
                else:
                    base_pos = (base_pos[0] + w_nx * corr * sign, base_pos[1] + w_ny * corr * sign)
                    base_neg = (base_neg[0] - w_nx * corr * sign, base_neg[1] - w_ny * corr * sign)
            if end_deg == 1 and anchored_corner != 0 and single_end_dir is not None:
                # Enforce a short flat terminal base (trapezoid), using only the
                # two base corners and no extra midpoint to avoid collapses.
                nux, nuy = single_end_dir
                if anchored_corner > 0:
                    ax, ay = base_pos
                else:
                    ax, ay = base_neg
                current_span = math.hypot(base_pos[0] - base_neg[0], base_pos[1] - base_neg[1])
                cap_len = min(
                    stroke_px * max(0.0, WEDGE_SINGLE_END_TERMINAL_CLIP_STROKE_MULT),
                    width * max(0.0, WEDGE_SINGLE_END_TERMINAL_CLIP_WIDTH_MULT),
                    current_span * 0.85,
                )
                if cap_len > 0.0:
                    bx = ax - nux * cap_len
                    by = ay - nuy * cap_len
                    if anchored_corner > 0:
                        base_neg = (bx, by)
                    else:
                        base_pos = (bx, by)
                cap_overlay = min(
                    stroke_px * max(0.0, WEDGE_SINGLE_END_CAP_OVERLAY_STROKE_MULT),
                    width * max(0.0, WEDGE_SINGLE_END_CAP_OVERLAY_WIDTH_MULT),
                    current_span * max(0.0, WEDGE_SINGLE_END_CAP_OVERLAY_SPAN_MULT),
                )
                if cap_overlay > 0.0:
                    half_overlay = cap_overlay * 0.5
                    base_pos = (base_pos[0] + nux * half_overlay, base_pos[1] + nuy * half_overlay)
                    base_neg = (base_neg[0] + nux * half_overlay, base_neg[1] + nuy * half_overlay)
            end_extension = 0.0
            if end_deg >= 2 and self._wedge_join_end:
                # Re-anchor each wide-end corner toward the bond on its side so
                # large terminal width still looks connected, not detached.
                pos_neighbor = None
                neg_neighbor = None
                pos_best = -1e9
                neg_best = -1e9
                any_best = -1e9
                any_neighbor = None
                for nux, nuy, nwidth in self._wedge_join_end:
                    cross = w_ux * nuy - w_uy * nux
                    dot = w_ux * nux + w_uy * nuy
                    score = abs(cross) + max(0.0, dot) * 0.25
                    if score > any_best:
                        any_best = score
                        any_neighbor = (nux, nuy, nwidth)
                    if cross > 0.03 and score > pos_best:
                        pos_best = score
                        pos_neighbor = (nux, nuy, nwidth)
                    if cross < -0.03 and score > neg_best:
                        neg_best = score
                        neg_neighbor = (nux, nuy, nwidth)
                if pos_neighbor is None:
                    pos_neighbor = any_neighbor
                if neg_neighbor is None:
                    neg_neighbor = any_neighbor
                if pos_neighbor is not None:
                    p_nux, p_nuy, _ = pos_neighbor
                    p_overlap = min(
                        stroke_px * WEDGE_MULTI_END_CORNER_JOIN_STROKE_MULT,
                        width * WEDGE_MULTI_END_CORNER_JOIN_WIDTH_MULT,
                    )
                    p_joint = (base_cx + p_nux * p_overlap, base_cy + p_nuy * p_overlap)
                    base_pos = (
                        base_pos[0] * (1.0 - WEDGE_MULTI_END_CORNER_BLEND)
                        + p_joint[0] * WEDGE_MULTI_END_CORNER_BLEND,
                        base_pos[1] * (1.0 - WEDGE_MULTI_END_CORNER_BLEND)
                        + p_joint[1] * WEDGE_MULTI_END_CORNER_BLEND,
                    )
                if neg_neighbor is not None:
                    n_nux, n_nuy, _ = neg_neighbor
                    n_overlap = min(
                        stroke_px * WEDGE_MULTI_END_CORNER_JOIN_STROKE_MULT,
                        width * WEDGE_MULTI_END_CORNER_JOIN_WIDTH_MULT,
                    )
                    n_joint = (base_cx + n_nux * n_overlap, base_cy + n_nuy * n_overlap)
                    base_neg = (
                        base_neg[0] * (1.0 - WEDGE_MULTI_END_CORNER_BLEND)
                        + n_joint[0] * WEDGE_MULTI_END_CORNER_BLEND,
                        base_neg[1] * (1.0 - WEDGE_MULTI_END_CORNER_BLEND)
                        + n_joint[1] * WEDGE_MULTI_END_CORNER_BLEND,
                    )
                # Re-apply terminal width floor after corner re-anchoring so the
                # wide end remains conical and does not collapse.
                sep = (base_pos[0] - base_neg[0]) * w_nx + (base_pos[1] - base_neg[1]) * w_ny
                min_sep = width * WEDGE_MULTI_END_MIN_SEP_MULT
                if abs(sep) < min_sep:
                    corr = (min_sep - abs(sep)) * 0.5
                    sign = 1.0 if sep >= 0.0 else -1.0
                    base_pos = (base_pos[0] + w_nx * corr * sign, base_pos[1] + w_ny * corr * sign)
                    base_neg = (base_neg[0] - w_nx * corr * sign, base_neg[1] - w_ny * corr * sign)

                # Keep existing fork shape, but extend each wide-end corner a
                # little along its neighbor direction to close tiny visual gaps.
                contact_push = min(
                    stroke_px * WEDGE_MULTI_END_CONTACT_STROKE_MULT,
                    width * WEDGE_MULTI_END_CONTACT_WIDTH_MULT,
                )
                if contact_push > 0.0:
                    if pos_neighbor is not None:
                        p_nux, p_nuy, _ = pos_neighbor
                        base_pos = (
                            base_pos[0] + p_nux * contact_push,
                            base_pos[1] + p_nuy * contact_push,
                        )
                    if neg_neighbor is not None:
                        n_nux, n_nuy, _ = neg_neighbor
                        base_neg = (
                            base_neg[0] + n_nux * contact_push,
                            base_neg[1] + n_nuy * contact_push,
                        )

                # True longitudinal extension: moves the wide terminal outward
                # along the wedge axis (does not change atom positions).
                end_extension = min(
                    stroke_px * WEDGE_MULTI_END_EXTENSION_STROKE_MULT,
                    width * WEDGE_MULTI_END_EXTENSION_WIDTH_MULT,
                )
                if end_extension > 0.0:
                    base_pos = (
                        base_pos[0] + w_ux * end_extension,
                        base_pos[1] + w_uy * end_extension,
                    )
                    base_neg = (
                        base_neg[0] + w_ux * end_extension,
                        base_neg[1] + w_uy * end_extension,
                    )

            fork_point: tuple[float, float] | None = None
            if end_deg >= 2:
                # Add a small inward "fork" at the wide end to mimic ChemDraw's
                # terminal bifurcation when wedge meets two simple bonds.
                depth = min(
                    stroke_px * WEDGE_MULTI_END_FORK_DEPTH_STROKE_MULT,
                    width * WEDGE_MULTI_END_FORK_DEPTH_WIDTH_MULT,
                )
                mid_x = (base_pos[0] + base_neg[0]) * 0.5
                mid_y = (base_pos[1] + base_neg[1]) * 0.5
                # Move fork inward (toward the narrow tip), not outward.
                target_x = base_cx - w_ux * depth
                target_y = base_cy - w_uy * depth
                fork_point = (
                    mid_x * (1.0 - WEDGE_MULTI_END_FORK_BLEND) + target_x * WEDGE_MULTI_END_FORK_BLEND,
                    mid_y * (1.0 - WEDGE_MULTI_END_FORK_BLEND) + target_y * WEDGE_MULTI_END_FORK_BLEND,
                )
                if end_extension > 0.0:
                    fork_point = (
                        fork_point[0] + w_ux * end_extension,
                        fork_point[1] + w_uy * end_extension,
                    )

            path.moveTo(tip_pos[0], tip_pos[1])
            path.lineTo(base_pos[0], base_pos[1])
            if fork_point is not None:
                path.lineTo(fork_point[0], fork_point[1])
            path.lineTo(base_neg[0], base_neg[1])
            path.lineTo(tip_neg[0], tip_neg[1])
            path.closeSubpath()
            # Fill-only wedge avoids a 1px outline that makes the base look bulky.
            self.setPen(QPen(Qt.PenStyle.NoPen))
            self.setBrush(QBrush(color))
            
        elif self.style == BondStyle.HASHED:
            steps = self._style.hash_count
            if self._bond_in_ring:
                use_trim_start = trim_start if trim_start > 0 else 0.0
                use_trim_end = trim_end if trim_end > 0 else 0.0
                if use_trim_start > 0.0 or use_trim_end > 0.0:
                    h1x, h1y, h2x, h2y = p1x, p1y, p2x, p2y
                else:
                    h1x, h1y, h2x, h2y = x1, y1, x2, y2
                h_trim_start = use_trim_start
                h_trim_end = use_trim_end
            else:
                h1x, h1y, h2x, h2y = p1x, p1y, p2x, p2y
                h_trim_start = trim_start
                h_trim_end = trim_end
            e1x, e1y, e2x, e2y = self._extend_line_endpoints(
                h1x, h1y, h2x, h2y, ux, uy, h_trim_start, h_trim_end, stroke_px
            )
            for i in range(1, steps + 1):
                t = i / steps
                px = e1x + (e2x - e1x) * t
                py = e1y + (e2y - e1y) * t
                width = (
                    self._style.hash_min_px
                    + (self._style.hash_max_px - self._style.hash_min_px) * t
                ) * stroke_scale
                path.moveTo(px + nx * width / 2, py + ny * width / 2)
                path.lineTo(px - nx * width / 2, py - ny * width / 2)
            hash_stroke = max(self._style.hash_stroke_px * stroke_scale, stroke_px * 0.85)
            pen = QPen(color, hash_stroke)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            
        elif self.style == BondStyle.WAVY:
            wavelength = max(8.0, stroke_px * 3.5)
            amplitude = max(stroke_px * 0.9, wavelength * 0.28)
            e1x, e1y, e2x, e2y = self._extend_line_endpoints(
                p1x, p1y, p2x, p2y, ux, uy, trim_start, trim_end, stroke_px
            )
            path = _build_wavy_path(
                QPointF(e1x, e1y),
                QPointF(e2x, e2y),
                amplitude,
                wavelength,
            )
            pen = QPen(color, stroke_px)
            pen.setCapStyle(self._style.cap_style)
            pen.setJoinStyle(self._style.join_style)
            self.setPen(pen)
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))

        self.setPath(path)

    def set_style(self, style: DrawingStyle, atom1: Atom, atom2: Atom) -> None:
        """Actualiza estilo.

        Args:
            style: Estilo de dibujo aplicado.
            atom1: Átomo inicial del enlace.
            atom2: Átomo final del enlace.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._style = style
        self.update_positions(atom1, atom2)


class WavyAnchorItem(QGraphicsPathItem):
    """Wavy anchor stub attached to an atom."""

    def __init__(
        self,
        start: QPointF,
        end: QPointF,
        style: DrawingStyle = CHEMDOODLE_LIKE,
    ) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Args:
            start: Punto inicial en escena.
            end: Punto final en escena.
            style: Estilo de dibujo aplicado.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__()
        self._style = style
        self._start = QPointF(start)
        self._end = QPointF(end)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setZValue(30)
        self.update_positions(start, end)

    def update_positions(self, start: QPointF, end: QPointF) -> None:
        """Actualiza posiciones.

        Args:
            start: Punto inicial en escena.
            end: Punto final en escena.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._start = QPointF(start)
        self._end = QPointF(end)

        dx = self._end.x() - self._start.x()
        dy = self._end.y() - self._start.y()
        length = math.hypot(dx, dy)
        if length <= 1e-6:
            length = 1.0
            dx = 1.0
            dy = 0.0

        wavelength = max(8.0, self._style.stroke_px * 3.5)
        amplitude = max(self._style.stroke_px * 0.9, wavelength * 0.28)
        path = _build_wavy_path(self._start, self._end, amplitude, wavelength)

        pen = QPen(QColor(self._style.bond_color), self._style.stroke_px)
        pen.setCapStyle(self._style.cap_style)
        pen.setJoinStyle(self._style.join_style)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setPath(path)

    def start_point(self) -> QPointF:
        """Devuelve el punto inicial.

        Returns:
            Punto inicial (QPointF).

        Side Effects:
            No tiene efectos laterales.
        """
        return QPointF(self._start)

    def end_point(self) -> QPointF:
        """Devuelve el punto final.

        Returns:
            Punto final (QPointF).

        Side Effects:
            No tiene efectos laterales.
        """
        return QPointF(self._end)


class ArrowItem(QGraphicsPathItem):
    """Elemento gráfico que representa una flecha de reacción/anotación."""

    def __init__(
        self,
        start: QPointF,
        end: QPointF,
        head_at_end: bool = True,
        kind: str | None = None,
        style: DrawingStyle = CHEMDOODLE_LIKE,
    ) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Args:
            start: Punto inicial en escena.
            end: Punto final en escena.
            head_at_end: Si la cabeza va en el extremo final.
            kind: Tipo de flecha/corchete.
            style: Estilo de dibujo aplicado.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__()
        self._style = style
        if kind is None:
            self._kind = "forward" if head_at_end else "retro"
        else:
            self._kind = kind
        self._start = QPointF(start)
        self._end = QPointF(end)
        pen = QPen(QColor(self._style.bond_color), self._style.stroke_px)
        pen.setCapStyle(self._style.cap_style)
        pen.setJoinStyle(self._style.join_style)
        self.setPen(pen)
        self.setBrush(QBrush(QColor(self._style.bond_color)))
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setZValue(5)
        self.update_positions(start, end)

    def paint(self, painter, option, widget=None) -> None:
        """Dibuja el elemento en la escena con el estilo actual.

        Args:
            painter: Pintor Qt usado para dibujar.
            option: Opciones de estilo del render.
            widget: Widget contenedor opcional.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        painter.setPen(self.pen())
        painter.setBrush(self.brush())
        painter.drawPath(self.path())

    def update_positions(self, start: QPointF, end: QPointF) -> None:
        """Actualiza posiciones.

        Args:
            start: Punto inicial en escena.
            end: Punto final en escena.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._start = QPointF(start)
        self._end = QPointF(end)
        dx = end.x() - start.x()
        dy = end.y() - start.y()
        length = math.hypot(dx, dy)
        if length < 1e-6:
            self.setPath(QPainterPath())
            return

        ux = dx / length
        uy = dy / length
        head_len = 12.0
        head_width = 6.0
        nx = -uy
        ny = ux

        open_head_kinds = {
            "forward_open",
            "retro_open",
            "both_open",
            "equilibrium_open",
            "retrosynthetic",
        }
        dashed_kinds = {
            "forward_dashed",
            "retro_dashed",
            "both_dashed",
            "equilibrium_dashed",
        }
        curved_kinds = {"curved", "curved_fishhook"}

        is_open = self._kind in open_head_kinds
        is_dashed = self._kind in dashed_kinds
        is_curved = self._kind in curved_kinds
        is_fishhook = self._kind == "curved_fishhook"

        self._apply_kind_style(is_open or is_fishhook, is_dashed)

        def add_head(
            path: QPainterPath,
            tip: QPointF,
            dir_x: float,
            dir_y: float,
            head_style: str,
        ) -> QPointF:
            """Dibuja la cabeza de la flecha.

            Args:
                path: Ruta QPainterPath donde se dibuja.
                tip: Punto de la punta de la flecha.
                dir_x: Componente X de la dirección.
                dir_y: Componente Y de la dirección.
                head_style: Estilo de cabeza (filled/open/half).

            Returns:
                None.

            Side Effects:
                Modifica el estado del item o la escena.
            """
            base_x = tip.x() - dir_x * head_len
            base_y = tip.y() - dir_y * head_len
            left = QPointF(base_x + nx * head_width, base_y + ny * head_width)
            right = QPointF(base_x - nx * head_width, base_y - ny * head_width)
            if head_style == "filled":
                path.moveTo(left)
                path.lineTo(tip)
                path.lineTo(right)
                path.closeSubpath()
            elif head_style == "open":
                path.moveTo(left)
                path.lineTo(tip)
                path.lineTo(right)
            elif head_style == "half":
                path.moveTo(left)
                path.lineTo(tip)
            return QPointF(base_x, base_y)

        path = QPainterPath()
        head_style = "half" if is_fishhook else ("open" if is_open else "filled")

        if self._kind in {"both", "both_open", "both_dashed"}:
            start_base = QPointF(start.x() + ux * head_len, start.y() + uy * head_len)
            end_base = QPointF(end.x() - ux * head_len, end.y() - uy * head_len)
            path.moveTo(start_base)
            path.lineTo(end_base)
            add_head(path, end, ux, uy, head_style)
            add_head(path, start, -ux, -uy, head_style)
        elif self._kind in {"equilibrium", "equilibrium_dashed"}:
            offset = 4.0
            top_start = QPointF(start.x() + nx * offset, start.y() + ny * offset)
            top_end = QPointF(end.x() + nx * offset, end.y() + ny * offset)
            bottom_start = QPointF(start.x() - nx * offset, start.y() - ny * offset)
            bottom_end = QPointF(end.x() - nx * offset, end.y() - ny * offset)

            top_base = QPointF(top_end.x() - ux * head_len, top_end.y() - uy * head_len)
            path.moveTo(top_start)
            path.lineTo(top_base)
            add_head(path, top_end, ux, uy, head_style)

            bottom_base = QPointF(bottom_start.x() + ux * head_len, bottom_start.y() + uy * head_len)
            path.moveTo(bottom_end)
            path.lineTo(bottom_base)
            add_head(path, bottom_start, -ux, -uy, head_style)
        elif self._kind == "retrosynthetic":
            offset = 2.5
            tip = end
            tail = start
            dir_x = ux
            dir_y = uy
            base = add_head(path, tip, dir_x, dir_y, "open")
            for sign in (-1, 1):
                path.moveTo(QPointF(tail.x() + nx * offset * sign, tail.y() + ny * offset * sign))
                path.lineTo(QPointF(base.x() + nx * offset * sign, base.y() + ny * offset * sign))
        elif is_curved:
            curve_offset = max(12.0, min(24.0, length * 0.25))
            control = QPointF(
                (start.x() + end.x()) * 0.5 + nx * curve_offset,
                (start.y() + end.y()) * 0.5 + ny * curve_offset,
            )
            tx = end.x() - control.x()
            ty = end.y() - control.y()
            tlen = math.hypot(tx, ty)
            if tlen < 1e-6:
                tx = ux
                ty = uy
            else:
                tx /= tlen
                ty /= tlen
            base = add_head(path, end, tx, ty, head_style)
            path.moveTo(start)
            path.quadTo(control, base)
        else:
            head_at_end = self._kind != "retro"
            if head_at_end:
                tip = end
                tail = start
                dir_x = ux
                dir_y = uy
            else:
                tip = start
                tail = end
                dir_x = -ux
                dir_y = -uy
            base = add_head(path, tip, dir_x, dir_y, head_style)
            path.moveTo(tail)
            path.lineTo(base)

        self.setPath(path)

    def kind(self) -> str:
        """Devuelve el tipo actual.

        Returns:
            Tipo de flecha/corchete.

        Side Effects:
            No tiene efectos laterales.
        """
        return self._kind

    def set_kind(self, kind: str) -> None:
        """Actualiza tipo.

        Args:
            kind: Tipo de flecha/corchete.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._kind = kind

    def start_point(self) -> QPointF:
        """Devuelve el punto inicial.

        Returns:
            Punto inicial (QPointF).

        Side Effects:
            No tiene efectos laterales.
        """
        return QPointF(self._start)

    def end_point(self) -> QPointF:
        """Devuelve el punto final.

        Returns:
            Punto final (QPointF).

        Side Effects:
            No tiene efectos laterales.
        """
        return QPointF(self._end)

    def set_style(self, style: DrawingStyle) -> None:
        """Actualiza estilo.

        Args:
            style: Estilo de dibujo aplicado.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._style = style
        self.update_positions(self._start, self._end)

    def _apply_kind_style(self, open_head: bool, dashed: bool) -> None:
        """Aplica el estilo según el tipo.

        Args:
            open_head: Si la cabeza es abierta.
            dashed: Si el trazo es discontinuo.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        pen = QPen(QColor(self._style.bond_color), self._style.stroke_px)
        pen.setCapStyle(self._style.cap_style)
        pen.setJoinStyle(self._style.join_style)
        if dashed:
            pen.setStyle(Qt.PenStyle.DashLine)
        self.setPen(pen)
        if open_head:
            self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        else:
            self.setBrush(QBrush(QColor(self._style.bond_color)))


class PreviewArrowItem(ArrowItem):
    """Flecha de previsualización para colocar anotaciones."""

    def __init__(self) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__(QPointF(0, 0), QPointF(1, 0), kind="forward")
        self._preview_pen = QPen(QColor("#4A90D9"), 1.5, Qt.PenStyle.DashLine)
        self._preview_pen.setCapStyle(self._style.cap_style)
        self._preview_pen.setJoinStyle(self._style.join_style)
        self.setPen(self._preview_pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable, False)
        self.setZValue(40)
        self.setVisible(False)

    def update_preview(self, start: QPointF, end: QPointF, kind: str) -> None:
        """Actualiza la previsualización.

        Args:
            start: Punto inicial en escena.
            end: Punto final en escena.
            kind: Tipo de flecha/corchete.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.set_kind(kind)
        self.update_positions(start, end)
        self.setPen(self._preview_pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        if not self.isVisible():
            self.setVisible(True)

    def hide_preview(self) -> None:
        """Oculta la previsualización.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.setVisible(False)

    def _apply_kind_style(self, open_head: bool, dashed: bool) -> None:
        """Aplica el estilo según el tipo.

        Args:
            open_head: Si la cabeza es abierta.
            dashed: Si el trazo es discontinuo.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        if not hasattr(self, "_preview_pen"):
            super()._apply_kind_style(open_head, dashed)
            return
        self.setPen(self._preview_pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))


class BracketItem(QGraphicsPathItem):
    """Elemento gráfico que representa corchetes alrededor de una región."""

    def __init__(
        self,
        rect: QRectF,
        kind: str = "[]",
        padding: float = 8.0,
        style: DrawingStyle = CHEMDOODLE_LIKE,
    ) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Args:
            rect: Rectángulo base de la anotación.
            kind: Tipo de flecha/corchete.
            padding: Margen interno alrededor del contenido.
            style: Estilo de dibujo aplicado.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__()
        self._padding = padding
        self._base_rect = QRectF(rect)
        self._rect = rect.adjusted(-padding, -padding, padding, padding)
        self._kind = kind
        self._style = style
        pen = QPen(QColor(self._style.bond_color), self._style.stroke_px)
        pen.setCapStyle(self._style.cap_style)
        pen.setJoinStyle(self._style.join_style)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(2)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self._update_path()

    def paint(self, painter, option, widget=None) -> None:
        """Dibuja el elemento en la escena con el estilo actual.

        Args:
            painter: Pintor Qt usado para dibujar.
            option: Opciones de estilo del render.
            widget: Widget contenedor opcional.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        painter.setPen(self.pen())
        painter.setBrush(self.brush())
        painter.drawPath(self.path())

    def base_rect(self) -> QRectF:
        """Devuelve el rectángulo base.

        Returns:
            Rectángulo base de la anotación.

        Side Effects:
            No tiene efectos laterales.
        """
        return QRectF(self._base_rect)

    def set_rect(self, rect: QRectF, padding: float | None = None) -> None:
        """Actualiza rectángulo.

        Args:
            rect: Rectángulo base de la anotación.
            padding: Margen interno alrededor del contenido.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        if padding is not None:
            self._padding = padding
        self._base_rect = QRectF(rect)
        self._rect = self._base_rect.adjusted(
            -self._padding, -self._padding, self._padding, self._padding
        )
        self._update_path()

    def _update_path(self) -> None:
        """Recalcula la geometría del trazado.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        rect = self._rect
        path = QPainterPath()
        height = rect.height()
        top = rect.top()
        bottom = rect.bottom()
        left = rect.left()
        right = rect.right()

        if self._kind in {"()", "(", ")"}:
            width = max(8.0, height * 0.12)
            mid = (top + bottom) / 2
            if self._kind in {"()", "("}:
                path.moveTo(left + width, top)
                path.quadTo(left, top + height * 0.25, left, mid)
                path.quadTo(left, bottom - height * 0.25, left + width, bottom)
            if self._kind in {"()", ")"}:
                path.moveTo(right - width, top)
                path.quadTo(right, top + height * 0.25, right, mid)
                path.quadTo(right, bottom - height * 0.25, right - width, bottom)
        elif self._kind in {"{}", "{", "}"}:
            mid = (top + bottom) / 2
            total_width = max(10.0, height * 0.18)
            spine_offset = total_width * 0.45
            corner = min(height * 0.08, total_width * 0.4)
            # Smaller neck => sharper middle cusp
            neck = max(4.0, height * 0.05)

            def draw_curly_side(tip_x: float, spine_x: float, waist_x: float) -> None:
                """Dibuja un lateral de llave curva.

                Args:
                    tip_x: Coordenada X de la punta.
                    spine_x: Coordenada X del lomo de la llave.
                    waist_x: Coordenada X de la cintura de la llave.

                Returns:
                    None.

                Side Effects:
                    Modifica el estado del item o la escena.
                """
                path.moveTo(tip_x, top)
                path.quadTo(spine_x, top, spine_x, top + corner)
                path.lineTo(spine_x, mid - neck)
                path.cubicTo(
                    spine_x,
                    mid - neck * 0.25,
                    waist_x,
                    mid - neck * 0.05,
                    waist_x,
                    mid,
                )
                path.cubicTo(
                    waist_x,
                    mid + neck * 0.05,
                    spine_x,
                    mid + neck * 0.25,
                    spine_x,
                    mid + neck,
                )
                path.lineTo(spine_x, bottom - corner)
                path.quadTo(spine_x, bottom, tip_x, bottom)

            if self._kind in {"{}", "{"}:
                draw_curly_side(
                    tip_x=left + total_width,
                    spine_x=left + spine_offset,
                    waist_x=left,
                )
            if self._kind in {"{}", "}"}:
                draw_curly_side(
                    tip_x=right - total_width,
                    spine_x=right - spine_offset,
                    waist_x=right,
                )
        else:  # "[]", "[", "]"
            arm = max(6.0, height * 0.08)
            if self._kind in {"[]", "["}:
                path.moveTo(left + arm, top)
                path.lineTo(left, top)
                path.lineTo(left, bottom)
                path.lineTo(left + arm, bottom)
            if self._kind in {"[]", "]"}:
                path.moveTo(right - arm, top)
                path.lineTo(right, top)
                path.lineTo(right, bottom)
                path.lineTo(right - arm, bottom)

        self.setPath(path)

    def set_style(self, style: DrawingStyle) -> None:
        """Actualiza estilo.

        Args:
            style: Estilo de dibujo aplicado.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._style = style
        pen = QPen(QColor(self._style.bond_color), self._style.stroke_px)
        pen.setCapStyle(self._style.cap_style)
        pen.setJoinStyle(self._style.join_style)
        self.setPen(pen)

class AromaticCircleItem(QGraphicsEllipseItem):
    """Círculo interior que indica aromaticidad en anillos."""
    
    def __init__(self, center_x: float, center_y: float, radius: float) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Args:
            center_x: Coordenada X del centro.
            center_y: Coordenada Y del centro.
            radius: Radio visual en píxeles.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__(
            center_x - radius, center_y - radius,
            radius * 2, radius * 2
        )
        self.setPen(QPen(QColor("#333333"), 1.5))
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(-10)  # Behind bonds and atoms


class HoverAtomIndicatorItem(QGraphicsEllipseItem):
    """Indicador circular para átomos en estado hover."""

    def __init__(self, radius: float = 10.0) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Args:
            radius: Radio visual en píxeles.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__(-radius, -radius, radius * 2, radius * 2)
        self._radius = radius
        pen = QPen(QColor("#E0A825"), 1.5)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(50)
        self.setVisible(False)

    def update_position(self, x: float, y: float) -> None:
        """Actualiza posición.

        Args:
            x: Coordenada X.
            y: Coordenada Y.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.setPos(x, y)
        if not self.isVisible():
            self.setVisible(True)

    def hide_indicator(self) -> None:
        """Oculta el indicador.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.setVisible(False)


class HoverBondIndicatorItem(QGraphicsPathItem):
    """Indicador de paréntesis para enlaces en hover."""

    def __init__(self, radius: float = 10.0, separation: float = 12.0) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Args:
            radius: Radio visual en píxeles.
            separation: Separación entre paréntesis del indicador.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__()
        self._radius = radius
        self._separation = separation
        pen = QPen(QColor("#E0A825"), 1.5)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(50)
        self.setVisible(False)

    def update_for_bond(self, p1: QPointF, p2: QPointF) -> None:
        """Actualiza estado para enlace.

        Args:
            p1: Punto final.
            p2: Punto final del segmento.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        mid = QPointF((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2)
        angle = math.degrees(math.atan2(p2.y() - p1.y(), p2.x() - p1.x()))

        r = self._radius
        sep = self._separation
        path = QPainterPath()

        left_rect = QRectF(-sep - r, -r, 2 * r, 2 * r)
        right_rect = QRectF(sep - r, -r, 2 * r, 2 * r)
        path.arcMoveTo(left_rect, 60)
        path.arcTo(left_rect, 60, 240)
        path.arcMoveTo(right_rect, 120)
        path.arcTo(right_rect, 120, 240)

        self.setPath(path)
        self.setPos(mid)
        self.setRotation(angle)
        if not self.isVisible():
            self.setVisible(True)

    def hide_indicator(self) -> None:
        """Oculta el indicador.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.setVisible(False)


class OptimizeZoneItem(QGraphicsEllipseItem):
    """Blue translucent optimize zone around sprout anchor."""

    def __init__(self, radius: float = 28.0) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Args:
            radius: Radio visual en píxeles.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__(-radius, -radius, radius * 2, radius * 2)
        self._radius = radius
        pen = QPen(QColor("#4A90D9"), 1.2)
        brush = QBrush(QColor(74, 144, 217, 40))
        self.setPen(pen)
        self.setBrush(brush)
        self.setZValue(30)
        self.setVisible(False)

    def update_center(self, x: float, y: float) -> None:
        """Actualiza centro.

        Args:
            x: Coordenada X.
            y: Coordenada Y.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.setPos(x, y)
        if not self.isVisible():
            self.setVisible(True)

    def set_radius(self, radius: float) -> None:
        """Actualiza radio.

        Args:
            radius: Radio visual en píxeles.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self._radius = radius
        self.setRect(-radius, -radius, radius * 2, radius * 2)

    def radius(self) -> float:
        """Devuelve el radio actual.

        Returns:
            Radio actual en píxeles.

        Side Effects:
            No tiene efectos laterales.
        """
        return self._radius

    def hide_zone(self) -> None:
        """Oculta la zona de optimización.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.setVisible(False)


class PreviewBondItem(QGraphicsPathItem):
    """Línea de previsualización para colocar enlaces."""

    def __init__(self) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__()
        pen = QPen(QColor("#4A90D9"), 1.5, Qt.PenStyle.DashLine)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(40)
        self.setVisible(False)

    def update_line(self, p1: QPointF, p2: QPointF) -> None:
        """Actualiza línea.

        Args:
            p1: Punto final.
            p2: Punto final del segmento.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        path = QPainterPath()
        path.moveTo(p1)
        path.lineTo(p2)
        self.setPath(path)
        if not self.isVisible():
            self.setVisible(True)

    def hide_preview(self) -> None:
        """Oculta la previsualización.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.setVisible(False)


class PreviewRingItem(QGraphicsPathItem):
    """Polígono de previsualización para colocar anillos."""

    def __init__(self) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__()
        pen = QPen(QColor("#4A90D9"), 1.5, Qt.PenStyle.DashLine)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(40)
        self.setVisible(False)

    def update_polygon(self, vertices: list[QPointF]) -> None:
        """Actualiza polígono.

        Args:
            vertices: Lista de vértices del polígono.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        if not vertices:
            self.setVisible(False)
            return
        path = QPainterPath()
        path.moveTo(vertices[0])
        for v in vertices[1:]:
            path.lineTo(v)
        path.closeSubpath()
        self.setPath(path)
        if not self.isVisible():
            self.setVisible(True)

    def hide_preview(self) -> None:
        """Oculta la previsualización.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.setVisible(False)


class PreviewChainItem(QGraphicsPathItem):
    """Polilínea de previsualización para colocar cadenas."""

    def __init__(self) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__()
        pen = QPen(QColor("#4A90D9"), 1.5, Qt.PenStyle.DashLine)
        self.setPen(pen)
        self.setBrush(QBrush(Qt.BrushStyle.NoBrush))
        self.setZValue(40)
        self.setVisible(False)

    def update_polyline(self, points: list[QPointF]) -> None:
        """Actualiza polilínea.

        Args:
            points: Lista de puntos del trazo.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        if len(points) < 2:
            self.setVisible(False)
            return
        path = QPainterPath()
        path.moveTo(points[0])
        for p in points[1:]:
            path.lineTo(p)
        self.setPath(path)
        if not self.isVisible():
            self.setVisible(True)

    def hide_preview(self) -> None:
        """Oculta la previsualización.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.setVisible(False)


class PreviewChainLabelItem(QGraphicsTextItem):
    """Etiqueta de previsualización con la longitud de cadena."""

    def __init__(self) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__()
        self.setDefaultTextColor(QColor("#4A90D9"))
        font = QFont("Arial", 10, QFont.Weight.Bold)
        self.setFont(font)
        self.setZValue(41)
        self.setVisible(False)

    def update_label(self, text: str, pos: QPointF) -> None:
        """Actualiza etiqueta.

        Args:
            text: Texto a mostrar.
            pos: Posición en escena.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.setPlainText(text)
        rect = self.boundingRect()
        self.setPos(pos.x() - rect.width() / 2, pos.y() - rect.height() / 2)
        if not self.isVisible():
            self.setVisible(True)

    def hide_label(self) -> None:
        """Oculta la etiqueta.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        self.setVisible(False)




class TextAnnotationItem(QGraphicsTextItem):
    """Elemento gráfico para anotaciones de texto libres."""
    
    def __init__(self, text: str = "", x: float = 0.0, y: float = 0.0) -> None:
        """Inicializa la instancia y configura el elemento gráfico.

        Args:
            text: Texto a mostrar.
            x: Coordenada X.
            y: Coordenada Y.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        super().__init__(text)
        self.setPos(x, y)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsFocusable)
        
        # Default style
        font = QFont("Arial", 12)
        self.setFont(font)
        self.setDefaultTextColor(QColor("black"))
        
        # Editor interaction
        # Start in "NoInteraction" to allow dragging by default
        self.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        
        self.setZValue(10)
        
        
        # Drag Handle logic removed - handled by Canvas selection overlay
        # self.handle = TextDragHandle(self)

    def paint(self, painter, option, widget=None) -> None:
        # Suppress default dashed line from QGraphicsTextItem if selected
        # because the canvas draws its own selection overlay.
        # But we must preserve other states like HasFocus for the cursor.
        """Dibuja el elemento en la escena con el estilo actual.

        Args:
            painter: Pintor Qt usado para dibujar.
            option: Opciones de estilo del render.
            widget: Widget contenedor opcional.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        option.state &= ~QStyle.StateFlag.State_Selected
        
        super().paint(painter, option, widget)
        
        # If we are in edit mode, maybe we want a subtle background or border?
        # For now, let's keep it clean since the selection overlay is handled by canvas.
        if (self.textInteractionFlags() & Qt.TextInteractionFlag.TextEditorInteraction):
             # When editing, show a very subtle dashed border to indicate focus area
             painter.setPen(QPen(QColor("#CCCCCC"), 1.0, Qt.PenStyle.DashLine))
             painter.drawRect(self.boundingRect())

    def mouseDoubleClickEvent(self, event) -> None:
        # Enable editing on double click
        """Maneja el doble clic del mouse.

        Args:
            event: Evento de Qt recibido.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        if self.textInteractionFlags() == Qt.TextInteractionFlag.NoTextInteraction:
            self.setTextInteractionFlags(Qt.TextInteractionFlag.TextEditorInteraction)
            self.setFocus()
        super().mouseDoubleClickEvent(event)

    def focusInEvent(self, event) -> None:
        # Remember last focused text item for restore on window activation
        """Maneja la ganancia de foco.

        Args:
            event: Evento de Qt recibido.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        try:
            scene = self.scene()
            if scene:
                for view in scene.views():
                    if hasattr(view, "remember_text_edit_item"):
                        view.remember_text_edit_item(self)
        except Exception:
            pass
        super().focusInEvent(event)
    
    def focusOutEvent(self, event) -> None:
        # If focus loss is due to window deactivation, keep edit mode.
        """Maneja la pérdida de foco.

        Args:
            event: Evento de Qt recibido.

        Returns:
            None.

        Side Effects:
            Modifica el estado del item o la escena.
        """
        if event.reason() == Qt.FocusReason.ActiveWindowFocusReason:
            super().focusOutEvent(event)
            return

        # When losing focus, go back to movable mode
        self.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        
        # Clear any internal text selection to avoid lingering highlights
        # (the green background for selection)
        cursor = self.textCursor()
        cursor.clearSelection()
        self.setTextCursor(cursor)
        
        # If empty, delete the item
        content = self.toPlainText()
        if not content or content.strip() == "":
             if self.scene():
                 self.scene().removeItem(self)

        try:
            scene = self.scene()
            if scene:
                for view in scene.views():
                    if hasattr(view, "remember_text_edit_item"):
                        view.remember_text_edit_item(None)
        except Exception:
            pass
             
        super().focusOutEvent(event)

    def shape(self) -> QPainterPath:
        # Ensure a minimum size for the hit-test even if empty
        """Devuelve la forma para hit-testing.

        Returns:
            QPainterPath usado para detección de colisiones.

        Side Effects:
            No tiene efectos laterales.
        """
        br = self.boundingRect()
        if br.width() < 10:
            br.setWidth(10)
        if br.height() < 10:
            br.setHeight(20) # Character-like height
            
        path = QPainterPath()
        path.addRect(br)
        return path
