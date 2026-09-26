# Convergencia visual: producción vs. spike PyQt6 aprobado

**Contrato visual**: [`../pyqt6-spike/`](../pyqt6-spike/) (spike aprobado, mockup de
referencia). En cualquier conflicto visual, **el spike gana**; la lógica funcional
histórica (drawing, acciones, atajos, menús, handlers) se conserva detrás de la
superficie nueva.

## Cómo se generó

- Entorno: el `.venv` del repositorio (PyQt6), `QT_QPA_PLATFORM=offscreen`.
- Script: [`make_captures.py`](make_captures.py) (1440×900 y 980×600, temas
  claro/oscuro, producción y spike).
- Salida numérica: [`captures/checks.json`](captures/checks.json).
- Capturas: [`captures/`](captures/) (ventana completa, 980×600, close-ups de
  rail y app bar, flyout de enlaces, contact sheets `contact-*.png`).

## Chequeos numéricos (13/13 OK)

| Chequeo | Esperado | Real |
|---|---|---|
| Altura de la app bar | 54 px | 54 px |
| Ancho del rail | 58 px | 58 px |
| Botones del rail | 42 px | 42 px |
| Barra de estado | 34 px | 34 px |
| `QMenuBar` visible | `False` | `False` |
| Text toolbar visible por defecto | `False` | `False` |
| Botones en el rail | 15 | 15 |
| Icono de hamburguesa | presente | presente |
| Tamaño mínimo de la ventana | 900×560 | 900×560 |
| El rail es un `QToolBar` | `False` | `False` |
| `QToolBar` visibles | `[]` | `[]` |
| Encaja en 980×600 | sí | sí |
| Botón 42 px a 980×600 | 42 px | 42 px |

## Comparación por regiones (media RGB, 1440×900)

| Región | Prod claro | Spike claro | Prod oscuro | Spike oscuro |
|---|---|---|---|---|
| App bar | (249, 250, 251) | (247, 248, 250) | (20, 29, 49) | (21, 31, 51) |
| Rail | (251, 252, 253) | (248, 250, 251) | (18, 27, 46) | (20, 30, 49) |
| Estado | (252, 253, 253) | (249, 250, 250) | (17, 26, 45) | (19, 29, 47) |

Los iconos del rail se midieron por *ink* de píxeles: primer botón (puntero)
en **x 23–35** tanto en producción como en el spike (antes del fix de
centrado, producción estaba en x 13–27: desplazado 8 px a la izquierda).
El fondo del rail ahora sigue el token `surface` en ambos temas (fix del
viewport del `QScrollArea`: Qt no estiliza el viewport con `> QWidget`; se
usa el patrón "nieto transparente + fondo del scroll").

## Matriz de aceptación visual

| # | Criterio | Veredicto | Evidencia |
|---|---|---|---|
| 1 | App bar de 54 px, sin `QMenuBar` visible | **SÍ** | `checks.json`, `prod-appbar-light.png` |
| 2 | Menú accesible (hamburguesa + tecla Alt) | **SÍ** | `test_alt_opens_menu_popup_and_hamburger_exists`, smoke funcional |
| 3 | Text toolbar oculta por defecto, contextual | **SÍ** | `test_text_toolbar_hidden_by_default_and_contextual` |
| 4 | Rail = `QWidget` de 58 px (sin `QToolBar` visible) | **SÍ** | `checks.json`, `prod-rail-*.png` |
| 5 | 15 botones de 42 px, iconos 21 px, **centrados**, sin recortes, sin badges kbd | **SÍ** | comparación de ink (x 23–35), `prod-rail-*.png`, `test_no_kbd_badges_in_rail` |
| 6 | Clean2D / Validar / Numerar fuera del rail, `QAction` vivos | **SÍ** | `test_clean2d_validate_numbering_remain_accessible_without_rail_buttons` |
| 7 | 1er clic activa · 2º clic (categoría activa) y clic derecho abren flyout | **SÍ** | `test_second_click_on_active_category_opens_flyout`, `prod-bond-flyout-light.png` |
| 8 | Barra de estado de 34 px (tool / IUPAC / carga) | **SÍ** | `checks.json`, `prod-full-*.png` |
| 9 | Encaja en 980×600 sin micro-iconos (scroll compacto invisible) | **SÍ** | `checks.json`, `prod-980x600-light.png`, `contact-980x600.png` |
| 10 | Temas claro/oscuro equivalentes al spike (franja shell) | **SÍ** | tabla de media RGB, `prod-full-dark.png` vs `spike-full-dark.png` |
| 11 | Acciones, atajos, menús, handlers intactos (funcional) | **SÍ** | suite UI dirigida (320+ tests), suite completa (solo 4 fallos de baseline preexistentes) |

### Residuos documentados (no bloquean la aceptación)

- **Lienzo**: la producción pinta la escena real (molécula + retícula); el
  spike pinta su demo. La diferencia de píxeles en la región del lienzo es
  esperada y no es un fallo de la shell.
- **Panel derecho**: el spike lo simula; por instrucción de esta tarea no se
  implementa un panel falso — la zona derecha sigue siendo lienzo.
- **Iconos de orbitales**: residuo SVG documentado en el OpenSpec de la
  Fase 4 (se reutilizan iconos `QPainter` existentes; migración a SVG
  diferida).

## Veredicto

**CONVERGENCIA ALCANZADA**: los 13 chequeos numéricos pasan, las regiones de
shell (app bar, rail, estado) coinciden con el spike en ambos temas dentro
del ruido de los iconos, y los 11 criterios de la matriz son SÍ. La función
se preserva íntegramente (ver `AGENT_REPORT.md`, sección "Convergencia
visual con el spike PyQt6 aprobado").
