# Convergencia visual: producción PyQt6 y diseño aprobado

**Fuente de tokens conservada:** [`../pyqt6-spike/theme.py`](../pyqt6-spike/theme.py), referencia normativa del contrato de tema. El demo del spike y las capturas intermedias de comparación se retiraron al cerrar la campaña; su salida numérica `captures/checks.json`, la aprobación real de KDE/Wayland y las imágenes finales sí se conservan.

## Medidas aceptadas

Los **13/13 chequeos numéricos** están en [`captures/checks.json`](captures/checks.json):

| Chequeo | Valor verificado |
|---|---:|
| App bar / rail / botón / estado | 54 / 58 / 42 / 34 px |
| Menú visible / toolbar de texto por defecto | `false` / `false` |
| Botones del rail | 15 |
| Icono hamburguesa | presente |
| Tamaño mínimo | 900×560 |
| Rail/QToolBars visibles | no es `QToolBar`; lista visible vacía |
| Encaje a 980×600 | sí; botón conserva 42 px |

En la comparación RGB a 1440×900, las medias de shell fueron: app bar claro `(249,250,251)` producción vs `(247,248,250)` referencia; oscuro `(20,29,49)` vs `(21,31,51)`; rail claro `(251,252,253)` vs `(248,250,251)`; oscuro `(18,27,46)` vs `(20,30,49)`; estado claro `(252,253,253)` vs `(249,250,250)`; oscuro `(17,26,45)` vs `(19,29,47)`. Medición de tinta del primer icono: x=23–35 en ambas superficies tras corregir el desplazamiento de 8 px.

## Matriz de aceptación

Los once criterios se aceptaron: app bar sin menubar visible; menú accesible por hamburguesa/Alt; toolbar de texto contextual; rail `QWidget`; botones/iconos centrados; acciones Clean2D/Validar/Numerar accesibles; flyout por segundo clic/clic derecho; barra de estado funcional; ajuste a 980×600; temas claro/oscuro; acciones, atajos, menús y handlers conservados. La evidencia vigente es el JSON numérico, tests UI archivados, capturas canónicas [`../after/`](../after/README.md) y la inspección KDE/Wayland.

La aprobación manual no se sustituye por offscreen: se conservan [`real-qt-wayland-light.png`](real-kde/real-qt-wayland-light.png), [`real-qt-wayland-dark.png`](real-kde/real-qt-wayland-dark.png), [`hidpi-icon-sheet.png`](real-kde/hidpi-icon-sheet.png), los diagnósticos de plataforma y `real-capture-comparison.json`.

## Residuos explícitos

- El canvas productivo dibuja la escena real; el prototipo mostraba una demo. No es diferencia de shell.
- No se creó un panel derecho falso; producción usa SidePanel/widgets reales.
- Los iconos orbitales `QPainter` quedan documentados por su contrato/OpenSpec; su migración no fue requisito para esta campaña.
- La comparación homogénea before/after no es posible: no hay PNG original de Fase 0 versionado.

La cronología, correcciones Qt y QA final están resumidos en [`docs/history/CAMPAIGNS.md`](../../history/CAMPAIGNS.md) y en el [reporte de Fase 8](../../../openspec/changes/ui-modernization-final-qa-release-2026-10-02/report.md).
