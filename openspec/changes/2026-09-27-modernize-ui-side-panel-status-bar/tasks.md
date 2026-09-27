# Tasks: Panel lateral moderno y barra de estado (Fase 5)

## 1. OpenSpec y baseline

- [x] 1.1 Confirmar rama, HEAD local/remoto y árbol limpio.
- [x] 1.2 Revisar plan, mockup, spike, implementación de docks y OpenSpec de Fase 4.
- [x] 1.3 Ejecutar compileall, colección de pytest, suite completa y Ruff; registrar resultados.
- [x] 1.4 Definir y validar este cambio con `openspec validate 2026-09-27-modernize-ui-side-panel-status-bar --strict` antes de tocar implementación.

## 2. Preferencias y contrato visual

- [x] 2.1 Añadir `SidePanelPreferences` y load/save normalizados para `active_tab` y `visible` en la infraestructura existente de settings.
- [x] 2.2 Añadir métrica del panel y QSS tokenizado para tabs, overflow, páginas QDockWidget embebidas y labels actuales de status.

## 3. SidePanel y montaje

- [x] 3.1 Implementar `SidePanel` con `SideTabRow + QStackedWidget`, cinco tabs principales, menú overflow y ancho objetivo de 324 px.
- [x] 3.2 Insertar como páginas las siete instancias históricas QDockWidget, quitando el registro/docking clásico y su title chrome sin cambiar sus contenidos.
- [x] 3.3 Montar el panel al lado del lienzo y restaurar defaults/preferencias sin estado duplicado.
- [x] 3.4 Asegurar que `show_page(key)` es la ruta única para menú Ver, overflow y navegación interna a Validación.
- [x] 3.5 Mantener el contrato existente de preferencias del AppBar y alcanzar Apariencia desde su menú hamburguesa → Ver.
- [x] 3.6 Adaptar solo el layout visual de los controles de Validación al ancho de 324 px, conservando los mismos widgets, señales y handlers.

## 4. Menú y barra de estado

- [x] 4.1 Sustituir los siete `toggleViewAction()` por acciones Ver que muestran el SidePanel y seleccionan la página pedida.
- [x] 4.2 Añadir/actualizar control de visibilidad del panel que persiste el valor sin convertir acciones de página en toggles.
- [x] 4.3 Refinar jerarquía QSS de la barra de estado, conservando QStatusBar, `showMessage()`, IUPAC y carga; no añadir datos sintéticos.

## 5. Tests y arquitectura

- [x] 5.1 Añadir tests de settings válidos, ausentes e inválidos.
- [x] 5.2 Añadir tests de identidad 1:1 dock→página, contenido reutilizado y ausencia de docks clásicos registrados/flotantes.
- [x] 5.3 Cubrir tabs principales, overflow, navegación Ver, visibilidad, persistencia, selección Inspector, rutas de Validación, tema light/dark y tamaño 980×600.
- [x] 5.4 Registrar `side_panel.py` y la decisión de reutilización en M08 de `architecture/modules.yml`.
- [x] 5.5 Comprobar que los widgets existentes de navegación/reporte/corrección de Validación caben dentro de la página embebida.

## 6. Evidencia y validación final

- [x] 6.1 Generar las capturas reproducibles pedidas en `docs/ui-modernization/side-panel-phase-shots/` sin declararlas aprobación visual.
- [x] 6.2 Ejecutar tests dirigidos, suite existente de docks, smoke Qt offscreen, compileall, Ruff scoped y `git diff --check`.
- [x] 6.3 Ejecutar suite completa y comparar con baseline (baseline: 1799 collected, 1744 passed, 55 skipped; final: 1811 collected, 1756 passed, 55 skipped, 0 failed; Ruff F401 preexistente).
- [x] 6.4 Revalidar OpenSpec strict, revisar alcance/diff y comprobar cero cambios en subsistemas protegidos.
- [x] 6.5 Crear el commit solicitado, hacer push sin force a `origin ui/modernization` y verificar la referencia remota.
