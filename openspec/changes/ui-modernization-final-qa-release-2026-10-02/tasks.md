# Tasks: Fase 8 — QA final y cierre de release

## 1. OpenSpec y baseline
- [x] Crear proposal/design/tasks/spec con el scope exclusivo de QA + release closure.
- [x] Registrar HEAD, estado Git, últimos commits, compileall, collect y suite de baseline en `baseline.md`.
- [x] Validar este cambio con `openspec validate ui-modernization-final-qa-release-2026-10-02 --strict` antes de continuar.

## 2. QA de regresión y estática
- [x] Ejecutar suite completa y comparar conteos/identidad de fallos con `baseline.md`: baseline y final tienen 1760 passed, 55 skipped y el mismo fallo histórico CompChem async.
- [x] Ejecutar `pytest tests/architecture -q`: 269 passed, sin violaciones nuevas.
- [x] Ejecutar Ruff scoped (`F401,F811,F821,E722,E741`), identificar y conservar el F401 histórico de un test Clean2D fuera de alcance.
- [x] Ejecutar compileall completo y `git diff --check`.

## 3. UI dirigida y contratos
- [x] Ejecutar los tests dirigidos de modernización UI: 304 passed; la suite completa también conserva el mismo único fallo preexistente.
- [x] Confirmar atajos `Ctrl+P`, `Ctrl+K`, `Ctrl+Shift+K`, `Ctrl+Alt+K` y shortcuts actuales de ToolRail con tests existentes y manual actualizado.
- [x] Verificar tabs/dirty, undo/redo, paneles/menú Ver, plantillas single-click, exports y documentos mediante la suite existente; smoke de paquete incluye round-trip `.cmsn`. Sin cambios de química.

## 4. Smoke Qt y empaquetado
- [x] Smoke offscreen light/dark en 1440×900 y 980×600, más `QT_SCALE_FACTOR=2`; AppBar 54 px, rail, canvas, tabs, SidePanel, palette, flyout y onboarding presentes.
- [x] Construir y probar los mecanismos oficiales disponibles fuera del checkout: Flatpak instalado en sandbox y ejecutable PyInstaller empaquetado por `build_appimage.sh`.
- [x] Verificar SVG, tokens/QSS, onboarding, SidePanel, CommandPalette y plantillas en el Flatpak; round-trip guardar/abrir `.cmsn`.
- [x] Registrar disponibilidad real y límites: `appimagetool` no está instalado; el script oficial produce un ejecutable portable PyInstaller con extensión `.AppImage`, no un contenedor Type 2.

## 5. Evidencia y documentación
- [x] Crear cinco capturas en `docs/ui-modernization/after/` más README de resolución/tema/vistas; no hay PNG baseline original para un montaje.
- [x] Actualizar PLAN a Fases 0–8 completas, reflejando los atajos y estado final; Clean2D y química de plantillas siguen fuera de alcance.
- [x] Auditar manual, notas de release y fuentes de versión; mantener `0.3.0-dev` y documentar recomendación `0.4.0-dev` sin bump.
- [x] Escribir el reporte final en `report.md` y observaciones conocidas en `docs/ui-modernization/KNOWN_ISSUES.md`.

## 6. Archivado y cierre
- [x] Confirmar los siete OpenSpec UI en strict válido antes de archivar.
- [x] Archivar Fases 1–7 con `openspec archive` oficial, sin mover carpetas manualmente; siete specs globales actualizadas.
- [x] Validar después del archivado: `openspec validate --all --strict` dio 42 passed, 0 failed; Fase 8 queda activa.
- [x] Revisar alcance final, `git diff --check`, estado y resumen de archivos.
- [x] No hubo regresiones nuevas ni pérdida de assets; no iniciar Fase 9.
