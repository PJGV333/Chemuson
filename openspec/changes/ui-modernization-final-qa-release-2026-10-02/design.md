# Design: QA final y cierre de release

## Baseline

Punto de partida inmutable: `140a080c336515650bbeea0a4e6ead67a9999b23`. Baseline registrada en `baseline.md` antes de crear artefactos de esta fase. En este entorno Python 3.14.7, compileall y colección pasan; suite: 1 fallo, 1760 passed, 55 skipped. El fallo es `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`, ya listado entre los fallos históricos en la baseline Fase 7. La comprobación final debe comparar identidad y conteos; no se oculta ni se corrige aquí.

## Estrategia de verificación

1. Completar compileall, arquitectura, Ruff scoped y todos los tests UI de contrato indicados en `tasks.md`; guardar logs largos fuera del repo y registrar resúmenes exactos.
2. Hacer smoke offscreen light/dark a 1440×900 y 980×600 y con `QT_SCALE_FACTOR=2`: inspeccionar ventana, AppBar (54 px), ToolRail, SidePanel, canvas, DocumentTabs, CommandPalette, thumbnails, onboarding e iconos SVG. Capturas automáticas son evidencia, no sustituyen la aprobación manual KDE/Wayland ya otorgada.
3. Auditar primero `pyproject.toml`, manifiestos y scripts existentes; usar sólo el mecanismo de empaquetado soportado. Instalar/probar fuera del checkout en un directorio temporal aislado para detectar assets que dependan accidentalmente del source tree. No se añade otro sistema de build.
4. Verificar flujos archivo/documentos/undo-redo/menú Ver/atajos/plantillas/exportaciones por tests o smoke con fixtures seguros; no generar cambios de química.
5. Crear hasta cinco capturas finales y README conciso. Actualizar PLAN, manual/release notes según convención existente. Auditar todas las fuentes de versión antes de proponer un bump; mantener versión si no existe contrato inequívoco.
6. Validar strict los siete OpenSpec UI antes de archivarlos. Usar `openspec archive <id> -y` para cada cambio; nunca mover carpetas manualmente. Mantener activo este cambio y volver a validar OpenSpec al finalizar.

## Decisiones y límites

- El único estado aceptable de tests es el baseline documentado o mejor, sin fallos nuevos. Un nuevo fallo funcional/arquitectónico o pérdida de assets detiene la fase.
- No se cambian código/tests para silenciar fallos históricos. Las excepciones conocidas de Ruff se registran sin corregir fuera de alcance.
- Evidencia final limitada a cinco vistas: principal light/dark, flyout, Plantillas en SidePanel y CommandPalette. Un montaje before/after es opcional si puede realizarse con capturas existentes.
- Archivar OpenSpec es una operación documental oficial que puede actualizar `openspec/specs/`; se inspeccionará el diff y se conservarán rastreables los cambios archivados.

## Riesgos

- Los tests asíncronos/RDKit pueden variar por entorno; la baseline identifica el fallo reproducible conocido.
- Packaging puede no ofrecer artefacto ejecutable en esta plataforma; si no hay mecanismo oficial instalable, se documenta y verifica el mecanismo realmente soportado, sin inventar uno.
- Onboarding aparece sólo con settings limpios; smoke/evidencia aislará `QSettings` en entorno temporal y restaurará el entorno.
- La campaña de plantillas tiene deuda química documentada que permanece explícitamente intacta.