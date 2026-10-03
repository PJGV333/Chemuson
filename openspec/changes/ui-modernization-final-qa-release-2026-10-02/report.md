# Reporte de cierre — Fase 8: QA de modernización UI

**Rama:** `release/ui-modernization-qa`<br>
**Punto de partida:** `origin/main` / `140a080c336515650bbeea0a4e6ead67a9999b23`<br>
**Versión de la aplicación:** `0.3.0-dev`
**Resultado:** sin regresiones nuevas frente a baseline. No se modificaron código, tests, `main`, Clean2D, ChemName, persistencia ni geometría de plantillas.

## Tests, arquitectura y estática

| Verificación | Resultado |
|---|---|
| `python -m compileall src tests tools packaging` | Pasa (exit 0). |
| `pytest --collect-only -q` | Pasa; 1816 tests. |
| `pytest -q` baseline | 1760 passed, 55 skipped, 1 failed, 1137.93 s. |
| `pytest -q` final | 1760 passed, 55 skipped, 1 failed, 1152.05 s. El único fallo conserva identidad: `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`; ya constaba en la baseline Fase 7. |
| `pytest tests/architecture -q` final | 269 passed en 9.42 s; sin nuevas violaciones observadas. |
| Tests UI dirigidos | 304 passed en 522.96 s. Incluyen contratos de AppBar/tabs, rail/flyouts, panel lateral, paleta, preferencias, onboarding, plantillas, iconos/HiDPI, atajos y documentos/exportación. |
| Ruff scoped (`F401,F811,F821,E722,E741`) | Un F401 preexistente: `import math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`. No se corrigió por estar fuera de alcance. |
| `git diff --check` | Pasa al cierre. |

## UI y smoke

- Smoke Qt offscreen con light/dark a 1440×900 y 980×600; AppBar de 54 px, rail, canvas, tabs, indicador dirty, flyout de enlaces, SidePanel/Plantillas, CommandPalette y onboarding presentes.
- HiDPI `QT_SCALE_FACTOR=2`: DPR 2.0; preview de 88×56 puntos lógicos renderizado a 176×112 píxeles.
- Cinco capturas finales PNG (1440×900): ventana light/dark, flyout, Plantillas y CommandPalette. No existe PNG original de Fase 0 versionado; por eso no se produjo un montaje before/after. KDE/Wayland ya había sido aprobado manualmente.
- Advertencia no fatal del plugin Qt `offscreen`: `This plugin does not support propagateSizeHints()`. No apareció excepción funcional en los smokes finales.

## Packaging y assets

### Flatpak

Se construyó un bundle con el manifiesto oficial existente:
`/tmp/f8-flatpak/Chemuson-v0.3.0-dev-linux-x86_64.flatpak` (aprox. 105 MB).
El smoke corrió el paquete en Flatpak, fuera del checkout: versión `0.3.0-dev`, módulos UI bajo `/app/lib/python3.13/site-packages/chemuson`, módulos de AppBar/rail/panel/paleta/onboarding presentes, 69 SVG `i-*.svg`, 7 plantillas incorporadas, preview 88×56 y round-trip guardar/abrir `.cmsn` (967 bytes). Resultado: `PACKAGE SMOKE PASS`.

### Ejecutable Linux / salida llamada AppImage por el proyecto

El flujo oficial CI usa `pyinstaller` + `packaging/linux/build_appimage.sh`; ambos se probaron. PyInstaller construyó el ejecutable y el script generó `/tmp/f8-appimage/Chemuson-v0.3.0-dev-linux-x86_64.AppImage` (aprox. 326 MiB), con sidecars de actualización. Se verificaron 4794 entradas del archivo —incluidos AppBar, CommandPalette, onboarding, SidePanel, tokens y SVG—; `--version` devolvió `0.3.0-dev`, y la GUI permaneció iniciada en offscreen hasta su terminación controlada por timeout.

**Límite del entorno:** no está instalado `appimagetool`. El script de packaging existente envuelve/nombra el ejecutable portable PyInstaller como `.AppImage`; no se afirmó que se haya construido un contenedor AppImage Type 2. PyInstaller también avisó de bibliotecas opcionales ausentes para plugins Qt3D/QML y drivers SQL no usados; el ejecutable sí arrancó y los assets UI requeridos están embebidos.

## Documentación, versión y OpenSpec

- `PLAN.md` marcado Fases 0–8 completas y actualizado con las desviaciones reales: paleta `Ctrl+P`, `Ctrl+K` reservado a Clean2D 1 paso, SidePanel integrado, onboarding, plantillas single-click y HiDPI.
- Manual revisado al RC `140a080`; notas curadas de release añadidas sin simular una publicación. La workflow existente de GitHub genera notas automáticas (`generate_release_notes: true`).
- Fuente canónica sigue en `0.3.0-dev`. Dado que la modernización integrada es un hito visible, `0.4.0-dev` es una propuesta coherente para la siguiente publicación de prueba, sujeta a aprobación; no hubo bump, tag, metadatos de release ni publicación.
- Los siete OpenSpec UI pasaron `--strict` (7/7) antes del archivado oficial. `openspec archive` creó las carpetas `2026-10-03-*` y actualizó sus siete specs globales. El primer archive reportó una tarea de commit sin marcar en Theme Foundation; se verificó el commit `f1264fdd5c288a576f5c72481dfa3d696171ea67` en `origin/ui/modernization` y se completó su checkbox en el archivo archivado.
- Después del archive, `openspec validate --all --strict`: **42 passed, 0 failed**. La Fase 8 permanece activa; Fase 9 no se inicia.

## Deudas y decisión final

Se mantienen sin cambios el fallo histórico CompChem async y el F401 de un test Clean2D. Clean2D y las mejoras de geometría/química de plantillas (Fischer, Haworth, silla y tetrandrina) siguen separadas y no bloquean el cierre de UI. No se hizo merge automático ni se modificó `main`.
