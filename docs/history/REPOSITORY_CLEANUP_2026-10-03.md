# Cierre de higiene del repositorio — 2026-10-03

## Alcance y resultado

Rama de trabajo `maintenance/repository-hygiene-closure`, iniciada en `origin/main` `1db4f63b52af79247745b3a8a220fb728348218c`. No se integra a `main`, no se reescribe Git history, no se cambia lógica química ni comportamiento del producto. El `git diff --name-status origin/main...HEAD` del commit de cierre es el manifiesto exacto de archivos; el resumen de grupos está debajo.

**Cierre:** rama de mantenimiento publicada con push normal; las 23 ramas remotas redundantes y tres ramas locales ancestro de `origin/main` fueron eliminadas después de registrar/verificar sus SHAs. `origin/main` sigue en `1db4f63b52af79247745b3a8a220fb728348218c`; ramas únicas/activas preservadas. No se hizo merge ni force-push.

## Cambios y evidencia

| Grupo | Acción y prueba de seguridad |
|---|---|
| `src/sys` | Retirado PostScript de ImageMagick, 11,708,416 bytes, sin imports, referencias, packaging, fixtures ni entrypoints consumidores. Su blob histórico permanece. |
| `src/repro_v2.png` | Retirado PNG 800×400 sin referencias ni uso de producto/test. |
| Parches Clean2D de raíz | Retirados `clean2d_failed_local_graph_attempt.patch` y `clean2d_failed_tetrandrine_selection_integrity.patch`, artefactos diff no importados. El contexto y commit de origen quedan en `CAMPAIGNS.md`; no se editó `src/chemuson/clean2d/` ni se aplicó un parche. |
| Archive orbital de tests | Retirado `tests/archive/orbitals_fit_report/`: 124 archivos, 246,667 bytes (119 PNG y 5 TXT), salida generada; ningún test/fixture/CI lo consume. Los fixtures actuales permanecen. Se resumieron métricas/limitaciones y los defaults de los generadores ahora escriben fuera del árbol fuente. |
| UI experimental e intermedia | Retirados demo ejecutable/widgets/iconos SVG del spike, HTML mockup, baseline viejo, cinco carpetas de capturas de fase, capturas de comparación redundantes y generadores one-shot F7. Se conservan la referencia normativa `pyqt6-spike/theme.py`, `checks.json`, README de decisiones, cinco PNG finales, tres capturas reales KDE/Wayland aprobadas y siete pruebas visuales F7 representativas (HiDPI, pequeña pantalla, onboarding y templates 200%). OpenSpec Markdown/specs siguen archivados. |
| Reportes duplicados | Retirados `docs/cleanup_dead_code_report.md`, `docs/cleanup_phase2_architecture_report.md`, `docs/archive/refactor_arquitectura_2026_04.md` y `docs/refactor_phase2.md`; las decisiones, fallos y cautelas se consolidaron en `CAMPAIGNS.md`. `AGENT_REPORT.md` se mantiene como registro corriente. |
| Documentación y tooling | Plan/known issues reducidos a estado vigente; enlaces `file://` del manual pasan a rutas relativas. Comentarios de tema ya no enlazan prototipos retirados. `orbital_fit_report.py` añade `--output-dir`; ambas herramientas de preview usan el temp del sistema por defecto y aceptan `CHEMUSON_ORBITAL_REPORT_DIR`. No se añadieron ignores globales. |

La referencia global de tema, assets de producción, tests actuales, fixtures, baselines vigentes y documentos Markdown OpenSpec se conservan. La implementación `local_graph_cleaner.py` y todo `src/chemuson/clean2d/` no se modifican.

## Verificación funcional y documental

| Comando/ejercicio | Resultado |
|---|---|
| `python -m compileall src tests tools packaging` | Pasa. |
| `pytest --collect-only -q` | 1816 tests. |
| `pytest -q tests/architecture` | 269 passed. |
| UI dirigida (13 contratos UI centrales más settings/templates/document/export) | 299 passed. |
| `pytest -q` | 1760 passed, 55 skipped y **1 fallo de baseline**: `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend`. Misma identidad/resultados en baseline y ambos runs de cierre; no se toca. |
| Ruff `F401,F811,F821,E722,E741` | Un único F401 preexistente (`math`) en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`; se deja intacto por ser test Clean2D fuera de alcance. |
| `openspec validate --all --strict` | 43 passed, 0 failed. |
| Smoke Qt fuente, 1440×900 y 980×600, claro/oscuro | PASS; app bar 54 px, sidepanel/rail/canvas, tabs dirty/new, flyout, CommandPalette y onboarding presentes. `XDG_CONFIG_HOME` y capturas aislados en `/tmp`. |
| Smoke HiDPI, `QT_SCALE_FACTOR=2` | PASS, DPR 2.0; preview 176×112 píxeles para 88×56 puntos lógicos. |
| Preview orbital (`render_orbital_family_preview.py`) | PASS; PNG en directorio temporal. |
| Fit report orbital opcional | Limitación preexistente observada: `_family_metric_strings` accede `PiBondingParams.ring`, atributo inexistente; termina con `AttributeError` tras escribir salidas parciales **solo** en `/tmp`. El fallo no está relacionado con la nueva ruta de salida; no se repara una herramienta orbital fuera de alcance. Queda `NEEDS_REVIEW`, no se declara smoke PASS. |
| `git diff --check origin/main...HEAD` | Pasa en el cierre; un primer control halló trailing spaces en tres metadatos del baseline OpenSpec, corregidos en commit separado antes del push. |
| Links Markdown relativos | 350 documentos examinados; 0 enlaces relativos rotos. Se cambiaron tres `file://` locales del manual por paths relativos. |

La primera invocación de smoke Qt sin `PYTHONPATH=src` falló por setup del comando (`ModuleNotFoundError: chemuson`); reintentada con `PYTHONPATH=src` y configuración Qt temporal, pasó a escala 1 y 2. No cambió preferencias del usuario.

No se reconstruyó Flatpak/PyInstaller: no se tocó código/assets incluidos en el paquete ni manifiestos; el smoke instalado Flatpak y arranque PyInstaller del cierre UI previo permanecen válidos. Las modificaciones Python aquí son defaults de salida de herramientas dev y docstrings, no runtime químico.

## Métricas antes/después

Métrica de baseline en `openspec/changes/repository-hygiene-closure-2026-10-03/baseline.md`. Medición del árbol tras la poda remota/local, antes del último commit de cierre documental:

| Medida | Antes (HEAD `1db4f63`) | Después |
|---|---:|---:|
| Árbol de archivos Git lógicos | 22,926,388 bytes / 1244 archivos | 7,922,534 bytes / 983 archivos (−15,003,854 bytes; −261 archivos). |
| PNG versionados | 222 / 4,597,068 bytes | 34 / 1,610,231 bytes (−188 PNG; −2,986,837 bytes). |
| `src/sys` en checkout | 11,708,416 bytes | 0; blob histórico conservado. |
| UI-modernization | 56 PNG; 132 archivos / 2,211,651 bytes | 8 PNG / 529,643 bytes; 18 archivos / 565,021 bytes. |
| Evidencia F7 OpenSpec | 27 PNG / 2,195,088 bytes | 7 PNG / 821,501 bytes. |
| `tests/archive/orbitals_fit_report` | 124 archivos / 246,667 bytes | 0. |
| Checkout sin `.git`/venv/cache/build | 25 MiB asignados | 11 MiB asignados. |
| `.git` y pack | 169 MiB; 8.05 MiB sueltos (784); 14,142 objetos en 2 packs / 160.27 MiB | ~170 MiB; 8.33 MiB sueltos (844 en la medición post-prune, previa al commit documental final); 14,142 objetos en 2 packs / 160.27 MiB. Sin `git gc`; las refs removidas no compactan objetos. |

Árbol lógico final por área: `src` 3,337,688 bytes/294 archivos; `tests` 1,677,814/250; `docs` 771,502/53; `openspec` 1,771,874/335; `packaging` 50,623/19. No quedan archivos trackeados en `tests/archive`.

No se ejecuta `git gc`/prune de objetos; la eliminación de refs no compacta el almacén local.

### Candidatos de reducción histórica futura

Los 6895 blobs `flatpak/beta/repo/objects/*` suman 159,512,305 bytes lógicos y 158,119,906 bytes empaquetados (~150.8 MiB). Están en la punta activa `origin/gh-pages`, que publica el remoto Flatpak: **no son recuperables de forma segura sin migrar primero ese servicio**. El blob histórico de `src/sys` empaqueta 430,236 bytes (~0.41 MiB); purgar historia por ese solo archivo no es proporcional. Por ello no se propone todavía `git filter-repo`, BFG ni migración de objetos. Si se decide mover el remoto Flatpak a Releases/object storage, debe ser otra propuesta con backup y plan de despliegue.

## Auditoría de ramas — estado previo a la poda

Inventario inicial: 32 ramas remotas reales (más `origin/HEAD` simbólico); 23 puntas candidatas fueron revalidadas por SHA y `git merge-base --is-ancestor <ref> origin/main` tras fetch: todas idénticas a la baseline y ancestros. Tras el push normal de esta rama y el gate, se eliminaron las 23 por nombre explícito. Fetch/prune final deja 10 ramas remotas reales: `main`, esta rama de mantenimiento y las ocho únicas/activas enumeradas abajo; `origin/HEAD` sigue a `origin/main`.

### Ramas remotas `SAFE_TO_DELETE` (SHA registrados)

| Ref bajo `origin/` | Punta verificada |
|---|---|
| `architecture/phase1-module-catalog-contracts` | `47f96a479829bf60034ecfb29db173e604ae8c5b` |
| `architecture/phase10-extract-canvas-selection-geometry` | `6aa3a5d37995fc8171c9cec93794e66e86193f3d` |
| `architecture/phase11-extract-canvas-selection-bounds` | `66d2b43c8b677f8d05477156ba453ec182e6e51a` |
| `architecture/phase2-decouple-utils-autosave` | `c329d0b7144b846a5c80c92b6d81886cffee0a49` |
| `architecture/phase3-extract-shared-molecular-view` | `05c0d8696d9dd0f75ac08af58f0081974c58d701` |
| `architecture/phase5-eliminate-persistence-gui-exception` | `7ea76b2983742bff862e4fea155b6b13b6ba1564` |
| `architecture/phase6-establish-application-composition-root` | `8f04c4adea9ec67fb248bda7f1161045553847cd` |
| `architecture/phase7-close-application-shell` | `22fc6860c1464b777135544f691274cff0e71499` |
| `architecture/phase7-extract-application-shell` | `3db53a10746d8065dfd4cef2305b4570cadf8d29` |
| `architecture/phase8-extract-main-window-background-workers` | `97a3f03261eeb2a8426400c503f58afa33638a84` |
| `architecture/phase9-extract-main-window-clean2d-geometry` | `a26fa638ccfa53c251c21e80763a38a0cac154fd` |
| `architecture/phases10-15-selection-modularization` | `cd0ee77e75465f083680619efe2576bcdefc2bd4` |
| `cleanup/architecture-map-and-boundaries` | `71798871e4298b24130a6502ed458dd512744c6c` |
| `cleanup/remove-dead-code-and-tests` | `933a27cb068b22f251fdad354e4c63cde8e875ee` |
| `codex/add-error-indication-for-multiple-links-in-atoms` | `3dc1257de890269f69d26a2142b049492fa056d8` |
| `codex/add-improvements-and-new-features-to-chemuson` | `7eb2e604156a3a5d192bd60c46b8b5fe2e2a7424` |
| `codex/refactor-architectural-structure-of-chemuson` | `8a031a6595dfc57869f9e6352a24a2327a352497` |
| `feature/arrow-pushing-snapping` | `9ad2b915ee7c7a5c943a06a56b23f04e02c7e808` |
| `feature/geometry3d-compchem` | `ff876bcab4eb6e3a2e1b9accd978d2a5c714728a` |
| `feature/next-level-workbench` | `15477451d55b1a73682bce6d728d53e46fcb5e86` |
| `fix/flatpak-reproducible-build` | `8a934785ecb970f8830dee953bd66b77290a4503` |
| `fix/ui-openspec-post-integration` | `140a080c336515650bbeea0a4e6ead67a9999b23` |
| `release/ui-modernization-qa` | `1db4f63b52af79247745b3a8a220fb728348218c` |

### Ramas locales eliminadas

Después del push se eliminaron con `git branch -d` solo estas ramas ya integradas: `fix/ui-openspec-post-integration` (`140a080c336515650bbeea0a4e6ead67a9999b23`), `integration/ui-modernization` (`1db4f63b52af79247745b3a8a220fb728348218c`) y `release/ui-modernization-qa` (`1db4f63b52af79247745b3a8a220fb728348218c`). Permanecen `main`, mantenimiento, artículo y `ui/modernization`.

### Preserve unique/active refs

| Ref | SHA | Ahead/behind `origin/main` | Motivo |
|---|---|---:|---|
| `clean2d/campaign-implementation` | `6b153e3699bf2b04e60189a4606040f766c8eb5c` | 18/31 | Campañas Clean2D no integradas. |
| `clean2d/campaign-policy` | `31c565c8dfe15f2b4c2287eb4ece56e2b67d9732` | 2/31 | Política única Clean2D. |
| `codex/implementar-sistema-de-apariencia-en-chemuson` | `2a40df84137ca07bc4fd2583fd3ec31ff9b0eb47` | 1/271 | Apariencia/tema pendiente de cotejo. |
| `docs/chemuson-interactive-article` | `cdcae3c2e37df6967c61c7e7b315db0e225cfb67` | 36/31 | Artículo interactivo y assets únicos. |
| `gh-pages` | `c856fc1445a11576a5a2718cd40380bc54056853` | 6/252 | Remoto Flatpak publicado. |
| `refactor/packaging-pyproject` | `d44115b349fdd4be3bd25599c89cf605cb33613b` | 2/337 | Runtime hooks únicos para PyInstaller. |
| `tema-chemuson` | `167ade37081140747bb014cc613646ffdcf3063d` | 4/271 | Decisiones visuales únicas sin revisión de equivalencia. |
| `ui/modernization` | `677b8294138357774e3937ceecbb97dd621b8f26` | 43/31 | Contiene UI ya parcialmente integrada y commits Clean2D únicos. |

## Cierre de ramas y medidas posteriores

- Push inicial y final de `maintenance/repository-hygiene-closure`: normal, sin `--force`; la rama sigue la remota. La última punta se registra en el `git rev-parse HEAD` del cierre.
- `git push origin --delete` retiró exactamente las 23 ramas de la tabla `SAFE_TO_DELETE`; `git fetch --prune` confirmó las refs remotas retiradas.
- SHAs de `origin/main`, `origin/gh-pages`, las dos Clean2D activas y las otras seis ramas no-ancestro verificados intactos. `origin/main` no cambió.
- Ramas finales: 10 remotas reales (más `origin/HEAD` simbólico) y 4 locales. `git status --short` limpio al cierre; `git diff --check` pasa.
- Se conserva el objeto Git histórico hasta una decisión futura de compactación; no se ejecutó ningún rewriter ni `git gc --prune`.
