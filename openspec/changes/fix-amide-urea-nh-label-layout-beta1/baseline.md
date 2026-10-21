# Baseline — layout visual de etiquetas N–H

## Rama y estado antes de cambios

El workspace indicado literalmente por el usuario (`/home/ccachyavgp/Documentos/ChemUSON`) no existe en este entorno. Se usó el checkout disponible en `/home/unison-pjgv/Documentos/GitHub/Chemuson`, que está en la rama base solicitada `release/v0.3.0-beta.1-prep` y con HEAD sincronizado con origin (`409581dc685fd166428bde47d8355633a8eb05da`). Antes de crear la rama, `git status --short` estaba vacío (worktree limpio). Se creó `fix/amide-urea-nh-label-layout-beta1`; logs: `baseline/00-branch-check.log`, `baseline/01-git-status.log` (salida vacía), salidas de comandos obligatorios en `baseline/02`–`05`.

## Comandos obligatorios (entorno `.venv`, `QT_QPA_PLATFORM=offscreen`)

| Comando | Exit | Resultado |
|---|---:|---|
| `git status --short` | 0 | Vacío antes de crear la rama y los documentos del OpenSpec. |
| `python -m compileall src tests tools packaging` | 0 | Compilación completa. Salida exacta: `baseline/02-compileall.log`. |
| `pytest --collect-only -q` | 0 | **2,592 tests** recolectados en 2.31 s. Salida exacta: `baseline/03-collect.log`. |
| `pytest -q` (límite operativo 600 s) | 124 | Interrumpido al llegar a 600 s, progreso ~30%; el log no contiene resumen ni fallo completo: `baseline/04-pytest.log`. Es un timeout basal, no una suite declarada verde. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | 1 | Único error conocido: F401 histórico `math` sin usar en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; salida exacta `baseline/05-ruff.log`. |

Los archivos `.exit` registran el código exacto. `03-collect.log` y `04-pytest.log` contienen la salida íntegra no comprimida. La línea de estado limpia se capturó antes del OpenSpec; el trabajo desde entonces está en la nueva rama.

## Diagnóstico visual previo

Capturas Qt post-fix de trim: `../improve-amide-urea-2d-continuity/evidence/after/urea-after-light-1x.png` y `.../urea-after-dark-1x.png`. El `N` del `NH` se desplaza con el texto completo y `H` solo aparece inline a un lado. En el fragmento doblado, ese renglón compite con el enlace vecinal; el problema no es el tramo centro-a-centro: los enlaces siguen en 40 unidades y el trim mejorado está activo. La corrección requerida actúa solo en el H visible, conservando los 40 y el cálculo actual de clearance.
