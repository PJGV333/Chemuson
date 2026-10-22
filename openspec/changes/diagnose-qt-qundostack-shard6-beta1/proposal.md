# Propuesta — diagnóstico SIGABRT Qt/QUndoStack en shard 6

## Síntoma

Tras publicar `aff565a064f17dccf98b0b4f23c7ddeba593ca31` en `stabilize/v0.3.0-beta.1`, el Actions run `39689871971` vuelve a abortar en `tests/test_ui_polish.py::test_onboarding_close_releases_overlay`. El stack nativo incluye `QUndoStack::~QUndoStack`/`cleanChanged` durante destrucción de widgets. El mismo shard abortó dos veces localmente; `tests/test_ui_polish.py` completo pasa aislado, lo que apunta a interacción con el orden/estado previo pero no demuestra causa.

## Objetivo

Comparar el comportamiento de la secuencia exacta afectada entre el commit base `61f7647b683db3eca19ea06e5ec82ab2897d640a` y el SHA actual, en snapshots de fuente aislados sin worktrees, y clasificar la causa como regresión, orden de pruebas, lifecycle Qt o entorno. Si se confirma raíz corregible, realizar el cambio mínimo y verificarlo con focales, ocho shards, Windows smoke y Flatpak smoke.

## Fuera de alcance

No tocar química, Clean2D, ChemName, persistencia `.cmsn`, representación gráfica, test acceptance manual o empaquetado. No añadir skips/xfails, no ampliar límites, no mover pruebas para esconder una falla, no crear worktrees permanentes, no mergear con la rama beta, ni crear tag/Release/publicación.
