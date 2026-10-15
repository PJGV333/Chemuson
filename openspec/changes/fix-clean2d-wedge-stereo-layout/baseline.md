# Baseline — Clean2D con cuñas

Captura previa a cualquier edición de código, rama `clean2d/wedge-bond-campaign` creada desde `origin/release/v0.3.0-beta.1-prep`.

## Git

- Rama base: `release/v0.3.0-beta.1-prep`, HEAD `9bc8192` (`docs(release): record beta preview acceptance gates`), coincidiendo con `origin/release/v0.3.0-beta.1-prep`.
- Rama de trabajo creada exactamente desde ese SHA. La Fase A (`fix/beta-save-shortcut-menu-access`, remoto `47c588d`) **no** es ancestro de esta rama: no se transportan cambios de GUI, `platform/settings.py`, `document_tabs.py`, `theme/qss.py` ni sus tests.
- `git status --short` antes de editar: sin salida.

## Comandos baseline (limitados externamente)

- `timeout 180 python -m compileall -q src tests tools packaging`: exit 0.
- `timeout 300 python -m pytest --collect-only -q`: exit 0; **2174 tests collected in 0.92s**.
- `timeout 180 ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 por un único hallazgo preexistente fuera de alcance: `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` (`import math`, F401). No se modificó.
- Grupo focalizado Clean2D sobre la base sin parche (14 archivos: `test_clean2d_engine_candidates`, `test_clean2d_invariants`, `test_clean2d_local_graph_cleaner`, `test_clean2d_safety`, `test_clean2d_complex_preserve`, `test_clean2d_tetrandrine_like_local`, `test_clean2d_deterministic_baselines`, `test_clean2d_regression_corpus`, `test_clean2d_smart_propose`, `test_clean2d_conformer_propose`, `test_clean2d_selection_bond_integrity`, `test_clean2d_contract_reporting`, `test_clean2d_quality_reporting`, `test_clean2d_debug_snapshots`): **240 passed in 38.70s**.
- Suite monolítica `python -m pytest -q`: no ejecutada en la captura de baseline; la validación completa se hizo por shards (ver `validation.md`).
- Entorno: Python local con RDKit y PyQt6 (`QT_QPA_PLATFORM=offscreen`).

## Reproducción del defecto en la base

Sonda desechable `run_clean2d_engine` (matriz de 7 casos) y sonda GUI con `ChemusonWindow` real en `offscreen`.

### Motor (base sin parche)

| caso | estado | fuente | motivo |
|---|---|---|---|
| 1 control acíclico sin estereo | `applied` | `rdkit_isolated` | — |
| 2 una cuña sólida | `failed-controlled` | ninguna | `geometria_actual_necesita_reconstruccion` / estable `worse-quality` |
| 3 una cuña rayada | `failed-controlled` | ninguna | igual que el caso 2 |
| 4 centro tetraédrico cuña + hash | `applied` | `rdkit_isolated` | — |
| 5 varios centros quirales | `applied` | `rdkit_isolated` | — |
| 6 selección con centro quiral | `failed-controlled` | ninguna | `geometria_actual_necesita_reconstruccion` |
| 7 selección con cuña en frontera | `failed-controlled` | ninguna | `geometria_actual_necesita_reconstruccion` |

Detalle del caso 2 (centro 2, vecinos 1/3/4, cuña 2–4), candidatos y motivo de descarte con la firma antigua:

- `current` → rechazado: `geometria_actual_necesita_reconstruccion`.
- `rdkit_isolated` → rechazado: **`cambio_firma_estereoquimica`**. Orden cíclico antes `(1,3,4)`; después `(3,4,1)` (rotación equivalente). Signos por par antes `+1,+1,+1`; después `+1,-1,+1`.
- `rdkit_direct` → rechazado: `rdkit_directo_deshabilitado:usar_rdkit_aislado`.
- `internal_templates` → rechazado: **`cambio_firma_estereoquimica`** (signos `+1,0,+1`).
- `clean2d_v2` → rechazado: `candidato_no_canonicaliza_geometria`.

### GUI (base sin parche, ventana real, `offscreen`)

Para una molécula con una cuña sólida y otra con una cuña rayada:

- `_on_clean_2d_full`, `_on_clean_2d` y `_on_clean_2d_publication`: **0 átomos movidos**, sin comando de undo, mensaje `Limpieza 2D omitida: geometria_actual_necesita_reconstruccion`.
- `_on_clean_2d_propose`: mensaje `La estructura debe optimizarse antes de proponer conformeros 2D` (requiere una base con calidad `good`, y la limpieza previa no existe).

## Estado de la estereoquímica en la base

Ninguna de las reproducciones invirtió la configuración: el fallo era por exceso de rechazo (nunca se aplicaba nada), no por corromper la molécula. La comprobación antigua, sin embargo, sí era ciega ante inversiones parciales que cruzan un vecino de frontera (ver `design.md` §2), riesgo que esta campaña cierra.
