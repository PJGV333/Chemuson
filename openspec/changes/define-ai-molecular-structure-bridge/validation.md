# Validación — AI Molecular Structure Bridge / Phase 1

Rama `ai/molecular-assistant-foundation`, iniciada desde `9f34dd184406dc2ca3d9d65b02db1aacd653de11`. No se modifican Clean2D, GUI, dependencias externas ni la rama `main`. La suite y el adaptador usan fakes; no se hicieron solicitudes a un proveedor real ni a Internet.

## Tests enfocados y arquitectura

```text
$ pytest -q tests/architecture tests/test_molecular_assistant.py
333 passed, 2 skipped in 11.35s

$ ./.venv/bin/python -m pytest -q tests/test_molecular_assistant.py
60 passed in 0.83s
```

El entorno predeterminado `/usr/bin/python` no dispone de RDKit y el worker aislado devuelve `rdkit_unavailable`; por eso se omiten sólo los dos tests de integración química que requieren el worker. Se ejecutaron los 60 tests con el Python del `.venv` que sí contiene RDKit, incluyendo aceptación/rechazo real y equivalencia de la estructura ChemIO con el importador ordinario.

## OpenSpec y estática

```text
$ openspec validate define-ai-molecular-structure-bridge --strict
Change 'define-ai-molecular-structure-bridge' is valid

$ openspec validate --all --strict
Totals: 44 passed, 0 failed (44 items)

$ python -m compileall src tests tools packaging
(exit 0)

$ pytest --collect-only -q
1882 tests collected in 0.83s
(exit 0)

$ git diff --check
(exit 0)
```

Ruff sobre todos los archivos del alcance:

```text
$ ruff check src/chemuson/molecular_assistant tests/test_molecular_assistant.py tests/architecture/test_molecular_assistant_boundary.py tests/architecture/test_module_catalog.py tests/architecture/test_editor2d_selection_ownership.py tests/architecture/test_update_boundary_audit.py tests/architecture/test_operational_resilience_audit.py --select F401,F811,F821,E722,E741
All checks passed!
```

El Ruff solicitado para todo el repositorio mantiene el único F401 preexistente, idéntico al baseline y fuera de alcance:

```text
$ ruff check src tests tools packaging --select F401,F811,F821,E722,E741
F401 [*] `math` imported but unused
 --> tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3:8
Found 1 error.
(exit 1; mismo hallazgo documentado en baseline.md)
```

## Suite completa

```text
$ pytest -q
FAILED tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend
1 failed, 1824 passed, 57 skipped in 1147.62s (0:19:07)
(exit 1)
```

El único fallo es el mismo test histórico que ya fallaba antes de implementar Phase 1 (`1 failed, 1760 passed, 55 skipped` en el baseline de `baseline.md`/historial). La identidad y el resultado del assert fallido coinciden; los 64 tests adicionales pasan y dos tests de RDKit quedan omitidos bajo Python de sistema. No se modificó el test ni se ocultó el fallo.

## Límites verificados

- M23 sólo depende de M00/M01; catálogo y AST prohíben las dependencias inversas desde M00/M01/M02 y Clean2D no importa M23.
- Resultados fallidos no contienen `MolGraph`; el request no acepta documento, canvas ni contexto de mutación.
- JSON exacto `{"smiles":"..."}`, errores ChemIO por identificador completo, límites, transporte inyectable y errores públicos sanitizados están cubiertos sin red.
- No se ejecutó la prueba de snapshot de canvas/undo/dirty-state: requiere la integración UI posterior y sigue fuera de Phase 1.
