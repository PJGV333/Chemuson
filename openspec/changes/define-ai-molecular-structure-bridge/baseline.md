# Baseline — define-ai-molecular-structure-bridge

Capturada el **2026-10-04**, antes de editar documentación o código. El punto de partida fue `main` limpio, actualizado contra `origin/main`; tras la captura se creó `ai/molecular-assistant-foundation` directamente desde ese `origin/main`.

## Git

```text
$ git fetch origin --prune
(exit 0)
$ git status --short --branch
## main...origin/main
$ git rev-parse HEAD
4068319e9a1deee7dbf239872191a8f509c8b641
$ git rev-parse origin/main
4068319e9a1deee7dbf239872191a8f509c8b641
$ git status --short
(sin salida; árbol limpio)
```

`origin/main` vigente y HEAD coincidían; el SHA de baseline es `4068319e9a1deee7dbf239872191a8f509c8b641` (`Clarify repository hygiene integration policy`). No se reutilizó ninguna rama prohibida.

## Baseline de calidad solicitado por `AGENTS.md`

```text
$ python -m compileall src tests tools packaging
(exit 0; sin errores)

$ pytest --collect-only -q
1816 tests collected in 0.84s
(exit 0)
```

`pytest -q`:

```text
FAILED tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend
1 failed, 1760 passed, 55 skipped in 1118.42s (0:18:38)
(exit 1)
```

Este fallo único ya está documentado como preexistente en `docs/history/CAMPAIGNS.md` para el mismo SHA; no se corrige ni se atribuye a esta campaña.

```text
$ ruff check src tests tools packaging --select F401,F811,F821,E722,E741
F401 [*] `math` imported but unused
 --> tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3:8
Found 1 error.
(exit 1)
```

F401 preexistente documentado por la campaña UI en `docs/history/CAMPAIGNS.md`; queda fuera de alcance y no se modifica.

## Inspección arquitectónica que fija el baseline

- M00 `core` define `MolGraph`; M01 `chemio` tiene `rdkit_io.smiles_to_molgraph`, que intenta primero `rdkit_safe`/worker RDKit aislado (8 s) y permite fallback existente a RDKit directo. También existe `rdkit_safe.smiles_to_molgraph_isolated`, sin fallback in-process; el diseño selecciona esta última para salida LLM no confiable y compara su grafo con el importador ordinario.
- M02 `clean2d` posee selección de depiction/candidatos SMILES y depende de M00/M01; no depende de IA.
- La importación visual actual de SMILES intenta primero M02 en `TemplateController._import_smiles_graph`; `on_import_smiles` inserta con `canvas._insert_molgraph` mediante comandos undoables.
- M16 `name2structure` resuelve nombres estáticos/PubChem y ya valida por parser aislado, pero no es una interfaz de LLM/proveedor.
- `architecture/modules.yml` contiene M00–M22; M23 es el siguiente ID disponible. No se modificó durante la planificación.

## Registro de comandos

Ejecutados antes de la primera escritura en el repositorio:

```text
git status --short
python -m compileall src tests tools packaging
pytest --collect-only -q
pytest -q
ruff check src tests tools packaging --select F401,F811,F821,E722,E741
```

La salida de colección enumera individualmente los 1816 tests; aquí se conserva el resumen exacto y el identificador/salida final de los fallos para mantener el baseline legible, según los baselines compactos previos del repositorio. No se modificaron archivos de baseline para ocultar fallos.