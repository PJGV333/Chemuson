# ChemName Acceptance Harness

Este documento describe el harness de pruebas de aceptación para el motor de nomenclatura (`ChemName`).

## Archivos

- Runner CLI/script: `tools/chemname_acceptance.py`
- Dataset de casos: `tests/data/chemname_acceptance_cases.yml`
- Fixtures `.mol` auxiliares: `tests/data/chemname_acceptance_mols/`
- Reporte JSON por defecto: `artifacts/chemname_acceptance_report.json`

## Formato de casos

Cada caso define:

- `id`: identificador único.
- `title`: descripción corta.
- `input_type`: `smiles`, `molblock` o `file`.
- `input`: SMILES, texto molblock o ruta de archivo.
- `options`: flags para `NameOptions`.
- `expect`: regex esperada o `ND`.
- `max_ms` (opcional): tiempo máximo permitido por caso.
- `skip_if_no_rdkit` (opcional): si `true`, el caso se marca como skip cuando el worker RDKit no está disponible.

Notas:

- El archivo `.yml` puede estar en sintaxis YAML o JSON (JSON es subconjunto de YAML).
- `expect: ND` pasa cuando el resultado es `N/D` o cuando se produce `ChemNameNotSupported`.

## Ejecución

Comando básico:

```bash
PYTHONPATH=src python tools/chemname_acceptance.py --cases tests/data/chemname_acceptance_cases.yml
```

Con salidas explícitas:

```bash
PYTHONPATH=src python tools/chemname_acceptance.py \
  --cases tests/data/chemname_acceptance_cases.yml \
  --report artifacts/chemname_acceptance_report.json \
  --junit artifacts/chemname_acceptance_report.xml
```

Opciones útiles:

- `--only id1,id2,id3`: ejecuta subset de casos por ID.
- `--limit N`: limita cantidad de casos.
- `--top-slowest N`: número de casos lentos a mostrar.
- `--quiet`: suprime tabla de stdout.

## Interpretación del reporte

El JSON contiene:

- `rdkit_worker_available`: disponibilidad del worker aislado.
- `total`, `passed`, `failed`, `skipped`, `errors`.
- `slowest`: lista resumida de los más lentos.
- `results`: detalle por caso (`status`, `name`, `reason`, `build_ms`, `name_ms`, `duration_ms`, etc.).

Estados:

- `pass`: expectativa cumplida.
- `fail`: expectativa no cumplida (regex/ND/max_ms).
- `skip`: caso omitido (por ejemplo, `skip_if_no_rdkit` con worker no disponible).
- `error`: problema de construcción/ejecución no mapeado a skip/fail.

## Flags de `NameOptions`

Casos pueden activar/desactivar banderas desde `options`, por ejemplo:

- `allow_coordination`
- `enable_exotic_hetero`
- `enable_special_templates`
- `enable_advanced_stereo`
- `enable_experimental`
- `rdkit_isolated`
- `strict`
- `return_nd_on_fail`

## Aislamiento RDKit (importante)

- Cuando `rdkit_isolated=true` y `input_type=smiles`, el runner usa worker aislado (`chemuson.chemio.rdkit_safe`) para convertir SMILES a molblock.
- Si el worker falla (timeout/crash/signal), el caso reporta `worker_error` y:
  - se marca `skip` si `skip_if_no_rdkit=true`,
  - o se evalúa como `ND` si `expect=ND`,
  - o se marca `error`/`fail` en los demás casos.

Esto evita que un crash nativo de RDKit tumbe la suite completa.

## Cómo agregar casos

1. Añade entrada en `tests/data/chemname_acceptance_cases.yml`.
2. Usa `input_type=file` con un `.mol` en `tests/data/chemname_acceptance_mols/` para casos estables sin dependencia RDKit.
3. Para SMILES dependientes de RDKit, define:
   - `options.rdkit_isolated: true`
   - `skip_if_no_rdkit: true`
4. Define `expect` como regex suficientemente robusta (o `ND` para negativos).
5. Ejecuta el runner y revisa `artifacts/chemname_acceptance_report.json`.

