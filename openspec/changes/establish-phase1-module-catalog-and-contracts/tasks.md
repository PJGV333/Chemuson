# Tasks: Establish Phase 1 Module Catalog and Contracts

## Fase 0: Baseline y Preparación

- [x] **0.1 Capturar baseline del repositorio**
    - **Descripción**: Ejecutar los comandos de baseline y guardar resultados en `baseline.md` dentro del directorio del cambio OpenSpec.
    - **Archivos**: `openspec/changes/establish-phase1-module-catalog-and-contracts/baseline.md`
    - **Comandos**:
        - `git status --short`
        - `python -m compileall src tests tools packaging`
        - `pytest --collect-only -q`
        - `pytest -q`
        - `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`
    - **Aceptación**: Todos los comandos se ejecutan. Se registra la salida exacta sin asumir fallos permitidos.

## Fase 1: Catálogo de Módulos

- [x] **1.1 Crear directorio `architecture/`**
    - **Descripción**: Crear el directorio que contendrá el catálogo.
    - **Archivos**: `architecture/`
    - **Aceptación**: El directorio existe.

- [x] **1.2 Escribir el esqueleto de `modules.yml` con los 20 módulos (M00-M19)**
    - **Descripción**: Crear el archivo YAML con las entradas para M00-M19 según la tabla del design.md. Cada entrada con `id`, `name`, `title`, `responsibility`, `paths`, `status`, `risk_level`, `current_dependencies`, `target_dependencies` inicializados. YAML parseable.
    - **Archivos**: `architecture/modules.yml`
    - **Aceptación**: El archivo contiene exactamente 20 entradas. Cada ID es único. Los paths existen en el repositorio. El YAML se parsea con `yaml.safe_load`.

- [x] **1.3 Completar `public_api` para cada módulo**
    - **Descripción**: `public_api` contiene contratos deliberadamente públicos: símbolos listados en `__all__`; reexportaciones deliberadas desde `__init__.py`; funciones asociadas a `[project.scripts]`; para un módulo de archivo único, una función deliberadamente catalogada, como `M19.main`. Un símbolo importado directamente desde un módulo interno por otro paquete NO se convierte automáticamente en API pública. La validación será mediante AST, sin importar módulos.
    - **Archivos**: `architecture/modules.yml`
    - **Referencia**: `src/chemuson/core/__init__.py`, `src/chemuson/clean2d/__init__.py`, `src/chemuson/chemcalc/__init__.py`, `src/chemuson/chemname/__init__.py`, `src/chemuson/geometry3d/__init__.py`, `src/chemuson/spectroscopy/__init__.py`, `src/chemuson/update/__init__.py`, `src/chemuson/name2structure/__init__.py`, `src/chemuson/markush/__init__.py`, `src/chemuson/gui/canvas/__init__.py`, `src/chemuson/gui/controllers/__init__.py`, `src/chemuson/gui/commands/__init__.py`.
    - **Aceptación**: Los símbolos declarados coinciden con la definición de API pública.

- [x] **1.4 Mapear dependencias actuales por análisis visual de imports**
    - **Descripción**: Para cada módulo, identificar qué otros módulos chemuson importa (lectura de statements `from chemuson.*`). Los imports dentro de funciones o lazy siguen siendo runtime. Los imports bajo `TYPE_CHECKING` no entran en `current_dependencies`. Completar `current_dependencies`.
    - **Archivos**: `architecture/modules.yml`
    - **Referencia**: `docs/architecture.md` (tabla de dependencias observadas) como guía inicial; verificar contra código.
    - **Aceptación**: Las dependencias coinciden con los imports reales del código, excluyendo TYPE_CHECKING.

- [x] **1.5 Definir dependencias objetivo y prohibidas**
    - **Descripción**: Para cada módulo, definir `target_dependencies` (arquitectura ideal) y `forbidden_dependencies` basado en las reglas de `docs/architecture.md`.
    - **Archivos**: `architecture/modules.yml`
    - **Reglas**: `core` no importa ningún otro paquete chemuson. `clean2d` solo importa `core` y `chemio`. `chemname` importa `core`, `chemcalc`, `chemio`, `utils`. La GUI puede importar todo. `update` es independiente. `utils` no importa dominio químico pesado.
    - **Aceptación**: Las reglas son consistentes con `docs/architecture.md`.

- [x] **1.6 Registrar excepciones temporales**
    - **Descripción**: Para cada dependencia en `current_dependencies` que viola `forbidden_dependencies` o no está en `target_dependencies`, crear una entrada en `temporary_exceptions` con todos los campos obligatorios (source_id, target_id, file, import_path, reason, debt_ref, elimination_condition, type_checking_only). Una excepción runtime representa una dependencia actual no permitida o prohibida; un cruce bajo TYPE_CHECKING puede documentarse como excepción con `type_checking_only: true`.
    - **Archivos**: `architecture/modules.yml`
    - **Casos conocidos actuales relevantes**: M01 → M02 runtime; M03 → M04 runtime; M15 → M01 runtime; M01 → M09 TYPE_CHECKING; M15 → M09 TYPE_CHECKING. No se fija un número mágico de excepciones.
    - **Aceptación**: Cada excepción tiene los 8 campos obligatorios. No existen comodines.

- [x] **1.7 Registrar dependencias circulares conocidas**
    - **Descripción**: Documentar los ciclos conocidos (`chemio`↔`clean2d`, `chemcalc`↔`chemname`) en `circular_dependencies`.
    - **Archivos**: `architecture/modules.yml`
    - **Aceptación**: Cada ciclo tiene `modules`, `edges`, `severity`, `resolution_plan`.

- [x] **1.8 Completar campos restantes del catálogo**
    - **Descripción**: Llenar `entry_points`, `tests`, `internal_api`, `notes` para cada módulo.
    - **Archivos**: `architecture/modules.yml`
    - **Aceptación**: Ningún campo obligatorio está vacío.

## Fase 2: Documentación por Módulo

- [x] **2.1 Crear directorio `docs/modules/`**
    - **Descripción**: Crear estructura de directorios para documentación detallada.
    - **Archivos**: `docs/modules/`
    - **Aceptación**: El directorio existe.

- [x] **2.2 Crear `docs/modules/README.md`**
    - **Descripción**: Índice que lista los módulos M00-M19 con título y referencia a documento individual.
    - **Archivos**: `docs/modules/README.md`
    - **Aceptación**: Cada módulo tiene una entrada con link al documento individual. Los módulos sin documentación existente se marcan como "pendiente" sin enlace roto. Fuente estructurada identificada como `architecture/modules.yml`. Los módulos sin documentación existente se marcan como "pendiente" sin enlace roto.

- [x] **2.3 Documentar el Canvas (M09)**
    - **Descripción**: Crear `docs/modules/M09-canvas.md` con: bases directas y composición de submixins, inventario de archivos (20 archivos verificados), estado compartido, flujo de eventos mouse/teclado, selección e hit testing, undo/redo con QUndoStack (propiedad de cada instancia), sincronización modelo-escena, clipboard, texto/edición, render/export, serialización, dependencias internas, zonas de riesgo, extracciones futuras seguras vs. peligrosas.
    - **Archivos**: `docs/modules/M09-canvas.md`
    - **Referencia**: Leer los 20 archivos de `src/chemuson/gui/canvas/`. La clase `ChemusonCanvas` en `canvas_view.py` hereda en orden directo: CanvasInputMixin, CanvasSelectionMixin, CanvasTextMixin, CanvasRenderMixin, CanvasStructureMixin, QGraphicsView. `CanvasInputMixin` agrega 6 sub-mixins. `CanvasToolsBondingMixin` en `canvas_tools_bonding.py` agrega 5 sub-mixins de bonding. Cada `ChemusonCanvas` crea y posee su propio `self.undo_stack = QUndoStack(self)`.
    - **Aceptación**: El documento cubre los 14 temas requeridos. No copia código ni explica cada método individual. Bases directas en orden, composición de submixins, conteo real de archivos, y propiedad del QUndoStack documentados correctamente.

- [x] **2.4 Documentar módulos críticos brevemente**
    - **Descripción**: Para M00 (core), M01 (chemio), M02 (clean2d), M04 (chemname), escribir 1 párrafo de contexto en `docs/modules/M##-nombre.md`. No requiere la profundidad del canvas.
    - **Archivos**: `docs/modules/M00-core.md`, `docs/modules/M01-chemio.md`, `docs/modules/M02-clean2d.md`, `docs/modules/M04-chemname.md`
    - **Aceptación**: Cada documento explica responsabilidad, APIs principales y riesgos conocidos.

- [x] **2.5 Actualizar `docs/architecture.md`**
    - **Descripción**: Añadir referencia al catálogo YAML y a `docs/modules/`. No duplicar tablas; indicar que el YAML es la fuente de verdad actualizada.
    - **Archivos**: `docs/architecture.md`
    - **Aceptación**: El documento existente se mantiene; se añade sección de referencia al catálogo.

## Fase 3: Contrato de Agentes

- [x] **3.1 Crear `AGENTS.md` raíz**
    - **Descripción**: Escribir el protocolo de agentes con: OpenSpec como fuente del alcance, baseline obligatorio, tests focalizados y completos, prohibición de refactor oportunista, prohibición de ocultar fallos, informe de desviaciones, lista de archivos modificados, `AGENT_REPORT.md` para trabajo autónomo, ramas/worktrees separados, un solo agente escritor por rama, criterios de parada, protección especial de GUI/Clean2D/ChemName/persistencia.
    - **Archivos**: `AGENTS.md`
    - **Aceptación**: Todas las reglas del design.md §6 están presentes.

## Fase 4: Tests de Límites Arquitectónicos

- [x] **4.1 Crear directorio `tests/architecture/`**
    - **Descripción**: Crear estructura para tests arquitectónicos.
    - **Archivos**: `tests/architecture/__init__.py`
    - **Aceptación**: El directorio y el init existen.

- [x] **4.2 Implementar `test_module_catalog.py`**
    - **Descripción**: Tests que validan: YAML parseable, todos los IDs M00-M19 presentes y únicos, cada path existe en disco, campos obligatorios no vacíos, `status` es valor válido, `risk_level` es valor válido. Futuras pruebas deben validar rutas exclusivas, self-dependencies y M19.
    - **Archivos**: `tests/architecture/test_module_catalog.py`
    - **Aceptación**: Los tests pasan contra el `modules.yml` creado.

- [x] **4.3 Implementar `test_import_boundaries.py`**
    - **Descripción**: Test que analiza todos los archivos `.py` de `src/chemuson/` con `ast`, extrae imports `from chemuson.*`, verifica contra `forbidden_dependencies` del catálogo, aplica `temporary_exceptions` como whitelist, distingue `TYPE_CHECKING`. Mensajes de fallo con archivo, línea, import, regla violada.
    - **Archivos**: `tests/architecture/test_import_boundaries.py`
    - **Aceptación**: Los tests pasan con las excepciones registradas. Fallo inmediato ante violación no documentada.

- [x] **4.4 Implementar `test_public_api_exists.py`**
    - **Descripción**: Para cada módulo con `public_api` no vacío, verificar con AST que los símbolos existen en el archivo del `__init__.py` correspondiente. No importar los módulos.
    - **Archivos**: `tests/architecture/test_public_api_exists.py`
    - **Aceptación**: Los tests pasan contra las APIs declaradas.

- [x] **4.5 Implementar `test_no_tools_in_src.py`**
    - **Descripción**: Verificar que ningún archivo en `src/chemuson/` importa `tools` o `chemuson.tools`.
    - **Archivos**: `tests/architecture/test_no_tools_in_src.py`
    - **Aceptación**: El test pasa.

- [ ] **4.6 Implementar `test_exceptions_no_growth.py`**
    - **Descripción**: Verificar que la lista de `temporary_exceptions` en el YAML no crezca accidentalmente. El test lee el YAML y cuenta excepciones por módulo; si aparece una excepción no en la lista baseline, falla.
    - **Archivos**: `tests/architecture/test_exceptions_no_growth.py`
    - **Aceptación**: El test pasa con las excepciones actuales. Bloquea adiciones futuras sin cambio explícito.

## Fase 5: Validación Final

- [ ] **5.1 Ejecutar tests arquitectónicos focalizados**
    - **Descripción**: Ejecutar `pytest tests/architecture/ -v`.
    - **Aceptación**: Todos los tests pasan.

- [ ] **5.2 Ejecutar suite completa de tests**
    - **Descripción**: Ejecutar `pytest -q` para verificar que los tests nuevos no rompen la suite existente.
    - **Aceptación**: Los resultados coinciden con el baseline (mismos tests passing/failing/skipping). No hay fallos nuevos.

- [ ] **5.3 Validar OpenSpec**
    - **Descripción**: Ejecutar `openspec validate establish-phase1-module-catalog-and-contracts --strict`.
    - **Aceptación**: Validación exitosa.

- [ ] **5.4 Verificar baseline final**
    - **Descripción**: Re-ejecutar los comandos de baseline y comparar con `baseline.md`.
    - **Comandos**: mismos de 0.1.
    - **Aceptación**: Los resultados son idénticos al baseline inicial. No hay cambios involuntarios.

## Criterios Negativos (Aplican a todas las fases)

Ninguna tarea de este cambio puede:
- Mover o renombrar archivos en `src/chemuson/`.
- Modificar código de producción.
- Cambiar comportamiento químico, visual, de persistencia o undo/redo.
- Introducir dependencias de paquete nuevas.
- Actualizar snapshots de tests existentes.
- Crear paquetes nuevos en `src/chemuson/`.
- Corregir bugs oportunistas.
