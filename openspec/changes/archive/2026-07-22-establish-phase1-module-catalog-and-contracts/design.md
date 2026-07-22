# Design: Modularization Phase 1 - Catalog and Contracts

## Overview

El diseño crea un "gemelo digital" de la arquitectura actual en formato estructurado y verificable. Sirve como plano para la modularización de fases posteriores.

## 1. Catálogo de Módulos (`architecture/modules.yml`)

### Esquema YAML

Cada entrada del catálogo es un documento YAML dentro de una lista `modules`. Los campos obligatorios y su significado son:

| Campo | Tipo | Obligatorio | Descripción |
| --- | --- | --- | --- |
| `id` | string | Sí | Identificador estable (formato `M\d\d`). No se reutiliza. |
| `name` | string | Sí | Nombre corto del paquete (p. ej. `core`, `gui.canvas`). |
| `title` | string | Sí | Título legible para humanos. |
| `responsibility` | string | Sí | Descripción concisa de la responsabilidad del módulo. |
| `paths` | list[string] | Sí | Rutas relativas desde raíz del repo que pertenecen al módulo. |
| `status` | enum | Sí | Uno de: `stable`, `evolving`, `legacy`, `empty`. Definiciones abajo. |
| `public_api` | list[string] | No | Símbolos reexportados en `__init__.py` o entry points instalados. |
| `internal_api` | list[string] | No | Módulos o clases usados solo internamente o por tests. |
| `current_dependencies` | list[string] | Sí | IDs de módulos que este módulo importa actualmente (runtime). |
| `target_dependencies` | list[string] | Sí | IDs de módulos que este módulo debería importar en la arquitectura objetivo. |
| `forbidden_dependencies` | list[string] | No | IDs de módulos que NUNCA deben importarse. |
| `temporary_exceptions` | list[object] | No | Excepciones documentadas. Esquema abajo. |
| `circular_dependencies` | list[object] | No | Dependencias circulares conocidas. Esquema abajo. |
| `entry_points` | list[string] | No | Funciones/clases que actúan como puertas de entrada. |
| `tests` | list[string] | No | Rutas a suites de tests relevantes. |
| `risk_level` | enum | Sí | Uno de: `low`, `medium`, `high`. |
| `notes` | string | No | Contexto adicional. |

**Campos eliminados intencionalmente**:

- `maintainers_or_agents`: información administrativa no verificable. El catálogo es técnico.
- `openspec_capabilities`: meta-campo sin significado verificable.
- `allowed_dependencies`: eliminado para evitar ambigüedad. Se usa `target_dependencies` para lo deseado y `temporary_exceptions` para lo tolerado temporalmente.
- `internal_api` es opcional: se llena solo cuando hay interfaces internas que otro agente debería conocer pero no invocar directamente.

### Estado (`status`)

- `stable`: módulo con comportamiento estable, tests passing, sin deuda conocida.
- `evolving`: módulo con cambios activos o tests known-failing documentados.
- `legacy`: módulo con deuda técnica reconocida, candidato a refactor futuro.
- `empty`: directorio existe pero sin código funcional (p. ej. `templates/`). No se asigna ID.

### Regla de Propiedad de Rutas

Cada archivo Python de `src/chemuson/` pertenece como máximo a un módulo. Las rutas declaradas no pueden superponerse. Toda ruta debe existir. Ningún módulo puede depender de sí mismo. Estas reglas se verifican en los tests de catálogo.

### Esquema de `temporary_exceptions`

Cada excepción es un objeto con campos obligatorios:

| Campo | Tipo | Descripción |
| --- | --- | --- |
| `source_id` | string | ID del módulo que importa. |
| `target_id` | string | ID del módulo importado (o patrón de import). |
| `file` | string | Ruta exacta del archivo que contiene la violación. |
| `import_path` | string | Cadena exacta del import ofensivo (p. ej. `from chemuson.gui.canvas import ...`). |
| `reason` | string | Justificación de por qué la dependencia persiste. |
| `debt_ref` | string | Referencia a ticket, cambio OpenSpec o especificación que resolverá la deuda. |
| `elimination_condition` | string | Condición verificable para eliminar la excepción. |
| `type_checking_only` | bool | Si el import está guardado por `TYPE_CHECKING`, no es dependencia runtime. |

Los tests de límites verificarán que:
- La lista de excepciones no crezca sin añadir entradas nuevas explícitamente.
- Cada entrada tenga todos los campos obligatorios.
- No existan comodines como `gui/*` o `utils/*` sin especificar archivo e import.

### Esquema de `circular_dependencies`

Cada dependencia circular conocida se registra como:

| Campo | Tipo | Descripción |
| --- | --- | --- |
| `modules` | list[string] | IDs de los módulos involucrados en el ciclo. |
| `edges` | list[object] | Pares `{source, target, file, import_path}` que forman el ciclo. |
| `severity` | enum | `high`, `medium`, `low`. |
| `resolution_plan` | string | Plan para resolver el ciclo en fase futura. |

### Dependencias: distinción de claves

- **`current_dependencies`**: imports runtime reales; imports dentro de funciones o lazy siguen siendo runtime; imports bajo `TYPE_CHECKING` no entran en current_dependencies.
- **`target_dependencies`**: arquitectura objetivo.
- **`forbidden_dependencies`**: dependencias prohibidas.
- **`temporary_exceptions`**: una excepción runtime representa una dependencia actual no permitida o prohibida; un cruce bajo TYPE_CHECKING puede documentarse como excepción con `type_checking_only: true`.

## 2. Identificadores Estables (M00-M23)

### Lista de módulos

Los IDs se asignan a los siguientes módulos basados en el código actual:

| ID | name | title |
| --- | --- | --- |
| M00 | core | Modelo molecular fundamental |
| M01 | chemio | Import/export químico y persistencia |
| M02 | clean2d | Limpieza/depiction 2D |
| M03 | chemcalc | Cálculos químicos auxiliares |
| M04 | chemname | Nomenclatura IUPAC-lite |
| M05 | geometry3d | Servicios 3D |
| M06 | compchem | Exportación química computacional |
| M07 | spectroscopy | Predicción espectral |
| M08 | gui | Orquestación de interfaz PyQt6 |
| M09 | gui.canvas | Canvas de edición molecular |
| M10 | gui.controllers | Controllers de la GUI |
| M11 | gui.commands | Comandos undo/redo |
| M12 | gui.dialogs | Diálogos de la GUI |
| M13 | gui.items | Items gráficos (átomos, enlaces, etc.) |
| M14 | update | Subsistema de auto-actualización |
| M15 | utils | Utilidades compartidas |
| M16 | name2structure | Resolución nombre → estructura |
| M17 | markush | Estructuras Markush y polímeros |
| M18 | version | Gestión de versión |
| M19 | bootstrap | Arranque y composición de la aplicación |

Los IDs M20-M23 quedan reservados para módulos futuros.

### Reglas de estabilidad de IDs

- Un ID asignado NUNCA se reutiliza, incluso si el módulo desaparece.
- Módulos nuevos reciben IDs nuevos del rango disponible.
- Renombrar un módulo no cambia su ID.
- Dividir un módulo: el módulo original conserva su ID; los nuevos reciben IDs nuevos. Se registra la transición en `notes`.
- Fusionar módulos: el ID del módulo superviviente se conserva; los IDs eliminados se marcan como `retired` y se documenta la fusión.
- El orden en el archivo YAML no define los IDs.

## 3. Fuentes de Verdad

La jerarquía de información es:

1. **`architecture/modules.yml`**: fuente de verdad estructurada para topología, APIs y dependencias. Es la referencia para tests automáticos.
2. **`docs/architecture.md`**: visión general para humanos. Se actualiza para reflejar cambios del catálogo. No duplica tablas completas; referencia el YAML.
3. **`docs/modules/M##-nombre.md`**: documentación detallada por módulo. Explica responsabilidades, decisiones de diseño y riesgos. Es complementaria, no sustitutiva del YAML.

No se implementará generación automática de documentación durante esta fase. La consistencia se verifica manualmente y con tests de integridad del catálogo.

## 4. APIs Públicas

### Definición

`public_api` contiene contratos deliberadamente públicos:

- Símbolos listados en `__all__`;
- Reexportaciones deliberadas desde `__init__.py`;
- Funciones asociadas a `[project.scripts]`;
- Para un módulo de archivo único, una función deliberadamente catalogada, como `M19.main`.

Un símbolo importado directamente desde un módulo interno por otro paquete NO se convierte automáticamente en API pública.

La validación futura será mediante AST, sin importar los módulos.

## 5. Dependencias Circulares Conocidas

El código actual tiene dos dependencias circulares que deben registrarse:

### 5a. `chemio` ↔ `clean2d` (Severidad: alta)

- `chemio/depiction_candidates.py` importa de `clean2d/geometry`, `clean2d/safety`.
- `chemio/rdkit_io.py` importa de `clean2d/geometry`, `clean2d/scaffold_depiction`, `clean2d/block_unwrap`.
- `clean2d/scaffold_depiction.py` importa de `chemio/depiction_candidates`.
- `clean2d/block_unwrap.py` importa de `chemio/depiction_candidates`.
- `clean2d/engine.py` importa de `chemio/rdkit_safe` (lazy, dentro de funciones).

Plan de resolución: mover `depiction_candidates.py` a `clean2d` o a un paquete compartido en fase futura.

### 5b. `chemcalc` ↔ `chemname` (Severidad: media)

- `chemcalc/formula.py` y `chemcalc/valence.py` importan de `chemname.molview` (`MolView`).
- `chemname/engine.py` y múltiples módulos importan de `chemcalc.valence` (`implicit_h_count`).

Plan de resolución: extraer `MolView` a `core/model.py` o mover `implicit_h_count` a `core` en fase futura.

### 5c. `chemio` → `gui` y `utils` → `gui` (TYPE_CHECKING, Severidad: baja)

- `chemio/persistence.py`: `from chemuson.gui.canvas import ChemusonCanvas` bajo `TYPE_CHECKING`.
- `utils/autosave.py`: `from chemuson.gui.canvas import ChemusonCanvas` bajo `TYPE_CHECKING`.

Estas NO son dependencias runtime. Se marcan como `type_checking_only: true`.

## 6. Protocolo de Agentes (`AGENTS.md`)

El `AGENTS.md` raíz establece reglas vinculantes para cualquier agente IA que trabaje en el repositorio:

### Reglas obligatorias

- **OpenSpec como fuente del alcance**: todo trabajo comienza leyendo `proposal.md`, `design.md`, specs y `tasks.md` del cambio OpenSpec activo.
- **Baseline antes de editar**: ejecutar y registrar `git status --short`, `compileall`, `pytest`, `ruff` antes de cualquier modificación.
- **Tests focalizados y completos**: ejecutar tests relevantes al cambio Y la suite completa después.
- **Prohibido refactor oportunista**: no corregir bugs, reorganizar código o mejorar estilo fuera del alcance del cambio.
- **Prohibido ocultar fallos**: no modificar snapshots, baselines o excepciones para hacer pasar tests.
- **Informar desviaciones**: cualquier desviación del plan se documenta en `AGENT_REPORT.md`.
- **Listar archivos modificados**: todo commit o reporte debe incluir lista completa de archivos tocados.
- **`AGENT_REPORT.md` para trabajo autónomo**: registro de acciones, razonamiento y resultados.
- **Ramas/worktrees separados**: un solo agente escritor por rama.
- **Criterios de parada**: detenerse si tests existentes fallan, si el alcance requiere mover código de producción, o si aparecen dependencias no catalogadas.

### Protección especial

Los siguientes subsistemas requieren precaución adicional:
- GUI (`gui/`): no alterar señales, eventos, orden de inicialización.
- Clean2D (`clean2d/`): no modificar heurísticas sin resolver fallos baseline primero.
- ChemName (`chemname/`): no cambiar reglas de nomenclatura sin tests acceptance.
- Persistencia (`chemio/persistence.py`): no cambiar formatos de serialización.

### AGENTS.md anidados

No se crearán AGENTS.md anidados en esta fase. El contrato raíz cubre todas las reglas. Si futuras fases requieren reglas específicas por directorio (p. ej. reglas especiales para `clean2d/` cuando se aborde la deuda), se documentará la necesidad y se creará entonces.

## 7. Tests de Límites Arquitectónicos

### Ubicación

`tests/architecture/` con los siguientes archivos:

- `test_module_catalog.py`: valida esquema YAML, IDs únicos, paths existentes.
- `test_import_boundaries.py`: verifica que no existan imports prohibidos usando AST.
- `test_public_api_exists.py`: verifica que símbolos declarados existen en los archivos.
- `test_no_tools_in_src.py`: verifica que `src/` no importa `tools/`.

### Enfoque técnico

- Los tests analizan imports mediante módulo `ast` de Python estándar.
- Resuelven imports absolutos y relativos.
- Distinguen `TYPE_CHECKING` guards (no cuentan como dependencia runtime).
- No importan los módulos inspeccionados (evita cargar PyQt6, RDKit).
- Producen mensajes de fallo con: archivo, línea, import ofensivo, regla violada.
- Las excepciones se leen desde `architecture/modules.yml` y se aplican como whitelist.

### Baseline de excepciones

La primera ejecución de tests producirá violaciones reales. Estas se registrarán como `temporary_exceptions` en el YAML. Los tests posteriores fallarán si aparecen violaciones nuevas no documentadas, impidiendo crecimiento accidental de deuda.

## 8. Documentación del Canvas

El canvas se documenta en `docs/modules/M09-canvas.md` con:

- **Composición**: `ChemusonCanvas` extiende `QGraphicsView` mediante 5 mixins (`CanvasInputMixin`, `CanvasSelectionMixin`, `CanvasTextMixin`, `CanvasRenderMixin`, `CanvasStructureMixin`).
- **Meta-mixin**: `CanvasInputMixin` agrega 6 sub-mixins (`CanvasSelectionInputMixin`, `CanvasKeyboardMixin`, `CanvasContextMenuMixin`, `CanvasToolsAnnotationsMixin`, `CanvasToolsBondingMixin`, `CanvasToolsRingsChainsMixin`).
- **`CanvasToolsBondingMixin`** agrega 5 sub-mixins de bonding.
- **Archivos**: inventario de los ~20 archivos con responsabilidad principal de cada uno.
- **Estado compartido**: variables de instancia clave (`model`, `scene`, `state`, `drawing_style`, variables de drag/selection).
- **Flujo de eventos**: mouse/keyboard a través de mixins, orden de dispatch.
- **Selección y hit testing**: lógica de selección, fragment pivot, bond hit testing.
- **Undo/redo**: integración con `QUndoStack`, comandos desde `gui/commands/`.
- **Sincronización modelo-escena**: cómo `CanvasStructureMixin` mantiene consistencia entre `MolGraph` y `QGraphicsScene`.
- **Clipboard**: copiar/pegar estructuras, imágenes, anotaciones.
- **Texto y edición**: `CanvasTextMixin`, edición rica, anotaciones.
- **Render/export**: `CanvasRenderMixin`, PNG/SVG/PDF, bounds, overlays.
- **Serialización**: molfile, SMILES, CMSN desde el canvas.
- **Dependencias internas**: imports desde `core`, `chemio`, `chemname`, `gui/items`, `gui/commands`, `gui/dialogs`.
- **Zonas de riesgo**: archivos grandes (`items.py` 5704 líneas, `canvas_selection.py` 5154 líneas, `canvas_structure.py` 3875 líneas).
- **Extracciones futuras**: qué es seguro extraer (helpers puros de cálculo) vs. qué es peligroso (orden de eventos Qt, geometría visual).

## 9. Baseline

La propuesta exige registrar resultados reales de:

- `git status --short`
- `python -m compileall src tests tools packaging`
- `pytest --collect-only -q`
- `pytest -q`
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`

El baseline se guarda en `openspec/changes/establish-phase1-module-catalog-and-contracts/baseline.md`.

El baseline debe reflejar el estado real del repositorio en el momento de captura, sin asumir fallos permitidos.

## 10. Criterios Negativos

Esta fase NO puede:

- Mover código de producción.
- Renombrar paquetes.
- Modificar comportamiento químico.
- Modificar geometría visual.
- Alterar señales o eventos Qt.
- Cambiar undo/redo.
- Cambiar persistencia.
- Actualizar snapshots.
- Corregir Clean2D.
- Corregir ChemName.
- Crear frameworks de inyección de dependencias.
- Eliminar workers por ausencia de imports tradicionales.
- Crear paquetes nuevos en `src/chemuson/`.
