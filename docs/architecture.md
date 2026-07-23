# Arquitectura de Chemuson

Este documento describe los límites actuales del código para que futuras limpiezas y refactors mantengan Chemuson entendible sin cambiar comportamiento químico o de UI por accidente.

## 1. Fuentes de Verdad

- **Catálogo de Módulos**: La fuente de verdad estructurada para la topología de módulos, APIs, rutas y dependencias es `architecture/modules.yml`. Para detalles individuales de cada módulo, consulte el directorio `docs/modules/`.
- **Reglas de Arquitectura**: Este documento (`docs/architecture.md`).

## 2. Entry points

| Entrada | Ubicación | Responsabilidad |
| --- | --- | --- |
| `chemuson.__main__:main` | `src/chemuson/__main__.py` | Parser CLI mínimo, `--version` y arranque diferido de la GUI. |
| Composition root gráfico | `src/chemuson/app/bootstrap.py` | Construcción de `QApplication`, crash reporter, ventana y ciclo de vida. |
| Script de consola `chemuson` | `pyproject.toml`, `[project.scripts]` | Entry point oficial instalado por packaging: `chemuson = "chemuson.__main__:main"`. |

El composition root se importa de forma diferida desde `main()` para que
importar el entry point, solicitar ayuda o ejecutar `chemuson --version` no
cargue PyQt6 ni subsistemas pesados. `chemuson.app` pertenece a M19 y es una
capa terminal: ningún módulo M00–M18 puede depender del bootstrap.

## Paquetes principales

### `chemuson.core`

| Aspecto | Detalle |
| --- | --- |
| Responsabilidad | Modelo molecular independiente de UI, estado químico, validación básica, análisis elemental y grafo químico multicapa. |
| Archivos principales | `model.py`, `layers.py`, `elemental_analysis.py`, `molecule.py`, `__init__.py`. |
| Puede importar | Biblioteca estándar y utilidades puras. Idealmente no depende de otros paquetes `chemuson`. |
| No debería importar | `chemuson.gui`, PyQt6, RDKit directo, `chemuson.clean2d`, `chemuson.chemio`. |
| API pública | `chemuson.core.__init__`, `MolGraph`, `Atom`, `Bond`, `ChemState`, `build_multilayer_chemical_graph`, funciones de análisis elemental. |
| Internos/privados | Helpers dentro de `model.py` y `layers.py` que no se reexportan. |

### `chemuson.chemio`

| Aspecto | Detalle |
| --- | --- |
| Responsabilidad | Import/export químico y persistencia: SMILES, Molfile, CML, CMSN y acceso seguro a RDKit. |
| Archivos principales | `rdkit_io.py`, `rdkit_safe.py`, `_rdkit_worker.py`, `persistence.py`, `cml_io.py`. |
| Puede importar | `chemuson.core`, workers/subprocess propios. |
| No debería importar | `chemuson.gui`, PyQt6 o `tools`, incluso bajo `TYPE_CHECKING`; la persistencia usa un contrato estructural local. |
| API pública | `rdkit_io.py`, `rdkit_safe.py`, `persistence.py`, `cml_io.py`. |
| Internos/privados | `_rdkit_worker.py` es privado aunque vivo: lo invocan workers/subprocess y no debe borrarse por no tener imports tradicionales. `PersistenceDocument` es el Protocol interno que desacopla CMSN de cualquier clase GUI concreta. |

### `chemuson.clean2d`

| Aspecto | Detalle |
| --- | --- |
| Responsabilidad | Limpieza/depiction 2D desacoplada de la GUI: candidatos, ranking, invariantes, seguridad geométrica, reparación local y políticas de moléculas complejas. |
| Archivos principales | `engine.py`, `safety.py`, `length_only.py`, `local_graph_cleaner.py`, `complex_policy.py`, `block_unwrap.py`, `scaffold_depiction.py`, `geometry.py`, `v2.py`. |
| Puede importar | `chemuson.core`, `chemuson.chemio` para conversiones/candidatos cuando sea necesario. |
| No debería importar | PyQt6, `chemuson.gui`, controllers, dialogs, actions. |
| API pública | Reexports en `chemuson.clean2d.__init__`, especialmente `run_clean2d_engine`, `generate_clean2d_candidates`, `rank_clean2d_candidates`, `evaluate_clean2d_layout`. |
| Internos/privados | Funciones con `_` en `engine.py`/`safety.py`; módulos de estrategia específica cuando no están reexportados. |

### `chemuson.chemcalc`

| Aspecto | Detalle |
| --- | --- |
| Responsabilidad | Fórmula, masa molecular, valencias típicas e hidrógenos implícitos auxiliares. |
| Archivos principales | `formula.py`, `mass.py`, `valence.py`, `__init__.py`. |
| Puede importar | `chemuson.core` y, con cuidado, utilidades de nomenclatura si no crea ciclos funcionales. |
| No debería importar | `chemuson.gui`, PyQt6, `chemuson.clean2d`, RDKit directo. |
| API pública | `chemuson.chemcalc.__init__`: `molecular_formula`, `format_formula`, `molecular_weight`, `implicit_h_count`. |
| Internos/privados | Helpers de conteo y normalización no reexportados. |

### `chemuson.chemname`

| Aspecto | Detalle |
| --- | --- |
| Responsabilidad | IUPAC-lite y nomenclatura semisistemática: selección de cadenas/anillos, plantillas, grupos funcionales, estereoquímica y coordinación. |
| Archivos principales | `engine.py`, `molview.py`, `parent_chain.py`, `rings.py`, `template.py`, `template_match.py`, `functional_groups.py`, `substituents.py`, `stereo.py`, `coordination.py`. |
| Puede importar | `chemuson.core`, `chemuson.chemcalc`, `chemuson.chemio` para fallback/plantillas cuando sea necesario, `chemuson.utils.resources`. |
| No debería importar | `chemuson.gui`, PyQt6, Clean2D, controllers. |
| API pública | `chemuson.chemname.__init__`: `iupac_name`, `NameOptions`. |
| Internos/privados | La mayoría de módulos de reglas (`rings.py`, `substituents.py`, `template_match.py`) son internos salvo que tests de dominio los ejerciten explícitamente. |

### `chemuson.geometry3d`

| Aspecto | Detalle |
| --- | --- |
| Responsabilidad | Modelos y servicios 3D, generación/proyección de conformeros, cache y export XYZ. |
| Archivos principales | `service.py`, `model.py`, `coordinates.py`, `rdkit_backend.py`, `obabel_backend.py`, `cache.py`, `export_xyz.py`. |
| Puede importar | `chemuson.core`, `chemuson.chemio` para conversión segura. |
| No debería importar | Widgets GUI, controllers, Clean2D heurístico salvo una integración explícita documentada. |
| API pública | Reexports en `chemuson.geometry3d.__init__`: `conformer_3d_for_graph`, `project_conformer_to_2d`, modelos 3D. |
| Internos/privados | Backends específicos (`rdkit_backend.py`, `obabel_backend.py`) detrás del servicio. |

### `chemuson.compchem`

| Aspecto | Detalle |
| --- | --- |
| Responsabilidad | Exportadores de química computacional y especificaciones de entrada para Gaussian/ORCA/NWChem. |
| Archivos principales | `exporters/common.py`, `exporters/gaussian.py`, `exporters/orca.py`, `exporters/nwchem.py`. |
| Puede importar | `chemuson.core`, `chemuson.geometry3d`. |
| No debería importar | PyQt6, `chemuson.gui`; la UI debe vivir en controllers/docks. |
| API pública | `chemuson.compchem.exporters`. |
| Internos/privados | Helpers de formateo por backend. |

### `chemuson.spectroscopy`

| Aspecto | Detalle |
| --- | --- |
| Responsabilidad | Predicción MVP de NMR/masa y registro de predictores espectrales. |
| Archivos principales | `service.py`, `__init__.py`. |
| Puede importar | `chemuson.core`, `chemuson.chemio` para expansión/conversión si hace falta. |
| No debería importar | `chemuson.gui`, PyQt6, controllers. |
| API pública | `predict_spectra`, `register_predictor`, modelos de picos reexportados en `__init__.py`. |
| Internos/privados | Predictores concretos y helpers no reexportados en `service.py`. |

### `chemuson.gui`

| Aspecto | Detalle |
| --- | --- |
| Responsabilidad | UI PyQt6: ventana principal, canvas/editor 2D, acciones, docks, controllers, comandos undo/redo, rendering y herramientas visuales. |
| Archivos principales | `main_window.py`, `canvas/`, `controllers/`, `items.py`, `docks.py`, `toolbar.py`, `actions/`, `commands/`, `dialogs/`. |
| Puede importar | Servicios de dominio (`core`, `chemio`, `clean2d`, `chemname`, `geometry3d`, `compchem`, `spectroscopy`, `update`, `name2structure`) para orquestar UI. |
| No debería importar | `tools`; no debería implementar lógica química pesada que pueda vivir en paquetes de dominio. |
| API pública | `chemuson.gui.main_window.ChemusonWindow`, `chemuson.gui.canvas.ChemusonCanvas` y constantes reexportadas en `canvas/__init__.py`. |
| Internos/privados | Mixins de `canvas/`, controllers concretos, builders, comandos y widgets auxiliares. `CanvasStructureMixin.rebuild_persistence_view()` es el hook público que satisface estructuralmente el contrato de persistencia y delega en la reconstrucción visual existente. |

### `chemuson.update`

| Aspecto | Detalle |
| --- | --- |
| Responsabilidad | Política, proveedor, seguridad, rollback, portable/windows y telemetría del auto-update. |
| Archivos principales | `core.py`, `provider.py`, `security.py`, `policy.py`, `portable.py`, `windows.py`, `rollback.py`, `telemetry.py`, `types.py`. |
| Puede importar | Biblioteca estándar, red/SSL y tipos propios. |
| No debería importar | `chemuson.gui`; la presentación de mensajes debe quedar en GUI/controllers. |
| API pública | Reexports en `chemuson.update.__init__`. |
| Internos/privados | Helpers de proveedor, validación y serialización dentro de cada módulo. |

### `chemuson.utils`

| Aspecto | Detalle |
| --- | --- |
| Responsabilidad | Recursos empaquetados, autosave y crash reporting compartido. |
| Archivos principales | `resources.py`, `autosave.py`, `crash_reporter.py`. |
| Puede importar | Biblioteca estándar y dependencias externas ya permitidas. `autosave.py` recibe persistencia y timers desde `gui.tab_manager`; `crash_reporter.py` usa PyQt6 para notificación visual. |
| No debería importar | Otros paquetes ChemUSON, dominio químico pesado, RDKit directo ni `tools`. `autosave.py` no importa PyQt6, GUI ni ChemIO. |
| API pública | Actualmente no reexporta API en `__init__`; usar módulos específicos. |
| Internos/privados | Helpers de rutas, serialización de crash/autosave y contratos estructurales mínimos de autosave. |

### `chemuson.name2structure`

| Aspecto | Detalle |
| --- | --- |
| Responsabilidad | Resolver nombres químicos a estructuras con conectores estáticos/PubChem y fallback seguro. |
| Archivos principales | `service.py`, `__init__.py`. |
| Puede importar | `chemuson.core`, `chemuson.chemio`. |
| No debería importar | `chemuson.gui`, PyQt6, Clean2D salvo normalización explícita posterior a importación. |
| API pública | `resolve_name_to_structure`, `NameToStructureResult`, conectores reexportados. |
| Internos/privados | Helpers de conectores y parseo dentro de `service.py`. |

## Dependencias observadas entre paquetes

| Paquete | Importa actualmente paquetes Chemuson |
| --- | --- |
| `core` | Ninguno. |
| `chemio` | `core`. La persistencia recibe en runtime un objeto que satisface `PersistenceDocument`, sin importar GUI. |
| `clean2d` | `core`, `chemio`. |
| `chemcalc` | `chemname` en estado actual. Revisar si puede invertirse o aislarse. |
| `chemname` | `core`, `chemcalc`, `chemio`, `utils`. |
| `geometry3d` | `core`, `chemio`. |
| `compchem` | `core`, `geometry3d`. |
| `spectroscopy` | `core`, `chemio`. |
| `gui` | Orquesta casi todos los subsistemas. |
| `update` | Ninguno. |
| `utils` | Ninguno entre módulos ChemUSON. `autosave.py` recibe serialización y temporizadores inyectados desde `gui.tab_manager`; `crash_reporter.py` conserva PyQt6 como dependencia externa para notificación visual. |
| `name2structure` | `core`, `chemio`. |

Estas dependencias describen el estado actual, no siempre el ideal. Las reglas siguientes definen el objetivo de mantenibilidad.

El catálogo mantiene actualmente cero `temporary_exceptions` y cero dependencias
circulares. La inyección runtime de un `ChemusonCanvas` en `PersistenceManager`
no crea una dependencia de módulo M01→M09: ChemIO sólo conoce el Protocol
estructural `PersistenceDocument`, mientras M09 implementa sus operaciones por
compatibilidad estructural.

## Reglas de arquitectura para mantener Chemuson entendible

| Regla | Motivo |
| --- | --- |
| `core` no debe importar GUI, RDKit ni servicios de alto nivel. | Mantiene el modelo químico testeable y estable. |
| `clean2d` no debe depender de widgets Qt ni de controllers. | Clean2D debe poder validarse como motor puro sobre `MolGraph`/coordenadas. |
| La GUI puede orquestar, pero no debe contener lógica química pesada. | Facilita tests sin PyQt y evita mezclar presentación con dominio. |
| RDKit debe entrar preferentemente por `chemio/rdkit_safe.py`, workers o backends aislados. | Evita aborts del proceso principal por extensiones nativas incompatibles. |
| `_rdkit_worker.py` y workers similares no deben borrarse por falta de imports entrantes. | Se invocan por subprocess y tienen ciclo de vida distinto al import normal. |
| Tests no deben usar `sys.path.append` manual. | `pyproject.toml` ya declara `pythonpath = [".", "src"]`. |
| `tools` no debe ser importado por `src` salvo justificación explícita. | Mantiene herramientas dev fuera del producto. |
| `chemio` no debería depender de `gui` en runtime. | La persistencia/conversión debe poder usarse sin PyQt. |
| `utils` debe permanecer sin dependencias ChemUSON. | La ausencia de dependencias M15 se refiere a módulos ChemUSON: `autosave.py` se prueba sin PyQt6, GUI, ChemIO ni RDKit, mientras `crash_reporter.py` conserva PyQt6 externo para notificación visual. |
| `update` debe permanecer desacoplado de la GUI. | Permite testear política/proveedor/seguridad sin interfaz. |
| Los tests PR/regression pueden existir, pero deben migrar gradualmente a nombres semánticos o acceptance data. | Reduce deuda organizacional sin perder casos químicos. |

## Deuda temporal

El catálogo mantiene un baseline vacío: no existen `temporary_exceptions` ni
dependencias circulares activas. Cualquier nueva desviación deberá justificarse
en un OpenSpec, registrarse explícitamente en `architecture/modules.yml` y
quedar acompañada por una condición verificable de eliminación.
