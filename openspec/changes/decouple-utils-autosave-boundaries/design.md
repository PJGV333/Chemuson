# Design: Decouple Utils Autosave Boundaries

## Estado actual

`utils.autosave` importa exactamente `PersistenceManager` en runtime y
`ChemusonCanvas` bajo `TYPE_CHECKING`; además importa `QObject` y `QTimer`.
`CanvasTabManager.create_document_tab()` es el único constructor de
`AutosaveManager`. Autosave usa el canvas exclusivamente como argumento de
`PersistenceManager.save_to_dict`, el undo stack exclusivamente mediante
`isClean()`, y los timers exclusivamente para iniciar, detener y ejecutar el
callback de autosave.

Antes:

```text
gui.tab_manager -> utils.autosave -> chemio.persistence
                                -> gui.canvas (TYPE_CHECKING)
                                -> PyQt6
```

Después:

```text
gui.tab_manager -> utils.autosave (contratos locales)
       |                  ^
       +-> PersistenceManager.save_to_dict y factory QTimer
```

## Diseño elegido

Se usa el patrón de callbacks inyectados con `typing.Protocol` locales:
`AutosaveDocument` representa el documento opaco, `AutosaveUndoStack` expone
`isClean`, y `AutosaveTimer` expone iniciar/detener. El constructor recibe un
serializador tipado y un factory que ya entrega timers configurados. El módulo
de utils no importa PyQt6, GUI, ChemIO ni usa importación dinámica.

`CanvasTabManager` posee los adaptadores concretos: configura `QTimer` con los
intervalos, single-shot y callback recibido, e inyecta
`PersistenceManager.save_to_dict`. Es el composition root autorizado porque
ya depende de GUI y ChemIO.

## Ownership y compatibilidad

Autosave mantiene ownership de rutas, UUID, hash, metadata, JSON, rotación y
decisiones clean/dirty. Persistencia mantiene ownership de serialización y
recuperación permanece en `RecoveryController`. El constructor cambia de modo
explícito para recibir colaboradores: no se conserva un default que obligue a
utils a reimportar dependencias prohibidas. Se actualiza su único call site.

No cambia la política de errores: fallos de listar/eliminar siguen manejando
`OSError` localmente y un fallo del serializador o de escritura sigue
propagándose como antes.

## Migración, tipos y recursos

No se crean paquetes ni dependencias. Los Protocols se limitan a las operaciones
observadas; no representan toda la API del canvas. Los fakes de pruebas son
estructurales y no importan PyQt6. Rutas, nombres, timestamps, metadata,
intervalos y formatos se preservan sin cambios.

## Catálogo y baseline

M15 quedará con `current_dependencies`, `target_dependencies` y
`temporary_exceptions` vacíos. El baseline congelado mantendrá las ocho
identidades M01/M03, validará que el catálogo real coincide con ellas y añadirá
una prueba que rechaza la reaparición de cualquiera de las dos identidades M15.

## Pruebas

Se caracterizan escritura y rotación con fakes. Las pruebas nuevas verifican
inyección de serializador, argumentos, propagación de errores, AST sin imports
ChemUSON/PyQt6, aislamiento en subprocess y composición concreta en tab manager.
Las suites arquitectónica y completa verifican que no aparezcan nuevas
violaciones.

## Riesgos, alternativas y rollback

El riesgo principal es alterar timers Qt; el factory conserva su configuración
existente en el mismo call site. Se descarta un puerto de persistencia de clase
por añadir una abstracción mayor para una sola operación, y se descartan imports
lazy por preservar la deuda. El rollback es revertir los commits de esta fase,
sin migración de datos porque el formato no cambia.

## Corrective design

### Núcleo genérico

`AutosaveManager` usa `DocumentT = TypeVar("DocumentT")` y
`Generic[DocumentT]`. `AutosaveSerializer` conserva el mismo tipo mediante
`Callable[[DocumentT], dict[str, object]]`; no se mantiene un Protocol vacío
para representar documentos arbitrarios. `AutosaveUndoStack`, `AutosaveTimer`
y `AutosaveTimerFactory` conservan sus contratos mínimos actuales.

### Contrato del controlador

`CanvasTabManager` almacena un `AutosaveController` pequeño con únicamente los
métodos usados por sus call sites: `start`, `stop`, `set_original_path`,
`restart_debounce`, `cancel_debounce` y `cleanup_after_save`. No contiene
métodos privados ni pretende representar toda la implementación del core.

### Adaptador Qt

`_QtAutosaveManager` hereda de `QObject`, contiene un
`AutosaveManager[ChemusonCanvas]` y delega explícitamente los métodos del
controlador. No hereda del core ni usa `__getattr__`. Una propiedad de solo
lectura `core` permite introspección dirigida en pruebas sin exponer estado
mutable. El adaptador conserva el ownership de los `QTimer`, su parent, los
intervalos, `singleShot` y los callbacks.

### Factory tipado

`gui.tab_manager` declara `AutosaveFactory` con parent `QObject`,
`ChemusonCanvas`, `AutosaveUndoStack` y `AutosaveSerializer[ChemusonCanvas]`;
devuelve `AutosaveController`. El registro es
`dict[ChemusonCanvas, AutosaveController]` y `autosave_manager_for()` devuelve
`Optional[AutosaveController]`. No se usa `Callable[..., ...]` ni `object`
como sustituto del documento real.

### Documentación

La documentación distingue dependencias entre módulos ChemUSON de dependencias
externas. `autosave.py` se aísla de PyQt6, GUI, ChemIO y RDKit;
`crash_reporter.py` conserva PyQt6 para notificación visual, que no es una
dependencia entre módulos ChemUSON. `utils` no importa otros módulos ChemUSON.

### Compatibilidad

El constructor del composition root, timers, rutas, hashes, filenames, payload,
metadata, cleanup, debounce, rotación y propagación de errores mantienen su
comportamiento. La revisión añade pruebas directas de estas garantías, incluidas
las transiciones start/stop y la ejecución de callbacks de temporizador.
