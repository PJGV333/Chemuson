## Context

`chemuson.chemio.persistence` serializa el `MolGraph` y delega el payload
visual al documento concreto. Su anotación bajo `TYPE_CHECKING` hacia
`ChemusonCanvas` crea la última dependencia M01 a M09; la carga también llama
al método privado `_rebuild_items_from_model`. FileController, RecoveryController
y autosave ya inyectan un canvas concreto en runtime y no requieren que ChemIO
conozca su clase.

Antes:

```
chemio.persistence --TYPE_CHECKING--> gui.canvas.ChemusonCanvas
chemio.persistence --> canvas._rebuild_items_from_model()
```

Después:

```
chemio.persistence --> PersistenceDocument (Protocol local)
gui.canvas structurally satisfies PersistenceDocument
chemio.persistence --> document.rebuild_persistence_view()
```

La inyección de un objeto GUI en runtime no es una dependencia de módulo desde
ChemIO hacia GUI.

## Goals / Non-Goals

**Goals:**

- Eliminar todo import M01 a M09, incluido `TYPE_CHECKING`.
- Definir un contrato pequeño junto a `PersistenceManager`.
- Preservar exactamente el payload CMSN, `VERSION`, parsing legado, escritura
  temporal y `os.replace`.
- Añadir un hook público que delegue una sola vez a la reconstrucción existente.
- Reducir catálogo y baseline de excepciones de uno a cero, sin ciclos.

**Non-Goals:**

- Rediseñar CMSN, Canvas, FileController, RecoveryController o autosave.
- Mover persistencia a GUI o Core, introducir adaptadores, dependencias o dos
  contratos especulativos.
- Limpiar parsing legado, excepciones amplias históricas, Ruff histórico o
  placeholders fuera de alcance.

## Decisions

### One structural contract

`PersistenceDocument` vive en `chemuson.chemio.persistence` y declara sólo:
`model: MolGraph`, `get_persistence_data`, `load_persistence_data` y
`rebuild_persistence_view`. Un único Protocol representa la superficie real de
los métodos públicos de guardado y carga; separar save/load duplicaría nombres
sin ventaja concreta. No usa `@runtime_checkable`, `isinstance`, herencia
nominal, `Any`, `object` o adapters.

### Public reconstruction hook

`CanvasStructureMixin` posee la implementación privada, por lo que añade
`rebuild_persistence_view()` junto a la persistencia y delega directamente a
`_rebuild_items_from_model()`. `PersistenceManager` llama exclusivamente al
hook tras restaurar átomos, enlaces, IDs y payload de canvas. No cambia el
orden ni duplica la reconstrucción.

### Enforcement and compatibility

Pruebas puras ejercitan el Protocol con un fake de `MolGraph` sin Qt; pruebas
de canvas verifican la delegación. El analizador AST bloquea cualquier M01 a
M09 en todos los ámbitos y las pruebas de excepción mantienen vacía la base
congelada, incluida la identidad histórica normalizada.

## Risks / Trade-offs

- [Orden de carga sensible] → El fake registra que canvas payload e IDs existen
  antes de la única reconstrucción; regresiones reales siguen cubiertas.
- [Compatibilidad del serializer de autosave] → La firma sigue siendo callable
  para `ChemusonCanvas` por tipado estructural y el test de TabManager conserva
  la inyección existente.
- [Reaparición de la deuda] → Tests del catálogo, analyzer real/sintético y
  baseline vacío fallan ante la identidad histórica.

## Migration Plan

1. Caracterizar el comportamiento CMSN y aislamiento de imports.
2. Sustituir la anotación GUI por el Protocol local y añadir el hook público.
3. Vaciar catálogo y baseline, reforzar límites y documentación.
4. Validar persistencia, arquitectura y suite completa. El rollback consiste en
   revertir los commits de esta fase; no existe migración de datos.

## Open Questions

No hay preguntas abiertas: el reconocimiento confirmó que no existe otra
implementación del documento ni un hook público previo apropiado.
