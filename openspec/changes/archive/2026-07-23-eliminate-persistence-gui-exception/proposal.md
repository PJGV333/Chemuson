## Why

ChemIO (M01) está aislado conceptualmente de la GUI, pero `persistence.py`
mantiene un import de `ChemusonCanvas` bajo `TYPE_CHECKING` para anotar una
superficie que ya es estructural. Además, `PersistenceManager` invoca un método
privado del canvas. Esta es la última `temporary_exception` del catálogo e
impide alcanzar un baseline arquitectónico vacío.

## What Changes

- Definir un contrato estructural GUI-neutral para el documento de persistencia.
- Hacer que `PersistenceManager` dependa de ese contrato sin nombrar clases GUI.
- Exponer un hook público mínimo del canvas para terminar la reconstrucción visual.
- Eliminar la dependencia M01 a M09, incluso bajo `TYPE_CHECKING`.
- Congelar el baseline de excepciones en cero y evitar la reaparición de la
  identidad histórica.
- Mantener sin cambios el formato CMSN, `VERSION`, atomicidad, recuperación,
  autosave y flujos de archivo.

## Capabilities

### New Capabilities
- `persistence-contract`: Contrato estructural GUI-neutral para guardar y restaurar documentos CMSN.

### Modified Capabilities
- `architecture-boundaries`: ChemIO deja de conocer GUI y no quedan excepciones temporales.
- `module-catalog`: M01 y el catálogo completo reflejan una topología sin deuda temporal.

## Impact

- Afecta `chemuson.chemio.persistence` y el mixin de estructura del canvas.
- Conserva las rutas públicas y métodos de `PersistenceManager`.
- Actualiza pruebas de persistencia, límites de importación, catálogo y documentación arquitectónica.
- No incorpora dependencias externas ni modifica los controllers, autosave o el formato CMSN.
