# Refactor fase 3

## Objetivo
Seguir reduciendo `main_window.py` delegando flujo de plantillas y SMILES a un controlador específico.

## Cambios
- Se crea `TemplateController` para operaciones de plantillas y categorías.
- `main_window.py` delega inserción/guardado/import/export/renombrado/borrado de plantillas.
- `main_window.py` delega import/export de SMILES y cadena lineal.
- `tests/test_clean_2d_strategy.py` deja de importar `main_window` y valida la estrategia vía `Clean2DController`.

## Notas
- Se mantienen `gui/style.py` y `gui/styles.py` como módulos separados por responsabilidad.
