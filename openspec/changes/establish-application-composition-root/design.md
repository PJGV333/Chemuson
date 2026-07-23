## Context

`src/chemuson/__main__.py` implementa un CLI pequeño y sólo importa la GUI
después de procesar `--version`. Sin embargo, el objeto `QApplication`, la
ventana y el manejo de errores se coordinan mediante `run_app()` al final de
`src/chemuson/gui/main_window.py`.

Estado actual:

```text
chemuson.__main__.main
  └─ lazy import gui.main_window.run_app
       ├─ crash_reporter.install
       ├─ QApplication(sys.argv)
       ├─ application name/version
       ├─ ChemusonWindow
       ├─ check_autosaves
       ├─ show
       ├─ app.exec
       └─ crash reporting / sys.exit
```

Estado objetivo:

```text
chemuson.__main__.main
  └─ lazy import app.bootstrap.run_app
       ├─ crash_reporter.install
       ├─ QApplication(sys.argv)
       ├─ application name/version
       ├─ gui.main_window.ChemusonWindow
       ├─ check_autosaves
       ├─ show
       ├─ app.exec
       └─ crash reporting / sys.exit
```

El paquete `app/` es una capa superior de composición, no un nuevo dominio ni
un contenedor genérico de servicios.

## Goals / Non-Goals

**Goals:**

- Hacer que M19 posea explícitamente el arranque y el ciclo de vida.
- Mantener `__main__.py` importable sin PyQt6, GUI, RDKit u otros subsistemas
  pesados.
- Preservar exactamente el orden, mensajes y efectos observables de
  `run_app()`.
- Mantener `ChemusonWindow` en su ruta pública actual.
- Registrar ownership, API, dependencias y tests de M19 en el catálogo.
- Dejar una base pequeña para el cambio posterior
  `extract-application-shell`.

**Non-Goals:**

- Extraer menús, toolbars, docks, acciones o lifecycle interno de la ventana.
- Dividir `main_window.py`, `items.py` o cualquier mixin del canvas.
- Crear service locator, dependency-injection framework, plugin registry o
  abstracciones especulativas.
- Modificar autosave, crash reporting, packaging, formato CMSN, renderizado o
  lógica química.
- Cambiar el resultado visual, la secuencia de recuperación o los códigos de
  salida.

## Decisions

### M19 owns a small app package

Se añadirá `src/chemuson/app/__init__.py` y un módulo de bootstrap explícito.
El paquete completo pertenecerá a M19 junto con `src/chemuson/__main__.py`.
`app/__init__.py` no reexportará símbolos que obliguen a importar Qt al cargar
el paquete; el dispatcher importará el módulo de bootstrap sólo cuando deba
iniciar la GUI.

### The CLI remains the public entry point

`chemuson.__main__:main` seguirá siendo la API pública y el target de
`[project.scripts]`. `_build_parser()` permanece interno. El CLI procesa
`--version` antes de importar el bootstrap gráfico y después delega una sola
vez en `app.bootstrap.run_app()`.

`run_app()` será una operación interna de composición, no una nueva API pública
del producto. No se mantendrá un alias de compatibilidad en
`gui.main_window`: la búsqueda del repositorio deberá demostrar que no existen
consumidores deliberados fuera de `__main__.py` y pruebas.

### Preserve lifecycle literally before improving it

La implementación moverá la secuencia existente con cambios mínimos:

1. instalar crash reporter;
2. crear `QApplication(sys.argv)`;
3. establecer nombre y versión;
4. crear `ChemusonWindow`;
5. comprobar autosaves;
6. mostrar la ventana;
7. ejecutar el event loop;
8. conservar el reporte y mensaje de error actuales;
9. terminar con el mismo exit code.

No se corregirán patrones históricos (`sys.exit`, catch amplio o mensajes) en
este cambio. Cualquier mejora posterior requerirá otro OpenSpec.

### Import isolation is enforced without a display server

Las pruebas del CLI usarán subprocess o módulos sustitutos para verificar que
importar `chemuson.__main__`, pedir `--version` o pedir ayuda no carga
`chemuson.app.bootstrap`, `chemuson.gui` ni PyQt6. Pruebas AST verificarán
ownership, dirección de dependencias y ausencia de `run_app` en
`main_window.py`.

La caracterización del orden de lifecycle podrá usar colaboradores falsos o
monkeypatching del bootstrap, sin construir widgets reales. Las pruebas GUI
existentes seguirán siendo la evidencia de comportamiento visual en entornos
con Qt completo.

### Catalog dependencies follow direct imports

M19 declarará ownership sobre `src/chemuson/app/` y
`src/chemuson/__main__.py`. Sus dependencias actuales y objetivo serán:

- M08 para `ChemusonWindow`;
- M15 para `crash_reporter`;
- M18 para `get_app_version`.

M00-M18 continuarán sin referenciar M19. Después de retirar el bootstrap de
`main_window.py`, las dependencias de M08 sólo cambiarán si el análisis AST
demuestra que M15 o M18 ya no tienen otros consumidores directos en M08.

## Risks / Trade-offs

- [Import eager accidental] → pruebas de subprocess inspeccionan
  `sys.modules` y bloquean la carga de Qt/GUI en rutas CLI no gráficas.
- [Cambio de orden de startup] → caracterización con fakes registra cada paso y
  exige una única ejecución en el orden histórico.
- [Ownership solapado] → pruebas del catálogo verifican que `app/` pertenece
  sólo a M19 y que todos sus paths existen.
- [API fantasma de compatibilidad] → búsqueda y AST garantizan que
  `main_window.py` ya no define ni reexporta `run_app`.
- [Refactor se expande al shell] → el diff queda limitado al bootstrap,
  catálogo, documentación y pruebas de arranque.

## Migration Plan

1. Capturar baseline y consumidores reales de `run_app`.
2. Añadir caracterización no Qt del CLI y del lifecycle.
3. Crear el paquete M19 y mover literalmente el bootstrap.
4. Reducir `__main__.py` al dispatch diferido y retirar `run_app` de M08.
5. Actualizar catálogo, ficha M19 y arquitectura.
6. Ejecutar validación focalizada, arquitectura, OpenSpec y suite disponible.

Rollback: revertir el cambio restaura `run_app()` en `main_window.py`; no hay
migración de datos, configuración ni packaging.

## Open Questions

No hay preguntas de diseño abiertas. El nombre del módulo interno de bootstrap
podrá confirmarse durante el reconocimiento, pero no altera ownership, API ni
contratos descritos aquí.
