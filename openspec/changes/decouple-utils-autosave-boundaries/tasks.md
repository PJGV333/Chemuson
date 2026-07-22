# Tasks: Decouple Utils Autosave Boundaries

## 0. Baseline

- [x] **0.1 Registrar baseline y tests relacionados**
  - Archivos: baseline de ejecución y `tests/test_autosave_manager.py`.
  - Aceptación: arquitectura y suite completa pasan antes de producción.

## 1. Reconocimiento

- [x] **1.1 Confirmar collaborators y composition root**
  - Archivos: `src/chemuson/utils/autosave.py`, `src/chemuson/gui/tab_manager.py`.
  - Aceptación: canvas, persistencia, undo y timers usados quedan identificados.

## 2. Contratos

- [x] **2.1 Definir contratos estructurales de autosave**
  - Archivos: `src/chemuson/utils/autosave.py`.
  - Aceptación: Protocols/callbacks tipados, sin `Any`, PyQt6, ChemIO ni GUI.

## 3. Adaptación de call sites

- [x] **3.1 Inyectar serializador y timers concretos desde GUI**
  - Archivos: `src/chemuson/gui/tab_manager.py`.
  - Aceptación: el único constructor conserva intervalos y callbacks Qt existentes.

## 4. Tests de comportamiento

- [x] **4.1 Caracterizar escritura, rotación y errores con fakes**
  - Archivos: `tests/test_autosave_manager.py`.
  - Aceptación: serializador, ruta, metadata, rotación y errores quedan cubiertos.

## 5. Actualización del catálogo

- [x] **5.1 Eliminar dependencia y excepciones M15**
  - Archivos: `architecture/modules.yml`.
  - Aceptación: M15 tiene listas actuales, objetivo y excepciones vacías.

## 6. Reducción del baseline de excepciones

- [x] **6.1 Congelar las ocho excepciones restantes**
  - Archivos: `tests/architecture/test_exceptions_no_growth.py`.
  - Aceptación: total 8, runtime 7, TYPE_CHECKING 1 y reaparición M15 rechazada.

## 7. Validación final

- [x] **7.1 Verificar aislamiento, arquitectura, documentación y suite completa**
  - Archivos: pruebas arquitectónicas, `docs/architecture.md` y este cambio.
  - Aceptación: strict válido, tareas completas, suite y Ruff sin diagnósticos nuevos.

## 8. Corrective Review

- [x] **8.1 Replace canonical specification Purpose placeholders**
  - Archivos: los cuatro specs canónicos de Fase 1 afectados.
  - Aceptación: cada Purpose describe su alcance real y no queda el placeholder heredado.
- [x] **8.2 Make autosave core and Qt composition type-safe**
  - Archivos: `src/chemuson/utils/autosave.py`, `src/chemuson/gui/tab_manager.py`.
  - Aceptación: core genérico, controller explícito, factory y registro tipados, sin proxy dinámico.
- [x] **8.3 Complete RDKit import-isolation enforcement**
  - Archivos: `tests/test_autosave_manager.py`.
  - Aceptación: AST y subprocess bloquean RDKit y sus submódulos además de GUI, ChemIO y PyQt6.
- [x] **8.4 Complete autosave lifecycle compatibility tests**
  - Archivos: `tests/test_autosave_manager.py`, `tests/test_main_window_tabs.py`.
  - Aceptación: start/stop, debounce, callbacks, cleanup, hashes, nombres, metadata y composición quedan cubiertos.
- [x] **8.5 Correct utils/PyQt6 architecture documentation**
  - Archivos: `docs/architecture.md`, `architecture/modules.yml`.
  - Aceptación: diferencia dependencias ChemUSON de PyQt6 externo y documenta autosave y crash reporter correctamente.

## 9. Corrective Validation

- [x] **9.1 Run focused autosave and architecture suites**
  - Archivos: pruebas focalizadas y arquitectura.
  - Aceptación: compileall, suites focalizadas, arquitectura y comprobaciones estructurales terminan con código 0.
- [x] **9.2 Run full pytest and Ruff baseline**
  - Archivos: `src`, `tests`, `tools`, `packaging`.
  - Aceptación: colección y suite completas pasan; Ruff conserva únicamente el F401 histórico permitido.
- [x] **9.3 Validate OpenSpec globally and prepare archive**
  - Archivos: artefactos del cambio y specs canónicos.
  - Aceptación: validación strict local y global, ocho excepciones y M15 limpio antes de archive.
- [ ] **9.4 Archive decouple-utils-autosave-boundaries**
  - Archivos: cambio archivado y spec canónico actualizado por OpenSpec.
  - Aceptación: archive único posterior a todas las validaciones y sin cambios activos relacionados.
