# Tasks: Decouple Utils Autosave Boundaries

## 0. Baseline

- [ ] **0.1 Registrar baseline y tests relacionados**
  - Archivos: baseline de ejecución y `tests/test_autosave_manager.py`.
  - Aceptación: arquitectura y suite completa pasan antes de producción.

## 1. Reconocimiento

- [ ] **1.1 Confirmar collaborators y composition root**
  - Archivos: `src/chemuson/utils/autosave.py`, `src/chemuson/gui/tab_manager.py`.
  - Aceptación: canvas, persistencia, undo y timers usados quedan identificados.

## 2. Contratos

- [ ] **2.1 Definir contratos estructurales de autosave**
  - Archivos: `src/chemuson/utils/autosave.py`.
  - Aceptación: Protocols/callbacks tipados, sin `Any`, PyQt6, ChemIO ni GUI.

## 3. Adaptación de call sites

- [ ] **3.1 Inyectar serializador y timers concretos desde GUI**
  - Archivos: `src/chemuson/gui/tab_manager.py`.
  - Aceptación: el único constructor conserva intervalos y callbacks Qt existentes.

## 4. Tests de comportamiento

- [ ] **4.1 Caracterizar escritura, rotación y errores con fakes**
  - Archivos: `tests/test_autosave_manager.py`.
  - Aceptación: serializador, ruta, metadata, rotación y errores quedan cubiertos.

## 5. Actualización del catálogo

- [ ] **5.1 Eliminar dependencia y excepciones M15**
  - Archivos: `architecture/modules.yml`.
  - Aceptación: M15 tiene listas actuales, objetivo y excepciones vacías.

## 6. Reducción del baseline de excepciones

- [ ] **6.1 Congelar las ocho excepciones restantes**
  - Archivos: `tests/architecture/test_exceptions_no_growth.py`.
  - Aceptación: total 8, runtime 7, TYPE_CHECKING 1 y reaparición M15 rechazada.

## 7. Validación final

- [ ] **7.1 Verificar aislamiento, arquitectura, documentación y suite completa**
  - Archivos: pruebas arquitectónicas, `docs/architecture.md` y este cambio.
  - Aceptación: strict válido, tareas completas, suite y Ruff sin diagnósticos nuevos.
