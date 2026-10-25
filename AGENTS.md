# Protocolo de Agentes para ChemUSON

Este documento establece las reglas vinculantes para cualquier agente de IA (LLM) que trabaje en este repositorio. El objetivo es garantizar que la integridad arquitectónica se mantenga durante la evolución del software.

## 1. El Ciclo de Trabajo del Agente

Antes de realizar cualquier modificación en el código o la documentación, el agente debe seguir estrictamente este flujo:

### 1.1 Análisis de Alcance (OpenSpec)
Todo trabajo comienza leyendo la especificación OpenSpec activa. El agente **debe** leer:
- `proposal.md`: Para entender el propósito y contexto.
- `design.md`: Para entender las decisiones de diseño y restricciones.
- `tasks.md`: Para identificar las tareas específicas a realizar.
- `specs/*.md`: Para entender los contratos técnicos detallados.

### 1.2 Captura de baseline proporcional al alcance
Antes de modificar código, pruebas, empaquetado o documentación, registrar `git status --short`, rama, HEAD y las verificaciones relevantes **antes** del cambio para poder comparar resultados. Aplicar la [política de pruebas y compilación](docs/testing/AGENT_TESTING_POLICY.md), que forma parte de este protocolo.

- **Código Python y empaquetado:** baseline con `python -m compileall src tests tools packaging`, pruebas focalizadas del área modificada y `tests/architecture/` cuando cambie código de producción o dependencias; Ruff focal sobre los archivos pertinentes.
- **Colección de pruebas:** usar `pytest --collect-only -q` y verificar el manifiesto de node IDs cuando el cambio afecte tests, su colección, parametrización o distribución. No es necesario repetir esa colección para una corrección puramente documental.
- **Documentación solamente:** comprobar rutas, enlaces, coherencia de los contratos descritos y `git diff --check`; no ejecutar pytest/compileall de forma innecesaria.
- **Errores históricos:** registrar el baseline real e identificar las excepciones existentes; no declararlas PASS, ocultarlas ni alterar evidencias antiguas.
- **Registro:** guardar comandos, resultados, entorno y SHA exacto en `AGENT_REPORT.md` o en el `baseline.md` / `validation.md` del OpenSpec correspondiente. No atribuir a CI resultados obtenidos localmente.

**La suite completa `pytest -q` no forma parte del baseline obligatorio para cambios ordinarios.** Se conserva íntegra y está disponible para auditoría extraordinaria, justificada por el riesgo/alcance y autorizada expresamente por el propietario. El workflow `test.yml` permanece manual (`workflow_dispatch`), con sus validaciones fail-closed.

### 1.3 Ejecución y verificación
- **Focalización:** los cambios deben ser quirúrgicos y limitarse al alcance de la tarea y del OpenSpec cuando corresponda.
- **Pruebas:** después de cada modificación, ejecutar las pruebas relevantes; al finalizar, repetir los controles focalizados necesarios y `git diff --check`. No ejecutar automáticamente la suite completa de ~2.600 pruebas.
- **Integridad:** investigar los fallos nuevos; corregir regresiones o registrar una excepción legítima sin rebajar aserciones ni modificar catálogos para obtener un falso PASS.
- **Gates oficiales:** no desactivar controles de identidad, procedencia, seguridad, checksum ni smokes acotados de los workflows de publicación. Construir un paquete no requiere incluir `tests/` en el ejecutable ni ejecutar el suite completo durante su compilación.

## 2. Reglas de Integridad Arquitectónica

### 2.1 Prohibiciones (No permitir bajo ninguna circunstancia)
- **Refactor oportunista**: No corregir bugs, mejorar estilo de código, o reorganizar archivos que no estén explícitamente en el alcance del OpenSpec.
- **Ocultar fallos**: No modificar snapshots de tests, archivos de baseline o la lista de excepciones en `architecture/modules.yml` simplemente para que los tests pasen. Las excepciones deben ser documentadas y justificadas.
- **Modificación de comportamiento**: No alterar la lógica química (Clean2D, ChemName), el renderizado de la GUI o la serialización de datos sin un OpenSpec que lo autorice.
- **Creación de paquetes nuevos**: No añadir directorios como `app/`, `services/`, etc., en `src/chemuson/` sin aprobación.

### 2.2 Gestión de Dependencias
- **Catálogo de Módulos**: Cualquier nueva dependencia (import) o cambio en la estructura de paquetes debe registrarse en `architecture/modules.yml`.
- **Excepciones**: Si un import viola las reglas de arquitectura (ej. `core` importando `gui`), el agente debe crear una entrada en `temporary_exceptions` con: `source_id`, `target_id`, `file`, `import_path`, `reason`, `debt_ref`, `elimination_condition` y `type_checking_only`.
- **Circularidades**: Si se detecta una nueva dependencia circular, debe registrarse en `circular_dependencies` del catálogo.

## 3. Reporte de Trabajo

### 3.1 Registro de Cambios
Cada commit o reporte de tarea debe incluir:
- Una lista clara de todos los archivos modificados.
- Un resumen de las decisiones tomadas respecto a la arquitectura.
- El resultado de las ejecuciones de baseline y tests.

### 3.2 Informe de Desviaciones (`AGENT_REPORT.md`)
Si el agente detecta que para cumplir con el OpenSpec es necesario violar una regla arquitectónica o si encuentra un problema no previsto, debe crear/actualizar un `AGENT_REPORT.md` explicando la situación y solicitando una decisión de diseño.

## 4. Protección de Subsistemas Críticos

Se requiere precaución extrema al trabajar con:
- **GUI (`src/chemuson/gui/`)**: No alterar la jerarquía de mixins, el orden de despacho de eventos de Qt o la lógica de los controllers.
- **Clean2D (`src/chemuson/clean2d/`)**: No modificar heurísticas de dibujo o reglas de validación sin un proceso de baseline/regression completo.
- **ChemName (`src/chemuson/chemname/`)**: No alterar las reglas de nomenclatura ni las plantillas de sustituyentes.
- **Persistencia (`src/chemuson/chemio/persistence.py`)**: No cambiar el formato de los archivos `.cmsn`.

## 5. Criterio de Parada
El agente debe detener su ejecución y solicitar ayuda si:
1. Los tests existentes fallan de forma inesperada.
2. El alcance de la tarea requiere mover archivos de producción o renombrar paquetes.
3. Se detectan dependencias no catalogadas que bloquean el progreso.
4. El trabajo requiere introducir dependencias externas nuevas.
