# QT-TEARDOWN — cierre seguro del stack de undo

## Purpose

Definir el resultado observable y la evidencia causal requerida cuando una prueba del shard 6 aborta dentro de `QUndoStack::cleanChanged` al destruir widgets.

## ADDED Requirements

### Requirement: Deferred-delete probes are scoped to the object under test

Una prueba que fuerza la entrega de `DeferredDelete` para verificar la destrucción de un widget SHALL pasar ese QObject como receiver a `QCoreApplication.sendPostedEvents`; MUST NOT vaciar globalmente los deletes pendientes de objetos creados por otros tests.

#### Scenario: Onboarding overlay is explicitly released
- **WHEN** `test_onboarding_close_releases_overlay` solicita la eliminación diferida del overlay
- **THEN** sólo entrega el `DeferredDelete` de ese overlay antes de comprobar `sip.isdeleted(overlay)`
- **AND** deja intactos los eventos pendientes de otros canvases/undo stacks.

#### Scenario: Ordered shard contains prior Qt widgets with deferred deletion
- **WHEN** la prueba de onboarding sigue a tests que dejaron objetos Qt pendientes
- **THEN** su cleanup no fuerza la destrucción de esos objetos ajenos
- **AND** el shard alcanza JUnit/timings completos sin SIGABRT por `QUndoStack::~QUndoStack`.

#### Scenario: Molecular-assistant dialog lifecycle flushes its own deferred delete
- **WHEN** `test_molecular_assistant_lifecycle.py` verifica el cierre/destrucción de un diálogo o ventana
- **THEN** `_process_qt_deletes` recibe ese diálogo/ventana como receiver y despacha su `DeferredDelete` antes de bombear eventos ordinarios
- **AND** se conservan las aserciones de `sip.isdeleted`, registros vacíos, cierre de ventana y ausencia de workers residuales
- **AND** no se entregan deletes pendientes de QUndoStack u otros widgets creados por tests previos.

#### Scenario: Isolated Qt probes delete the object under test
- **WHEN** un probe de shutdown/interoperabilidad llama `deleteLater()` sobre un widget, canvas o ventana concreto
- **THEN** `sendPostedEvents` recibe ese mismo QObject
- **AND** la limpieza no consume deletes ajenos aun si el probe se ejecuta dentro de un proceso dedicado.

### Requirement: Ordered reproduction distinguishes source regressions from test-order effects

La suite SHALL comprobar el mismo orden de IDs en snapshots del SHA base y del SHA candidato, con igual plataforma/runtime y sin worktrees permanentes. El análisis MUST registrar ejecución parcial/crash y MUST NOT cambiar skips/timeouts para producir PASS.

#### Scenario: Current shard reproduces SIGABRT
- **WHEN** el proceso aborta antes de generar JUnit
- **THEN** se conservan el exit/signal, el test iniciado, los últimos timings y el C stack; el shard y verifier siguen en failure.

#### Scenario: Base and current behavior differ
- **WHEN** el set/orden común de IDs se ejecuta en `61f7647b` y el SHA candidato
- **THEN** la diferencia se atribuye a fuente sólo si importa la ruta/vida útil del objeto; runtime y entorno se mantienen identificados y comparables.

### Requirement: Final CI evidence is complete and fail-closed

El cambio SHALL considerarse verificado sólo cuando los ocho shards del SHA final ejecutan todos los IDs exactamente una vez, sin fallo, timeout o señal fatal, y los smoke Windows/Flatpak requeridos pasan en la misma referencia.

#### Scenario: A shard aborts or omits JUnit
- **WHEN** pytest aborta antes de emitir resultados completos
- **THEN** el dictamen permanece NO-GO aunque los otros siete shards pasen.

#### Scenario: Smoke jobs pass but tests fail
- **WHEN** Windows smoke y Flatpak smoke terminan en PASS pero un shard pytest falla
- **THEN** esos smoke se registran como PASS automatizado separado y no vuelven verde la suite.
