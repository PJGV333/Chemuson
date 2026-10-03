# UI Modernization Final QA & Release Closure Specification

## Purpose

Define Fase 8 como aseguramiento de calidad y cierre documental/de release de la modernización UI ya aprobada, sin añadir funcionalidad.

## ADDED Requirements

### Requirement: Fase 8 SHALL ser QA y release closure, no desarrollo funcional

Fase 8 SHALL verificar la integración UI como Release Candidate, comparar resultados contra una baseline registrada, auditar el empaquetado existente, generar evidencia final, actualizar documentación y archivar formalmente los OpenSpec de Fases 1–7. Fase 8 SHALL NOT añadir funcionalidad, corregir Clean2D o química/geometría de plantillas, iniciar Fase 9 ni ocultar fallos históricos. Una regresión nueva o pérdida de assets en la distribución SHALL detener el cierre hasta una decisión/corrección mínima dentro del alcance aprobado.

#### Scenario: Baseline y fallos históricos
- **GIVEN** los comandos y resultados registrados en `baseline.md`
- **WHEN** se ejecuta la validación final
- **THEN** los fallos se comparan por identidad con la baseline
- **AND** un fallo nuevo bloquea el cierre
- **AND** un fallo histórico se reporta sin ocultarse ni corregirse fuera de alcance.

#### Scenario: Empaquetado usa el mecanismo del proyecto
- **GIVEN** los manifiestos y scripts de distribución existentes
- **WHEN** se prueba instalación o artefacto fuera del checkout
- **THEN** sólo se usa el mecanismo soportado por el repositorio
- **AND** se comprueba que los recursos de UI se distribuyen
- **AND** no se añade un sistema de packaging nuevo.

#### Scenario: Cierre documental y archivado oficial
- **GIVEN** QA y los siete OpenSpec UI en strict válido
- **WHEN** se cierra Fase 8
- **THEN** la evidencia y documentación reflejan los resultados reales
- **AND** Fases 1–7 se archivan mediante el comando oficial OpenSpec
- **AND** este cambio permanece activo y validable hasta completar su propia fase.

## OUT OF SCOPE

- Nuevas capacidades de UI o química.
- Fase 9, algoritmo Clean2D, ChemName, geometría/hit-testing del canvas, molecular graph, molblocks y correcciones de plantillas (Fischer, Haworth, silla, tetrandrina).
- Bump de versión sin proceso contractual confirmado.