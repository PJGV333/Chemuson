# Proposal: Establish Phase 1 Module Catalog and Contracts

## Context

ChemUSON está evolucionando hacia una arquitectura más modular. Para que futuros refactors (extraer lógica de la GUI, aislar motores químicos) sean seguros y predecibles, necesitamos un "contrato" verificable del estado actual. Actualmente, las dependencias se gestionan por convención y documentación, pero carecen de una fuente de verdad única y de ejecución automática.

Esta fase crea los cimientos: catálogo estructurado, reglas de límites, protocolo de agentes y documentación del canvas.

## Objetivo

Establecer una base verificable y estable para la modularidad de ChemUSON mediante:

1. Catálogo de módulos en formato YAML (`architecture/modules.yml`) como fuente de verdad estructural.
2. Documentación de responsabilidades, APIs y rutas de todos los paquetes principales.
3. Reglas de límites arquitectónicos (imports permitidos/prohibidos) con excepciones explícitas.
4. Diseño de tests automáticos (basados en AST) para ejecutar estas reglas.
5. Protocolo estandarizado para agentes IA mediante `AGENTS.md`.

Esta fase es **estrictamente documental y de validación**. No se moverá ni renombrará código de producción, ni se cambiará comportamiento funcional.

## Alcance

### Incluido

- **Catálogo de módulos**: Identificación de todos los bloques funcionales principales con IDs estables. El inventario se basa en los paquetes reales de `src/chemuson/`:
  - Paquetes de dominio: `core`, `chemio`, `clean2d`, `chemcalc`, `chemname`, `geometry3d`, `compchem`, `spectroscopy`, `name2structure`, `markush`.
  - Subsistema GUI (desglosado): `gui` (orquestación general), `gui/canvas`, `gui/controllers`, `gui/commands`, `gui/dialogs`, `gui/items`.
  - Utilidades: `update`, `utils`.
  - Módulos raíz: `version`/`_version`.
  - El directorio `templates/` está vacío y NO se cataloga como módulo funcional.
- **Definición de contratos**: Mapeo de paths, APIs públicas/internas y reglas de dependencia (actual vs. objetivo).
- **Mapeo del Canvas**: Documentación detallada de `src/chemuson/gui/canvas/` para preparar futuras extracciones.
- **Protocolo de agentes**: `AGENTS.md` raíz para guiar agentes IA en el mantenimiento de integridad arquitectónica.
- **Diseño de tests de límites**: Tests basados en análisis AST para detectar fugas de dependencia.
- **Captura de baseline**: Estado verificable del repositorio antes de empezar.

### Excluido (Prohibido)

- Mover o renombrar cualquier archivo en `src/chemuson/`.
- Modificar algoritmos químicos (Clean2D, ChemName, etc.).
- Alterar señales, eventos o renderizado de Qt.
- Cambiar formatos de persistencia o lógica undo/redo.
- Introducir dependencias o frameworks nuevos.
- Crear paquetes nuevos (p. ej. `app/`, `depiction/`).
- Corregir bugs existentes u oportunistas.
- Actualizar snapshots de tests.

## Resultados Esperados

- `architecture/modules.yml`: fuente de verdad para la topología de módulos.
- `docs/modules/`: directorio con contexto detallado por módulo.
- `AGENTS.md`: protocolo para agentes IA.
- `tests/architecture/`: suite de tests arquitectónicos nuevos.
- Mapa documentado del subsistema Canvas.
- Reporte de baseline registrado.

## Evaluación de Riesgos

- **Riesgo bajo**: El cambio se limita a documentación, tests nuevos y archivos YAML.
- **Mitigación**: Todo el trabajo está contenido en `openspec/changes/...`, `architecture/`, `docs/modules/`, `tests/architecture/` y `AGENTS.md`. Ninguna lógica existente se toca.
- **Riesgo residual**: El catálogo puede contener errores de interpretación que se detectarán durante la ejecución de tests. Los tests pueden fallar inicialmente hasta registrar excepciones correctas.
