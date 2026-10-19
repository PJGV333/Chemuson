# Spec: Integración de ramas para Beta 1

## ADDED Requirements

### Requirement: Se preservan los dos historiales de implementación
La rama temporal de integración SHALL conservar como ancestros tanto el HEAD de la rama visual Beta 1 fijado en el baseline como el HEAD de la rama `feature/chemuson-3d`, mediante un merge normal, sin reescritura ni selección parcial de commits.

#### Scenario: Merge aislado
- **GIVEN** la rama local temporal está en el SHA fijado de `origin/feat/molecular-depiction-visual-lab`
- **WHEN** se integra la referencia 3D con un merge commit local
- **THEN** ambos SHA originales son ancestros del resultado y las ramas de origen/remotas no cambian

### Requirement: La integración conserva las capacidades de 2D y 3D
La integración SHALL preservar los contratos existentes de dibujo/depiction, Clean2D y proyecciones visuales de Beta 1, además de las capacidades 3D de viewport, edición, Live MM, ChemIO y persistencia documentadas por sus OpenSpecs. Los conflictos compartidos SHALL resolverse y probarse según ambos contratos; no basta con aceptar una estrategia global `ours` o `theirs`.

#### Scenario: Regresión cruzada
- **GIVEN** se han resuelto los conflictos de ChemIO, arquitectura y UI
- **WHEN** se ejecutan las regresiones focales de 2D, 3D, serialización/importación y arquitectura
- **THEN** las capacidades previas continúan pasando y las nuevas capacidades 3D funcionan sin cambios químicos o de formato no autorizados

### Requirement: Los gates de release siguen siendo explícitos
La rama de release SHALL permanecer intacta mientras cualquier validación crítica, shard requerido o aceptación manual del propietario esté pendiente o fallando. Una ejecución parcial SHALL reportarse como parcial y no como aceptación de Beta 1.

#### Scenario: Gate incompleto
- **GIVEN** la suite, una regresión crítica o la aceptación manual aún no tiene resultado satisfactorio
- **WHEN** se reporta el resultado del merge temporal
- **THEN** el estado queda como bloqueado/no listo y no se modifica `release/v0.3.0-beta.1-prep`
