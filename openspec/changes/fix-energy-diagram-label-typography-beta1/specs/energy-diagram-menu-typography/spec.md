# Especificación — tipografía de menús de diagramas electrónicos

## Propósito

Mostrar notación química y configuraciones científicas correctas en los controles de selección de diagramas electrónicos, manteniendo los contratos internos y datos de diagramas existentes.

## ADDED Requirements

### Requirement: Presentación tipográfica de presets
La UI SHALL mostrar las fórmulas diatómicas, cargas, hibridaciones y conteos d especificados mediante subíndices/superíndices Unicode: `B₂`, `C₂`, `N₂`, `O₂`, `F₂`, `O₂⁻`, `O₂²⁻`, `Híbrido sp`, `Híbrido sp²`, `Híbrido sp³`, `d⁴`, `d⁵`, `d⁶`, `d⁸`. SHALL preservar el resto de la descripción de geometría/espín y SHALL mantener `2s`/`2p` en línea base.

#### Scenario: Opciones moleculares y campo de ligandos
- **GIVEN** el menú de presets electrónicos está abierto
- **WHEN** se muestran las opciones MO y de campo de ligandos
- **THEN** sus títulos visibles contienen la notación tipográfica aprobada y el texto de geometría/espín intacto
- **AND** el formato no introduce HTML literal ni glifos ausentes/recortados

### Requirement: Identidad interna y datos científicos estables
El formato visible SHALL estar separado de los IDs y datos científicos. Los IDs originales SHALL seguir seleccionándose/emitiéndose; los presets SHALL construir los mismos diagramas semánticos con electrones, carga, ocupaciones, energías, degeneración, orden y geometría idénticos. La presentación SHALL NOT modificar títulos ni payloads persistidos.

#### Scenario: Selección de superóxido/peróxido e hibridación
- **GIVEN** una acción visible `O₂²⁻` o `Híbrido sp³`
- **WHEN** el usuario la activa
- **THEN** se emite exactamente `mo:O2^2-` o se selecciona `hybrid_sp3`, respectivamente
- **AND** su diagrama/datos siguen coincidiendo con el preset de la rama base

#### Scenario: Documento existente
- **GIVEN** un documento persistido con diagramas/títulos en la línea base
- **WHEN** se abre o guarda con esta versión
- **THEN** el formato de persistencia y los títulos existentes permanecen iguales
