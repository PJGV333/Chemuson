# M00: core

**Responsabilidad**: Modelo molecular fundamental. Define las estructuras base de `Atom`, `Bond`, `MolGraph` y las reglas de validación estructural (`ValidationIssue`, `ValidationCorrectionAction`).

**APIs Principales**:
- `chemuson.core.Atom`
- `chemuson.core.Bond`
- `chemuson.core.MolGraph`
- `chemuson.core.ValidationIssue`

**Riesgos**:
Es el módulo con mayor criticidad. Cualquier cambio en las estructuras base (`Atom`, `Bond`) tiene un impacto masivo en todos los subsistemas. Se requiere máxima estabilidad y cobertura de tests.
