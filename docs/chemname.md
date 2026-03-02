# Motor de Nomenclatura (ChemName)

Este documento resume el alcance práctico del motor `chemuson.chemname` (IUPAC-lite), con foco en estabilidad y degradación segura (`N/D`) en lugar de cobertura completa Blue Book.

## Cobertura principal

- Cadenas lineales y cicloalcanos simples.
- Aromáticos comunes: benceno, naftaleno, antraceno, fenantreno, pireno.
- Heteroaromáticos: furanos/tiofenos/pirroles, piridinas/diazinas, triazoles, tetrazoles, isoxazoles, oxazinas y variantes fusionadas seleccionadas.
- Sistemas fusionados seleccionados: benzotriazol, quinolinas/isoquinolinas/indol, quinonas aromáticas (`benzene-1,4-dione`, `naphthalene-1,2/1,4-dione`).

## Grupos funcionales soportados

- Alcoholes, aminas, carbonilos básicos, ácidos carboxílicos/ésteres/amidas/nitrilos.
- Tiol (`thiol`), tióteres/sulfuros (prefijo), sulfóxidos, sulfonas.
- Ácido sulfónico / sulfonato / sulfonamida.
- Azida (`azido`) y peróxido (`peroxy`).
- Detección iónica básica para carboxilatos/sales comunes según carga formal.

## Cargas, isótopos y radicales

- Modelo con `formal_charge`, `isotope`, `radical_electrons`, `oxidation_state`.
- Round-trip en `.cmsn` y, cuando aplica, en Molfile/RDKit.
- Render mínimo de isótopos y radicales (con fallback seguro).

## Plantillas especiales (fase A)

- Carbohidratos (exact/subestructura best-effort):  
  `alpha-d-glucopyranose`, `beta-d-glucopyranose`, `beta-d-fructofuranose`, `d-ribose`.
- Esteroides (MVP): `androstane`, `cholestane`.
- Sustitución simple sobre plantilla (MVP): `hydroxy`, `oxo`, `amino`; si no hay cobertura, se degrada a `"<base> (substituted)"`.

## Coordinación (experimental)

- Detección de metal central y ligandos mínimos:
  `carbonyl`, `ammine`, `aqua`, halo, `cyano`, `η5-cyclopentadienyl`.
- Estado de oxidación heurístico.
- Casos básicos `cis/trans` en complejos tipo Pt(II) cuadrado plano.

## Estereoquímica avanzada (fase B, best-effort)

- Soporte de descriptores: `M/P`, `R_a/S_a`, `endo/exo`, `si/re`.
- Prioridad de render: `M/P` > `R_a/S_a` > `R/S/E/Z` > `endo/exo/si/re`.
- Fuente:
  - RDKit aislado en subproceso si está disponible.
  - Anotaciones del modelo como fallback.
- Si no hay información confiable, se omite descriptor (sin excepción).

## Flags relevantes (`NameOptions`)

- `enable_experimental`
- `enable_special_templates`
- `enable_advanced_stereo`
- `allow_coordination`
- `enable_exotic_hetero`
- `rdkit_isolated`
- `return_nd_on_fail`

## Limitaciones conocidas

- No cubre toda la nomenclatura IUPAC oficial.
- Plantillas especiales y coordinación son MVP orientados a casos comunes.
- En ausencia de RDKit, funcionalidades avanzadas dependen de anotaciones/plantillas exactas.
- Ante ambigüedad o no-soporte, el motor prioriza estabilidad y retorna `N/D` (o excepción en modo `strict`).

