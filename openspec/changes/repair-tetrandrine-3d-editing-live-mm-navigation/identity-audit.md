# Auditoría de identidad: tetrandrina local

## Alcance y preservación

Se leyó la entrada exacta `Productos naturales / tetrandrina` de `~/.config/Chemuson/template_library.json` sin instanciar una biblioteca que pudiera escribir al cargar y sin modificar la configuración. El MolBlock coincide byte por byte con `tests/data/tetrandrina_local_template.mol` (SHA-256 `3cea74aee637ded3263673c3eeb7935436de2a935b6974cda3ad75af8f600b44`). La importación directa del MolBlock y `graph_from_template_payload` producen los mismos conteos y elementos; en esta entrada no se usa el SMILES de respaldo.

La referencia es tetrandrina neutra, PubChem CID 73078, fórmula `C38H42N2O6`. La comparación topológica descrita abajo ignora la estereoquímica de la referencia; no se asignó estereoquímica a la plantilla local.

## Resultados por etapa

| Etapa | Átomos / enlaces | Composición | Componentes / ciclos | Observación |
|---|---:|---|---|---|
| MolGraph local importado | 43 / 48 | `C:33`, `CH3:2`, `N:2`, `O:2`, `OCH3:4` | 1 componente; rango de ciclos 6 | Cero centros o enlaces estereoquímicos especificados. |
| Copia tras expansión de abreviaturas | 47 / 52 | `C:39`, `N:2`, `O:6` | 1 componente; rango de ciclos 6 | Cada `OCH3` se expande a O–C, agregando 1 vértice; cada `CH3` se sustituye por C en el mismo vértice. La expansión en sí conserva los enlaces de anclaje y explica exactamente el cambio de conteos. |
| Copia después de la normalización aromática actual | 47 pesados + 46 H = 93 / 98 | `C39H46N2O6` | 1 componente; rango de ciclos 6; RDKit informa 7 anillos | Sanitiza y genera coordenadas, pero no es la identidad de referencia. Presenta 18 átomos y 18 enlaces aromáticos. |
| Referencia PubChem 73078 | 46 pesados + 42 H = 88 / 52 | `C38H42N2O6` | 1 componente; rango de ciclos 7; RDKit informa 8 anillos | 24 átomos y 24 enlaces aromáticos. |

La inspección de la plantilla local localiza un defecto de **autoría en el MolGraph almacenado**, no de importación ni de expansión:

- El átomo documental `11` es un carbono terminal (grado 1) unido al átomo `5` por el enlace `9`, marcado aromático y con orden sencillo.
- Los átomos `6` y `11` tienen las mismas coordenadas 2D. El átomo `6` ya está enlazado a `1` y `7`; no existe el enlace `5–6` que cerraría el anillo esperado.
- La expansión convierte correctamente los cuatro `OCH3` en cuatro grupos metoxi y los dos `CH3` en carbonos. No añade el carbono discrepante: el átomo `11` ya está en el MolGraph original.

## Qué hace —y qué no hace— la normalización actual

En la copia de cálculo, desmarcar aromaticidad en el puente terminal `5–11` permite que RDKit sanee el grafo. El carbono terminal `11` pasa a interpretarse como metilo con tres H implícitos. El resultado es fórmula `C39H46N2O6`, una escena de 93 átomos y una malla finita. Las valencias típicas quedan satisfechas (C 4, N 3, O 2), por lo que la comprobación de valencia no descubre el defecto. Esa sanitización **enmascara**, no corrige, la conectividad incorrecta: el enlace `5–11` deja de ser aromático, pero no se crea el cierre `5–6` ni se elimina el duplicado.

Como prueba diagnóstica aislada, se ensayó sobre un objeto RDKit descartable eliminar el duplicado `11` y sustituir `5–11` por el enlace aromático `5–6`. Ese candidato da `C38H42N2O6`, 46 átomos pesados, 52 enlaces, una componente, 24 átomos/enlaces aromáticos y es isomorfo al grafo de referencia al ignorar estereoquímica. **No se aplicó esa transformación** al MolGraph, a la plantilla del usuario ni a la fixture: evidencia la causa probable, pero no autoriza una edición silenciosa del documento fuente.

## Diagnóstico y límite de aceptación

- Defecto: plantilla/molgraph local, átomo duplicado `11` y enlace de cierre mal asignado (`5–11` en vez de `5–6`).
- Importación: fiel al MolBlock local; fixture y MolBlock son idénticos.
- Expansión: conteos y grupos expandidos correctos; no explican el carbono adicional.
- Aromaticidad: el flag imposible del enlace 9 es una manifestación del defecto estructural, no un problema aislado de sanitización.
- Identidad: el conformero/escena de 93 átomos generado por el código actual corresponde al grafo local defectuoso, no a tetrandrina `C38H42N2O6`.
- No cambiar la entrada local sin autorización explícita. La aceptación química y gráfica de tetrandrina queda pendiente de una corrección aprobada del MolGraph fuente y de una nueva validación; no se afirma que OpenGL la dibuje correctamente.
