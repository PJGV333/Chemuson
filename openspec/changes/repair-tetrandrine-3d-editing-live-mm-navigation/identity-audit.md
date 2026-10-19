# Auditoría de identidad: tetrandrina local

## Alcance y preservación

Se leyó la entrada exacta `Productos naturales / tetrandrina` de `~/.config/Chemuson/template_library.json` sin instanciar una biblioteca que pudiera escribir al cargar. Antes de la reparación, el MolBlock coincidía byte por byte con `tests/data/tetrandrina_local_template.mol` (SHA-256 `3cea74aee637ded3263673c3eeb7935436de2a935b6974cda3ad75af8f600b44`). La importación directa y `graph_from_template_payload` eran fieles; la entrada no usaba el SMILES de respaldo.

Tras la autorización del propietario, se actualizó exclusivamente la entrada ID `tpl_f79a6425241340fab69f51e0b9d65870`: se eliminó el registro del átomo duplicado 11, se sustituyó el enlace 5–11 por el cierre aromático 5–6 y se regeneró el SMILES canónico no estereoespecífico desde el grafo expandido. Los índices de los demás extremos se ajustaron tras quitar el átomo; las coordenadas 2D y atributos de los otros 42 átomos y los 47 enlaces no reemplazados se conservaron. La fixture original se mantiene intacta como evidencia histórica; `tetrandrina_local_template_repaired.mol` es la fixture independiente corregida.

Antes del reemplazo atómico se creó y verificó el respaldo completo `~/.config/Chemuson/template_library.json.bak-20261017T193637027947-0700` (10,844 bytes; SHA-256 `75051c9bc73e92d34bc53f109a123aa909cbf1e47186f8e9f6dd7fafdab5c68d`). El JSON pasó a 10,778 bytes (SHA-256 `0f34a831b22577c1d954788638cf6950dd53cb2960cbd2c533b31446c474eb14`). Solo cambiaron los valores de las propiedades `molblock` y `smiles` del objeto autorizado; el resto de bytes y los datos parseados de toda la biblioteca quedaron iguales. Ni el respaldo ni la configuración local se añaden a Git.

La referencia es tetrandrina neutra, PubChem CID 73078, fórmula `C38H42N2O6`. La comparación topológica descrita abajo ignora la estereoquímica de la referencia; no se asignó estereoquímica a la plantilla local.

## Resultados por etapa

| Etapa | Átomos / enlaces | Composición | Componentes / ciclos | Observación |
|---|---:|---|---|---|
| MolGraph local antes de reparar | 43 / 48 | `C:33`, `CH3:2`, `N:2`, `O:2`, `OCH3:4` | 1 componente; rango de ciclos 6 | Sin centros/enlaces estereoquímicos especificados; incluía el defecto 5–11. |
| Copia de cálculo previa tras expansión | 47 / 52 | `C:39`, `N:2`, `O:6` | 1 componente; rango de ciclos 6 | Expansión correcta; aún contenía el carbono terminal duplicado. |
| Escena histórica tras normalización anterior | 47 pesados + 46 H = 93 / 98 | `C39H46N2O6` | 1 componente; rango de ciclos 6; RDKit informa 7 anillos | Generaba coordenadas al demarcar el puente; identidad incorrecta y 18 átomos/enlaces aromáticos. |
| MolGraph local reparado | 42 / 48 | `C:32`, `CH3:2`, `N:2`, `O:2`, `OCH3:4` | 1 componente; rango de ciclos 7 | Átomo 11 eliminado; cierre aromático 5–6 presente; restantes coordenadas conservadas. |
| Copia de cálculo reparada tras expansión | 46 / 52 | `C:38`, `N:2`, `O:6` | 1 componente; rango de ciclos 7 | Fórmula `C38H42N2O6`, 24 átomos y 24 enlaces aromáticos; conectividad coincide con CID 73078 sin comparar estereoquímica. |
| Escena derivada del conformero reparado | 88 / 94 | `C38H42N2O6` | 1 componente | 46 átomos pesados + 42 H; la malla CPU finita no demuestra presentación OpenGL. |
| Referencia PubChem 73078 | 46 pesados + 42 H = 88 / 52 | `C38H42N2O6` | 1 componente; rango de ciclos 7; RDKit informa 8 anillos | 24 átomos y 24 enlaces aromáticos; comparación estereoagnóstica. |

La auditoría pre-reparación localizó un defecto de **autoría en el MolGraph almacenado**, no de importación ni de expansión:

- El átomo documental `11` es un carbono terminal (grado 1) unido al átomo `5` por el enlace `9`, marcado aromático y con orden sencillo.
- Los átomos `6` y `11` tienen las mismas coordenadas 2D. El átomo `6` ya está enlazado a `1` y `7`; no existe el enlace `5–6` que cerraría el anillo esperado.
- La expansión convierte correctamente los cuatro `OCH3` en cuatro grupos metoxi y los dos `CH3` en carbonos. No añade el carbono discrepante: el átomo `11` ya está en el MolGraph original.

## Qué hacía —y qué no hacía— la normalización anterior

La implementación anterior desmarcaba aromaticidad en el puente terminal `5–11` de una copia de cálculo, lo que permitía a RDKit sanear el grafo. El carbono terminal `11` pasaba a interpretarse como metilo con tres H implícitos. El resultado era `C39H46N2O6`, una escena de 93 átomos y una malla finita. Las valencias típicas quedaban satisfechas (C 4, N 3, O 2), por lo que la comprobación de valencia no descubría el defecto. Esa sanitización **enmascaraba**, no corregía, la conectividad incorrecta: el enlace `5–11` dejaba de ser aromático, pero no se creaba el cierre `5–6` ni se eliminaba el duplicado. La normalización silenciosa ya fue retirada.

Antes de la autorización, se ensayó la reparación solo sobre un objeto RDKit descartable: eliminar el duplicado `11` y sustituir `5–11` por el enlace aromático `5–6` daba `C38H42N2O6`, 46 átomos pesados, 52 enlaces, una componente, 24 átomos/enlaces aromáticos e isomorfismo con la referencia al ignorar estereoquímica. Esa prueba **no modificó entonces** el documento fuente. Después de la autorización, la transformación se aplicó únicamente a la entrada local indicada y se validó de nuevo con el InChI y el SMILES regenerado; no se añadió estereoquímica.

## Diagnóstico y límite de aceptación

- Defecto original: el MolGraph previo tenía el átomo terminal duplicado `11` y el enlace aromático mal asignado (`5–11` en vez de `5–6`); la entrada local fue reparada dentro del alcance autorizado.
- Importación y expansión: la fixture histórica reproduce el estado original; la fixture corregida y la configuración local ahora expanden a 46 átomos pesados y 52 enlaces, sin mutar el documento ni inferir estereoquímica.
- Aromaticidad: el enlace 5–6 cierra el anillo. El backend ya no elimina flags aromáticos para forzar sanitización; detecta y rechaza puentes aromáticos antes del worker con `error_code=invalid_aromatic_bridge`, `stage=calculation_graph_validation` y los IDs afectados en las rutas de conformero, optimización y streaming.
- Identidad: fórmula calculada `C38H42N2O6`; la conectividad coincide con PubChem CID 73078 ignorando estereoquímica. El SMILES canónico almacenado vuelve a producir el mismo InChI.
- Geometría: con semilla 173 y MMFF94 se obtuvieron 46 coordenadas finitas; la minimización no convergió (`method=unoptimized`, `energy=None`). La escena CPU tiene 88 átomos y 94 enlaces y una malla finita; esto no demuestra píxeles ni aceptación OpenGL.
- Estado de aceptación: no se asignó estereoquímica. La verificación gráfica en una sesión real de CachyOS/OpenGL permanece pendiente y se documenta en `manual-test.md`.
