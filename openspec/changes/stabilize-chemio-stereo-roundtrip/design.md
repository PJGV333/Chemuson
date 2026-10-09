# Diseño — preservar la identidad estereoquímica en ChemIO

## Diagnóstico previo a la implementación

Baseline en `baseline.md`: ambos tests P0 fallan por ausencia de `Atom.stereo_cip`, aunque el grafo sí contiene cuña/discontinua. Para `C[C@H](O)F`, el parser MOL fallback ordena el par de átomos del enlace con cuña, cambiando el extremo de referencia sin invertir la cuña. En `N[C@@H](C)C(=O)O` el orden original coincide con el orden normalizado, por eso la cuña sobrevive visualmente. `molgraph_to_smiles` manda al worker átomos/enlaces sin reconstruir `stereo_cip` ni `BondStyle.WEDGE/HASHED`; la salida observada es aquiral (`CC(O)F` / `CC(N)C(=O)O`). El serializador MOL fallback siempre escribe código de estereoquímica 0. La ruta RDKit tampoco transfiere estilo/paridad en la conversión MolGraph→RDKit. El worker no recibe coordenadas, por lo que no puede apoyarse en geometría de un bloque MOL para reconstruir dobles enlaces.

El SMILES histórico de tetrandrina produce cero centros especificados y dos centros potenciales sin asignar con RDKit 2026.03.6; no incluye E/Z. El test que espera cuñas es científicamente incorrecto para esa entrada.

## Decisiones

1. Preservar el orden de extremos de enlaces dirigidos al parsear CTAB; ordenar pares sólo para la clave de deduplicación y mantener el sentido asociado a wedge/hash.
2. Convertir cuña/hash y `Atom.stereo_cip` explícitos en datos estereoquímicos RDKit antes de SMILES/MOL, preservando el orden de vecinos; el worker aislado aplicará el mismo contrato. No asignar tags a centros sin anotación.
3. Serializar códigos wedge/hash/either en MOL fallback. Para E/Z, transferir `Bond.stereo_ez` y sus átomos de referencia al backend; la ruta de importación conservará la información estéreo que el formato/backend realmente proporcione. Ante estereoquímica especificada no representable, fallar explícitamente en lugar de degradar en silencio.
4. Mantener fórmula, conectividad, identidad/orden de átomos y enlaces, cargas, isótopos, órdenes y coordenadas del grafo sin normalizaciones químicas oportunistas. Las representaciones SMILES pueden diferir textualmente; se compara equivalencia química estereoquímica con RDKit.
5. No cambiar geometría Clean2D. Las cuñas son una representación 2D de una asignación química existente; no se agregan a centros no especificados.
6. Dejar sin cambios interfaces públicas y dependencias; las nuevas pruebas usan RDKit ya instalado y workers aislados para reproducción.

## Riesgos

- Las cuñas dependen del átomo inicial del enlace; normalizar sus extremos sin invertir la dirección puede cambiar o borrar el centro.
- `stereo_cip` es una etiqueta derivada, no sustituye a la paridad/química efectiva ni se puede copiar como texto SMILES `@`/`@@` sin considerar el orden de vecinos.
- Los molfiles pueden codificar estereoquímica mediante dirección de enlace, paridad, o geometría según versión/backend. Los casos no fielmente representables deben fallar o advertir; una salida no estereoespecífica no pasa como equivalencia.
- Las pruebas se ejecutarán con `timeout` externo; no se ejecuta el pytest monolítico por la restricción P0 de duración/Qt.

## Validación

Ver `validation.md`. La aceptación exige enantiómeros opuestos no equivalentes, equivalencia de representación alternativa, centros no especificados que siguen sin asignar, controles moleculares y MOL/SDF focalizados. No se ejecutará la suite completa de Clean2D ni se afirmará validación integral de otros subsistemas; sólo se incluyen regresiones focalizadas de consumidores ChemIO para comprobar que no cambió su contrato.
