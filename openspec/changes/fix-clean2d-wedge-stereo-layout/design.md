# Diseño — Clean2D con enlaces en cuña sólida y rayada

## 1. Qué comprueba realmente la protección estereoquímica

`run_clean2d_engine` descarta un candidato si `stereo_layout_signature(graph, before, selected)` difiere de `stereo_layout_signature(graph, after, selected)` (`src/chemuson/clean2d/engine.py`, motivo `cambio_firma_estereoquimica` → `stereo-risk`). La firma debe responder a una sola pregunta: **¿cambió la configuración estereoquímica que el dibujo comunica?**

En un dibujo 2D esa configuración la determina el **orden cíclico** de los vecinos alrededor de cada centro estereogénico junto con el estilo de cada enlace (cuña sólida, rayada o normal). No la determinan los ángulos concretos ni el orden en que el modelo almacena los vecinos.

## 2. Por qué fallaba el control anterior

La implementación anterior guardaba, para cada centro, los signos del producto vectorial normalizado de cada **par** de vecinos seleccionados:

```
signed = (lx·ry − ly·rx) / (llen·rlen)   →   sign ∈ {−1, 0, +1}
```

Ese signo cambia cuando el hueco angular entre dos vecinos cruza 180°, sin que el orden cíclico cambie. Caso real medido sobre el caso `single-wedge` (centro 2; vecinos 1, 3 y 4; cuña 2–4):

| disposición | ángulos de vecinos respecto al centro | orden cíclico | signos por par (1,3) (1,4) (3,4) |
|---|---|---|---|
| dibujo original | 1 = −173.4°, 3 = −56.3°, 4 = −10.0° | (1, 3, 4) | +1, +1, +1 |
| candidato `rdkit_isolated` | 1 = 154.8°, 3 = −85.2°, 4 = 34.7° | (3, 4, 1) ≡ (1, 3, 4) | +1, **−1**, +1 |

El candidato deja el orden cíclico intacto, pero el par (1, 4) pasa de un hueco antihorario de 163.4° a 239.9° y su signo se invierte. La firma cambia, el motor etiqueta `cambio_firma_estereoquimica` y descarta el candidato. Como los re-encuadres razonables siempre mueven huecos alrededor de los 180°, **todos** los candidatos geométricos se descartaban y la molécula con una sola cuña se quedaba sin candidato aplicable. El candidato `internal_templates` además producía signos `0` (vecinos casi colineales respecto al centro), otro artefacto del mismo criterio.

Segundo defecto: la firma sólo consideraba vecinos **dentro** de la selección (`neigh in coords` y adyacencia construida sólo con enlaces seleccionados). Un centro parcialmente seleccionado podía reordenarse alrededor de un vecino de frontera sin que la firma cambiara: riesgo real de inversión silenciosa en limpiezas parciales.

## 3. Decisión

Reescribir `stereo_layout_signature()` sobre el orden cíclico canónico:

1. **Enlaces considerados**: los estructurales (`bond_affects_valence`) con **al menos un extremo** en la selección, es decir, incluidos los enlaces de frontera (`_structural_bonds_touching_selection`).
2. **Centro**: se mantiene `_stereo_centers(...) & selected`, de modo que sólo se firman centros dentro de la selección.
3. **Orden cíclico** (`_cyclic_neighbour_order`): para cada centro se ordenan antihorariamente **todos sus vecinos estructurales reales** (dentro o fuera de la selección) por `atan2`, se desempata por ID de átomo y se canoniciza por rotación mínima. Cada elemento del ciclo lleva el estilo de enlace normalizado, así que una cuña que se convierte en rayada cambia la firma aunque la geometría no se mueva.
4. **Coordenadas de vecinos fuera de la selección** (`_layout_point`): si la propuesta no las incluye, se usan las del grafo, que es donde realmente están durante una limpieza parcial.
5. **Enlaces con estereoquímica explícita**: se conserva el bloque existente (ID, extremos, estilo y `stereo`), ahora también para enlaces de frontera.

Consecuencia exacta: una rotación rígida, un cambio de escala o un re-encuadre que conserva el orden cíclico producen la **misma** firma; una imagen especular, un intercambio de dos sustituyentes, un cambio de sentido de cuña o el cruce del rayo de un vecino de frontera producen una firma **distinta**.

## 4. Alternativas descartadas

- **Relajar o desactivar el control para estructuras con cuñas.** Rechazado: ocultaría inversiones reales y violaría el requisito existente `Stereo risk rejected`.
- **Comparar descriptores CIP (R/S) de RDKit por candidato.** Rechazado: exigiría una ida y vuelta completa a RDKit por cada candidato (coste y aislamiento de worker ya acotados en otras campañas), y el CIP depende de valencias e hidrógenos implícitos, no de lo que el dibujo comunica; además no distingue cuña sólida de rayada cuando el centro no es especificado.
- **Suma de áreas con signo del abanico de vecinos.** Equivalente en poder de detección pero sin canonicización por rotación: quedaría sensible al orden de almacenamiento de los vecinos. Se prefiere el mínimo canónico por rotación, que es determinista y no depende del orden interno del grafo.

## 5. Invariantes que se preservan

Identidad de átomos y enlaces, conectividad, órdenes, aromaticidad, cargas formales, isótopos, metadatos de estereoquímica de átomo y enlace, sentido y procedencia de cada cuña, estereoquímica E/Z explícita, fronteras de selección (no se mueven átomos fuera de la selección), determinismo (misma entrada, misma firma, mismo `geometry_hash`) y Undo/Redo (la aplicación de coordenadas sigue pasando por el mismo comando de undo).

## 6. Límites conocidos y admitidos

- Con dos o más centros quirales, `complex_policy` sigue devolviendo `multiple_stereo_centers` y el resultado es `preserve_only`: no se redibuja, y no se toca esa política en esta campaña.
- Una selección parcial que deja estereoquímica en la frontera puede seguir terminando sin candidato aplicable (`failed-controlled` / no-op). El contrato que se fija aquí es que **nunca** se aplica una inversión, no que siempre exista un re-encuadre seguro.
- `propose` sigue exigiendo una base con calidad `good`; sobre un dibujo sucio informa `La estructura debe optimizarse antes de proponer conformeros 2D`. Tras una limpieza correcta, `propose` sí propone sobre la molécula con cuña.
- No se añade dependencia nueva ni se modifica el vocabulario estable de motivos de rechazo.
