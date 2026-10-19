# Diseño — geometría local y recorte de etiquetas N–H

## Rama y baseline

- Rama solicitada/verificada: `fix/initial-2d-geometry-beta1`, HEAD inicial `ef241d4a06434fe719a1d36752138eb2806cd53c` (tres commits locales conservados sobre `3e0ef606`). No hay cambios tracked/staged; permanece únicamente el prompt previo sin seguimiento y fuera de alcance.
- Baseline completa y logs exactos: `baseline.md` y `baseline/`. Colección: 2,578 IDs; compileall pasa; `pytest -q` excedió 600 s. `pytest -x` repite el fallo previo del menú de visor 3D tras 366 pasadas. Ruff conserva el F401 histórico `math` en el test Clean2D.
- Familias relevantes previas: **90 passed** (hidrógenos/labels, geometría sp3, dibujo, cuñas/dobles, persistencia y Visual Lab).

## Reproducción real del fragmento referencia

`evidence/before/reproduce_urea_baseline.py` genera un benceno con el gesto de anillo Qt, añade `Ar–NH–C(=O)–NH–CH3` mediante clics cortos y arrastres reales, toma preview antes del release y captura label paths/trazos. El archivo JSON y las imágenes anteriores conservan el estado sin parche.

| Enlace (referencia 40) | Centro-a-centro | Ángulo local relevante | Trim actual | Trazo principal visible |
|---|---:|---:|---:|---:|
| Ar–N1H | 40.000000 | N1: 120° | 5.631875 en N | 34.368125 |
| N1–C(=O) | 39.999999999999986 | N1: 120° | 15.620745 en N | 24.379255 |
| C=O | 39.999999999999964 | C: 120.514° con drag (120° en clic-only) | 8.204943 en O | 31.795057 por cada línea sigma/pi |
| C(=O)–N2H | 39.999999999999986 | C: 119.486° con drag; N2: 119.596° | 15.743313 en N | 24.256687 |
| N2–CH3 | 39.999999999999986 | N2: 119.596° | 5.679473 en N | 34.320527 |

Con sólo clics cortos, los ángulos locales de ambos N y el carbonilo son 120°; las coordenadas no explican el segmento corto. El drag produce hasta `120.514°`, acorde al redondeo de eventos Qt. Los extremos de etiqueta son el único trim; `_endpoint_trim_start/end` valen cero en estos enlaces. Los dos segmentos NH–carbonilo pierden alrededor de 8 unidades adicionales cada uno frente a lo que la forma real del texto requiere.

## Glifos y causa confirmada

- N1 centro `(479.0, 499.0)`. El texto es `NH` a 11 pt: glyph N bbox `(476.9025, 497.4637, 8.375, 10.3125)`; H bbox `(487.7306, 497.4637, 8.375, 10.3125)`. El N se ancla ópticamente al átomo más el offset local `(2.09, 3.62)`; H queda a su derecha.
- N2 centro `(479.0, 429.7180)`. N bbox `(476.8897, 420.9344, 8.375, 10.3125)`; H bbox `(487.7179, 420.9344, 8.375, 10.3125)`. Offset óptico aproximado `(2.08, -3.63)`.
- En N1→carbonilo, la proyección actual devuelve `13.320745` de extent + `2.30` de clearance = `15.620745` de recorte. Al dilatar el path real `QGraphicsTextItem.shape()` por el clearance de `2.30` y cruzarlo con la línea central, el extremo requerido es `8.203442`; el recorte vigente excede `7.417304` unidades.
- En carbonilo→N2, el extent proyectado es `13.443313` + `2.30` = `15.743313`; la dilatación/intersección requiere `8.375382`. Exceso: `7.367931` unidades. Para Ar–N1H y N2H–CH3 el nuevo extremo coincide con el recorte actual con diferencias `<0.04`, controlando que la corrección no acorte enlaces por sistema.

La proyección de puntos de glyph es una cota de soporte sobre el eje, no la longitud de la línea necesaria para despejar su forma. El NH está desplazado ópticamente hacia el sector libre; puntos laterales/lejanos inflan el recorte aun cuando la trayectoria centralizada, una vez extendida por stroke y clearance, deja de tocar la forma. La dilatación explícita del path seguida de intersección con la línea central da el extremo geométrico requerido. La reproducción confirma un defecto de trimming visual; la geometría del grafo permanece en 40 unidades.

## Decisiones de diseño por validar

1. La geometría automática consulta topología/órdenes reales: carbonilo `C=O` unido por enlace simple al N indica amida; dos N unidos al mismo carbonilo indican urea. Enlace/aromaticidad del propio N y coordinaciones se evalúan primero/por separado. No inspeccionar `item.label`, `NH` ni `NH2` para clasificar química.
2. Para N amida/urea y el carbonilo, elegir candidatos locales trigonal-planos ~120° que mantengan las coordenadas previas, las colisiones evitables y la longitud de referencia (40). No reubicar la estructura para mejorar texto. Amina primaria/secundaria/terciaria, N aromático y coordinado conservan sus reglas salvo caso medido.
3. Corregir sólo el enlace simple, plano y no aromático que termina en `NH`/`NH2` de una amida/urea neutra: dilatar la forma visible label/charge por el clearance vigente e intersectarla con la línea central. Usar ese extremo geométrico en vez del soporte proyectado remoto. Dobles/triples, enlaces aromáticos, otros labels/N, wedges, hashes y estilos Bold conservan el trimming previo; así se evita alterar junctions y renderer editorial.
4. `AtomItem` conserva el anclaje óptico del token del elemento; sólo tocarlo si mediciones de glyph N/H prueban un desplazamiento distinto de su offset local aprobado.
5. Preview y commit consultan el mismo candidato. El arrastre explícito mantiene la dirección del cursor; el trim no altera `MolGraph`.

## Riesgos y pruebas

- El cambio queda fuera de estilos estereo/multienlace; una regresión de coverage asegura que el gating continúe acotado.
- Márgenes cero/inadecuados pueden solapar el glyph; márgenes grandes vuelven a acortar el trazo. Medir gap y recorte por separado.
- Contexto incompleto antes de que se dibuje el carbonilo puede clasificar un N de forma prematura; probar el urea/amida por eventos y verificar actualización al completar C=O.
- Matriz: `Ar–NH–C(=O)–NH–CH3`, acetanilida, aminas primarias/secundarias/terciarias, N aromáticos, coordinación/carga/valencia alta y labels NH/NH2 en direcciones múltiples; preview/commit, undo/redo, CMSN, zoom/DPR y temas. Mantener controles de anillo aromático.
