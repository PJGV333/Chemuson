# Diseño — política única de geometría esqueletal inicial

## Base y auditoría

- Rama aislada: `fix/initial-2d-geometry-beta1`.
- Base verificada tras `git fetch origin`: `origin/release/v0.3.0-beta.1-prep` en `3e0ef606edb90929c1fb32dddcd5385fabcc5d82`, que contiene el merge `3e0ef606`.
- Traza: `mousePressEvent` resuelve átomo/enlace y comienza la herramienta; `mouseMoveEvent` actualiza preview; `CanvasBondGeometryMixin`/`CanvasToolsRingsChainsMixin` determinan puntos; el release construye los comandos `AddBondCommand`, `AddChainCommand` o `AddRingCommand`; los comandos escriben coordenadas en `MolGraph`; `AtomItem` toma esas coordenadas y `BondItem.update_positions()` dibuja trazos, recortes y offsets sin recolocar los átomos.

### Causas clasificadas por separado

1. **Modelo, cadena acíclica:** `SP3_BOND_ANGLE_DEG=109.5` se usaba como ángulo de dibujo en `candidate_directions_deg`, `_sp3_display_angle_deg` y `zigzag_turn = 180 - sp3_angle`. Se midieron ángulos internos 109.5° en la cadena libre y en la sugerencia sp3 con un vecino. Son coordenadas realmente distintas de la convención editorial 2D de 120°.
2. **Modelo, sustituyente rotado:** con `fixed_angles=True`, `angle_step_deg=30` y benceno girado 17°, el candidato local radial estaba a 17°, pero el snapping absoluto lo movió a 0°; los ángulos con las dos aristas quedaron en 137° y 103°, no en 120°. El cursor no expresaba una dirección en un clic corto, por lo que la rejilla global no debía sustituir la geometría local.
3. **Flujo de clic corto:** en un clic iniciado en un vértice aromático, el release también encontraba una arista adyacente bajo el átomo. Se reprodujo que el gesto dejaba el grafo en 6 átomos/6 enlaces en vez de crear el séptimo átomo.
4. **Benceno/render:** la herramienta genera seis lados `40.0` (error máximo observado < `4e-14`) y ángulos internos `120.0°` (error < `1.2e-13°`). Las capturas 1× y 1.25× del renderer de producción muestran el ciclo cerrado y visualmente equilibrado. Los trazos pi se retranquean/offsetean intencionalmente; esos píxeles no son coordenadas de vértices. No se confirmó defecto de clipping, junction, antialias o estilo.
5. **Caso completo:** la secuencia Ph–CH₂–C(=O)–CH₂–CH₃ dibujada con arrastres dirigidos produce 11 átomos, 11 enlaces, anillo regular y longitudes de 40 unidades; los ángulos locales quedan entre 119.69° y 120.32° debido a cuantización de coordenadas de evento Qt. Esta secuencia no reproduce el defecto de orientación automática; la reproducción numérica separada sí lo hace.

### Política

| Operación | Regla de geometría 2D | Interacción preservada |
|---|---|---|
| Enlace sp3 nuevo por sugerencia/clic corto | Direcciones vecinas a 120° respecto de enlaces existentes; seleccionar la alternativa libre con menor colisión | No interpreta el punto de liberación cercano como una orientación explícita |
| Arrastre dirigido | El cursor expresa orientación; `fixed_angles`, `angle_step_deg`, `fixed_lengths`, `Alt` y `Shift` conservan su flujo actual | No autoajustar globalmente una parte ya dibujada |
| Cadena libre/asistida | Longitudes uniformes y zigzag con giro de trayectoria ±60° (ángulo interno 120°); la primera orientación sigue al cursor/snap existente | `Alt`/`fixed_angles` mantienen el control actual del zigzag |
| Centro sp2/carbonilo | Distribución trigonal de 120° según enlaces existentes y espacio libre | Conserva orden de enlace y conectividad |
| Centro sp | Direcciones lineales de 180° | Sólo cuando la inferencia actual lo clasifica sp |
| Anillo regular n=3/4/5/6 | Polígono regular con ángulos internos 60°/90°/108°/120° y lados iguales | En anclaje a enlace, la arista existente fija el lado y no se mueven sus extremos |
| Benceno exocíclico | La sugerencia corta usa la bisectriz radial exterior del vértice sp2 (120° local) | El arrastre explícito sigue pudiendo elegir otra orientación |

La constante 109.5° permanecerá nombrada como referencia tetraédrica y no será la regla de colocación esqueletal plana. Se centraliza el valor 120° de dibujo 2D en `geom.py` para evitar duplicar heurísticas. El paso de rejilla no se cambia. Para el clic corto automático se omite sólo el snap a ejes globales, ya que el gesto carece de dirección expresada; el arrastre y demás call-sites conservan el snap.

## Implementación mínima

- Añadir una constante explícita para ángulo esqueletal sp3 2D; usarla en generación de candidatos, estado de zigzag y cadena libre.
- Hacer explícito en la selección de candidatos si hay intención angular del cursor. El endpoint automático usa ángulos relativos locales sin rejilla absoluta; el arrastre no cambia.
- En el release, no ejecutar el fallback de hit-test de enlace sobre el mismo átomo que inició un clic corto. El hit-test directo de un clic en enlace permanece intacto.
- No modificar `_compute_ring_vertices` ni `BondItem` salvo que una nueva prueba numérica/pixel demuestre una falla separada. La regularidad del polígono se verifica antes y después.
- Preview y commit usan el mismo helper; los comandos existentes permanecen undoables/deterministas. Sólo se agregan coordenadas nuevas, no se realinean vecinos.

## Riesgos y cobertura

El riesgo principal es cambiar una expectativa histórica de `fixed_angles` en clics cortos. Se acota a sugerencias sin dirección de cursor y se prueba por separado que el arrastre mantiene el snapping. Otro riesgo es dibujar un enlace coincidente con un vecino existente; se cubren colisión/snap y el click sobre enlace. Las pruebas numéricas separan coordenadas del modelo de trazos renderizados, incluyen anillos girados/trasladados/zoom/DPR y la estructura completa, y contrastan undo/redo/persistencia sin tocar `.cmsn`.
