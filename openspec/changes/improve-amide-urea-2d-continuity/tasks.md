# Tareas — continuidad 2D de amidas y ureas

## A. Alcance y baseline
- [x] A.1 Confirmar `fix/initial-2d-geometry-beta1`, HEAD anterior, ausencia de cambios tracked/staged y preservación del prompt sin seguimiento.
- [x] A.2 Leer OpenSpec activo, `canvas_text.py`, geometría de enlaces, `AtomItem`/`BondItem` y contratos de Visual Lab/uniones.
- [x] A.3 Ejecutar baseline obligatorio: status, compileall, colección, pytest acotado y Ruff global; guardar la salida y el fallo basal de visor 3D.
- [x] A.4 Construir `Ar–NH–C(=O)–NH–CH3` con eventos Qt reales; registrar grafo, longitudes, ángulos, glifos, trim y segmentos visibles.
- [x] A.5 Confirmar por intersección geométrica que proyección vs forma real causa el recorte; validar el contrato OpenSpec estricto antes de producción.

## B. Entorno químico y colocación
- [ ] B.1 Añadir clasificación pura basada en conectividad/órdenes: amida, urea, aminas, N aromático y coordinación/otras valencias; sin parsear labels.
- [ ] B.2 Aplicar candidatos locales ~120° a N amida/urea cuando factible, conservando átomos previos, longitud 40, colisiones y continuidad.
- [ ] B.3 Mantener las reglas de amina/N aromático/centros no relacionados; el drag explícito conserva orientación y controles existentes.

## C. Trimming y anclaje de glifos
- [ ] C.1 Corregir `_label_projection_extent_for_atom()`/`_update_bond_label_shrinks()` para medir intersección real de la envolvente del trazo con el path visible de label/charge.
- [ ] C.2 Elegir margen óptico controlado por grosor/perfil, conservar estilos de Visual Lab, sigma/pi, wedges, hashes, junctions y anclaje óptico de `AtomItem`.
- [ ] C.3 Probar enlace lineal/oblicuo a N/NH/NH2, varios órdenes y direcciones: sin solape, hueco acotado, centro-a-centro intacto.

## D. Regresiones y evidencia
- [ ] D.1 QTests de urea, N,N′-dimetilurea, N-metil-N′-(4-nitrofenil)urea, acetanilida y amida secundaria sencilla.
- [ ] D.2 Controles de amina primaria/secundaria/terciaria, N aromático, coordinación, benceno sustituido y sustituyentes/labels N, NH, NH2 en orientaciones múltiples.
- [ ] D.3 Cubrir preview/commit, drag explícito, colisiones, undo/redo, CMSN round-trip, zoom, DPR, claro/oscuro.
- [ ] D.4 Renderizar before/after comparable con el renderer actual; guardar métricas por enlace, snapshots y hashes.

## E. Validación y entrega
- [ ] E.1 Ejecutar focales de dibujo/labels/Visual Lab/proyecciones/persistencia/arquitectura, compileall, Ruff scoped, OpenSpec estricto y `git diff --check`.
- [ ] E.2 Registrar fallos basales aparte; no declarar completa la monolítica si mantiene timeout/fallos preexistentes.
- [ ] E.3 Crear commits descriptivos en la rama actual; worktree limpio excepto el prompt previo preservado sin seguimiento; no push/merge/PR/tag.
