# Tareas — continuidad 2D de amidas y ureas

## A. Alcance y baseline
- [x] A.1 Confirmar `fix/initial-2d-geometry-beta1`, HEAD anterior, ausencia de cambios tracked/staged y preservación del prompt sin seguimiento.
- [x] A.2 Leer OpenSpec activo, `canvas_text.py`, geometría de enlaces, `AtomItem`/`BondItem` y contratos de Visual Lab/uniones.
- [x] A.3 Ejecutar baseline obligatorio: status, compileall, colección, pytest acotado y Ruff global; guardar la salida y el fallo basal de visor 3D.
- [x] A.4 Construir `Ar–NH–C(=O)–NH–CH3` con eventos Qt reales; registrar grafo, longitudes, ángulos, glifos, trim y segmentos visibles.
- [x] A.5 Confirmar por intersección geométrica que proyección vs forma real causa el recorte; validar el contrato OpenSpec estricto antes de producción.

## B. Entorno químico y colocación
- [x] B.1 Añadir clasificación pura basada en conectividad/órdenes: amida, urea, aminas, N aromático y coordinación/otras valencias; sin parsear labels.
- [x] B.2 Aplicar geometría trigonal local ~120° a N amida/urea cuando corresponda, conservando átomos previos, longitud 40 y controles de colisión existentes.
- [x] B.3 Mantener las reglas de amina/N aromático/centros no relacionados; el drag explícito conserva orientación y controles existentes.

## C. Trimming y anclaje de glifos
- [x] C.1 Añadir a `_update_bond_label_shrinks()` la intersección con el path real para NH de amida/urea; conservar el helper de proyección para los estilos/casos no afectados.
- [x] C.2 Usar el margen óptico vigente para el path; conservar estilos de Visual Lab, sigma/pi, wedges, hashes, junctions y el anclaje óptico de `AtomItem`.
- [x] C.3 Probar enlaces N/NH/NH2 en gestos/etiquetas relevantes: no solape, hueco acotado y centro-a-centro intacto.

## D. Regresiones y evidencia
- [x] D.1 QTests del urea asimétrico de referencia `Ar–NH–C(=O)–NH–CH3` y acetanilida (amida secundaria).
- [x] D.2 Controles de amina primaria/secundaria/terciaria, N aromático, coordinación/carga/valencia alta, benceno sustituido y labels N, NH, NH2.
- [x] D.3 Cubrir preview/commit, drag explícito, geometría sin mover centros, undo/redo, CMSN round-trip, zoom, DPR y claro/oscuro.
- [x] D.4 Renderizar before/after comparable con el renderer actual; guardar métricas por enlace, snapshots y hashes.

## E. Validación y entrega
- [x] E.1 Ejecutar focales de dibujo/labels/Visual Lab/proyecciones/persistencia/arquitectura, compileall, Ruff scoped, OpenSpec estricto y `git diff --check`.
- [x] E.2 Registrar fallos basales aparte; no declarar completa la monolítica si mantiene timeout/fallos preexistentes.
- [x] E.3 Crear commits descriptivos en la rama actual; worktree limpio excepto el prompt previo preservado sin seguimiento; no push/merge/PR/tag.
