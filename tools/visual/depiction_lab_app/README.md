# ChemUSON Molecular Depiction Visual Lab

Desarrollo/calibración visual aislada: los dos paneles ejecutan `ChemusonCanvas`, `AtomItem` y `BondItem` de producción. Current es de solo lectura; Working cambia exclusivamente el estilo de su canvas. No se editan documentos, preferencias normales, geometría química ni Undo/Redo. Los perfiles no se aplican al editor principal.

## Inicio

Desde la raíz del repositorio y el entorno virtual de desarrollo:

```bash
.venv/bin/python tools/visual/depiction_lab.py
```

La galería A1–C4 incluye 12 IDs deterministas; B4 está marcado explícitamente como fixture geométrico no molecular. Selecciona un caso, tema, zoom y rotación visual en la barra superior. Los controles se agrupan en trazos/multienlaces, cuñas y etiquetas; cada valor actualiza los `BondItem`/`AtomItem` Working sin reconstruir el canvas. El zoom y los scrollbars se sincronizan. “Diagnóstico” añade centros, endpoints, paths/contornos, bounding boxes y marcas objetivas de intersección; no calcula calidad estética. “Ocultar controles” maximiza los paneles.

Los 15 parámetros expuestos son: grosor de trazo; offset y trim interno de dobles; offset bold-double y triple opcionales; anchura máxima y punta opcional de cuña sólida; cantidad, anchos mínimo/máximo, grosor y desplazamiento terminal opcional de hashes; clearance común opcional de etiqueta; terminación y join Qt. Los valores geométricos son dimensionless ratios: offsets/anchuras principales usan la longitud de enlace de referencia (40 scene units); el ancho de punta se refiere al máximo de cuña, el shift hash a su intervalo y el clearance a la longitud del enlace local. El conteo y las opciones son discretos. `null` en los cinco overrides opcionales significa “Auto: cálculo actual”, no cero.

No se exponen factores internos de miter/fork/bifurcación, ancho estrecho contextual de cuña, óptica tipográfica por elemento/dirección, trims de anillos ni orientación sigma/pi: no son controles aislables seguros en esta iteración. La UI documenta estas limitaciones. Los valores Current igualan exactamente `CHEMDOODLE_LIKE`; toda extensión de `DrawingStyle` tiene defaults `None` que conservan la implementación previa.

## Perfiles

- “Restablecer todo” vuelve a Current; cada parámetro tiene su propio botón de reset.
- Escribe un nombre y usa “Guardar perfil” para guardar JSON en la carpeta de configuración del laboratorio: `QStandardPaths.AppConfigLocation/depiction-lab/profiles/` (la aplicación se identifica como `ChemUSON/DepictionVisualLab`).
- “Abrir JSON…” valida schema, versión, tipos, rangos y claves desconocidas. La migración v0→v1 está explícita. Abrir/guardar usa un diálogo sólo en esa operación; calibrar controles no abre diálogos.
- Los perfiles y QSettings normales de ChemUSON permanecen separados.

## Exportación reproducible

Desde la ventana exporta el caso actual o toda la galería; sin interfaz, desde la raíz:

```bash
.venv/bin/python tools/visual/depiction_lab.py --export-gallery /tmp/chemuson-depiction-lab --theme light --profile /ruta/candidate.json
.venv/bin/python tools/visual/depiction_lab.py --export-case B1 --output /tmp/chemuson-depiction-case --theme dark --profile /ruta/candidate.json
```

La hoja de contacto muestra Current/Working de los 12 casos; el bundle individual incluye PNGs, comparación, perfil JSON y `report.md` con versión, SHA de origen, estado clean/modified del worktree, parámetros, plataforma/Qt y checksums. SVG se emite vía `QSvgGenerator` de la misma escena Qt si el módulo `PyQt6.QtSvg` está disponible. Overlays diagnósticos se excluyen del export regular. Resolución de panel PNG: 1200×820; hoja de contacto conserva escala igual entre ambos renders.

Prueba de arranque acotada:

```bash
.venv/bin/python tools/visual/depiction_lab.py --smoke
```

## Invariantes

La rotación es una transformación de `QGraphicsView`; no cambia coordenadas ni `BondStereo`. Los controles no cambian IDs, conectividad, órdenes, aromaticidad, cargas, isótopos, coordenadas o datos estéreo. Los gráficos se crean desde fixtures; no se ejecuta Clean2D. La aceptación estética y cualquier preset editorial futuro pertenecen al propietario. Para promover ajustes se requiere otra propuesta OpenSpec; el laboratorio nunca cambia defaults de producción.
