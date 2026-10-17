# Evidencia visual reproducible

## Haworth C3 antes/después

- Baseline: `87e298c46535e1dffb0477f36412b3dedc87bd98`, exportada desde worktree limpio.
- Geometría refinada: `0f501cc70388b429d34c7bf7132ad227849a8aad`, exportada desde worktree limpio.
- Mismo caso `C3`, renderer de producción, tema claro, estilo `Current`, anti-aliasing activado y paneles idénticos de 1200×820 px. `haworth/before-after-C3.png` los presenta a escala 1:1 en una lámina 2400×868; no hay resize independiente.
- Bundles originales con reportes, perfiles, PNG y SVG: `haworth/baseline/` y `haworth/refined/`.
- SHA-256 de la lámina: `806eb83f2e18ed7dc6a5d4599e78e358bc27c2deea41821b984b3d117a497a5b`.

C3 está excluido de la cohorte editorial AB3. Se usa sólo como geometría Haworth de referencia para revisar desplazamiento del grupo C3 OH y el espacio con su vecino; no calibra ni promueve un preset.

## Viewport del laboratorio

`viewport/depiction-lab-C3-synced.png` muestra Current y Working en la misma ventana real a 100%/0°, con geometría no vacía y centrada. La probe Qt generó la captura durante un ciclo de 12 comprobaciones: A1/C3, dos tamaños de ventana, tres niveles/ángulos, conteos iguales de átomos/enlaces visibles, error de centro ≤10 px y cierre limpio. SHA-256: `1994ec67608aaef7ac4587dde7b37a49a8f5a834a29da8f52b1c2674d2630706`.

## Haworth, Fischer y silla con semántica tipada

`projections/projection-geometry-gallery.png` (1360×1350) contiene cinco paneles renderizados mediante `ProductionCanvasAdapter` y `ChemusonCanvas` reales: β-D-glucopiranosa, furanosa genérica, Fischer con cruces alternadas y una silla A/B con el mismo OH up pasando de axial a ecuatorial. `projection-geometry-gallery.md` registra el SHA de fuente, geometrías/semánticas y que todos los enlaces gráficos conservan `BondStereo.NONE`. SHA-256 de la imagen: `681a35cba206250b9135d3a1341998efea14a7c10f26baa5d7a3b1deb5a95642`.

Las capturas y pruebas verifican geometría y ausencia de cruces/colisiones bajo las configuraciones explícitas cubiertas; **no son aceptación estética del propietario** ni modifican AB3, defaults de producción o `.cmsn`.
