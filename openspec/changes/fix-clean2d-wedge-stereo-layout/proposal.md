# Propuesta — Clean2D con enlaces en cuña sólida y rayada

## Por qué

Clean2D no redibuja una molécula que contenga una sola cuña. En `release/v0.3.0-beta.1-prep` (`9bc8192`), una cadena de cuatro átomos con un enlace `WEDGE` (o `HASHED`) produce en los tres accesos —limpieza rápida, limpieza completa y publicación— el estado `failed-controlled` con mensaje `Limpieza 2D omitida: geometria_actual_necesita_reconstruccion`: cero átomos movidos y ninguna propuesta aplicable. Con dos o más centros quirales el motor además marca `multiple_stereo_centers` y queda en `preserve_only`.

La causa no es que el motor rechace deliberadamente las cuñas ni que no genere candidatos: genera candidatos RDKit e internos correctos. El fallo está en la comprobación de preservación estereoquímica `stereo_layout_signature()` (`src/chemuson/clean2d/local_graph_cleaner.py`), que compara los **signos** del producto vectorial entre pares de vecinos elegidos del centro estereogénico. Ese signo cambia cuando el hueco angular entre dos vecinos cruza 180°, aunque el orden cíclico —y por tanto la configuración dibujada— sea idéntico. Con tres vecinos el hueco medio es de 120°, así que prácticamente cualquier reconstrucción geométrica invierte algún signo y el motor descarta **todos** los candidatos con `cambio_firma_estereoquimica`.

El mismo cálculo sólo mira vecinos **dentro** de la selección. En una selección parcial que deja un vecino fuera, el orden cíclico respecto de ese vecino de frontera puede invertirse sin que la comprobación lo detecte.

## Alcance

- Reemplazar la firma estereoquímica de layout por una basada en el **orden cíclico canónico de los vecinos estructurales reales** de cada centro estereogénico, incluyendo los vecinos que quedan fuera de la selección usando las coordenadas del grafo.
- Conservar en la firma el estilo (`WEDGE`/`HASHED`) y la estereoquímica explícita de los enlaces, de modo que cambiar el sentido de una cuña siga siendo una diferencia.
- Añadir la matriz de reproducción pedida como tests permanentes: control acíclico, una cuña sólida, una cuña rayada, centro tetraédrico con ambas, varios centros quirales, selección con centro quiral y selección con enlace estereoquímico en la frontera.
- Fijar como contrato que un re-encuadre que conserva el orden cíclico no es riesgo estereoquímico, y que una inversión real (imagen especular, intercambio de vecinos, cambio de cuña, cruce de un vecino de frontera) sí lo es.

## Fuera de alcance

No se modifican ChemName, ChemIO, los modelos químicos (`core`), la GUI general, los menús, la persistencia `.cmsn`, CompChem, packaging ni workflows. No se cambian las heurísticas de dibujo de Clean2D, la política de complejidad (`multiple_stereo_centers` → `preserve_only`), el vocabulario estable de motivos de rechazo ni el ranking de candidatos. No se promite que toda molécula con cuñas se redibuje: una selección parcial con estereoquímica en la frontera puede seguir terminando en `preserve-only`/no-op.

## Módulos probables

`src/chemuson/clean2d/local_graph_cleaner.py` (único módulo de producción tocado) y `tests/test_clean2d_wedge_stereo.py`.
