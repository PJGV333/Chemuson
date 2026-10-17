# Aceptación manual ChemUSON-3D en CachyOS

**Estado:** pendiente de ejecución por el usuario. No se inició ninguna GUI durante esta campaña. El cierre puede terminar en SIGSEGV según las observaciones previas; guarda el `.cmsn` antes de cada cierre.

## Preparación

1. En una terminal normal (sin `sudo`), desde la raíz del repositorio y con `.venv` existente, iniciar con:

   ```bash
   PYTHONPATH=src ./.venv/bin/python -m chemuson
   ```

2. Anotar fecha, `git rev-parse --short HEAD`, resolución y tema claro/oscuro. Usar un `.cmsn` temporal y no documentos de trabajo.
3. Si el proceso requiere aprobación para abrir el display, detenerse y solicitarla; no cambiar permisos, configuración del sistema ni buscar otro método de lanzamiento.

## Recorrido A — Canvas 2D → conformero 3D → persistencia

1. Crear un documento vacío y dibujar una cadena de butano (`C–C–C–C`) en el canvas 2D.
2. Abrir CompChem y ejecutar **Generar conformero 3D**. En **Ver → Visor molecular 3D**, verificar que aparecen los cuatro carbonos y los H visuales derivados; anotar el estado/número de átomos y enlaces.
3. Verificar cámara: arrastre izquierdo para orbitar, derecho para pan, rueda para zoom; la molécula debe permanecer encuadrada y sin clipping. Seleccionar átomos y comprobar que el resaltado se distingue del fondo oscuro.
4. En el dock, activar **Editar geometría 3D** y abrir **Controles de edición**. Cambiar una distancia o el dihedro de cuatro átomos consecutivos usando sus IDs reales en el orden del enlace (no inferirlos por la posición en pantalla). Comprobar visualmente que los H del fragmento acompañan el cambio. Ejecutar undo y redo; cancelar un drag con Escape y comprobar que vuelve exactamente a la geometría anterior.
5. Activar **Live MM (RDKit)**. Verificar que CompChem recibe varios pasos/energías y que la geometría final se confirma; mover de nuevo un átomo y comprobar que el último trabajo sustituye al anterior, sin aplicar resultados obsoletos.
6. Guardar como archivo `.cmsn` temporal, cerrar el documento, reabrirlo y confirmar: canvas 2D intacto, topología/document IDs estables, geometría 3D presente, H todavía derivados (no documentales) y undo disponible para la edición final.

## Recorrido B — Construcción desde el visor 3D

1. Crear otro documento vacío, abrir **Ver → Visor molecular 3D**, activar **Editar geometría 3D**, y en **Controles de edición** elegir **Construir molécula**, elemento C y enlace **Sencillo**.
2. Hacer clic en espacio vacío para crear el primer átomo. Hacer clic en ese átomo para seleccionarlo como ancla; hacer clic en espacio vacío para crear el siguiente. Repetir ancla→espacio para formar una cadena de cuatro carbonos. Verificar que cada enlace nuevo queda a distancia de enlace razonable y que el estado del dock se actualiza.
3. Añadir dos enlaces adicionales al mismo carbono central hasta completar cuatro enlaces sencillos. Intentar añadir un quinto enlace: la operación debe rechazarse por valencia sin añadir átomo/enlace ni crear una entrada de undo parcial.
4. Activar Live MM y verificar actualización progresiva. Probar undo/redo sobre una construcción válida.
5. Guardar/reabrir el `.cmsn`; comprobar que conectividad y coordenadas corresponden al documento construido.
6. Crear una tercera pestaña/documento, generar otra estructura y alternar entre pestañas. Confirmar que escena, geometría, trabajo Live MM y selección no se filtran a los otros documentos.

## Matriz química y visual

Repetir generación/inspección de conformero para agua (`O`), etano (`CC`), etanol (`CCO`), butano (`CCCC`) y benceno (`c1ccccc1`) mediante las rutas de entrada existentes. Confirmar que los H derivados aparecen en el visor, pero no como átomos documentales. En benceno comprobar aromaticidad y visibilidad del anillo; en caso de dibujar Kekulé, comprobar que no se convierte arbitrariamente ciclohexano en aromático.

Para estéreo explícito, usar la importación existente con estos pares, si está disponible en la interfaz:

- Tetraédrico: `F[C@](Cl)(Br)I` y `F[C@@](Cl)(Br)I`.
- Alqueno: `F/C=C/Cl` y `F/C=C\Cl`.

Comprobar que los dos miembros de cada par no se convierten en la misma configuración y que la geometría optimizada no invierte silenciosamente R/S o E/Z. Para marcas ausentes/ambiguas, verificar que la aplicación informa que no puede garantizar la configuración, en lugar de inventarla.

## Usabilidad y estabilidad

- Redimensionar el dock y el panel CompChem: los controles principales de Live MM deben seguir accesibles, sin tapar el viewport ni deteriorar el canvas 2D.
- Probar tema oscuro y claro: los carbonos se distinguen del fondo; enlaces simples/dobles/triples/aromáticos, selección y profundidad son legibles.
- Registrar el cierre por separado: (a) sin abrir el visor, (b) visor abierto sin trabajo, (c) tras terminar Live MM y (d) en un intento independiente con trabajo activo. Para cada caso anotar código de salida y si aparece SIGSEGV. No combinar estos resultados con el SIGSEGV de pytest.
- Si ocurre el SIGSEGV, conservar el registro/backtrace y el core disponible; no usar `os._exit`, demoras arbitrarias ni limpieza forzada. La guía del caso histórico está en `validation.md`.

## Resultado a devolver

Por recorrido y molécula, anotar `PASS`/`FAIL`/`NOT TESTED`, síntoma, archivo temporal utilizado y código de salida. Adjuntar capturas del visor sólo si se autoriza y no contienen documentos privados. La aprobación visual y el cierre limpio sólo se consideran verificados después de esta prueba real; los tests offscreen no la sustituyen.
