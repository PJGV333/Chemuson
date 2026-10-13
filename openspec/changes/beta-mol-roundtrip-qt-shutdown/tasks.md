# Tareas — microcampaña P0

## 1. Alcance y baseline
- [x] Verificar que la rama base continúa en `c75026d5656972e0ee3f002f615d80060113bc68` y crear la rama independiente solicitada.
- [x] Leer OpenSpec beta y contratos existentes ChemIO/Qt; registrar baseline acotada sin ejecutar pytest monolítico.
- [x] Reproducir ambos enantiómeros, rutas de guardado y apertura; fijar excepciones y escalas previas.
- [x] Reproducir cierre en subprocesos normal/modificado/exportado/recuperado; localizar callback `QUndoStack` con GDB; C superpuestas no reproducibles en recuperación controlada.

## 2. MOL y flujo de archivos
- [x] Implementar normalización MOL compartida, positiva y uniforme, sin mutar MolGraph fuente; aplicarla a RDKit/worker/fallback.
- [x] Añadir extensión MOL por filtro, actualizar ruta sólo tras escritura exitosa y conservar guardado `.cmsn`.
- [x] Reimportar fixture S y control R desde archivos MOL propios y RDKit mediante FileController.
- [x] Validar elementos, conteos, conectividad, órdenes, cargas, isótopos, configuración CIP, extensión y coordenadas.
- [x] Cubrir fallback cuando representa fielmente el caso y fallo explícito cuando no puede preservar stereo.
- [x] Retirar anuncio SDF y rechazar `.sdf` con mensaje; documentar implementación multi-registro fuera de alcance.
- [x] Cubrir errores de lectura/escritura visibles y ausencia de falso estado guardado.

## 3. Cierre Qt y recuperación
- [x] Cubrir cierre normal, posterior a exportar, con documento modificado, destrucción de canvas y señales async en subprocesos acotados.
- [x] Probar recuperación de autosave conservando grafo y representación; investigar C superpuestas sin reescribir persistencia.
- [x] Reproducir el fallo de teardown de canvas con SIGABRT/GDB, corregir la señal activa retenida y repetir el caso; SIGSEGV exacto del AppImage queda sujeto a retest.

## 4. Verificación y entrega
- [x] Ejecutar tests focales ≤5 min y cada comando ≤10 min; matar/limpiar timeout; no repetir idéntica estrategia tras bloqueo.
- [x] Ejecutar arquitectura, OpenSpec estricto, compileall, Ruff focal y `git diff --check`.
- [x] Registrar reportes en OpenSpec; ejecutar test de `.cmsn`/autosave sin alterar formato.
- [x] Crear dos commits funcionales independientes: ChemIO/MOL y teardown Qt demostrado; OpenSpec/reporte en commit documental aparte.
- [x] Consultar autorización normal y publicar sólo la rama solicitada; corregir el manifiesto explícito de node IDs (+15) y confirmar CI por shards para el SHA publicado.
- [x] No declarar DATA-11 completado. Build Preview sólo admite `release/**-prep`, por lo que no se generó paquete desde la rama fix; aceptación manual del propietario queda pendiente.
