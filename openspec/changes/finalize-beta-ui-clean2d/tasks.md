# Tareas — cierre conjunto UI + Clean2D

## 1. Consolidación y baseline
- [x] 1.1 Verificar HEAD remoto/local, árbol limpio y base beta; crear rama conjunta desde `9bc8192` y preservar ramas fuente.
- [x] 1.2 Leer los OpenSpec de beta, guardado/menús/pestañas y Clean2D; registrar baseline antes de nuevos cambios.
- [x] 1.3 Incorporar ambas ramas con trazabilidad y reconciliar node IDs contra la colecta real sin perder tests vigentes.
- [x] 1.4 Ejecutar en un mismo proceso los tests UI y Clean2D; reproducir además la aplicación real para asegurar que el problema del propietario no está oculto por pruebas sintéticas.

## 2. Cabecera de dos filas
- [x] 2.1 Retirar marca/matraz/versión de AppBar, conservar branding/versionado en los lugares oficiales existentes.
- [x] 2.2 Integrar Nuevo/Abrir/Guardar/Pegar como controles compactos con las QAction originales y sin botón Nuevo duplicado; corregir duplicado actual de Pegar en Editar.
- [x] 2.3 Ocultar la tercera fila de acciones en defaults nuevos con migración de preferencias, mantener los menús tradicionales visibles y conservar Undo/Redo solo a la derecha.
- [x] 2.4 Añadir pruebas UI de distribución, accesibilidad, shortcut Ctrl+S, pestañas limpias/sucias/X, Save/Discard/Cancel, temas, panel oculto y tamaños compactos.
- [x] 2.5 Generar captura de ventana reproducible y guardar métricas/evidencia de la cabecera resultante.

## 3. Diagnóstico y corrección estereoquímica
- [x] 3.1 Registrar la molécula de referencia `C[C@H](F)CO`, conectividad y cuñas reales; cubrir documento multicomponente con selección de una molécula, todo y selección parcial con el centro.
- [x] 3.2 Registrar clasificación/política, IDs seleccionados, fuentes/candidatos/rechazos, candidato, operación, coordenadas y resultado de quick/Ctrl+K, publicación/Ctrl+Shift+K y propuesta/Ctrl+Alt+K en `validation.md`.
- [x] 3.3 Corregir únicamente la creación duplicada de bloques de centro estereoquímico por ambos extremos de un enlace dirigido, manteniendo intacta la política para centros distintos.
- [x] 3.4 Corregir `wedge_signature(graph) == wedge_signature(graph)` para usar una referencia capturada antes del cambio/candidato.
- [x] 3.5 Añadir asignación CIP R/S independiente + round-trip para dos enantiómeros; comprobar también que una cuña aislada sobre un centro no especificado/tied no inventa asignación química.
- [x] 3.6 Verificar invariantes de conectividad/identidad/cuñas, selección, determinismo, candidatos rechazados y Undo; no aceptar ni forzar candidato inseguro.
- [x] 3.7 Añadir un gate GUI independiente que rechaza un fallback cacheado con R/S invertido aunque la firma cíclica no cambie.

## 4. Validación acotada y entrega
- [x] 4.1 Ejecutar tests focalizados UI y Clean2D; arquitectura, OpenSpec estricto, Ruff focal, compileall y `git diff --check`.
- [x] 4.2 Ejecutar la suite completa por shards de CI con límite externo ≤300 s por shard; reconciliar node IDs y reportar CI exacto en `validation.md`.
- [x] 4.3 Registrar archivos separados por UI/Clean2D, commits, SHAs, CI y limitaciones; mantener pendiente la aceptación manual del propietario.
- [x] 4.4 Publicar por push normal solo `fix/beta-final-ui-clean2d`; nunca actualizar beta/main ni publicar tags, releases o paquetes oficiales.
