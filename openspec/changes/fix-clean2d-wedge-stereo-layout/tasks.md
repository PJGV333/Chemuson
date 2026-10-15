# Tareas — Clean2D con cuñas

## 1. Alcance, Git y baseline
- [x] 1.1 Verificar ramas/upstream reales y árbol limpio; crear `clean2d/wedge-bond-campaign` desde `origin/release/v0.3.0-beta.1-prep` en `9bc8192`, sin transportar la Fase A.
- [x] 1.2 Registrar baseline acotado (compileall, colecta, Ruff focal, grupo Clean2D de 240 casos) en `baseline.md`; no ejecutar suite monolítica en la captura.
- [x] 1.3 Reproducir el defecto con sonda de motor y sonda GUI sobre la base sin parche; documentar estados, fuentes y motivos.

## 2. Diagnóstico
- [x] 2.1 Confirmar que el motor **sí** genera candidatos para moléculas con cuñas y que el rechazo proviene de `cambio_firma_estereoquimica`, no de la ausencia de candidatos.
- [x] 2.2 Medir el artefacto de la firma antigua: orden cíclico conservado pero signo de un par invertido al cruzar el hueco angular los 180°.
- [x] 2.3 Detectar la ceguera ante vecinos de frontera en selecciones parciales.
- [x] 2.4 Descartar causas de importación/representación: el estilo `WEDGE`/`HASHED` llega intacto al grafo y al motor.

## 3. Mecanismo
- [x] 3.1 Reescribir `stereo_layout_signature()` sobre el orden cíclico canónico de los vecinos estructurales reales de cada centro estereogénico.
- [x] 3.2 Incluir los vecinos y enlaces de frontera usando las coordenadas del grafo (`_structural_bonds_touching_selection`, `_layout_point`).
- [x] 3.3 Conservar estilo de enlace y estereoquímica explícita en la firma, para que cambiar el sentido de una cuña siga detectándose.
- [x] 3.4 Mantener sin tocar la política de complejidad, el ranking, el vocabulario estable de motivos, los modelos químicos, ChemName, ChemIO, la GUI y la persistencia.

## 4. Matriz de regresión permanente
- [x] 4.1 Añadir `tests/test_clean2d_wedge_stereo.py` con los siete casos pedidos: control acíclico, cuña sólida, cuña rayada, centro tetraédrico con ambas, varios centros quirales, selección con centro quiral y selección con cuña en frontera.
- [x] 4.2 Fijar los contratos de la firma: rotación y re-encuadres que conservan el orden cíclico no cambian la firma; imagen especular, intercambio de vecinos, cambio de cuña y cruce de un vecino de frontera sí la cambian.
- [x] 4.3 Verificar invariantes químicas duras con `assert_clean2d_invariants` y comparación independiente de la configuración dibujada (orden cíclico + estilos).
- [x] 4.4 Añadir determinismo (dos pasadas, mismo `geometry_hash` y mismas coordenadas) y camino GUI con Undo.

## 5. Validación y entrega
- [x] 5.1 Ejecutar red/green: la matriz nueva falla en la base sin parche y pasa con el parche.
- [x] 5.2 Ejecutar compileall, colecta, Ruff focal y global, suite de arquitectura y grupo Clean2D completo, cada comando con timeout externo.
- [x] 5.3 Ejecutar la suite completa por shards con `tools/ci_pytest_shards.py` (plan/run/verify) y registrar conteos, omitidos y fallos.
- [x] 5.4 Validar OpenSpec en modo estricto y `git diff --check`.
- [x] 5.5 Registrar archivos, decisiones, resultados y estado Git en `validation.md`.
- [x] 5.6 Crear commits pequeños y publicar sólo la rama `clean2d/wedge-bond-campaign`; sin merge a beta ni `main`, sin tags ni releases.
