# Tareas — disposición visual N–H

## A. Diagnóstico y baseline
- [x] A.1 Confirmar workspace disponible, rama base, estado Git limpio y crear rama aislada solicitada.
- [x] A.2 Leer OpenSpec de continuidad previa y diagnosticar visualmente sus capturas después del fix de trim.
- [x] A.3 Ejecutar y registrar baseline obligatorio (status, compileall, collect-only, pytest completo y Ruff).
- [x] A.4 Capturar el estado visual anterior con gestos reales del editor para urea, acetanilida y aminas.

## B. Layout de glifos
- [x] B.1 Mantener N centrado sobre su átomo y renderizar H/H₂ como glifo independiente, sin rotación.
- [x] B.2 Elegir lado/arriba/abajo por enlaces vecinos y disponibilidad, evitando trazos, carbonilos, cargas y etiquetas próximas.
- [x] B.3 Incluir el H separado en shape/bounds para que el recorte existente siga midiendo la forma visible completa.
- [x] B.4 Mantener sin cambios etiquetas de otros heteroátomos, estilos, geometría del grafo y Clean2D.

## C. Regresiones y evidencia
- [x] C.1 QTests de `Ar–NH–C(=O)–NH–R`, acetanilida y aminas primaria/secundaria mediante eventos Qt reales.
- [x] C.2 Verificar N anclado, orientación H, separación/no solape de glifos y enlaces, longitudes/ángulos, conectividad, preview/commit, undo/redo y CMSN.
- [x] C.3 Capturar comparaciones antes/después en temas claro y oscuro; registrar métricas/reproducibilidad.

## D. Verificación y entrega
- [x] D.1 Correr tests focalizados, compileall, Ruff scoped y `git diff --check`.
- [x] D.2 Ejecutar la suite completa final por shards y reportar los fallos basales sin ocultarlos (`validation.md`).
- [ ] D.3 Crear commit local solo si todos los gates pasan; no push, merge ni PR. Bloqueado por los fallos basales documentados.
