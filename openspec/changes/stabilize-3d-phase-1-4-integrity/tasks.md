# Tasks — estabilización 3D de fases 1–4

## 1. Preparación y baseline

- [x] 1.1 Confirmar rama, upstream, HEAD, worktree y commits recientes.
- [x] 1.2 Leer los OpenSpecs activos de viewport y editor 3D y la evidencia previa sin reemplazarla.
- [x] 1.3 Capturar compilación, colección, Ruff y pruebas focalizadas del HEAD inicial; documentar por qué no se repite `pytest -q`.
- [x] 1.4 Validar este cambio con OpenSpec strict antes de editar producción.

## 2. Hidrógenos derivados

- [ ] 2.1 Añadir primero una regresión RED de rotación rígida que verifique vector C–H global, longitudes y ángulos.
- [ ] 2.2 Transformar H metadata con el mismo traslado/rotación que el fragmento en helpers geométricos, preview y commit.
- [ ] 2.3 Cubrir etano, etanol y butano: rotación de componente, distancia/ángulo/dihedro, cancelación, undo/redo y Live MM posterior.
- [ ] 2.4 Confirmar que H derivados no entran en MolGraph ni `.cmsn` y que una topología incompatible los invalida/regenera.

## 3. Estereoquímica RDKit

- [ ] 3.1 Trazar campos disponibles en `Atom`, `Bond`, importación/persistencia y payload worker; comprobar la asociación de etiquetas por índice/ID estable.
- [ ] 3.2 Añadir tests RED para importación R/S y E/Z, pares R/S y E/Z, controles no especificados, y verificación CIP/bond stereo y geometría.
- [ ] 3.3 Implementar mapeo sólo para encoding inequívoco; rechazar explícitamente etiquetas ambiguas/no soportadas.
- [ ] 3.4 Verificar configuraciones tras Embed y tras cada resultado final de minimización; cubrir que Live MM no invierta configuración sin aviso.

## 4. Aromaticidad y enlaces múltiples

- [ ] 4.1 Añadir pruebas RDKit para benceno explícito/aromático y Kekulé, ciclohexano, eteno y etino.
- [ ] 4.2 Preservar aromaticidad y orden en escena y render; no aromatizar anillos por heurística visual.
- [ ] 4.3 Probar mallas para enlaces dobles, triples y aromáticos sin regresión de picking.

## 5. Estabilidad y usabilidad

- [ ] 5.1 Documentar análisis del core y formular un repro mínimo seguro para el SIGSEGV Qt/SIP; mantener separado el fallo de suite.
- [ ] 5.2 Añadir regresiones deterministas permitidas para cierre/concurrencia; no usar offscreen para declarar resuelto el fallo gráfico.
- [ ] 5.3 Compactar controles del dock y mejorar contraste/rendimiento sólo con evidencia; conservar Live MM accesible.
- [ ] 5.4 No reintentar GUI si vuelve a requerir aprobación adicional; entregar pasos manuales y dejar aceptación visual pendiente.

## 6. Validación y entrega

- [ ] 6.1 Ejecutar pruebas focalizadas por cambio, regresiones 2D/arquitectura, compileall, Ruff scoped, OpenSpec strict y `git diff --check`.
- [ ] 6.2 No repetir la suite completa SIGSEGV conocida; documentar la excepción y resultados previos sin reescribirlos.
- [ ] 6.3 Añadir evidencia nueva como addendum, separada de informes históricos.
- [ ] 6.4 Commits locales descriptivos por área real; rama feature solamente, sin push/merge/PR; verificar worktree limpio.