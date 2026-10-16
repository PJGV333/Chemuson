# Validación de consolidación — viewport molecular 3D

**Addendum:** 2026-10-16. Complementa, no reemplaza, `baseline.md` ni la evidencia histórica de la demostración. Rama verificada: `feature/chemuson-3d`, upstream `origin/feature/chemuson-3d`; HEAD de inicio de esta continuación: `a7dc0dcadf3d6a5fe735ed556f23f1ab5373cf6e`.

## Evidencia focalizada

- `pytest` focalizado para viewport, geometría, persistencia, comandos, Live MM, pestañas y shutdown: **100 passed in 51.15s**, Qt offscreen.
- Contratos arquitectónicos: **280 passed in 22.55s**.
- `python -m compileall -q src tests tools packaging`: exit 0.
- Ruff scoped sobre los archivos Python modificados (`F401,F811,F821,E722,E741`): All checks passed.
- `openspec validate add-interactive-3d-molecular-viewport --strict` y `openspec validate add-editable-3d-molecular-editor --strict`: ambos `Change is valid`; el CLI avisó que ignoraba el formato de `rules` en artefactos.
- `git diff --check`: exit 0.

## Aceptación química observada

Se generaron conformeros RDKit desde los grafos de agua, etano, etanol, butano y benceno; la escena mantuvo los IDs documentales, todos los XYZ fueron finitos y los hidrógenos fueron derivados. Los rangos Z medidos fueron respectivamente **0.348, 2.579, 2.316, 3.395 y 4.710 Å**, por lo que no se devolvió una simple proyección con Z=0. Agua obtuvo un ángulo H–O–H de **103.978°**. Benceno se representó con seis enlaces en ciclo y una forma Kekulé de tres simples/tres dobles; los flags `is_aromatic` explícitos no aparecen en esta ruta de importación.

El worker Live MM de butano produjo **40 snapshots progresivos por campo** para MMFF94, MMFF94s y UFF (iteraciones 2–80), manteniendo fijo el átomo restringido; los tests validan coordenadas y energía cambiantes.

En una corrida pequeña de un solo proceso, la generación RDKit tardó **0.164–0.239 s** por conformero y ensamblar `SceneMolecule3D`, **0.000067–0.000106 s**. No se midió FPS ni interacción GPU; esos valores no justifican optimizar shaders o geometría de navegación.

## Límites no resueltos

- No se realizó una nueva demostración gráfica real ni se inspeccionó la captura/OpenGL; los cambios visuales, el encuadre, el resize y la compresión del dock quedan sin validar.
- El shard automatizado de cierre Qt offscreen terminó correctamente (2 tests, exit 0), pero no reproduce el cierre del demo con contexto OpenGL. No se atribuye causa común a ese caso y al SIGSEGV del worker en la suite baseline; ninguno se declara resuelto.
- La reconstrucción de quiralidad desde `BondStereo` wedge/hash o `stereo_cip` en el worker no está verificada; no afirmar preservación estereoquímica 3D. El ciclo de benceno tampoco conserva flags aromáticos explícitos en la ruta de importación probada.
- No se reejecutó la suite completa por el SIGSEGV de baseline ya documentado.

Commits locales de código/tests: `58383965` (`fix(3d): reject invalid molecular builder edits`) y `1b7fb1ec` (`test(3d): add chemistry and live-mm acceptance coverage`). El informe de viewport forma parte del commit documental separado; no se hizo push.

Prueba manual pendiente en CachyOS desde la raíz del repo: `QT_QPA_PLATFORM=xcb PYTHONPATH=src ./.venv/bin/python -m chemuson`. Ejecutar ambos recorridos de aceptación y confirmar visualmente el renderizado y el cierre limpio; esta instrucción no constituye evidencia de que se haya ejecutado aquí.
