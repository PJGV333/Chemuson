# Política de higiene del repositorio ChemUSON

Esta política define cómo reducir deuda histórica sin convertir una limpieza en refactor funcional ni borrar contexto útil. Aplicar junto con `AGENTS.md`, el OpenSpec activo y `docs/architecture.md`; si hay conflicto, se detiene el cambio y se solicita decisión.

## 1. Fuente de alcance y registro

1. Toda campaña comienza con OpenSpec activo: propuesta, diseño, tareas y contratos.
2. Registrar `HEAD`, estado del árbol, inventario/medidas y baseline de pruebas antes del cambio. Guardar comandos, salidas, identidades de fallos y limitaciones en el OpenSpec o el reporte de cierre; no sustituir resultados exactos por “todo pasa”.
3. Mantener una memoria de campaña compacta en `docs/history/CAMPAIGNS.md`: objetivo, decisiones, intentos útiles/no adoptados, fallos, ramas/SHAs relevantes y estado final. Los contratos detallados permanecen en OpenSpec; no duplicarlos íntegramente en informes.
4. El reporte de una limpieza incluye archivos eliminados/mantenidos, motivo, evidencias, métricas before/after, pruebas y ramas revisadas. No se archivan logs temporales si una tabla/resultado verificable basta.

## 2. Política de ramas e historia Git

- Hacer fetch/prune antes de inventariar; registrar nombre, SHA, fecha, ahead/behind, asunto y resultado de `merge-base --is-ancestor <tip> origin/main`.
- Clasificar cada rama como `SAFE_TO_DELETE`, `PRESERVE_UNIQUE`, `ACTIVE/PROTECTED` o `NEEDS_OWNER_REVIEW`.
- `SAFE_TO_DELETE` exige que la punta remota sea ancestro de `origin/main`, que su SHA quede registrado y que la campaña confirme que el trabajo está representado en la rama principal. Borrar solo después de los gates y de un push normal de la rama de mantenimiento.
- Una rama no ancestro no se borra por su nombre, edad, similitud visual o porque parte de su contenido parezca integrado. Revisar commits/árbol; documentar y conservar trabajo único hasta decisión explícita. Proteger ramas de publicación (`gh-pages`) y experimentos químicos activos.
- No borrar `main`, no integrar la rama de higiene en `main`, no usar force-push, rebase/cherry-pick/merge commit ni reescribir historia durante mantenimiento.
- Eliminar refs no equivale a compactar objetos: no ejecutar `git gc --prune`, BFG, `filter-repo` o `filter-branch` como parte de la limpieza. Una propuesta futura para cambiar historia necesita otro alcance, lista de ramas afectadas, backup, estimación real de pack y decisión del propietario.

## 3. Auditoría y eliminación de código

Antes de borrar una API, función, módulo o script:

1. Buscar imports, llamadas, `entry_points`, scripts CLI, plugins, `importlib`, referencias de packaging, herramientas de release, CI, documentación operativa y consumidores por nombre/ruta.
2. Revisar configuración de descubrimiento de pytest, fixtures, datos, snapshots, tests parametrizados/indirectos y usos dinámicos; ausencia de import directo no basta.
3. Comprobar pertenencia al paquete con `git ls-files`, configuración de build y manifests; diferenciar un script en `src/` de uno dev en `tools/`.
4. Clasificar como `KEEP`, `SAFE_DELETE` o `NEEDS_REVIEW`, y anotar evidencia. Borrar únicamente `SAFE_DELETE` dentro del alcance OpenSpec.
5. Si se elimina un test, debe existir cobertura semántica equivalente identificable caso por caso. Nunca quitar/renombrar tests químicos, fixtures, snapshots o baselines para obtener verde.

No cambiar resultados, geometría, serialización, APIs de producto, dependencias, jerarquía/eventos Qt o heurísticas de Clean2D/ChemName por conveniencia de higiene. Una importación/dependencia nueva requiere actualizar `architecture/modules.yml` y cumplir sus límites.

## 4. Artefactos, binarios y documentación

- Mantener como fuente canónica el código, tests actuales, assets de runtime, fixtures/base de comparación vigentes, especificaciones globales y Markdown de OpenSpec archivado.
- Mantener una cantidad pequeña de evidencia visual final, referenciada desde su reporte. Para una misma UI, conservar una captura canonical antes/después o salida representativa; retirar capturas intermedias, contact sheets redundantes y outputs que se puedan regenerar, tras resumir su resultado/decisión en Markdown.
- Mantener evidencia manual única si no la sustituye una prueba reproducible (p. ej. captura real de plataforma aprobada). Mantener JSON de baseline solo si describe un contrato/resultado reproducible, tiene dueño y no es secreto/cache.
- Los generadores de pruebas visuales deben escribir a una ruta explícita temporal o entregar su salida fuera del árbol fuente por defecto. Nunca añadir `__pycache__`, builds, caches, dumps, logs, PDFs/PS de depuración o imágenes sin consumidor a `src/`, `tests/` o un archive de producto.
- No crear ignores globales para tapar candidatos: añadir reglas exactas solo para outputs generados recurrentes, justificados y nunca requeridos como fixtures.
- Una documentación histórica obsoleta se resume en `CAMPAIGNS.md` y se retira solo después de comprobar enlaces entrantes y preservar decisiones/no-resultados. Mantener `AGENT_REPORT.md` como registro actual breve y no duplicar todo OpenSpec.
- Conservar outputs de infraestructura publicada mientras sirvan tráfico/instalación. En particular, el repositorio Flatpak servido por `gh-pages` no es basura local; migrarlo requiere plan de publicación separado.

## 5. Verificación de cierre

Como mínimo, ejecutar lo aplicable desde el OpenSpec y comparar con baseline:

```text
python -m compileall src tests tools packaging
pytest --collect-only -q
pytest -q
pytest -q tests/architecture
ruff check src tests tools packaging --select F401,F811,F821,E722,E741
openspec validate --all --strict
git diff --check
```

Además, validar smoke de UI/packaging/assets cuando el alcance borre artefactos o tooling relacionado. Investigar un fallo nuevo y detenerse ante regresión, pérdida de fixture/recurso, dependencia no catalogada o cambio de comportamiento no autorizado. Reportar con precisión los fallos heredados, herramientas no disponibles y verificaciones no ejecutadas.
