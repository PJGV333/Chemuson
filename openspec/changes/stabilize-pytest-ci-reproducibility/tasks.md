# Tareas — pytest y CI reproducibles

## 1. Alcance, base y evidencia
- [x] 1.1 Verificar branch, SHA exacto, upstream y árbol limpio; crear `fix/ci-pytest-stabilization` desde el upstream solicitado.
- [x] 1.2 Leer OpenSpec aplicable y registrar baseline antes de tocar código/documentación.
- [x] 1.3 Consultar runs #342/#343, reproducir cada nodo nombrado por separado y formar los grupos mínimos de contaminación.
- [ ] 1.4 Resolver en Actions si el `Inspector x=-3` reportado es clipping estable o layout transitorio. Localmente, en geometría estable, se observa `x=2`; la prueba sigue exigiendo ausencia de clipping. Pendiente de un run/log remoto porque no hay autenticación.

## 2. Regresiones focalizadas
- [x] 2.1 Verificar con spy la llamada real a `clean2d_isolated` y comparar el candidato eliminado con la geometría deduplicada equivalente.
- [x] 2.2 Aislar/restaurar `QSettings.NativeFormat/UserScope` en el fixture del diálogo y cubrir los valores predeterminados.
- [x] 2.3 Ejecutar los dos contratos de smoke de GUI en procesos independientes y conservar íntegra la detección de módulos/imports.
- [x] 2.4 Instrumentar `worker.finished`, `QThread.finished` y `job_finished` en la regresión CompChem sin ampliar el timeout.
- [x] 2.5 Esperar geometría de tabs estable y mantener las cotas originales de layout; cambiar producción sólo si el clipping estable se reproduce.

## 3. Plan de cobertura y GitHub Actions
- [x] 3.1 Crear y revisar el manifiesto ordenado de 2.130 IDs recolectados; el plan verifica coincidencia exacta.
- [x] 3.2 Implementar partición determinista, comprobación de cobertura sin huecos/duplicados, supervisor externo de 300 s y resultados/JUnit/log por shard.
- [x] 3.3 Reemplazar el job pytest monolítico por plan + ocho shards + resumen fallido ante resultados faltantes; mantener `windows-smoke` y `flatpak-smoke`.
- [x] 3.4 Actualizar contratos estáticos para probar límites, matriz, manifiesto y ausencia de exclusiones/errores enmascarados.

## 4. Validación y entrega
- [x] 4.1 Ejecutar cada nodo/grupo relevante en proceso nuevo, cada comando con timeout externo ≤300 s; no repetir bloqueos sin cambiar el diagnóstico.
- [x] 4.2 Ejecutar una vez todos los tests mediante los ocho shards locales; cobertura JUnit exacta: 2.130/2.130, 2.111 pasaron y 19 skips preexistentes, cero fallos. Máximo shard 180.759 s; no ejecutar pytest monolítico.
- [x] 4.3 Ejecutar compileall, colección acotada, Ruff focal, OpenSpec estricto, contratos de workflow, arquitectura y `git diff --check`.
- [x] 4.4 Registrar archivos y resultados en `validation.md`; crear commits pequeños sin merge/PR/tag/release.
- [ ] 4.5 Hacer push normal sólo a `fix/ci-pytest-stabilization` y comprobar CI real. Bloqueado: `gh auth status` no detecta sesión y el push HTTPS no pudo solicitar usuario; informar `NOT READY`.
