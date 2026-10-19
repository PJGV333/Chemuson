# Tareas — integración aislada 3D + visual lab Beta 1

## 1. Alcance y baseline
- [x] 1.1 Confirmar `AGENTS.md`, los OpenSpecs de integración/visualización y el OpenSpec 3D de origen.
- [x] 1.2 Verificar rama temporal, base, checkpoint, refs y árbol limpio; registrar baseline en `baseline.md` antes de editar.
- [x] 1.3 Hacer una previsualización de conflictos sin cambiar el árbol y registrar la lista.

## 2. Integración con historial
- [x] 2.1 Guardar localmente la propuesta/diseño/tareas/spec/baseline y validar OpenSpec estricto.
- [x] 2.2 Iniciar el merge normal `--no-ff --no-commit` de la referencia 3D; conservar ambos historiales.
- [x] 2.3 Resolver los seis conflictos manualmente, con atención a catálogo, worker/ChemIO y builder GUI.
- [x] 2.4 Revisar los merges automáticos de canvas, shell, docks y pruebas 3D; confirmar que se preservan las capacidades 2D/visuales probadas.
- [x] 2.5 Inspeccionar el diff y las invariantes químicas, de persistencia, Clean2D, proyecciones y arquitectura; sin cambios fuera del alcance 3D y la corrección de fidelidad de marca MOL ambigua.

## 3. Verificación y entrega
- [x] 3.1 Ejecutar pruebas focales 2D y 3D de ambos lados, más regresiones de conflictos, ChemIO/`.cmsn`, arquitectura y workflow de pruebas.
- [x] 3.2 Ejecutar compilación, colección, Ruff scoped, validación OpenSpec y `git diff --check`; documentar Ruff global baseline.
- [x] 3.3 Documentar el workflow de shards: plan de ocho shards válido; ejecución completa y `pytest -q` omitidos por la restricción Qt/SIP registrada.
- [x] 3.4 Actualizar `validation.md` con resultados reales, capturas/aceptación manual pendiente y límites de readiness.
- [x] 3.5 Crear sólo el merge commit local después de superar las regresiones críticas; verificar ancestros, rama, diff y worktree limpio. Sin push, PR, tag ni merge a release.
- [x] 3.6 Mantener `release/v0.3.0-beta.1-prep` sin cambios mientras haya pruebas, aceptación visual o gates de propietario pendientes.