# Tareas — integración aislada 3D + visual lab Beta 1

## 1. Alcance y baseline
- [x] 1.1 Confirmar `AGENTS.md`, los OpenSpecs de integración/visualización y el OpenSpec 3D de origen.
- [x] 1.2 Verificar rama temporal, base, checkpoint, refs y árbol limpio; registrar baseline en `baseline.md` antes de editar.
- [x] 1.3 Hacer una previsualización de conflictos sin cambiar el árbol y registrar la lista.

## 2. Integración con historial
- [ ] 2.1 Guardar localmente la propuesta/diseño/tareas/spec/baseline y validar OpenSpec estricto.
- [ ] 2.2 Iniciar el merge normal `--no-ff --no-commit` de la referencia 3D; conservar ambos historiales.
- [ ] 2.3 Resolver los seis conflictos manualmente, con atención a catálogo, worker/ChemIO y builder GUI.
- [ ] 2.4 Revisar los merges automáticos de canvas, shell, docks y pruebas 3D; confirmar que se preservan las capacidades 2D/visuales.
- [ ] 2.5 Inspeccionar el diff completo y las invariantes químicas, de persistencia, Clean2D, proyecciones y arquitectura antes de concluir el merge.

## 3. Verificación y entrega
- [ ] 3.1 Ejecutar pruebas focales 2D y 3D de ambos lados, más regresiones de conflictos, ChemIO/`.cmsn`, arquitectura y workflow de pruebas.
- [ ] 3.2 Ejecutar compilación, colección, Ruff scoped, validación OpenSpec y `git diff --check`; explicar cada hallazgo Ruff baseline.
- [ ] 3.3 Ejecutar o documentar los shards/validaciones integradas requeridas; no declarar éxito si falta un gate crítico.
- [ ] 3.4 Actualizar `validation.md` con resultados reales, capturas/aceptación manual pendiente y decisión de readiness.
- [ ] 3.5 Si los gates críticos automatizados permiten integrar, hacer sólo un commit de merge local; verificar ancestros, rama, diff y worktree limpio. No push, PR, tag ni merge a release.
- [ ] 3.6 Mantener `release/v0.3.0-beta.1-prep` sin cambios mientras haya pruebas, aceptación visual o gates de propietario pendientes.