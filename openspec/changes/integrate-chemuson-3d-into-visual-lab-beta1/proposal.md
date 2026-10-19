# Propuesta — integración aislada de ChemUSON 3D en Beta 1

## Por qué

La implementación de 3D está conservada en `origin/feature/chemuson-3d` (HEAD `ec486626333ccd2e0ae86671861708b30c508c5b`, 36 commits propios desde `origin/main`). La rama visual de Beta 1, `origin/feat/molecular-depiction-visual-lab` (HEAD `43a8a27c45cdabfd157bf3f2a52c9216b24a1552`), contiene trabajo posterior de Clean2D, depiction, proyecciones y calibración visual que debe preservarse. Se requiere una integración con historial completo en un worktree local separado, no una copia de archivos ni cherry-picks.

## Alcance

- Integrar mediante merge normal los historiales de 3D y de la rama visual en `integration/chemuson-3d-beta1`.
- Resolver y revisar explícitamente los seis conflictos de fusión; auditar además los archivos compartidos fusionados automáticamente.
- Ejecutar pruebas focales de 2D/3D/ChemIO, arquitectura, serialización y flujos de GUI; verificar colección, compilación, Ruff y OpenSpec.
- Documentar resultados, bloqueos, SHA y estado final del worktree.

## Fuera de alcance e invariantes

No reescribir, squash, cherry-pick, publicar ni alterar los commits de origen; no cambiar Clean2D, ChemName, identidad química, estereoquímica, formatos `.cmsn`, proyecciones editoriales ni la salida/calibración visual de la rama objetivo. No cambiar `main`, `release/v0.3.0-beta.1-prep` ni ramas remotas; no push, PR, tag o release. La rama de integración temporal es local. Si un conflicto no puede resolverse preservando ambos contratos sin cambio conductual no aprobado, detenerse y documentarlo.

## Módulos probables afectados

`src/chemuson/chemio/` (worker RDKit y límites seguros), `src/chemuson/gui/` (shell, canvas y docks), `architecture/modules.yml`, `AGENT_REPORT.md`, pruebas GUI/ChemIO y artefactos OpenSpec de ambas líneas.