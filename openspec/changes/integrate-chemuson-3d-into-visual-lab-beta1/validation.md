# Validación — integración 3D en visual lab Beta 1

## Estado

**EN PREPARACIÓN; todavía no se inició el merge de código.** El baseline quedó registrado en `baseline.md`. La base visual está fijada en `43a8a27c45cdabfd157bf3f2a52c9216b24a1552`; la referencia 3D, en `ec486626333ccd2e0ae86671861708b30c508c5b`. `git merge-tree` anticipa seis conflictos y varios merges automáticos que requieren inspección.

## Baseline confirmado

- Compilación: exit 0 con Python 3.13.14 del entorno del repositorio.
- Colección: 2358 tests.
- Regresión focal 2D/depiction/arquitectura: 450 passed, 3 skipped.
- Ruff global selectivo: un único F401 histórico (`math`, `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`); no corregido por estar fuera de alcance.
- Suite monolítica: no ejecutada por la restricción vigente sobre el SIGSEGV Qt/SIP histórico.

## Gates pendientes

No hay todavía resultados posteriores al merge. No declarar integración aprobada ni lista para release hasta resolver conflictos, auditar merges automáticos, ejecutar pruebas cruzadas 2D/3D/ChemIO y arquitectura, validar todos los cambios y conservar reportados los resultados fallidos, incompletos o dependientes de Qt. La aceptación manual integrada OpenGL/GUI y la revisión visual del propietario no se ejecutarán en este entorno.
