# Validación — dependencias y sandbox Flatpak beta.1

## Implementación actual

La eliminación de Pillow y el estrechamiento del permiso Flatpak están implementados en la rama aislada `release/v0.3.0-beta.1-security-prep`. El SHA final todavía no existe (cambios no comprometidos). `release/v0.3.0-beta.1-prep` y `main` no se han movido.

Pillow fue eliminado de `pyproject.toml`, `requirements.txt` y el módulo wheel Flatpak. No se fijó wheel sustituto: las rutas de anotaciones usan `QImage`/`QPixmap`, y el módulo nativo `rdMolDraw2D` no requiere PIL. Se corrigió la descripción de dependencias del manual. El contrato de dependencias/Flatpak y la aceptación manual `DIST-07` cubren la ausencia del paquete y el acceso por portal.

Flatpak ya no concede el home completo. Conserva sólo `home/.chemuson:create` y `home/.config/Chemuson:create`; `--share=network` permanece para el Assistant configurado y las consultas externas explícitas. Las notas públicas advierten que el uso de esas funciones puede enviar datos a servicios elegidos y que la interacción de archivos depende del portal. La prueba de portal instalado no está realizada.

## Verificación local focalizada

- Auditoría sin PIL: import blocker rechaza todo `PIL`; `QImage` se crea y `rdkit.Chem.Draw.rdMolDraw2D` importa. `tests/test_image_clipboard_persistence.py`: **9 passed in 0.43s** bajo ese blocker (antes de los cambios de dependencias).
- Después del cambio: `chem/bin/python -m pytest -q tests/test_ai_assistant_integration.py tests/test_chemname_packaged_smoke.py tests/test_image_clipboard_persistence.py tests/test_linux_distribution_scripts.py tests/test_preview_workflow_contract.py tests/test_appimage_validation.py`: **48 passed in 2.67s**.
- `chem/bin/python -m pytest tests/architecture -q`: **283 passed in 10.51s**.
- Ruff focal en `tests/test_ai_assistant_integration.py tests/test_chemname_packaged_smoke.py`: **All checks passed**.
- `flatpak-builder --show-manifest packaging/flatpak/io.github.PJGV333.Chemuson.yml`: parse PASS; lista permisos contiene `--share=network`, las dos rutas scoped y no hay módulo Pillow.
- `geometry3d.flow_smoke.run_smoke(timeout_s=30.0)`: `ok=true`, origen RDKit, 3/3 coordenadas 3D, MMFF94 convergió; elapsed 0.164 s. Es smoke de fuente con worker real, no evidencia de paquete.
- OpenSpec estricto y `git diff --check`: PASS.
- `python -m compileall src tests tools packaging`: PASS. `chem/bin/python -m pytest --collect-only -q`: **2.624 tests collected in 0.68s**; el manifiesto no cambió. Ruff global conserva exactamente el hallazgo basal F401 `math` sin usar en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no se tocó. `git diff --check`: PASS. Suite pytest completa y `test.yml` no se ejecutan por instrucción explícita; se conserva la colección completa.

La consulta OSV del 2026-10-23 encontró Pillow 12.0.0 con 77 registros (incluye CVE-2026-42311 y CVE-2026-54058), Pillow 12.3.0 con 39 registros (incluye CVE-2026-96637, CVE-2026-96643, CVE-2026-96644 y CVE-2026-96778), y cero registros al consultar Pillow 13.0.0. No se incorpora ninguna de esas versiones porque Pillow es innecesario.

## Build Preview y artefactos exactos

Pendiente. Debe existir un único run manual de Build Preview para el SHA final, con Windows portable, Windows installer, AppImage Type 2 y Flatpak, `publication=false`, hashes/provenance concordantes y validación RDKit/ChemName. Verificar el smoke 2D→3D congelado del AppImage y los controles empaquetados disponibles; no atribuirlo a aceptación manual.

## Riesgos, aceptación y decisión

- **P1 pendiente:** el propietario debe instalar y ejecutar los cuatro artefactos exactos. Para Flatpak, especialmente probar abrir/guardar un `.cmsn` fuera de los directorios de estado y reabrirlo mediante el portal sin permiso global al home.
- El permiso de red Flatpak sigue siendo general; habilita conectividad para endpoints configurados y resolución externa opt-in. Usar esas funciones transmite datos al proveedor/destino elegido; la política de red no limita destinos ni sustituye consentimiento.
- P0 de persistencia/estereoquímica/cierre y otros P1 de la matriz manual siguen sin aceptación exacta.
- Protección `v*` y environments `beta`/`flatpak-beta` continúan pendientes según la auditoría previa; revalidar antes de un tag.
- Sin tag, GitHub Release ni despliegue público. Dictamen: **NO-GO** hasta completar los retests y controles externos.
