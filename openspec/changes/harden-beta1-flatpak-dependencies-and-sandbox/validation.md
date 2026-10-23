# Validación — dependencias y sandbox Flatpak beta.1

## Implementación actual

La eliminación de Pillow y el estrechamiento del permiso Flatpak están implementados en la rama aislada `release/v0.3.0-beta.1-security-prep`. El SHA final todavía no existe (cambios no comprometidos). `release/v0.3.0-beta.1-prep` y `main` no se han movido.

Pillow fue eliminado como dependencia directa de `pyproject.toml`/`requirements.txt` y del módulo wheel Flatpak. La inspección del preview #39899351520 detectó que la metadata de RDKit aún lo resolvía transitivamente en el entorno PyInstaller y que sus CArchive contenían recursos Pillow; por eso este follow-up añade `excludes=["PIL"]` a `chemuson.spec` y una comprobación fail-closed de recursos `PIL`/`Pillow` en `validate_packaged_rdkit_worker.py`. No se fija wheel sustituto: el runtime de ChemUSON no usa PIL (QImage/QPixmap; `rdMolDraw2D` importa con PIL bloqueado). El paquete Flatpak usa RDKit `--no-deps` y su árbol no contiene la distribución top-level Pillow. La dependencia transitoria de build se distingue explícitamente de lo distribuido.

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

La consulta OSV del 2026-10-23 encontró Pillow 12.0.0 con 77 registros (incluye CVE-2026-42311 y CVE-2026-54058), Pillow 12.3.0 con 39 registros (incluye CVE-2026-96637, CVE-2026-96643, CVE-2026-96644 y CVE-2026-96778), y cero registros para Pillow 13.0.0. Ninguna versión se añade como dependencia directa ni al Flatpak; 13.0.0 sólo fue resuelto transitivamente por metadata upstream de RDKit en el entorno de build del preview inicial y se excluye del paquete congelado.

## Follow-up de auditoría de CArchive y verificación local (2026-10-23)

El preview inicial de seguridad #39899351520 (SHA `35f5888cd92f0048b0e2c15c97534952682d830d`) no cumple la ausencia de recursos Pillow en Windows/AppImage: Windows portable contiene cinco entradas `PIL/*`; AppImage contiene las mismas cinco más 13 `pillow.libs/*`. Su log registra Pillow 13.0.0 transitivo de RDKit; OSV no devolvió advisories para 13.0.0 en la consulta del 2026-10-23. El bundle Flatpak se importó y extrajo en OSTree temporal: no hay `site-packages/PIL` ni metadata Pillow; hay recursos internos RDKit `rdkit/sping/PIL`, no equivalentes a la distribución Pillow. El run previo queda preservado como histórico y no sirve para retest del candidato corregido.

Para cerrar el hallazgo, el nuevo CArchive validator inspecciona cada ejecutable PyInstaller y rechaza entradas top-level `PIL`/`Pillow`; la spec excluye `PIL`. La cobertura incluye rutas Windows con `\\`, módulos Linux, metadatos Pillow y falsos positivos de `rdkit/sping/PIL`/PyQt.

- Baseline follow-up previo a edición: worktree limpio en `35f5888c`; tests RDKit/paquete ChemName/release contract: **20 passed in 2.35s**.
- Tras implementar la exclusión: **73 passed in 3.80s** en tests de worker/RDKit, ChemName package, dependencia, imagen, distribución Linux, AppImage, preview/release workflows y manifest node-ID. No se ejecutó la suite completa ni los ocho shards.
- `compileall src tests tools packaging`: PASS; colección **2.630 node IDs**, planificador de 8 shards válido (estima ~208.7 s por shard). Ruff focal con las reglas F401/F811/F821/E722/E741: PASS. OpenSpec estricto y `git diff --check`: PASS.
- PyInstaller local (CachyOS, Python 3.14, PyInstaller 6.22.3) produjo un ejecutable temporal con **0** entradas top-level `PIL`/`Pillow`; `validate_packaged_rdkit_worker.py` pasó sobre ese binario congelado (worker nativo, descriptors, SMILES y MMFF94 3D). El build local emitió advertencias de Qt/ldd del host; el paquete no sustituye Build Preview CI.

## Build Preview y artefactos exactos

Pendiente: ejecutar el Build Preview del SHA con la exclusión PyInstaller. El run debe entregar Windows portable, Windows installer, AppImage Type 2 y Flatpak, `publication=false`, hashes/provenance concordantes y validación RDKit/ChemName, con **cero** entradas Pillow en CArchive y sin el paquete top-level en Flatpak. Verificar el smoke 2D→3D congelado de AppImage y controles empaquetados; no atribuirlos a aceptación manual.

## Riesgos, aceptación y decisión

- **P1 pendiente:** el propietario debe instalar y ejecutar los cuatro artefactos exactos. Para Flatpak, especialmente probar abrir/guardar un `.cmsn` fuera de los directorios de estado y reabrirlo mediante el portal sin permiso global al home.
- El permiso de red Flatpak sigue siendo general; habilita conectividad para endpoints configurados y resolución externa opt-in. Usar esas funciones transmite datos al proveedor/destino elegido; la política de red no limita destinos ni sustituye consentimiento.
- P0 de persistencia/estereoquímica/cierre y otros P1 de la matriz manual siguen sin aceptación exacta.
- Protección `v*` y environments `beta`/`flatpak-beta` continúan pendientes según la auditoría previa; revalidar antes de un tag.
- Sin tag, GitHub Release ni despliegue público. Dictamen: **NO-GO** hasta completar los retests y controles externos.
