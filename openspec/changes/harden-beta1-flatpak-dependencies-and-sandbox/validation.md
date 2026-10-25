# Validación — dependencias y sandbox Flatpak beta.1

## Implementación actual

La eliminación de Pillow como dependencia directa, exclusión de recursos PyInstaller y estrechamiento del permiso Flatpak están implementados en la rama aislada `release/v0.3.0-beta.1-security-prep`, commit `e4710b938f0096ff2d02f35010117bb9a6ed8286`. El push fue normal y limitado a esa rama; `release/v0.3.0-beta.1-prep` (`08a90c199f526639453a8ae8b82da957a1458d08`) y `main` (`f0fde603371255902bf0e630ca1ce032b8f16ad3`) no se movieron.

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

## Build Preview y artefactos exactos — #39900808044

Build Preview [#39900808044](https://github.com/PJGV333/Chemuson/actions/runs/39900808044) se ejecutó por push del commit e4710b93 y terminó `success` en source validation, Windows portable/installer, AppImage Type 2, Flatpak y reporte. `preview-build-report.json` (artifact [12478957557](https://github.com/PJGV333/Chemuson/actions/runs/39900808044/artifacts/12478957557)) declara `complete=true`, `publication=false`, versión `0.3.0-beta.1`, rama `release/v0.3.0-beta.1-security-prep` y SHA completo `e4710b938f0096ff2d02f35010117bb9a6ed8286`.

| Formato | Artifact ID | Archivo | Tamaño (bytes) | SHA-256 |
|---|---:|---|---:|---|
| Windows portable | [12479081901](https://github.com/PJGV333/Chemuson/actions/runs/39900808044/artifacts/12479081901) | `Chemuson-v0.3.0-beta.1-preview-e4710b93-windows-x86_64-portable.exe` | 127027888 | `d107a43dadb64c137dafb6e0278cc38118e2b911559c411654034ca458846754` |
| Windows installer | [12478962316](https://github.com/PJGV333/Chemuson/actions/runs/39900808044/artifacts/12478962316) | `Chemuson-v0.3.0-beta.1-preview-e4710b93-windows-x86_64-setup.exe` | 128133031 | `9859c922b2be0707a38afdb3c2d316b8d7125a0c516c6d4dc25e57d9e3c58f64` |
| AppImage Type 2 | [12478253483](https://github.com/PJGV333/Chemuson/actions/runs/39900808044/artifacts/12478253483) | `Chemuson-v0.3.0-beta.1-preview-e4710b93-linux-x86_64.AppImage` | 169608384 | `4abb85948b3c4eccd41d51340cded2469cf7e8b5a3d04ce5eca1bbfc53009169` |
| Flatpak | [12478758199](https://github.com/PJGV333/Chemuson/actions/runs/39900808044/artifacts/12478758199) | `Chemuson-v0.3.0-beta.1-preview-e4710b93-linux-x86_64.flatpak` | 104897448 | `828faf67f74496040594fa396530be0c1e2af92271b9eb9f6ab226e45d2556c3` |

Los cuatro artifacts se descargaron; `sha256sum -c checksums.sha256` pasó para cada paquete y su `preview-provenance.json`. Hash, tamaño, versión, source SHA, rama y `publication=false` coinciden entre archivo y provenance. El run registra Pillow 13.0.0 transitivo de RDKit en los entornos Windows y Linux PyInstaller; la consulta OSV del 2026-10-23 no reportó advisories para esa versión.

**Auditoría de recursos descargados:** el CArchive Windows portable y el ejecutable `usr/bin/Chemuson` extraído del AppImage Type 2 tienen **0 entradas** top-level `PIL`/`Pillow` al pasar por `_pillow_archive_entries`; los logs Windows/AppImage también reportan `pillow_archive_entries: []`. El bundle Flatpak se importó a un repo OSTree temporal y su árbol carece de `site-packages/PIL` y `Pillow-*.dist-info`; permanece `rdkit/sping/PIL`, que son recursos legacy de RDKit. El manifiesto conserva sólo permisos scoped de home más `--share=network`; no se concedió permiso total. La build CI valida 9/9 plantillas y 6/6 casos ChemName en Flatpak, y 9/9 + 6/6 en AppImage; Windows valida ChemName y worker RDKit. AppImage también pasó los smokes de iconos, RDKit worker, flujo 2D→3D congelado y checks Type 2/ausencia de sidecars públicos.

Estas comprobaciones acreditan empaquetado/procedencia, no instalación, portal ni aceptación GUI manual. Los casos P0/P1 del propietario siguen `NOT TESTED`; no se eleva ninguna fila por los smokes automatizados.

Tras registrar run/artifacts en OpenSpec y notas de aceptación: OpenSpec estricto PASS; `python -m compileall -q src tests tools packaging` PASS; pruebas focalizadas de contratos Flatpak/distribución/preview/AppImage/RDKit/ChemName: **58 passed in 3.23s**; `git diff --check` PASS. No se ejecutó la suite completa ni `test.yml`.

## Riesgos, aceptación y decisión

- **P1 pendiente:** el propietario debe instalar y ejecutar los cuatro artefactos exactos. Para Flatpak, especialmente probar abrir/guardar un `.cmsn` fuera de los directorios de estado y reabrirlo mediante el portal sin permiso global al home.
- El permiso de red Flatpak sigue siendo general; habilita conectividad para endpoints configurados y resolución externa opt-in. Usar esas funciones transmite datos al proveedor/destino elegido; la política de red no limita destinos ni sustituye consentimiento.
- P0 de persistencia/estereoquímica/cierre y otros P1 de la matriz manual siguen sin aceptación exacta.
- Protección `v*` y environments `beta`/`flatpak-beta` continúan pendientes según la auditoría previa; revalidar antes de un tag.
- Sin tag, GitHub Release ni despliegue público al 2026-10-23. Dictamen histórico: **NO-GO** hasta completar los retests y controles externos.

## Addendum — integración y autorización de publicación posterior (2026-10-24)

La rama de seguridad se integró localmente a `release/v0.3.0-beta.1-prep` por fast-forward-only después de comprobar cuatro commits ahead y cero behind. El propietario autorizó la publicación de esta versión como prerelease beta y aportó observaciones manuales limitadas por formato; la matriz contiene alcance exacto y falta de artifact/SHA de esas observaciones. No se convierten otras funciones ni artifacts anteriores en PASS.

El Build Preview #39900808044 (`e4710b938f0096ff2d02f35010117bb9a6ed8286`) es sólo evidencia histórica de ese SHA. El preview final #40144590463 (`85fa4ed2bd0afe1d31953440228fb73a69b54e7e`) pasó en los cuatro formatos y sus provenance/checksums y ausencia de Pillow se verificaron; los artifacts continúan temporales. El tag `v0.3.0-beta.1` se empujó al SHA validado, pero el subject `[skip ci]` impidió disparar `release.yml`; el intento de workflow_dispatch dio HTTP 422 (trigger ausente). No hay GitHub Release, assets oficiales/sidecars ni publicación Flatpak. No mover/recrear tag; se solicita decisión del propietario. Evidencia completa en `openspec/changes/final-stabilization-chemuson-0-3-0-beta-1/validation.md`.
