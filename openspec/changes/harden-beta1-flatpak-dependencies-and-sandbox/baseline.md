# Baseline — 2026-10-23

## Identidad y árbol

- Inicio: rama `release/v0.3.0-beta.1-prep`, HEAD y remoto `08a90c199f526639453a8ae8b82da957a1458d08`; `origin/main=f0fde603371255902bf0e630ca1ce032b8f16ad3`.
- `git status --short`: salida vacía (árbol limpio).
- Rama aislada creada antes de editar: `release/v0.3.0-beta.1-security-prep`, desde el SHA exacto `08a90c199f526639453a8ae8b82da957a1458d08`.
- Python: `/home/ccachyavgp/Documentos/ChemUSON/.venv/bin/python`, Python 3.14.8; Flatpak 1.18.4, flatpak-builder 1.4.12.

## Comandos de baseline previos a cambios

Comandos ejecutados antes de crear/editar esta campaña; salidas completas guardadas temporalmente en `/tmp/chemuson-beta1-pillow-fix-baseline/`:

```text
python -m compileall src tests tools packaging
exit: 0

chem/bin/python -m pytest --collect-only -q
exit: 0
2624 tests collected in 0.67s

ruff check src tests tools packaging --select F401,F811,F821,E722,E741
exit: 1
```

Salida Ruff exacta:

```text
F401 [*] `math` imported but unused
 --> tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3:8
  |
1 | from __future__ import annotations
2 |
3 | import math
  |        ^^^^
4 |
5 | from chemuson.clean2d import (
  |
help: Remove unused import: `math`
  |
2 |
  - import math
3 |
  |
Found 1 error.
[*] 1 fixable with `--fix`.
```

`pytest -q` completo **NO EJECUTADO**: lo prohíben la solicitud actual y la política proporcional vigente; se preservan todos los tests y gates. No se invocó `test.yml`.

## Auditoría de dependencias e imágenes

- Búsqueda inicial: sólo `pyproject.toml`, `requirements.txt`, el módulo `python3-pillow` del manifiesto Flatpak, una aserción de test y una descripción del manual nombraban Pillow/PIL; no hay import PIL en el runtime.
- QImage/clipboard/persistencia de anotaciones usan PyQt6. Con import blocker de PIL: `QImage` creó una imagen 16×16; `rdkit.Chem.Draw.rdMolDraw2D` importó correctamente.
- Con el mismo blocker: `PYTHONPATH=/tmp/chemuson-no-pillow-probe:src QT_QPA_PLATFORM=offscreen chem/bin/python -m pytest tests/test_image_clipboard_persistence.py -q` → `9 passed in 0.43s`.
- Consultas OSV API por paquete PyPI/Pillow: 12.0.0 → 77 registros, incluye CVE-2026-42311 y CVE-2026-54058; 12.3.0 → 39 registros, incluye CVE-2026-96637, CVE-2026-96643, CVE-2026-96644 y CVE-2026-96778; 13.0.0 → 0 registros a la fecha de consulta. Decisión: retirar dependencia no utilizada, no sustituir su wheel.

## Permisos observados antes de cambios

El manifiesto Flatpak concede `--share=network` y `--filesystem=home`. Las rutas directas encontradas para estado son `~/.chemuson` (autosave, crash/update logs/caches, name2structure cache, rollback/downloads) y `~/.config/Chemuson` (biblioteca de plantillas y fallback de configuración). Se usa `QFileDialog` para flujos de selección/guardado; retest portal sin permiso general de home queda pendiente.

## Baseline de follow-up: auditoría de paquetes congelados — 2026-10-23

- En `35f5888cd92f0048b0e2c15c97534952682d830d`, rama `release/v0.3.0-beta.1-security-prep`, `git status --short` estaba vacío antes de extender la exclusión. Focus baseline: `chem/bin/python -m pytest -q tests/test_rdkit_packaged_worker.py tests/test_chemname_packaged_smoke.py tests/test_release_workflow_contract.py` → **20 passed in 2.35s**. Pytest completo no se ejecutó.
- Se descargaron los artifacts del Build Preview #39899351520 (source SHA `35f5888c…`) e inspeccionaron los CArchive PyInstaller reales: Windows portable tenía 5 entradas `PIL/*`; AppImage tenía las mismas 5 extensiones más 13 bibliotecas `pillow.libs/*`. El install log de Windows confirma resolución transitiva Pillow 13.0.0 desde RDKit.
- El Flatpak se importó a un OSTree temporal y se inspeccionó el árbol `/app`: no existe el directorio top-level `site-packages/PIL` ni metadata `Pillow-*.dist-info`; sí existen stubs/recursos legacy `rdkit/sping/PIL`, que no son la distribución Pillow.
- El artifact preview previo pasa sus smokes funcionales, pero no satisface la ausencia de recursos Pillow en los paquetes congelados; se conserva como histórico. La solución de follow-up agrega exclusión PyInstaller y lectura fail-closed del CArchive.
