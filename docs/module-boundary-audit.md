# Auditoría de fronteras de módulos de dominio existentes

Fecha: 2026-09-19

Esta auditoría usa `architecture/modules.yml` como fuente de ownership. Los
módulos ya cohesivos se conservan sin extracción estructural:

| ID | Módulo | Paths auditados | Dependencias actuales / objetivo | Excepciones | Ciclos | Estado |
| --- | --- | --- | --- | --- | --- | --- |
| M05 | `geometry3d` | `src/chemuson/geometry3d/` | M00, M01 / M00, M01 | 0 | 0 | audited / no structural change required |
| M06 | `compchem` | `src/chemuson/compchem/exporters/` | M00, M05 / M00, M05 | 0 | 0 | audited / no structural change required |
| M07 | `spectroscopy` | `src/chemuson/spectroscopy/` | M00, M01 / M00, M01 | 0 | 0 | audited / no structural change required |
| M14 | `update` | `src/chemuson/update/` | ninguno / ninguno | 0 | 0 | audited / no structural change required |
| M16 | `name2structure` | `src/chemuson/name2structure/` | M00, M01 / M00, M01 | 0 | 0 | audited / no structural change required |
| M17 | `markush` | `src/chemuson/markush/` | M00 / M00 | 0 | 0 | audited / no structural change required |
| M18 | `version` | `src/chemuson/__init__.py`, `_version.py`, `version.py` | ninguno / ninguno | 0 | 0 | audited / no structural change required |
| M19 | `bootstrap` | `src/chemuson/__main__.py`, `src/chemuson/app/` | M18, M08, M15 / M18, M08, M15 | 0 | 0 | audited / no structural change required |
| M20 | `gui.editor2d.selection` | `src/chemuson/gui/editor2d/selection/` | ninguno / ninguno | 0 | 0 | audited / no structural change required |

## Audit decisions

- M05: audited / no structural change required
- M06: audited / no structural change required
- M07: audited / no structural change required
- M14: audited / no structural change required
- M16: audited / no structural change required
- M17: audited / no structural change required
- M18: audited / no structural change required
- M19: audited / no structural change required
- M20: audited / no structural change required

## APIs y ownership

- M05 conserva sus modelos, servicios 3D, backends, cache y exportación XYZ;
  no importa GUI.
- M06 conserva únicamente exportadores Gaussian/ORCA/NWChem y depende de M00 y
  M05; no importa GUI.
- M07 conserva predicción NMR/masa y registro de predictores; no importa GUI.
- M14 ya contiene política, provider, seguridad, rollback, portable, Windows y
  telemetría de actualización; no se crea un segundo módulo de update.
- M16 conserva conectores de resolución nombre→estructura y su fallback seguro;
  no importa GUI.
- M17 conserva los modelos y servicios Markush/polímeros; no requiere una
  frontera adicional.
- M18 posee la metadata de versión y mantiene `chemuson.__version__` y
  `get_app_version`.
- M19 conserva `main` como API pública y el composition root terminal.
- M20 conserva exclusivamente las cinco políticas deterministas de selección;
  el namespace `gui.editor2d` sigue disponible para hermanos futuros.

## Resultado

No se detectó ownership incorrecto, dependencia catalogada incorrectamente,
excepción temporal, ciclo ni documentación objetiva que requiera una
reorganización de estos módulos. Las siguientes fases deben auditar settings,
resilience y shell sin duplicar estas responsabilidades.
