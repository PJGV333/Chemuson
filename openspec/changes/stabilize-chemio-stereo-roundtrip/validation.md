# Validación — ChemIO stereo round-trip

## Identidad y alcance

- Rama local: `fix/chemio-stereo-roundtrip`, creada directamente desde `release/v0.3.0-beta.1-prep` en `6aeef19028ffd047f8a39bc8f6063ea0b57210bf`.
- Cambio OpenSpec independiente: `stabilize-chemio-stereo-roundtrip`; commit documental inicial `a5acf76198088cac27b5b2826bfaef0407f19076`.
- Commit de implementación local: `8cf4cbc` (`fix(chemio): preserve stereo through round-trips`). No contiene cambios de ChemName, Clean2D de producción, GUI, Persistencia, dependencias, ni workflows.
- El alcance acepta sólo cambios en las tres rutas ChemIO y su matriz focalizada de tests, además del OpenSpec independiente.

### Archivos del cambio

- `src/chemuson/chemio/rdkit_io.py`
- `src/chemuson/chemio/rdkit_safe.py`
- `src/chemuson/chemio/_rdkit_worker.py`
- `tests/test_smiles_stereo_import.py`
- `openspec/changes/stabilize-chemio-stereo-roundtrip/.openspec.yaml`
- `openspec/changes/stabilize-chemio-stereo-roundtrip/README.md`
- `openspec/changes/stabilize-chemio-stereo-roundtrip/baseline.md`
- `openspec/changes/stabilize-chemio-stereo-roundtrip/design.md`
- `openspec/changes/stabilize-chemio-stereo-roundtrip/proposal.md`
- `openspec/changes/stabilize-chemio-stereo-roundtrip/specs/chemio-stereo-roundtrip/spec.md`
- `openspec/changes/stabilize-chemio-stereo-roundtrip/tasks.md`
- `openspec/changes/stabilize-chemio-stereo-roundtrip/validation.md`

## Decisiones y resultado

- El parser MOL usa pares de átomos normalizados sólo como claves de deduplicación; guarda el orden CTAB original para conservar la dirección wedge/hash.
- La conversión RDKit directa y el worker reconstruyen tags tetraédricos desde `stereo_cip` y dirección de enlace, y transfieren E/Z explícito usando prioridades CIP. El payload lleva coordenadas para MolBlock; dobles no asignados se marcan como desconocidos en vez de inferir configuración por el dibujo.
- El importador transfiere metadatos RDKit a los átomos/enlaces ChemIO. El exportador MOL fallback escribe códigos wedge/hash/either y las rutas que no pueden representar anotaciones explícitas fallan en lugar de degradarlas silenciosamente.
- La matriz compara SMILES canónico isomérico, fórmula, centros, E/Z, conectividad, orden/enlace, carga e isótopo con RDKit independiente. Incluye ambos P0, enantiómeros opuestos, ordenamientos de vecinos alternativos, centro multicentro, E/Z y FC=CF no asignado, MOL/SDF, tetrandrina y vancomicina.
- Para la SMILES histórica de tetrandrina, RDKit observa dos centros potenciales no asignados, cero centros especificados y ningún E/Z asignado. El control ya no exige cuñas artificiales.
- No se añadió dependencia. No se alteró código de producción Clean2D ni se afirma validación integral de Clean2D; el grupo incluye sólo regresiones focalizadas de consumidores ChemIO.

## Evidencia de validación (comandos limitados externamente)

| Verificación | Resultado |
| --- | --- |
| Cada P0, en proceso pytest independiente, con `timeout 120s` | `test_chiral_smiles_import_creates_wedge_or_hash`: 1 passed; `test_amino_acid_chiral_smiles_import_creates_wedge_or_hash`: 1 passed |
| Grupo focal ChemIO/RDKit/MOL y consumidores, `timeout 300s` | 87 passed, 9 skipped en 24.20s |
| `tests/architecture`, `timeout 180s` | 280 passed en 13.04s |
| `python -m compileall -q src tests tools packaging`, `timeout 120s` | PASS |
| `pytest --collect-only -q`, `timeout 180s` | 2128 tests collected en 0.68s |
| Ruff focal de los tres módulos ChemIO y test estéreo (`F401,F811,F821,E722,E741`) | PASS |
| Ruff global `src tests tools packaging` con las mismas reglas | Único fallo baseline: `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`, igual que antes del cambio; fuera de alcance y sin modificar |
| `openspec validate stabilize-chemio-stereo-roundtrip --strict` | válido |
| `git diff --check` | PASS |
| Suite completa `pytest -q` | NO ejecutada; la campaña limita la suite monolítica por duración/Qt |

Los 9 skips del grupo focal son los skips condicionales preexistentes del modo de importación directa RDKit; la ruta directa sí se ejercita en un subproceso aislado habilitando explícitamente `CHEMUSON_ENABLE_DIRECT_RDKIT=1` dentro de `test_rdkit_direct_and_isolated_conversion_paths_preserve_stereo`.

## Git y publicación

- Commit de implementación: `8cf4cbc` sobre el commit OpenSpec `a5acf76`; ambos son locales en `fix/chemio-stereo-roundtrip`.
- Se intentó únicamente el push normal previsto: `git push -u origin fix/chemio-stereo-roundtrip`. GitHub rechazó el intento antes de transmitir commits: `fatal: could not read Username for 'https://github.com': No existe el dispositivo o la dirección`.
- `gh auth status` confirma que no hay sesión iniciada en GitHub. No se solicitaron/guardaron credenciales, no se cambió la configuración global de Git y no se intentaron otros mecanismos de publicación.
- Estado histórico al cerrar la campaña original: el push separado quedaba pendiente de autenticación. No hubo merge, PR, tag, release ni publicación de artefactos desde esa rama.

## Publicación descendiente e integración beta autorizada

- El commit de implementación `8cf4cbc` quedó publicado como ancestro de `fix/ci-pytest-stabilization` en SHA `571e012ad45a4973e18d42d9d8943ae204cfd9b3`. El run Actions [#37999001943](https://github.com/PJGV333/Chemuson/actions/runs/37999001943) aprobó el plan, ocho shards y smokes Windows/Flatpak: **2.110 passed, 20 skipped, 0 failed**. Esto verifica la integración del código en una rama descendiente; no equivale al push independiente de `fix/chemio-stereo-roundtrip`.
- Por solicitud explícita de integración en `release/v0.3.0-beta.1-prep`, el cambio se incorporó mediante fast-forward-only desde dicho SHA. Los tests focalizados post-FF pasan (**12** tests ChemIO stereo) y la campaña beta registra refs, alcance y baselines.
- La aceptación manual de estereoquímica en los cuatro paquetes continúa `NOT TESTED`; un smoke CI y la comparación RDKit no la convierten en aprobación manual. No se han creado tags/releases ni publicado canales.
