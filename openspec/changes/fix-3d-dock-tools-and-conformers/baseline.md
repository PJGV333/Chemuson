# Baseline — fix-3d-dock-tools-and-conformers

- Fecha local: 2026-10-16 21:15 MST (UTC−07:00).
- Rama: `feature/chemuson-3d`.
- HEAD: `875ede535fc5d8237b0e7cbfae54c52bfc39a43e` (`feat(viewport3d): expand molecular editing tools`).
- `git status --porcelain`: solo `?? PROMPT_LUNA_TEMP.txt`; prompt temporal no se añade a ningún commit.
- `./.venv/bin/python -m compileall -q src tests tools packaging` → exit 0.
- `./.venv/bin/pytest --collect-only -q` → exit 0; 2201 tests collected in 0.45s.
- `./.venv/bin/ruff check src tests tools packaging --select F401,F811,F821,E722,E741` → exit 1 por un F401 preexistente: `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`, `import math` sin uso. No relacionado con este cambio.
- `pytest -q` completo no se repite: ya se conoce SIGSEGV Qt/SIP de baseline y la instrucción activa pide priorizar regresiones focalizadas.
- RDKit y PyQt6 disponibles en `./.venv`.

## Diagnóstico reproducible previo a cambios

Se construyó `MolGraph` desde el isomeric SMILES auténtico de tetrandrina (PubChem CID 73078; la cadena usada y la diferencia respecto al documento local se registran en validation.md). El worker aislado devolvió 46 átomos, 52 enlaces y 46 posiciones documentales; dos centros R/S quedaron especificados/preservados; etapa completa ~0.33 s, método MMFF, energía 243.577994 kcal/mol, error nulo. Esto demuestra que el embedding para esa referencia funciona en el entorno; no demuestra la lectura del dibujo `.cmsn` del usuario ni su ruta UI.

## Logs completos del baseline

Los comandos con stdout/stderr exactos se guardaron fuera del repositorio para no introducir artefactos temporales:

- `/home/ccachyavgp/.hermes/cache/scratch/chemuson-repair-baseline-compile.log`
- `/home/ccachyavgp/.hermes/cache/scratch/chemuson-repair-baseline-collect.log`
- `/home/ccachyavgp/.hermes/cache/scratch/chemuson-repair-baseline-ruff.log`
