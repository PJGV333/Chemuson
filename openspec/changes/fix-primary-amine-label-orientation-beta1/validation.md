# Validación — aminas ordinarias horizontales

## Gates del cambio

- Reproducción y post-fix con script diagnóstico (`evidence/diagnose_label_orientation.py`,
  Qt offscreen, posiciones reales de glifos, no solo texto):
  - antes: etilamina `above`, anilina `above`, dietilamina `below`,
    1,1-diamino-3-fenilpropan-2-ona `above`+`below` (defecto exacto reportado),
    acetamida `below` (contrato amida), urea `left`+`below`, amonio sin glifo H.
  - después: etilamina/anilina/dietilamina/ambos amino del dicetona **horizontales**
    (`right`, sin apilado, sin solapar enlaces); acetamida y urea **sin cambios**.
- Tests dirigidos de etiquetas/H/enlaces (subconjuntos `grep`-dirigidos por
  `implicit_hydrogen|atom_label|separated_nh|...` y
  `add_bond_item|remove_bond_item|AddBondCommand|bond_order`):
  `tests/test_amine_label_orientation.py tests/test_amide_urea_2d_continuity.py
  tests/test_functional_group_hydrogen_adjustment.py
  tests/test_aromatic_nitrogen_regression.py tests/test_implicit_h_substitution.py
  tests/test_drag_move_undo.py tests/test_selection_tightness.py` → **54 passed**;
  batería adicional de enlaces/modelo/3D (`test_core_model.py`,
  `test_initial_2d_drawing_geometry.py`, `test_sp3_bond_placement.py`, etc.) →
  **195 passed, 1 skipped**. En total 249 tests dirigidos en verde, sin ejecutar la
  suite completa (≈2.610 IDs) por instrucción explícita del propietario.
- `chem/bin/python -m compileall -q src tests` — PASS.
- Ruff focal (`F401,F811,F821,E722,E741`) sobre los tres módulos y los dos archivos de
  test tocados — PASS.
- `openspec validate fix-primary-amine-label-orientation-beta1 --strict` — PASS.
- `git diff --check` — PASS.
- Capturas comparativas claro/oscuro en `evidence/before`, `evidence/after_light`,
  `evidence/after_dark` (no equivalen a aceptación visual en el AppImage nativo; el
  propietario revisará visualmente antes del merge).

## Clasificación del baseline monolítico

Los 63 fallos de `baseline.md` son previos a este cambio y ajenos a él: los mismos
archivos (`test_aromatic_nitrogen_regression.py`,
`test_functional_group_hydrogen_adjustment.py`, `test_implicit_h_substitution.py`)
pasan de forma aislada en HEAD (`git stash`) y con el cambio aplicado — 15 passed en
HEAD aislado, 18 passed con el cambio. Es interferencia entre tests con Qt en la
corrida monolítica local (Python 3.14.8/rdkit 2026.03.6 vs CI 3.11, donde la suite
completa estaba verde). Ninguno de esos tests fue modificado para este cambio.

## Invariantes verificadas

El nuevo test afirma longitudes de enlace 40 ± 1e-7 en la reproducción completa
(con anillo aromático), conteos de H implícitos por `implicit_h_count`, ausencia de
glifo H para amonio cargado, y horizontalidad estable tras añadir enlace, undo, redo,
arrastre con undo y cambio de tema. El recorte de enlaces se ejercita indirectamente
en las aserciones de no-solape entre forma del glifo H y `shape()` de los enlaces.

## Validación de la Revisión 2

- Reproducción cuantitativa (`evidence/diagnose_chains.py`, Qt offscreen, carbonos
  implícitos como en el dibujo normal):
  - antes: zigzag simétrico con N central → `hpos=right`, trazos N–C 31.5 px vs
    22.1 px (trim 8.5 vs 17.9); poliamina lineal y N-metilanilina con el H del
    secundario sobre el sector de un enlace;
  - después: `hpos=above`/`below` según sector abierto; trazos N–C 31.5 px vs
    31.5 px (trim 8.5/8.5); terminales `H2N`/`NH2` sin cambios.
- Tests dirigidos: `tests/test_amine_label_orientation.py` (14 node IDs),
  `test_amide_urea_2d_continuity.py`, `test_functional_group_hydrogen_adjustment.py`,
  `test_aromatic_nitrogen_regression.py`, `test_implicit_h_substitution.py`,
  `test_drag_move_undo.py`, `test_selection_tightness.py` → **57 passed**;
  batería de geometría/enlaces/modelo (`test_initial_2d_drawing_geometry.py`,
  `test_sp3_bond_placement.py`, `test_core_model.py`,
  `test_multilayer_chemical_graph.py`) → **44 passed**. No se ejecutó la suite
  completa (≈2.6 k node IDs) por instrucción del propietario.
- Manifiesto CI: `update-expected --confirm` → 2624 node IDs;
  `pytest_file_durations.json` +1 archivo (254); contador fijado 253→254 en
  `tests/test_ci_pytest_shards.py` → 6 passed;
  `tools/ci_pytest_shards.py plan` → exit 0, 8 shards, distribución 2624.
- `compileall` src+tests PASS; Ruff focal sobre los archivos tocados PASS;
  `openspec validate fix-primary-amine-label-orientation-beta1 --strict` PASS;
  `git diff --check` limpio.
- Capturas `evidence/chain_light|chain_dark` (antes/después) de zigzag, V de
  dietales, segmento horizontal, poliamina y N-metilanilina. No sustituyen la
  revisión visual en el AppImage nativo.
