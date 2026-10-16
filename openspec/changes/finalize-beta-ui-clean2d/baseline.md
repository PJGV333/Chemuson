# Baseline — consolidación UI + Clean2D

Fecha: 2026-10-15. Captura anterior a los cambios de producción de esta tarea.

## Refs, ramas y árbol

- `origin/release/v0.3.0-beta.1-prep`: `9bc8192db3835262191ddb5247a8ff87f4d9169f`.
- `origin/fix/beta-save-shortcut-menu-access`: `47c588db286c785012da28fc6c70c9e106f7f872`.
- `origin/clean2d/wedge-bond-campaign`: `a3b804fc20565fcc1041839335b904cf6b093e26`.
- Después de `git fetch origin --prune`, los HEAD locales coincidían con los tres HEAD remotos y `git status --short` estaba vacío.
- Las dos ramas de campaña comparten como base `9bc8192`; UI aporta 6 commits y Clean2D 3. Las ramas originales no se modifican.
- Rama nueva `fix/beta-final-ui-clean2d`: merge UI `932605ce369b9c476261a451c7b34fb65f9acacf`; merge Clean2D `81b36e44f93166dd9b58f9b3056f3dc287a536f2`.

## Captura de herramientas y tests

Antes del merge, en `clean2d/wedge-bond-campaign`:

- `python -m compileall src tests tools packaging`: **exit 0**.
- `pytest --collect-only -q`: **2174 tests collected in 0.86s**, exit 0.
- `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: **exit 1**, único hallazgo baseline conocido `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no relacionado.
- No se ejecutó `pytest -q` monolítico: la solicitud exige suite por shards. Los informes fuente registran las suites por shards ya pasadas; CI de origen Clean2D `38819097932` fue success y UI `38769058477` fue success en sus respectivos SHAs.

Tras ambos merges, antes de editar producto:

- `pytest --collect-only -q`: **2220 tests collected in 0.79s**.
- IDs colectados únicos: **2220**; manifiesto `tests/ci/expected_pytest_nodeids.txt`: **2220**, sin IDs ausentes, extra ni duplicados. Las cuatro IDs de tests shell reemplazados por la fase A no eran pruebas vigentes; sus escenarios están cubiertos por las versiones actualizadas. Las 20 pruebas nuevas Clean2D y las 45 nuevas de UI están registradas.
- Ejecución combinada `pytest -q tests/test_beta_compact_shell_and_tabs.py tests/test_beta_save_shortcut_menu_access.py tests/test_clean2d_wedge_stereo.py` (timeout externo 300 s): **65 passed in 73.28s**.
- `git merge-tree` reportó merges limpios; Git auto-fusionó el manifiesto. La comprobación contra la colecta real fue satisfactoria; no se actualizó el manifiesto para silenciar fallos.

## Sonda de aplicación real antes de corregir

Se cargó `C[C@H](F)CO` por la ruta de representación SMILES de ChemUSON en una ventana real `ChemusonWindow`; se transfirió el grafo importado al canvas y se activó la QAction existente de limpieza completa (atajo `Ctrl+K`). La primera sonda cerró el proceso con `QThread: Destroyed while thread '' is still running` (exit 134): el `QTimer.singleShot(1200, _maybe_check_updates_startup)` pudo disparar el worker de actualización mientras una sonda offscreen cerraba inmediatamente. Se aisló la sonda deshabilitando sólo ese callback externo a la prueba y cerrando la ventana normalmente; no se cambió producción.

Salida reproducible de la sonda aislada:

```text
UI {'menubar_visible': True, 'main_toolbar_visible': True,
    'main_toolbar_actions': ['Nuevo', 'Abrir...', 'Guardar'],
    'appbar_brand': 'ChemUSON', 'appbar_version': '0.3.0-beta.1',
    'appbar_flask_icon_null': False, 'side_panel_hidden': True}
CLEAN2D {'smiles': 'C[C@H](F)CO', 'atom_count': 5,
         'bonds': [(2, 1, 'hashed', 'down')], 'moved_atoms': [],
         'status': 'Estructura 2D compleja preservada; no se redibujó globalmente'}
PROBE_COMPLETE
EXIT=0
```

La misma estructura en la sonda de motor produjo `stereo_center_count=2`, `preserve_only=true`, razón `multiple_stereo_centers`; quick/publication devolvieron sólo candidato actual y no movieron coordenadas. El controlador retornó antes de generar candidatos geométricos. Esto no procede del backend, del filtro de calidad ni de la firma cíclica de Clean2D; la representación de un bond estereoquímico estaba creando un bloque por cada extremo. La versión del menú/toolbar UI sí coexiste en el proceso, pero el nombre/matraz/versión y la tercera fila visible contradicen la nueva aceptación del propietario.

La estructura exacta del archivo privado del propietario no está incluida en el repositorio ni en el mensaje. `C[C@H](F)CO` se usa como referencia química real ya registrada en la aceptación beta; la aceptación de sus archivos originales queda pendiente.
