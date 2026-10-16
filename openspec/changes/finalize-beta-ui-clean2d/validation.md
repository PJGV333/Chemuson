# Validación conjunta — UI + Clean2D

Estado: verificación local acotada; la aceptación manual del propietario sigue pendiente. No se ha declarado beta lista ni se han publicado artefactos oficiales.

## Procedencia y archivos

- Base beta prep: `9bc8192db3835262191ddb5247a8ff87f4d9169f`.
- UI source branch `fix/beta-save-shortcut-menu-access`: `47c588db286c785012da28fc6c70c9e106f7f872`; Clean2D source branch `clean2d/wedge-bond-campaign`: `a3b804fc20565fcc1041839335b904cf6b093e26`. Ambas ramas siguen apuntando a esos mismos SHAs.
- Integración de sólo lectura de fuentes: merge UI `932605ce369b9c476261a451c7b34fb65f9acacf`; merge Clean2D `81b36e44f93166dd9b58f9b3056f3dc287a536f2`.
- Commits de implementación del cierre: Clean2D `f7e703c` (`fix(clean2d): preserve wedge stereochemistry in GUI`); UI `f3cfc79` (`fix(ui): consolidate beta header actions`). La evidencia OpenSpec/manifiesto CI se registra en el commit documental posterior.
- Producción UI: `app_bar.py`, `actions/view_actions.py`, `command_registry.py`, `document_tabs.py`, `main_window.py`, `main_window_ui_builder.py`, `shell/assembly.py`, `theme/qss.py` y `platform/settings.py`.
- Producción Clean2D/estéreo: `core/layers.py` (conteo de centro desde origen de enlace dirigido) y `gui/controllers/clean2d_controller.py` (gate CIP aislado sólo para candidatos con cuñas/hash).
- Tests UI/settings: `test_beta_compact_shell_and_tabs.py`, `test_beta_save_shortcut_menu_access.py`, `test_platform_settings.py`, `test_ui_app_bar_tabs.py`, `test_ui_polish.py`, `test_ui_tool_rail.py`. Tests estéreo: `test_clean2d_wedge_stereo.py`. El manifiesto CI recoge la suma de los node IDs vigentes.
- Documentación: `docs/MANUAL_USUARIO.md` y este cambio OpenSpec (`proposal`, `design`, `spec`, `tasks`, baseline, validación y captura).

## Cabecera reproducible

- Evidencia: [`evidence/compact-header-light-1440x900.png`](evidence/compact-header-light-1440x900.png), captura offscreen 1440×900, PNG 30 KiB, SHA-256 `6814cbfb571eae8339f77b47fd6e8d1832dc83f01f21b1ce6a46b50a59ea2dba`.
- `tests/test_beta_compact_shell_and_tabs.py::test_default_header_has_one_compact_action_row_and_shared_file_actions` verifica los tamaños 1440, 980 y 900; menú tradicional visible, AppBar de 54 px, toolbar principal oculta, ausencia de marca/matraz/versión en la AppBar, botón `+` oculto en la ventana real, cuatro QAction compartidas y Pegar una sola vez en Editar. En las tres anchuras los controles de archivo caben en su grupo.
- La preferencia de toolbar migra a defaults versión 2: menús visibles y toolbar opcional oculta, incluso desde configuraciones con la toolbar anterior marcada visible. Los tests de tema, geometría dirty/X, accesibilidad, Ctrl+S y Save/Discard/Cancel permanecen activos.

## Diagnóstico Clean2D con estructura química importada

La referencia `C[C@H](F)CO` se importó por la ruta SMILES a un grafo de cinco átomos y cuatro enlaces, con un único enlace `2→1 HASHED/DOWN` y asignación RDKit independiente S. Antes del fix el perfil devolvía `stereo_center_count=2`, `multiple_stereo_centers`, preserve-only, y no generaba candidatos geométricos. Ahora cuenta un centro, no activa preserve-only y ejecuta los candidatos.

Resultado de `run_clean2d_engine` sobre grafo nuevo; átomos movidos son IDs respecto de la entrada. `rdkit_directo_deshabilitado:usar_rdkit_aislado` es el rechazo esperado de la ruta directa en el proceso principal; los demás candidatos seguros/no seguros no se fuerzan.

| Caso/selección | Modo | Perfil centros | Elegido / estado | Átomos movidos | Candidatos aceptados | Rechazos relevantes |
|---|---|---:|---|---|---|---|
| Referencia completa | quick | 1 | `current` / applied | 1,2,3,4,5 | `current` | `rdkit_direct`: backend-failure (`rdkit_directo_deshabilitado:usar_rdkit_aislado`); `internal_templates`: worse-quality (`candidato_no_canonicaliza_geometria`) |
| Referencia completa | publication | 1 | `current` / applied | 1,2,3,4,5 | `current` | Mismos dos rechazos que quick |
| Referencia completa | propose | 1 | `propose_reflection` / applied | 1,2,3,4,5 | reflexión, 8 rotaciones, proyección 3D, plantillas | 7 rotaciones `collision-risk` (`colision_no_enlazada`); `current` no es geometría alternativa |
| Desconectada + selección de sólo la referencia | quick/publication | 1 | `current` / applied | 1,2,3,4,5 | `current` | `rdkit_direct` deshabilitado; plantillas peores. Componente CCO intacto |
| Desconectada + selección de sólo la referencia | propose | 1 | `propose_reflection` / applied | 1,2,3,4,5 | reflexión, rotaciones, proyección 3D, plantillas | 7 rotaciones con colisión; `current` no alternativa |
| Wedge+hash con mismo origen, completa | quick/publication | 1 | `current` / applied | 1,3,4,5 | `current`, plantillas | `rdkit_direct` deshabilitado |
| Wedge+hash con mismo origen, completa | propose | 1 | `internal_templates` / applied | 1,3,4,5 | plantillas y alternativas | 8 rotaciones rechazadas por `invariant-violation` (`cambio_firma_estereoquimica`) |
| Wedge+hash, parcial `{1,2}` | quick/publication | 1 | `current` / applied | 1,2 | `current` | `rdkit_direct` deshabilitado |
| Wedge+hash, parcial `{1,2}` | propose | 1 | `propose_rotation` / applied | 1,2 | 2 rotaciones | `current` no alternativa; otras 2 rotaciones `cambio_firma_estereoquimica` |
| Centro real, parcial `{centro, vecino wedge}` | quick/publication | 1 | `current` / applied | sólo 1,2 | `current` | `rdkit_direct` deshabilitado |
| Centro real, parcial `{centro, vecino wedge}` | propose | 1 | `propose_rotation` / applied | sólo 1,2 | 4 rotaciones | `current` no alternativa |

El grafo combinado asigna IDs 1–5 a la referencia y 6–8 a CCO; ningún candidato contiene IDs fuera de la selección. La selección parcial conserva las coordenadas de los átomos exteriores y las cuñas.

## Verificación R/S independiente y guard GUI

- La prueba elimina `stereo_cip` cacheado de las copias de test, serializa Molfile y pide asignación en un proceso RDKit nuevo. `C[C@H](F)CO` se mantiene S y `C[C@@H](F)CO` se mantiene R y con SMILES isomérico distinto en quick, publication y propose.
- El control de cuña dibujada sin sustituyentes distinguibles no inventa R/S. Una molécula con dos centros reales sigue contando dos y conserva la política `multiple_stereo_centers → preserve_only`.
- Se reprodujo un caso no cubierto por la firma cíclica existente: tras recordar la primera propuesta segura (`propose_reflection`, S), el motor selecciona en la segunda propuesta una `propose_rotation` cuyo orden cíclico sigue idéntico, pero la asignación independiente pasa S→R. El gate GUI sobre la geometría completa, calculado por worker aislado sin CIP cacheado, rechaza ese candidato con `stable_rejection_reason=stereo-risk`, `expected_cip_by_atom={2:S}`, `observed_cip_by_atom={2:R}`. No se relajan filtros ni se aplica el candidato.
- Acciones de ventana con `C[C@H](F)CO` cargada: Ctrl+K/quick termina con `Estructura 2D alternativa propuesta`, mueve IDs 1–5 y conserva S; después, publication y propuesta no encuentran alternativa segura, mueven cero átomos y conservan S. No vuelven a informar el falso `multiple_stereo_centers`.

## Validación final

- CI local: `tools/ci_pytest_shards.py plan --shards 8` recolectó 2236 node IDs, distribución `[280,280,280,280,279,279,279,279]`, SHA-256 `ee85cc882d953b21fba5e8c93606c5ad804a3c6f5ab413ca34b2d5475391015e`. Se regeneró el manifiesto sólo por los node IDs de tests añadidos/reparametrizados; `test_test_workflow_contract.py` confirma cobertura exacta.
- Resultado final de ocho procesos bounded (≤270 s internos, timeout externo 290 s): **2216 passed, 20 skipped, 0 failed**, 2236/2236 IDs verificados por `tools/ci_pytest_shards.py verify`; ningún shard tuvo timeout ni excedió 137 s.
- En la primera pasada completa, shard 6 reveló que el nuevo gate inspeccionaba coordenadas de stubs sin cuñas en tests históricos. Se añadió el short-circuit por ausencia de enlace wedge/hash; `tests/test_clean2d_safety.py` pasó 47/47 y se repitieron los ocho shards, todos verdes.
- `pytest -q tests/architecture`: 280 passed. `pytest -q tests/test_test_workflow_contract.py`: 4 passed. Tests focalizados UI/Clean2D/settings/cabecera también pasaron.
- `openspec validate finalize-beta-ui-clean2d --strict --no-interactive`: válido.
- `python -m compileall src tests tools packaging`: pasó. Ruff focal (`F401,F811,F821,E722,E741`): permanece sólo el hallazgo baseline ajeno `F401 math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`; no hay hallazgos nuevos. `git diff --check`: pasó.
- Esta validación es local; la corrida remota asociada al push se registra por separado. La aceptación manual del propietario continúa pendiente.
