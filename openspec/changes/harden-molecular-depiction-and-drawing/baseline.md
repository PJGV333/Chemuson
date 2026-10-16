# Baseline — robustez de depiction y dibujo

## Identidad Git/entorno

- Rama: `feat/molecular-depiction-excellence`; HEAD/upstream inicial: `c38d954b9b7fb88a1aba04c06fc9376d5494ecc4`.
- `origin/feat/molecular-depiction-excellence` coincide con HEAD. `main` y `release/v0.3.0-beta.1-prep` se inspeccionaron y no se modificaron.
- `git status --short`: salida vacía (árbol limpio).
- Qt: `QT_QPA_PLATFORM=offscreen`; baseline realizado antes de cambios de producción/documentación de esta campaña. Logs detallados se conservaron temporalmente bajo `/tmp/chemuson-depiction-campaign-baseline/`.

## Comandos de baseline

| Comando | Resultado |
|---|---|
| `python -m compileall src tests tools packaging` | Exit 0; compiló todos los árboles solicitados, sin errores. |
| `pytest --collect-only -q` | Exit 0; `2267 tests collected in 0.71s`. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | Exit 1 por un hallazgo histórico único, no corregido: `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3:8`, `F401 math` sin uso. |
| `QT_QPA_PLATFORM=offscreen timeout 480 pytest -q` | Exit 124 por timeout externo a los 480 s; el progreso llegó aproximadamente al 34%; no se mostró fallo de aserción. El monolítico no se declara aprobado. |
| `tools/ci_pytest_shards.py plan --shards 8 --expected tests/ci/expected_pytest_nodeids.txt` | Exit 0; 2267 IDs; distribución `[284, 284, 284, 283, 283, 283, 283, 283]`; SHA256 de IDs `9670b28e58c90ad47fccb1539399f8c7317f11f9b04efa1c8604f950812772d4`. |

## Suite por shards y cobertura recuperada

Ejecuciones correctas: shard 0 `284 passed`; 1 `284 passed`; 2 `280 passed, 4 skipped`; 3 `278 passed, 5 skipped`; 5 `283 passed`; 6 `274 passed, 9 skipped`; 7 `282 passed, 1 skipped`.

El shard 4 expiró a los 300 s (SIGTERM, tras ejecutar hasta `tests/test_ui_tool_rail.py::test_group_for_tool_normalization[-None]`; sin error de aserción reportado). No se repitió el comando combinado. Se cubrieron sus 283 IDs en particiones: `tests/test_ui_tool_rail.py` 64 passed; 80 IDs Clean2D 80 passed; 139 IDs restantes 138 passed, 1 skipped. Así, los 2267 IDs quedaron cubiertos por ejecuciones que sí terminaron: **2247 passed, 20 skipped, 0 fallos de aserción**.

La cobertura completa por particiones no convierte el shard 4 oficial en success: el verificador oficial de ocho shards queda **no verificado** debido al timeout agregado. La suite completa monolítica también queda no verificada. No hubo cambios de código entre esas ejecuciones.

## Hallazgo Ruff preexistente

El único hallazgo global es el `F401 math` señalado arriba. Se conserva sin modificación; Ruff focalizado sobre futuros archivos cambiados se ejecutará después de cada fase.

## Alcance del baseline visual

La galería Qt de uniones existente (`tools/visual/render_bond_junction_gallery.py`) y la evidencia previa `docs/visual/bond-junctions-before-after.png` sirven como control de la fase anterior. Las capturas exploratorias nuevas se guardarán en `/tmp`; antes de editar geometría se actualizará una única lámina comparativa versionada.
