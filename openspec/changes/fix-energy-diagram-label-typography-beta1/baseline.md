# Baseline — tipografía de diagramas de energía

- Rama creada desde `release/v0.3.0-beta.1-prep` en `582fd0e2f80131a2fc9ae4e5a1a187c565bcb470`; `git fetch origin` confirmó el mismo SHA remoto. El checkout autorizado estaba limpio antes de crear esta carpeta OpenSpec.
- Worktrees inventariados: se usa solo `/home/unison-pjgv/Documentos/GitHub/Chemuson`, limpio y sin operaciones Git en curso. Los otros worktrees no se usaron; dos conservan archivos no rastreados ajenos a esta campaña (ver el informe final).
- La suite completa no se repitió: el OpenSpec de beta registra el timeout monolítico de 600 s y tres fallos Qt/SIP ya reproducidos en baseline. La solicitud permite omitir esa repetición; no se alteran esos fallos.

Salidas exactas (incluyendo códigos de salida) en `baseline/`:

| Comando | Resultado |
|---|---|
| `git status --short` | Solo la nueva carpeta OpenSpec sin seguimiento; ningún cambio de producto al inicio. |
| `python -m compileall src tests tools packaging` | PASS, exit 0. |
| `.venv/bin/pytest --collect-only -q` | 2.597 tests collected, exit 0. |
| `.venv/bin/pytest -q` monolítico | Omitido con autorización de alcance; razón y referencia al baseline conocido en `04-pytest.log`. |
| `.venv/bin/ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | Exit 1 por el único F401 preexistente `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py`; no corregido aquí. |
| Suite focal de presets, builders y `tests/test_ui_tool_rail.py` | 86 passed. |

## Evidencia anterior

- Capturas en `evidence/before/`: menú MO, campo de ligandos y paleta; temas claro/oscuro.
- Capturas tomadas con Qt 6 / QPA `xcb` en sesión KDE, pantalla `HDMI-A-1`, 1600×830, DPR 1.0; una segunda ejecución usa `QT_SCALE_FACTOR=1.5` (DPR 1.5; área disponible 1067×553). Se usó el menú real de `SymbolPaletteToolbar` y la captura Qt del popup, no una composición HTML. No equivale a revisión/aceptación manual del propietario.
- Fuente expuesta por Qt: `Sans Serif`; cobertura de glifos Unicode consultados (`₂`, `⁻`, `²`, `⁴`–`⁸`) reportada como presente. Las métricas y bounds están en `menu-metrics-*.json`.
- `evidence/semantic-preset-fingerprints.json` guarda SHA-256 del `to_json_dict()` canónico de los 25 presets de la rama base para comparar que la notación solo afecta la UI.
