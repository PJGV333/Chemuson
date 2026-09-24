# Fondos de referencia — fase 1 (fundación de temas)

Capturas offscreen (`QT_QPA_PLATFORM=offscreen`) de la **ventana real**
(`ChemusonWindow`) con la fundación de temas aplicada, para comparación
visual contra el spike aprobado
(`docs/ui-modernization/pyqt6-spike/`, commit `59e977d`).

- `foundation-light.png` — tema `light`
- `foundation-dark.png` — tema `dark`

Generadas durante la validación del OpenSpec
`2026-09-24-modernize-ui-theme-foundation` (Fase 1 del PLAN de
modernización de la UI).

Para regenerar (entorno efímero uv sobre el venv del proyecto, sin
instalaciones nuevas):

```bash
QT_QPA_PLATFORM=offscreen HOME=/tmp/chemuson-shots-home UV_CACHE_DIR=~/.cache/uv \
  uv run --no-project --offline \
    --python /home/unison-pjgv/Documentos/GitHub/Chemuson/.venv/bin/python \
    --with pytest --with ruff --with PyQt6 --with numpy --with Pillow \
    --with rdkit --with certifi --with PyYAML \
    -- python /tmp/chemuson_shots.py docs/ui-modernization/foundation-shots
```

(`chemuson_shots.py` crea la ventana, aplica `light`/`dark` vía
`ChemusonWindow._apply_theme` y captura `grab()`.)

Observaciones (conservador por diseño, Fase 1):
- La estructura de la ventana NO cambia (app bar/rail/side panel son fases
  posteriores); lo que cambia es el sistema de colores/superficies/bordes/
  botones/tipografía, ahora gobernado por tokens
  (`src/chemuson/gui/theme/tokens.py`).
- El canvas permanece como hoja blanca en ambos temas (consistente con la
  exportación PNG; ver PLAN.md Fase 7.3).
- Los iconos siguen siendo los de `gui/icons.py` (QPainter); su migración a
  SVG es la Fase 2 (el `IconProvider` ya existe en
  `src/chemuson/gui/theme/icon_provider.py`).
