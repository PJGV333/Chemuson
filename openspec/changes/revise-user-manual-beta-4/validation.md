# Validación — revisión del manual de usuario 0.3.0-beta.4

Cambio puramente documental: conforme a `docs/testing/AGENT_TESTING_POLICY.md` y al protocolo de `AGENTS.md`, se comprueban rutas, enlaces, coherencia interna y `git diff --check`. No se ejecuta `pytest` ni `compileall`, no se compilan instaladores y no se generan PDF.

## Entorno

- Rama: `docs/manual-usuario-beta4-campana` sobre `main` @ `bc39ef35c5b7baac4bbbc4a7e9a395528f747f99`.
- Referencia funcional del contenido: tag `v0.3.0-beta.4` (sin diferencias en `src/` respecto a `main`).
- Las capturas se generaron ejecutando la aplicación con `QT_QPA_PLATFORM=offscreen` y PyQt6 sobre Python 3.14.

## Comprobaciones y resultados

| Comprobación | Resultado |
|---|---|
| Equivalencia de código con la versión publicada: `git diff --stat v0.3.0-beta.4..main -- src` | Sin diferencias. El manual describe el código de `0.3.0-beta.4`. |
| `git diff --check` | Sin avisos de espacios finales ni conflictos de marcador. |
| Validador de enlaces e imágenes del manual | **37 imágenes, 37 leyendas `*Figura N.`**, **64 enlaces internos**. Sin rutas de archivo inexistentes, sin anclas rotas, sin imágenes sin leyenda y sin figuras huérfanas en disco. |
| Consecutividad de la numeración de figuras | 1…37 sin saltos ni repetidos. |
| `openspec validate revise-user-manual-beta-4 --strict` | `Change 'revise-user-manual-beta-4' is valid`. |
| Enlaces de entrada al manual | `README.md` y `README.en.md` mantienen 3 referencias cada uno a `docs/MANUAL_USUARIO.md`. |
| Manual histórico `docs/manual_usuario.md` | Sigue enlazando al manual principal; solo se corrigió su aviso de cabecera, que citaba `0.3.0-beta.3` al enlazar la release de `0.3.0-beta.4`. |
| Política de PDF | `find . -name '*.pdf'` no devuelve archivos en el repositorio y no se añadió, reemplazó ni subió ningún PDF. El manual enlaza Markdown. |
| Alcance del árbol modificado | Solo `docs/MANUAL_USUARIO.md`, `docs/manual_usuario.md`, `docs/manual-figures/` y `openspec/changes/revise-user-manual-beta-4/`. Sin cambios en `src/`, `tests/`, `packaging/`, `.github/`, `architecture/modules.yml` ni tags. |

## Pasadas cruzadas adicionales (segunda revisión)

| Comprobación | Resultado |
|---|---|
| Las **55 rutas de menú** mencionadas en el manual, contrastadas con el volcado de la interfaz (`menus.txt`, `rail.txt`, `shortcuts.txt`) | Todas resueltas; no queda ninguna ruta inventada. |
| Barrido inverso: **229 elementos** de la interfaz enumerados en los volcados y buscados en el manual | 45 sin mención literal, revisados uno a uno: variantes agrupadas (Anillo 4–12, glifos de orbital, tamaños de papel) o nombres ingleses ya citados. |
| Huecos reales detectados por ese barrido y **documentados en esta revisión** | Barra de formato de texto contextual (fuente, tamaño 6–144, B/I/U, sub/superíndice, alineación, color, `Op %`); **Color de etiquetas ▸ Por elemento / Negro**; las seis entradas de **Ver → Numeración**; **Grosor de enlace/flecha/corchete** y **Longitud de lineas** con sus etiquetas literales; **Atomic Presets**, **Diatomic MO Presets** y **Ligand Field Presets**. |
| Etiquetas literales del menú contextual verificadas en `canvas_context_menu.py` | **Incrementar (+1)** / **Disminuir (-1)**, **Diagrama de energia**, **Diagrama electrónico**, **Editar lóbulo ▸**, **Tamano** / **Pequena**, **Distribuir ligandos alrededor**. |
| Inventario de figuras de la sección 22 regenerado desde las leyendas del cuerpo | 37 filas frente a 37 figuras, **0 desajustes** de número, archivo o descripción. |
| Atajos del manual contrastados con el volcado de `QShortcut`/`QAction` | Coinciden, incluidas las excepciones declaradas (`Ctrl++` compartido con superíndice, `Ctrl+Q` sin atajo, `Ctrl+A`/`Ctrl+E` gestionados en el lienzo). |

## Ensayo de las capturas

`docs/manual-figures/capturar_figuras.py` vuelve a generar el conjunto completo desde el repositorio:

```
QT_QPA_PLATFORM=offscreen python docs/manual-figures/capturar_figuras.py ventana menu flyout contenido contextual paneles
```

Resultado: 31 capturas `OK` y ninguna `FALLO`. Las dos entradas duplicadas que produjo una primera versión del script (`33-contextual-analisis.png` y una segunda copia de la paleta de comandos) se eliminaron del script y del árbol.

## Huecos declarados

- **Visor molecular 3D sin figura.** Se dibuja con OpenGL y no pudo capturarse de forma fiable en modo sin pantalla. El capítulo 15 lo dice explícitamente en lugar de publicar una imagen no auténtica.
- **Menús con submenús profundos** (p. ej. los presets de diagramas electrónicos) se describen en texto; no se capturó cada nivel para no saturar el manual.
- **Comprobación visual en escritorio real** (KDE/Wayland, HiDPI) no repetida: corresponde a evidencia de campañas anteriores de modernización de interfaz.
