# Validación — revisión del manual de usuario 0.3.0-beta.4

Cambio puramente documental: conforme a `docs/testing/AGENT_TESTING_POLICY.md` y al protocolo de `AGENTS.md`, se comprueban rutas, enlaces, coherencia interna y `git diff --check`. No se ejecuta `pytest` ni `compileall`, no se compilan instaladores y no se generan PDF.

## Entorno

- Rama: `docs/manual-usuario-beta4-campana` sobre `main` @ `bc39ef35c5b7baac4bbbc4a7e9a395528f747f99`.
- Referencia funcional del contenido: tag `v0.3.0-beta.4` (sin diferencias en `src/` respecto a `main`).
- Las capturas se generaron ejecutando la aplicación con `QT_QPA_PLATFORM=offscreen` y PyQt6 sobre Python 3.14.

## Baseline de microcampaña (antes de modificar)

- Rama `docs/manual-usuario-beta4-campana`, HEAD y PR #13 `4699bea0ccb4a8d296e9afdfa719262706548627`; `git status --short` vacío. Referencias remotas actualizadas con `git fetch origin --prune`.
- Linux, Python 3.14.8, PyQt6 disponible; tag `v0.3.0-beta.4` = `514b2013ba7406d3a03b6e613d0ad471f2684a39`; `git diff v0.3.0-beta.4 HEAD -- src` vacío.
- `python -m compileall -q src tests tools packaging`: éxito; `ruff check docs/manual-figures/capturar_figuras.py --select F401,F811,F821,E722,E741`: éxito; `openspec validate revise-user-manual-beta-4 --strict`: éxito; `git diff --check`: limpio.
- Comando documentado de captura con `CHEMUSON_OUT` temporal: salida 0 pero solo **29/37** PNG. Faltan 20–26 por selección de grupo y 33 por no solicitar el submenú; no se informó error. Evidencia local: `/tmp/chemuson-manual-baseline.log`.
- Los cinco hilos de Codex estaban abiertos; nombre del AppImage y actualización por plataforma ya se corrigieron en commits previos. No se ejecutó pytest masivo.

## Validación de la microcampaña

- **Alcance:** manual, generador documental, dos PNG y registros OpenSpec; ningún cambio en producción o dependencias. El HEAD base de esta validación es `4699bea0ccb4a8d296e9afdfa719262706548627`, con los cambios de microcampaña en el árbol de trabajo.
- **Cinco hallazgos:** privacidad y generador corregidos aquí; el ejercicio previo de 2-propanol se había sustituido por 1-propanol en `8f534285` pero aún omitía el cuarto átomo, corregido aquí; AppImage coincide con la lista de assets de la release y fue corregido en `89c93660`; actualizaciones por plataforma se corrigieron en `71e113f2`.
- **Interfaz:** prueba Qt offscreen con tres arrastres y selección del oxígeno mediante `set_active_element`: cuatro átomos C–C–C–O y tres enlaces simples; `chem/bin/python` (RDKit 2026.03.6) devolvió SMILES `CCCO`. El atajo `Ctrl+K` se contrastó en `src/chemuson/gui/actions/structure_actions.py`. Una primera prueba con Python del sistema obtuvo la misma conectividad pero el worker informó `rdkit_unavailable` porque ese intérprete no tiene RDKit; no se presenta como PASS.
- **Generador:** `QT_QPA_PLATFORM=offscreen CHEMUSON_OUT=/tmp/chemuson-verify-all python docs/manual-figures/capturar_figuras.py ventana menu flyout contenido contextual paneles` terminó con salida 0, **37/37 PNG**, sin `FALLO`. Los dos PNG modificados coinciden con la regeneración temporal. `fig flechas` genera siete PNG con salida 0; nombres de grupo o figura desconocidos devuelven 2. El script previo también reproducía SIGSEGV al finalizar `fig flechas` en este entorno; la nueva ventana de cierre evita esa destrucción tardía de Qt sin tocar la aplicación.
- **Revisión visual:** hojas de contacto de los 37 PNG y lectura a tamaño real de figuras 17 y 25: flechas completas; manchas TLC 0,25/0,55/0,80 y cinco bandas reales visibles. Las escenas se generan con objetos de beta.4; el gel es ilustrativo. No se reemplazaron otros PNG. Las figuras de paneles mantienen espacio vacío y texto pequeño inherentes a las capturas originales; conviene revisarlas en la maquetación institucional posterior.
- **Validación documental local:** 22 capítulos, 37 imágenes/leyendas consecutivas, 37 filas de inventario correspondientes, 82 enlaces revisados para rutas y anclas locales, ambos README enlazan al manual; `openspec validate revise-user-manual-beta-4 --strict`, Ruff focal y `git diff --check` sin errores. El nombre del AppImage se cotejó con los assets publicados de `v0.3.0-beta.4` mediante `gh api`.
- **Limitaciones editoriales:** las descripciones de varias filas de los dos índices de figuras heredados están truncadas; las leyendas completas junto a cada imagen son la referencia. No se inspeccionó el renderizado en escritorio/HiDPI ni se generó PDF. El propietario debe revisar y aprobar el Markdown antes del merge.

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

Historial de la primera revisión: se documentaron 31 capturas `OK` y ninguna `FALLO`, antes de la microcampaña. Aquella versión todavía omitía imágenes de contenido y del submenú Análisis. Las entradas duplicadas (`33-contextual-analisis.png` y una segunda copia de la paleta de comandos) se habían eliminado.

**Validación vigente tras la microcampaña:** la ejecución con los seis grupos en un directorio temporal genera **37 de 37 PNG** con código de salida 0 y sin `FALLO`; la variante `fig flechas` genera los siete archivos de contenido. Se cotejó cada imagen del inventario con las leyendas del manual. El índice de figuras del capítulo 22 y el inventario de `docs/manual-figures/README.md` se completaron editorialmente para eliminar descripciones truncadas, sin modificar PNG adicionales.

## Huecos declarados

- **Visor molecular 3D sin figura.** Se dibuja con OpenGL y no pudo capturarse de forma fiable en modo sin pantalla. El capítulo 15 lo dice explícitamente en lugar de publicar una imagen no auténtica.
- **Menús con submenús profundos** (p. ej. los presets de diagramas electrónicos) se describen en texto; no se capturó cada nivel para no saturar el manual.
- **Comprobación visual en escritorio real** (KDE/Wayland, HiDPI) no repetida: corresponde a evidencia de campañas anteriores de modernización de interfaz.
