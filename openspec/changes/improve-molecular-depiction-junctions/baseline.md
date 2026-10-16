# Baseline — uniones de enlaces moleculares

## Identidad del repositorio

- Repositorio: `PJGV333/Chemuson`; remoto `origin=https://github.com/PJGV333/Chemuson.git`.
- Rama de trabajo: `feat/molecular-depiction-excellence`, publicada y con upstream en `origin/feat/molecular-depiction-excellence`.
- SHA inicial: `9a70ee1c908a39e700098bcc0a54d464a606e24c`.
- `origin/release/v0.3.0-beta.1-prep` y la rama local prep estaban exactamente en el mismo SHA inicial. La rama experimental se creó directamente desde ese ref; no se modificaron prep, `main`, tags ni releases.
- `git status --short` al iniciar: salida vacía (árbol limpio).
- Entorno de pruebas: Python 3.14.8, pytest 9.1.1, PyQt6 y RDKit; Qt offscreen.

## Captura de baseline (antes de cambios en producción)

| Comando | Resultado exacto relevante |
|---|---|
| `python -m compileall src tests tools packaging` | exit 0. Recorrió `src`, `tests`, `tools` y `packaging`; sin errores de compilación. |
| `pytest --collect-only -q` | exit 0; `2238 tests collected in 0.59s`. |
| `pytest -q tests/test_wedge_geometry.py` | exit 0; `11 passed in 0.22s`. |
| `ruff check src tests tools packaging --select F401,F811,F821,E722,E741` | exit 1 por un único hallazgo histórico fuera de alcance: `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3:8`, F401 `math` importado y no usado. No se corrigió. |
| `tools/ci_pytest_shards.py plan --shards 8 --expected tests/ci/expected_pytest_nodeids.txt` | exit 0; distribución `[280, 280, 280, 280, 280, 280, 279, 279]`; digest de asignación `d013b567a97a9f162c44598b632b417f066a004737d7d0161d8f32974bc281c4`. |
| Ocho ejecuciones `tools/ci_pytest_shards.py run --shard-index 0..7 --timeout-seconds 300` | Las ocho terminaron con exit 0 y `timed_out=false`; duraciones (s): `[95.288, 36.868, 85.252, 121.651, 66.383, 92.295, 48.441, 37.854]`. |
| `tools/ci_pytest_shards.py verify --shard-count 8 --plan-status success --shards-status success` | `pytest shard coverage verified: 2238 unique node IDs, 8 complete successful shards`. JUnit agregado: 2238 tests, 2218 passed, 20 skipped, 0 failures, 0 errors. |

La primera invocación de verificación no encontró el plan porque se guardó como `plan.json`, no con el nombre que requiere el verificador (`pytest-plan.json`). No se reejecutaron tests; se copió el plan idéntico al nombre contractual y `verify` pasó sobre los ocho artefactos existentes. `pytest -q` monolítico no se ejecutó: la validación completa de esta rama se hará por shards, como solicita esta campaña y exige el flujo CI acotado del repositorio.

## Referencia visual de base y diagnóstico inicial

- Capturas exploratorias renderizadas con `ChemusonCanvas`/Qt: `/tmp/depiction-before-light.png` y `/tmp/depiction-before-dark.png` (1680×1160). Los crops de inspección están bajo `/tmp/crop_*.png`.
- Galería exploratoria: 20 casos; simple/wedge, simple/hash, doble/wedge, doble/hash, bold/wedge, bold/hash, solid/hash, wedge/wedge, bifurcaciones, casi colineal, N/O/S/cargas, triple y escala, anillo aromático y fusión.
- Defecto visual confirmado para investigar: en algunas terminaciones hash que coinciden con un enlace simple/bold, el trazo transversal más ancho deja pequeños salientes fuera del trazo vecino. No se ha decidido aún si recortar, retraer o conservar el trazo terminal; eliminar una raya puede crear un hueco perceptible.
- Observación de implementación confirmada: el adaptador de unión de cuña sólida sólo recibe un ancho escalar por enlace vecino y no transmite explícitamente cada línea de un doble enlace. El doble renderizado separa sigma/pi y retranquea la línea pi en ciertos estilos. Que esta diferencia produzca un artefacto visual en el caso observado sigue siendo hipótesis, no se autoriza a conectar la pi al átomo.
- Cuñas en bifurcaciones/anillos pueden parecer densas, pero la forma observada no basta para calificarlas como defectos: se conservan como controles sin uniformar sus siluetas.
- El primer render exploratorio oscuro tenía un encabezado blanco sobre la hoja blanca por un error del generador temporal; no es un defecto del producto. El generador reproducible deberá corregir la composición del encabezado antes de guardar evidencia final.

## Alcance de la captura

No se modificó código de producción, Clean2D, modelo químico, catálogo de módulos, serialización ni tests existentes antes de completar esta captura. No se declara aceptación visual del producto por estas mediciones.

## Evidencia visual reproducible de la campaña

- Generador versionado: `tools/visual/render_bond_junction_gallery.py`; usa `ChemusonCanvas`/`BondItem`, soporta claro/oscuro, atom numbering real y composición de una lámina 2×2. Las capturas exploratorias se escriben en `/tmp` por defecto.
- Capturas base (antes del cambio de hash geometry): `/tmp/chemuson-bond-junction-baseline/bond-junctions-light.png` y `bond-junctions-dark.png` (1680×1160 cada una).
- Capturas posteriores: `/tmp/chemuson-bond-junction-after/bond-junctions-light.png` y `bond-junctions-dark.png` (1680×1160 cada una); numeración opcional inspeccionada en `/tmp/chemuson-bond-junction-numbered/`.
- Lámina única versionada: [`docs/visual/bond-junctions-before-after.png`](../../../docs/visual/bond-junctions-before-after.png), 3360×2408; columnas antes/después y filas tema claro/oscuro.
- La comparación del agente observa menor saliente terminal en hash→simple/bold y un perfil de marcas más simétrico. La separación al átomo se debe revisar en el tamaño/zoom de uso real; esta comparación no constituye aceptación visual del propietario.
