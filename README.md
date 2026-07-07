# ChemUSON

**ChemUSON** es un editor molecular 2D libre y de código abierto, orientado a la docencia, la investigación y la comunicación científica en química. El proyecto se desarrolla activamente con el propósito de ofrecer una herramienta accesible, moderna y útil para el dibujo, análisis y exportación de estructuras químicas en entornos académicos y de investigación.

**Autores:** Pedro Jancarlo Gomez Vega*, José Octavio Juárez Sánchez, Alondra Karolina Montaño Preciado, Ramón Moreno Corral, Felipe Medrano Valenzuela, David Octavio Corona Martínez y Karen L. Ochoa Lara*

**Afiliación:** Grupo de Química Supramolecular, Universidad de Sonora

\* Autores de correspondencia

## Estado del proyecto

> **Proyecto en construcción.**
>
> Chemuson está en desarrollo activo. Las funciones, formatos y resultados pueden cambiar mientras el proyecto madura.
>
> Estado de junio de 2026: la rama de trabajo usa SemVer de aplicación `0.3.0-dev` en `src/chemuson/_version.py` como fuente canónica. Esto prepara Workbench v0.3; no implica release publicada ni tag.

## ¿Qué hace Chemuson?

Chemuson permite dibujar estructuras químicas 2D, analizarlas y exportarlas en formatos comunes para trabajo académico y de investigación.

Capacidades actuales (resumen):
- Editor gráfico 2D con herramientas para átomos, enlaces, anillos, cadenas y anotaciones.
- Soporte de estilos de enlace (simple, doble, triple, aromático, cuña, punteado, coordinativo, etc.).
- Limpieza geométrica (`Limpiar 2D`) y opciones visuales (hidrógenos, carbonos, aromáticos como círculo, rejilla/reglas).
- Orientación inteligente de dobles enlaces en anillos (línea pi hacia dentro por defecto) con inversión manual mediante `Ctrl+Alt+Rueda`.
- Numeración visual no destructiva (átomos, estructuras o ambos), configurable para exportación.
- Importación y exportación química: `SMILES`, `Molfile` y archivos propios `.cmsn`.
- Autosave automático con backups rotativos por documento y recuperación de sesiones al iniciar.
- Exportación gráfica: `PNG`, `SVG`, `PDF`.
- Herramientas de análisis integradas: nombre, fórmula, masa exacta, peso molecular, `m/z` y análisis elemental.
- Dock de propiedades químicas con fórmula, masas, análisis elemental, issues de valencia y descriptores RDKit aislados (`logP`, TPSA, HBD/HBA, enlaces rotables y alertas Lipinski) sin bloquear la UI.
- Validación química interactiva con reporte exportable, navegación de issues y correcciones undoables cuando la acción es determinista.
- Dock **3D / CompChem** para generar conformeros, optimizar con RDKit aislado u Open Babel opcional, proyectar a 2D con confirmación/undo y exportar XYZ/inputs ORCA, Gaussian y NWChem.
- Dock de espectros MVP para predicciones heurísticas de ^1H NMR, ^13C NMR y MS con fuente, confianza, copia/exportación de tablas y selección pico↔átomo.
- Motor de nomenclatura **IUPAC-lite** para un subconjunto de estructuras orgánicas.
- Campo de estado en UI: **Nombre IUPAC** del documento activo (degradación a `N/D` sin crash).
- Preferencias de nomenclatura: `Nombre avanzado (fase 4/6)` y `Usar RDKit aislado`.
- Conversión Nombre→Estructura con fuentes offline/PubChem, caché y ejecución en worker.
- Biblioteca de plantillas (incluye base predefinida + categorías de usuario importables/exportables).

## Autosave y recuperación

- Chemuson guarda autosaves automáticamente cada ~2 minutos y también después de unos segundos de inactividad cuando hay cambios.
- Al iniciar, si existen autosaves pendientes, aparece un diálogo para **recuperar** o **descartar** cada sesión.
- Al guardar manualmente, se limpian autosaves obsoletos del documento y se conserva un respaldo rotativo reciente.

Rutas locales de trabajo:
- Autosaves: `~/.chemuson/autosave/`
- Autosaves archivados tras recuperar/descartar: `~/.chemuson/autosave/old/`
- Logs de crash: `~/.chemuson/crash_logs/crash_YYYYMMDD_HHMMSS.txt`

## Fortalezas

- Arquitectura modular (núcleo químico, GUI, nomenclatura, IO y persistencia separados).
- Integración con RDKit mediante workers/subprocesos para tareas pesadas o propensas a fallos nativos, además de rutas de respaldo internas.
- Cobertura de pruebas amplia para una base en evolución (800+ pruebas y skips controlados para backends opcionales).
- Orientado a uso docente y de laboratorio con enfoque práctico de dibujo y análisis rápido.

## Debilidades y límites actuales

- El motor de nombres es **IUPAC-lite**: no cubre toda la nomenclatura IUPAC ni todos los sistemas cíclicos/fusionados.
- Algunas rutas de cálculo/exportación usan subconjuntos de elementos o aproximaciones, según el caso.
- Las predicciones espectrales son heurísticas/estimadas y no sustituyen datos experimentales.
- La generación/optimización 3D depende de RDKit aislado; Open Babel solo se usa si el ejecutable `obabel` está disponible.
- No hay soporte CDXML, MRV ni PDB todavía.
- Al ser un proyecto en construcción, puede haber cambios de comportamiento entre versiones.
- No sustituye herramientas de validación química regulatoria o flujos de producción validados.

## Alcance de IUPAC-lite (actual)

Incluye soporte para:
- Cadenas lineales y algunas insaturaciones.
- Sustituyentes comunes (halógenos, alquilos simples, algunos grupos funcionales).
- Cicloalcanos simples, benceno y heteroaromáticos ampliados (incluye triazoles/tetrazol, oxazoles/isoxazoles, diazinas, oxazina).
- Algunos sistemas aromáticos fusionados (por ejemplo: naftaleno, antraceno, fenantreno, pireno; y casos concretos heterofusionados como benzotriazol).
- Carbonilos aromáticos tipo quinona (`benzene-1,4-dione`, `naphthalene-1,4-dione`, `naphthalene-1,2-dione`).
- Grupos funcionales adicionales (thiol, sulfoxide/sulfone/sulfonic acid/sulfonate, azido, peroxy).
- Cargas/isótopos/radicales básicos en modelo y round-trip (`formal_charge`, `isotope`, `radical_electrons`).
- Coordinación (experimental, MVP): carbonyl, ammine, aqua, halo, cyano, η5-cyclopentadienyl, con estado de oxidación y `cis/trans` básico.
- Plantillas especiales (MVP): carbohidratos seleccionados (`alpha/beta-d-glucopyranose`, `beta-d-fructofuranose`, `d-ribose`) y esteroides (`androstane`, `cholestane`), con sustitución simple (`hydroxy/oxo/amino`).
- Estereoquímica avanzada best-effort: `M/P`, `R_a/S_a`, `endo/exo`, `si/re` (si hay anotación confiable o soporte RDKit).

Limitaciones relevantes:
- Cobertura parcial en anillos fusionados complejos.
- Cobertura parcial en heterociclos fuera del subconjunto soportado.
- Casos con múltiples grupos funcionales y escenarios avanzados aún no completamente soportados.
- No cubre todo Blue Book; se prioriza estabilidad y degradación segura a `N/D` cuando procede.
- Para estereoquímica vía RDKit se recomienda modo aislado en subproceso (ver [docs/rdkit.md](docs/rdkit.md)).
- Descripción detallada del motor: [docs/chemname.md](docs/chemname.md).

## Stack técnico

- **GUI**: PyQt6
- **Química**: RDKit aislado para tareas de cálculo/conversión, Open Babel opcional para optimización 3D, motores internos para validación, Clean2D y heurísticas.
- **Lenguaje**: Python 3.10+

### Clean2D Quality Reporting

Clean2D exposes lightweight reporting helpers in `chemuson.clean2d.quality_reporting` for stable diagnostics in tests and controller/reporting surfaces. The normalized `score` in this contract is reporting-only; it must not change Clean2D layout algorithms, candidate ranking, candidate selection, geometry, RDKit/CoordGen behavior, MolGraph, or canvas behavior.

## Distribución por plataforma

Formatos soportados (estrategia híbrida):

- **Windows**
  - Portable: `Chemuson-vX.Y.Z-windows-x86_64-portable.exe`
  - Instalador: `Chemuson-vX.Y.Z-windows-x86_64-setup.exe`
- **Linux**
  - Portable: `Chemuson-vX.Y.Z-linux-x86_64.AppImage`
  - Instalable recomendado: **Flatpak** (canal oficial con remoto actualizable)

Convención de nombres en release:

- `Chemuson-v<version>-windows-x86_64-portable.exe`
- `Chemuson-v<version>-windows-x86_64-setup.exe`
- `Chemuson-v<version>-linux-x86_64.AppImage`

### Windows: instalador oficial + portable

- Se mantiene el `.exe` portable actual.
- Instalador oficial elegido: **Inno Setup** (balancea simplicidad operativa, soporte de modo silencioso para updater y convivencia con portable).
- Documento técnico y comandos completos: [docs/windows-distribution.md](docs/windows-distribution.md)

Build local rápido en Windows (PowerShell):

```powershell
python -m pip install --upgrade pip
pip install -r requirements.txt
pip install pyinstaller
pyinstaller --clean --noconfirm chemuson.spec
$env:CHEMUSON_VERSION = "<version>"
& "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" "packaging\windows\Chemuson.iss"
```

### Linux: Flatpak principal + AppImage portable

- Canal instalable principal: **Flatpak** (`.flatpak`).
- Canal portable: **AppImage** (`.AppImage`) sin instalación.
- Documento técnico y comandos completos: [docs/linux-distribution.md](docs/linux-distribution.md)

Instalacion rapida Flatpak estable (Arch/CachyOS, usuario final):

```bash
sudo pacman -S flatpak
flatpak --user remote-add --if-not-exists flathub https://dl.flathub.org/repo/flathub.flatpakrepo
flatpak install --user https://pjgv333.github.io/Chemuson/flatpak/stable/Chemuson-stable.flatpakref
flatpak run io.github.PJGV333.Chemuson
```

Fallback con `.flatpakrepo`:

```bash
flatpak remote-add --user --if-not-exists --from chemuson-stable https://pjgv333.github.io/Chemuson/flatpak/stable/Chemuson-stable.flatpakrepo
flatpak install --user chemuson-stable io.github.PJGV333.Chemuson//stable
flatpak run io.github.PJGV333.Chemuson
```

Canal beta:

```bash
flatpak install --user https://pjgv333.github.io/Chemuson/flatpak/beta/Chemuson-beta.flatpakref
```

Notas:
- El indice real de canales publicados vive en `https://pjgv333.github.io/Chemuson/` y solo enlaza canales que ya existen en `gh-pages`.
- Si `stable` todavia no aparece ahi, el remoto estable aun no ha sido sembrado; usa `beta` temporalmente o publica primero una release estable.
- La ruta recomendada es `.flatpakref` porque registra el remoto oficial para que `flatpak update` reciba futuras versiones.

Ejecución AppImage:

```bash
chmod +x Chemuson-vX.Y.Z-linux-x86_64.AppImage
./Chemuson-vX.Y.Z-linux-x86_64.AppImage
```

## Actualización automática (MVP)

- Se añadió núcleo de update en `src/chemuson/update/`.
- Fuente de versión única: `src/chemuson/_version.py`.
- Política de versión: SemVer de app. Las ramas de desarrollo pueden usar sufijo prerelease/dev como `0.3.0-dev`; releases futuras deben sincronizar CI, empaquetado y updater con esa fuente canónica antes de publicar.
- Canales soportados: `stable` y `beta`.
- Verificación de integridad:
  - hash SHA-256 (`.sha256`),
  - firma (`.sig`, MVP con HMAC-SHA256; preparado para ampliar a Ed25519).
- Rollback básico: si falla reemplazo de binario, se restaura backup local.
- En Windows instalado, el updater prioriza el asset `setup.exe` y lo aplica al cierre de la app.
- En Flatpak se deshabilita update in-app y se delega al remoto oficial publicado para `flatpak update`.
- Fallback seguro: si GitHub no responde, se usa caché local reciente de releases (si existe).
- Telemetría local mínima del updater (sin datos sensibles): `~/.chemuson/update_logs/events.jsonl`.

Preferencias persistentes de actualización (en `QSettings`):

- `update/enabled`
- `update/channel`
- `update/mode` (`notify` / `silent`)
- `update/check_interval_hours`
- `update/last_check_iso`

En GUI: `Editar -> Preferencias -> Actualizaciones`.
Chequeo manual: `Ayuda -> Buscar actualizaciones...`.

En CLI:

```bash
chemuson --version
```

## Migración sin romper portable actual

- El flujo portable actual **se mantiene**.
- Puedes seguir abriendo Chemuson con el ejecutable/AppImage sin instalar.
- Si migras a instalador en Windows, conserva tus archivos de trabajo (`.cmsn`) y configuración local.
- Los mecanismos nuevos de update no eliminan compatibilidad con releases portables existentes.

## CI/CD Windows (resumen)

- `test.yml` incluye smoke tests de updater/build de instalador en Windows.
- `release.yml` publica portable + setup y deja firma Authenticode preparada con secretos:
  - `WINDOWS_CODESIGN_CERT_BASE64`
  - `WINDOWS_CODESIGN_CERT_PASSWORD`
  - `WINDOWS_CODESIGN_TIMESTAMP_URL` (opcional)
- Limitación conocida: CI no ejecuta instalación real sobre el runner; se simulan pasos críticos y se valida compilación de setup.

## CI/CD Linux (resumen)

- `release.yml` construye y publica:
  - `Chemuson-vX.Y.Z-linux-x86_64.flatpak`
  - `Chemuson-vX.Y.Z-linux-x86_64.AppImage`
- `test.yml` valida manifiesto Flatpak con job `flatpak-smoke`.
- Para AppImage se publican sidecars de metadata de update (`.updateinfo` y `.update.json`) y `.zsync` cuando esté disponible.

## Instalación para desarrollo

```bash
python3 -m venv chemuson
./chemuson/bin/pip install -r requirements.txt
```

Opcional (modo editable):

```bash
./chemuson/bin/pip install -e .
```

## Ejecución

Desde el entorno virtual, con el paquete instalado en modo editable:

```bash
./chemuson/bin/chemuson
```
