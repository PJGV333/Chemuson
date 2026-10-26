# Diseño: beta.4 con identidad visual del proyecto y preview no publicador

## Base y aislamiento

- Base: `main`/`origin/main` en `9c171a8267d632e9bd914392680c316c28fe5e87`; `git fetch origin main` confirmó divergencia `0 0`.
- Rama: `release/v0.3.0-beta.4-prep`.
- `v0.3.0-beta.1`, `.2` y `.3` son tags anotados existentes y se preservan. No crear `v0.3.0-beta.4`.
- `gh auth status` indica que no hay sesión autenticada. La escritura remota, PR y Actions oficial quedan condicionados a credenciales accesibles; no se simulan.

## Fuente visual y variantes

- El PNG aprobado es una imagen RGB 1254×1254, PNG íntegro, SHA-256 `aab2600d2d46d87dba9dcea4266ae89ee4bbceafef0b05892b10713c310e29e6`.
- Se conserva como maestro byte a byte en `assets/branding/logo_ChemUSON.png`.
- Las variantes raster se generan por remuestreo de alta calidad sin recorte, recoloración ni cambios de proporción. Un `.ico` contiene tamaños de Windows; el SVG de distribución es un contenedor XML con el PNG original codificado en base64, no una vectorización.
- El paquete PyInstaller incluye la variante de ventana en su CArchive y usa el `.ico` como recurso PE en Windows. La aplicación carga el PNG desde la raíz de recursos de PyInstaller o desde el checkout fuente. Inno Setup utiliza el mismo `.ico` para el instalador.
- El SVG de marca sustituye el icono previo de AppImage/Flatpak, se instala bajo el ID estable existente y se declara como icono stock en AppStream. Los IDs, nombres, remotos y canales no cambian.
- README y manual presentan el emblema como identidad del proyecto ChemUSON; no como escudo oficial de Universidad de Sonora. El contraste claro/oscuro se verifica mediante rasterizado Qt de la imagen original opaca en ambos estilos y queda también sujeto a inspección visual del propietario.

## Empaquetado y limpieza

- El baseline clasifica `dist/`, `dist-appimage/`, `dist-appimage-final/`, `dist-preview/`, `dist-preview-final/` como paquetes antiguos ignorados y `build/` como intermedios PyInstaller; suman 2,499,985,408 bytes de espacio asignado. `dist-appimage-diagnostic/` contiene dos manifiestos trackeados y se conserva completo.
- No usar `git clean -fdX`: el dry-run también lista `.venv`, cachés y pycache que no se deben borrar en esta limpieza.
- Los scripts Linux existentes ya respetan `CHEMUSON_OUTPUT_DIR`. `packaging/windows/Chemuson.iss` será el único ajuste local si el test del preprocesador confirma que su `OutputDir` fijo ignora esa variable. Sin variable, se mantiene `dist-installer`, que es el contrato de los workflows Windows.
- No cambiar los comandos o rutas de los workflows oficiales por una preferencia de compilación local.

## Versión, documentación y actualizadores

- `_version.py` sigue siendo la fuente canónica. `set_version.py` actualiza versión/entrada AppStream para `0.3.0-beta.4`; las entradas históricas se conservan.
- Las notas beta.4 son bilingües y separan mantenimiento/identidad de funciones ya presentes en beta.3. Los README mantienen enlaces a la última Release que realmente existe; la vista previa no se presenta como descarga oficial.
- La API pública de GitHub confirmó en 2026-10-26 que `v0.3.0-beta.3` tiene `prerelease=false` y aparece en `/releases/latest`. Documentar como decisión excepcional del propietario. El workflow oficial todavía deriva `prerelease=true` para versiones del canal beta; no se modifica.
- El updater integrado usa SemVer de `tag_name`, no orden de fecha/título/flag, y `channel_accepts_version` excluye prereleases por versión en canal stable. La prueba focal confirma selección por canal y orden `beta.N`.
- El AppImage beta.3 publicado conserva `gh-releases-zsync|PJGV333|Chemuson|prerelease|...zsync`. La inspección estática del upstream AppImageUpdate (`AppImageCommunity/AppImageUpdate`, commit `a211784dfc746fdb6d8d32d6bb39add451c1ddeb`) demuestra que `prerelease` no es un selector: se interpreta como tag literal y consulta `/releases/tags/prerelease`. Los selectores especiales son `latest` (`/releases/latest`), `latest-pre` (primer prerelease real) y `latest-all` (primera release de cualquier tipo). La beta usa `latest-all` para permitir beta/RC y posterior graduación a estable según la política existente; stable mantiene `latest`. Beta.3 es inmutable y seguirá requiriendo actualización manual; beta.4 corrige el metadato futuro. No se cambian preferencias del usuario. `appimageupdatetool` no está instalado, por lo que no se afirma prueba dinámica del ejecutable.
- Flatpak actualiza por el remoto/rama instalado (`flatpak update`), no por el updater integrado; beta.4 conserva la rama beta y el bundle preview sigue aislado.
- Proteger el downgrade y exigir checksum/provenance se mantiene sin cambios.

## Revisión del gate de arranque no Qt

El propietario reportó que el Build Preview #40371397035, desde SHA `32ea9f05045a771874ec5e525e301bf855e8e9f6`, falló `tests/test_application_startup_non_qt.py::test_startup_success_order`: el test extrae `run_app()` mediante AST y su namespace simulado no incluía `_application_icon_path`/`QIcon`; la captura secundaria falló también porque faltaba `crash_reporter.write_crash_log`. La API pública del run devolvió HTTP 504 en el entorno, por lo que el diagnóstico del run se registra como aportado por el propietario; el nodo sí se reprodujo localmente en el SHA base.

La corrección es exclusivamente del harness: inyectar un path temporal existente, fake `QIcon` no nulo, y mocks de `setWindowIcon` que comprueban el mismo icon object en QApplication y ventana. Preservar la lista exacta de eventos de arranque, añadir caso de recurso ausente (`FileNotFoundError` → crash log → QMessageBox, sin cargar icono/crear ventana) y mantener intactos los tests de errores genéricos de arranque. Al agregar un node ID nuevo, regenerar y revisar explícitamente `tests/ci/expected_pytest_nodeids.txt` y verificar la colección completa; no ejecutar la suite completa.

## Validación y publicación

- Tests focales existentes para branding, iconos, packaging, release/version y updater; compileall, Ruff focal, arquitectura, contratos estrictos OpenSpec y `git diff --check`. Colección/manifiesto node IDs sólo si la lista de tests cambia.
- Validar recursos en source y PyInstaller, AppStream, scripts, metadatos y smoke esencial de apertura/editor 2D/flujo 2D→3D en los entornos disponibles. Las pruebas visuales de Windows y temas en paquetes exactos corresponden al propietario.
- Build Preview oficial requerido: Windows portable, Windows installer, Linux AppImage Type 2 y Flatpak, mismo SHA con checksums/provenance. No llamar “preview exitoso” sin run oficial completo y revisar los cuatro artifacts.
- Detener antes del tag. No hay publicación accidental ni modificación de distribución.

## Reparación de la preparación del entorno Flatpak (bloqueo de Build Preview)

Build Preview #40375931525 (SHA `670a7670235ec1f8833cb03c7d2475bbf5a2589f`) canceló
`Linux Flatpak bundle` por `timeout-minutes: 45`. La API pública de Actions registró los
tiempos por paso: `Install Flatpak build tools` 43 min 59 s y `Build local Flatpak bundle
only` sólo 2 min 24 s antes de la cancelación. El mismo par de pasos en los runs correctos
de beta.1–beta.3 tardó 0,3–0,7 min de preparación y 3,7–5,0 min de empaquetado, y en el
run inmediatamente anterior (#40371397035) la preparación ya había subido a 37 min 24 s
con el job terminando en 41 min 28 s. La conclusión es que el empaquetado no cambió; la
preparación del anfitrión sí.

Causa raíz externa, no de ChemUSON: el cierre `flatpak flatpak-builder ostree` baja ~99 MB
en ~123 paquetes. En los runners de GitHub el origen apt se resuelve por
`mirror+file:/etc/apt/apt-mirrors.txt` con `azure.archive.ubuntu.com` en prioridad 1; si el
nodo asignado estrangula el enlace, apt sigue en ese origen y la descarga avanza a decenas
de kB/s sin error. Está documentado upstream (actions/runner-images #14823, #14594,
#12949) y ya ocurrió en el job Flatpak de `test.yml` (123 paquetes/99 MB; ~38 min). El
`apt-get` del job AppImage del mismo run tardó 22 s, lo que confirma que no fue una caída
general de apt sino el enlazado lento del origen asignado a ese runner.

Decisión de diseño (respeta el orden de preferencia del propietario):

1. Entorno determinista: los jobs Flatpak quedan en `runs-on: ubuntu-24.04` en lugar de
   `ubuntu-latest`. Hoy apuntan a la misma imagen, pero la preparación queda anclada a la
   base sobre la que está medida; una migración futura de `-latest` no cambia la release.
2. Herramientas: se conserva apt como única fuente, con los paquetes del cierre y sus
   dependencias funcionales intactas. No se usa una imagen externa ni acciones de terceros
   de procedencia dudosa; la única acción nueva es `actions/cache`, oficial de GitHub, fijada
   al commit `55cc8345863c7cc4c66a329aec7e433d2d1c52a9` (etiqueta `v6.1.0`).
3. Caché segura: `packaging/flatpak/toolchain.lock` fija paquetes, versiones mínimas,
   mirrors oficiales permitidos y límites; su hash alimenta la clave de la caché de `.deb`.
   La instalación posterior usa `--no-download`, de modo que apt verifica cada archivo
   contra el índice firmado antes de instalarlo; un archivo obsoleto o alterado se vuelve a
   descargar. No se cachan índices ni repositorios de terceros.
4. Origen y reintentos: `packaging/linux/install_flatpak_tools.sh` sondea los cuatro
   orígenes oficiales de Ubuntu y reordena el `mirror+file` del runner (o la línea `URIs:`
   del stanza principal, preservando `Signed-By` y el stanza `noble-security`), con rebote
   al siguiente origen si el elegido falla. Se fija que la reescritura del `mirror+file`
   tenga efecto real: si `/etc/apt/apt-mirrors.txt` del runner no es escribible, el stanza
   principal se apunta directamente al origen oficial sondeado y, si tampoco se puede
   escribir la lista de fuentes, el job falla cerrado en vez de quedarse en el nodo asignado.
   Se fijan `Acquire::Retries`, timeouts HTTP/HTTPS
   y `ForceIPv4`, y cada fase corre bajo `timeout`.
5. Timeouts: no se aumenta ningún `timeout-minutes` (45 min en los jobs de empaquetado,
   25 min en la publicación Flatpak). El control queda dentro de las fases: 300 s para
   `apt-get update`, 600 s para la descarga acotada, 300 s para la instalación offline y,
   en `build_flatpak.sh`, 1500/300/300 s para builder/repo/bundle con salida 124 traducida a
   un diagnóstico explícito. Todas las fases del preparador comparten además un plazo global
   de 900 s (`prep_deadline_seconds`): cada límite efectivo es el mínimo entre el de la fase
   y lo que resta del plazo, así que la preparación no puede superar la mitad del presupuesto
   del job ni aunque los dos rebotes tarden. Si aun así se excede el límite, el job falla
   cerrado indicando la fase y el reporte del preview sigue marcando el Flatpak como no exitoso.
   Se fija además `DPkg::Lock::Timeout` (60 s) para el caso real de un runner actualizando
   en segundo plano, y el resumen publica `cached .deb archives published for reuse: N/M`, con
   aviso cuando ningún `.deb` pudo publicarse, para que la evidencia de la caché nunca sea
   optimista.

Riesgos de regresión y mitigación: el script escribe en `/etc/apt` del runner efímero, no en
el repositorio; se probó contra un arnés local con apt/curl/dpkg simulados en once escenarios
(espejos sanos, nodo estrangulado, rebote de update, rebote de descarga, caché cálida,
descarga fallida, paquetes de publicación, lista de mirrors no escribible, fallo cerrado,
plazo global agotado, update colgado cortado a los 2 s y caché de `.deb` no escribible) y contra los cuatro orígenes
reales para verificar la forma de la sonda. `build_flatpak.sh` sólo gana guardas de tiempo y
registro de fases: el manifiesto, los `sha256` pinneados, el smoke ChemName empaquetado, el
branch `preview-<sha>`, el nombre del bundle y las validaciones OSTree/GPG no cambian. Los
contratos estáticos nuevos (`tests/test_flatpak_toolchain_preparation.py`, 13 node IDs)
verifican que no se debilita la firma, que las fases están acotadas, que el empaquetado
científico sigue pinneado, que la caché se restaura antes de instalar y que los workflows no
usan contextos de paso (`runner`, `steps`) donde el parser de Actions los rechaza.
