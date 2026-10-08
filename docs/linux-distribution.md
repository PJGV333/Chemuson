# Distribución Linux dual: Flatpak + ejecutable portable

Fecha: 2026-10-08

## Objetivo

- **Canal principal instalable:** Flatpak.
- **Canal portable:** AppImage Type 2 auténtico, generado desde el ejecutable PyInstaller dentro de un AppDir validado. El AppImage sólo se empaqueta si el smoke del ejecutable congelado comprueba los 69 SVG y rasteriza los iconos esenciales mediante QtSvg en ambos temas.

## Estrategia Flatpak reproducible

- El manifiesto Flatpak declara dependencias Python como modulos explicitos con `url + sha256` pinneados.
- El modulo de Chemuson instala el paquete con:
  - `pip3 install --prefix=/app --no-build-isolation --no-deps .`
- Resultado:
  - evita resolucion dinamica de PyPI en el paso de instalacion de Chemuson,
  - hace el build mas reproducible (mismas fuentes y checksums),
  - previene fallos tipo `No matching distribution found for PyQt6` durante `pip install .`.
- Runtime KDE actualizado:
  - de `6.6` (EOL) a `6.10` (rama soportada actual para Qt6/KDE6).

## Artefactos publicados

- Flatpak bundle: `Chemuson-vX.Y.Z-linux-x86_64.flatpak`
- AppImage Type 2: `Chemuson-vX.Y.Z-linux-x86_64.AppImage`.
- Sidecars del updater oficial (no se generan para previews):
  - `Chemuson-vX.Y.Z-linux-x86_64.AppImage.updateinfo`
  - `Chemuson-vX.Y.Z-linux-x86_64.AppImage.update.json`
  - `Chemuson-vX.Y.Z-linux-x86_64.AppImage.zsync` (requerido para el canal AppImageUpdate oficial)

## Uso para usuarios finales

### Flatpak (principal instalable)

Pagina indice de canales realmente publicados:

```text
https://pjgv333.github.io/Chemuson/
```

Instalar desde el canal oficial estable:

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

Desinstalar:

```bash
flatpak uninstall io.github.PJGV333.Chemuson
```

Nota:
- El remoto oficial se publica en GitHub Pages bajo `flatpak/<canal>/repo/`.
- El indice de GitHub Pages solo enlaza canales ya publicados para evitar 404.
- Si `flatpak/stable/...` aun no existe, el remoto estable todavia no fue publicado; usa `beta` o ejecuta primero una release estable.
- Si instalas desde `.flatpakref` o desde un bundle generado con `CHEMUSON_FLATPAK_REPO_URL`, `flatpak update` encuentra futuras versiones automaticamente.
- Si solo instalaste un bundle local sin remote persistente, deberas reinstalar manualmente.

### AppImage Type 2 Linux

```bash
chmod +x Chemuson-vX.Y.Z-linux-x86_64.AppImage
./Chemuson-vX.Y.Z-linux-x86_64.AppImage
```

## Uso para mantenedores

### Build local AppImage Type 2

Requiere `curl`, `desktop-file-utils`, `appstream` (`appstreamcli`) y las dependencias Qt del runner. `zsyncmake` se reutiliza desde el AppImageKit `appimagetool` oficial pinneado. El helper verifica el SHA-256 del asset, su commit/versión y falla si la descarga cambia. No requiere ni instala `linuxdeploy`: PyInstaller ya contiene el bundle de PyQt6/RDKit y `appimagetool` empaqueta ese AppDir.

```bash
pyinstaller --clean --noconfirm chemuson.spec
python packaging/release/validate_packaged_icons.py --executable dist/Chemuson
bash packaging/linux/build_appimage.sh \
  "<version>" \
  "dist" \
  "dist-appimage" \
  "PJGV333" \
  "Chemuson" \
  "stable" \
  "v<version>" \
  "$(git rev-parse HEAD)"
python packaging/release/validate_appimage.py \
  --appimage "dist-appimage/Chemuson-v<version>-linux-x86_64.AppImage" \
  --version "<version>" --source-sha "$(git rev-parse HEAD)" \
  --build-type release --channel stable --tag "v<version>"
```

### Build local Flatpak

Requiere `flatpak` y `flatpak-builder`.

```bash
bash packaging/linux/build_flatpak.sh \
  "<version>" \
  "stable" \
  "dist-flatpak" \
  "packaging/flatpak/io.github.PJGV333.Chemuson.yml"
```

Opcional:
- definir `CHEMUSON_FLATPAK_REPO_URL` para que el script enlace el bundle a un remoto oficial y emita `.flatpakrepo` + `.flatpakref`.
- definir `CHEMUSON_FLATPAK_GPG_KEY_ID`, `CHEMUSON_FLATPAK_GPG_HOMEDIR` y `CHEMUSON_FLATPAK_PUBLIC_KEY_FILE` para firmar el repo y publicar `GPGKey`.
- definir `CHEMUSON_FLATPAK_REMOTE_NAME` si quieres sobreescribir el nombre sugerido del remote (`chemuson-<canal>`).
- exportar `ARCH` si se requiere override de arquitectura (default: `x86_64`).

## Troubleshooting Flatpak

- Error de ruta de manifiesto:
  - mensaje esperado: `Flatpak manifest not found: ...`
  - accion: validar el cuarto argumento de `build_flatpak.sh` o `FLATPAK_MANIFEST`.
- Error de red / DNS / sandbox:
  - mensaje esperado: `Unable to configure flathub remote` o `network/DNS or sandbox egress restrictions`.
  - accion: verificar conectividad a `https://dl.flathub.org`, DNS y politicas de sandbox.
- Error de dependencias Python en `flatpak-builder`:
  - mensaje esperado: `failed while installing Python dependencies inside flatpak-builder`.
  - accion: revisar wheels pinneados, `sha256` y compatibilidad ABI de Python con el runtime seleccionado.

## Pipeline CI/CD

- `release.yml`
  - `build_linux`: AppImage Type 2 real desde AppDir, validado por firma/extracción/resources/version/arranque headless y smoke real de iconos del ejecutable extraído, con el contrato existente de AppImageUpdate.
- `build-preview.yml`
  - compila desde ramas `release/**-prep` y adjunta portable Windows, installer Inno, portable Linux y Flatpak como artifacts Actions con SHA/procedencia; no publica Release, tags, remoto, `gh-pages` ni manifests públicos.
  - `build_flatpak`: build de bundle Flatpak + repo OSTree + validacion explicita de `.flatpakrepo/.flatpakref` y `repo/summary`.
  - `publish_flatpak_remote`: valida el payload antes y despues de moverlo entre jobs, publica el repo oficial por canal en `gh-pages` y verifica las URLs publicadas.
  - `release`: agrega checksums, firma opcional HMAC y publica assets en GitHub Releases.
- Requisito de plataforma:
  - habilitar GitHub Pages apuntando a la rama `gh-pages` para exponer el remoto oficial.
  - opcional pero recomendado: configurar secrets `FLATPAK_GPG_PRIVATE_KEY_BASE64` y `FLATPAK_GPG_KEY_ID` para firmar el remoto oficial.
  - si el repo solo tiene releases beta hasta ahora, ejecutar una vez una release estable para sembrar `flatpak/stable/...`.

- `test.yml`
  - `flatpak-smoke`: validacion sintactica del manifiesto Flatpak.

## Previews y checksums

Para ejecutar una compilación previa a publicación y descargar sus cuatro artifacts, sigue [docs/release/PREVIEW_BUILDS.md](release/PREVIEW_BUILDS.md). El preview Linux omite metadata pública de update y su bundle Flatpak no enlaza un remoto Chemuson.

## Checksums y firma

En release, todos los artifacts se procesan con:

- `packaging/release/generate_checksums.py` -> `checksums.txt`
- `packaging/release/sign_hmac.py` -> `*.sig` (si `CHEMUSON_SIGN_KEY` existe)

Esto aplica tambien a `.flatpak` y `.AppImage`.
