# Distribucion Linux dual: Flatpak + AppImage

Fecha: 2026-02-26

## Objetivo

- **Canal principal instalable:** Flatpak.
- **Canal portable:** AppImage (se mantiene).

## Artefactos publicados

- Flatpak bundle: `Chemuson-vX.Y.Z-linux-x86_64.flatpak`
- AppImage portable: `Chemuson-vX.Y.Z-linux-x86_64.AppImage`
- Sidecars AppImage (best-effort):
  - `Chemuson-vX.Y.Z-linux-x86_64.AppImage.updateinfo`
  - `Chemuson-vX.Y.Z-linux-x86_64.AppImage.update.json`
  - `Chemuson-vX.Y.Z-linux-x86_64.AppImage.zsync` (si `zsyncmake` esta disponible)

## Uso para usuarios finales

### Flatpak (principal instalable)

Instalar desde bundle descargado:

```bash
flatpak install --user ./Chemuson-vX.Y.Z-linux-x86_64.flatpak
flatpak run io.github.PJGV333.Chemuson
```

Desinstalar:

```bash
flatpak uninstall io.github.PJGV333.Chemuson
```

Nota:
- En el MVP actual, el update de Flatpak depende del origen remoto configurado.
- Si solo instalas desde bundle local y sin remote persistente, el update tipico es reinstalar el bundle de una version mas nueva.

### AppImage (portable)

```bash
chmod +x Chemuson-vX.Y.Z-linux-x86_64.AppImage
./Chemuson-vX.Y.Z-linux-x86_64.AppImage
```

## Uso para mantenedores

### Build local AppImage

```bash
pyinstaller --clean --noconfirm chemuson.spec
bash packaging/linux/build_appimage.sh \
  "0.2.1" \
  "dist" \
  "dist-appimage" \
  "PJGV333" \
  "Chemuson" \
  "stable" \
  "v0.2.1"
```

### Build local Flatpak

Requiere `flatpak` y `flatpak-builder`.

```bash
bash packaging/linux/build_flatpak.sh \
  "0.2.1" \
  "stable" \
  "dist-flatpak" \
  "packaging/flatpak/io.github.PJGV333.Chemuson.yml"
```

Opcional:
- definir `CHEMUSON_FLATPAK_REPO_URL` para que el script emita un `.flatpakref`.

## Pipeline CI/CD

- `release.yml`
  - `build_linux`: AppImage + metadata de update.
  - `build_flatpak`: build de bundle Flatpak.
  - `release`: agrega checksums, firma opcional HMAC y publica assets en GitHub Releases.

- `test.yml`
  - `flatpak-smoke`: validacion sintactica del manifiesto Flatpak.

## Checksums y firma

En release, todos los artifacts se procesan con:

- `packaging/release/generate_checksums.py` -> `checksums.txt`
- `packaging/release/sign_hmac.py` -> `*.sig` (si `CHEMUSON_SIGN_KEY` existe)

Esto aplica tambien a `.flatpak` y `.AppImage`.
