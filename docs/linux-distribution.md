# Distribucion Linux dual: Flatpak + AppImage

Fecha: 2026-02-26

## Objetivo

- **Canal principal instalable:** Flatpak.
- **Canal portable:** AppImage (se mantiene).

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
