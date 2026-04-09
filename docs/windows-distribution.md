# Distribucion Windows: portable + instalador oficial

Fecha: 2026-02-26

## Decision de instalador

Se usa **Inno Setup** para el paquete instalable oficial de Windows.

Justificacion tecnica:
- Ya estaba alineado con el RFC y con el flujo de empaquetado actual.
- Permite convivir sin friccion con el binario portable (`.exe`).
- Soporta ejecucion silenciosa (`/VERYSILENT`) para integracion con updater.
- Se integra bien con firma Authenticode en pipeline.

## Artefactos publicados

- Portable: `Chemuson-vX.Y.Z-windows-x86_64-portable.exe`
- Instalador: `Chemuson-vX.Y.Z-windows-x86_64-setup.exe`

Ambos se publican en paralelo en cada release.

## Updater para instalaciones

- El instalador crea el marcador `.chemuson-installed` dentro del directorio de la app.
- En runtime, Chemuson detecta si corre como instalado o portable.
- Si esta instalado en Windows, el updater prioriza assets `installer`.
- En modo `notify`, se ofrece descargar el setup y aplicarlo al cerrar.
- En modo `silent`, prepara el setup en segundo plano y lo aplica al cierre.

## Build local (Windows)

PowerShell:

```powershell
python -m pip install --upgrade pip
pip install -r requirements.txt
pip install pyinstaller
pyinstaller --clean --noconfirm chemuson.spec
```

Generar instalador:

```powershell
$env:CHEMUSON_VERSION = "0.2.3-beta.3"
& "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" "packaging\windows\Chemuson.iss"
```

Salida esperada:
- `dist\Chemuson.exe` (portable)
- `dist-installer\Chemuson-v0.2.3-beta.3-windows-x86_64-setup.exe` (instalador)

## Firma Authenticode (pipeline preparada)

Script local/pipeline:
- `packaging/windows/sign_artifacts.ps1`

Secretos esperados:
- `WINDOWS_CODESIGN_CERT_BASE64` (PFX en base64)
- `WINDOWS_CODESIGN_CERT_PASSWORD`
- opcional `WINDOWS_CODESIGN_TIMESTAMP_URL`

Ejemplo local:

```powershell
$env:WINDOWS_CODESIGN_CERT_BASE64 = "<pfx_base64>"
$env:WINDOWS_CODESIGN_CERT_PASSWORD = "<password>"
& "packaging\windows\sign_artifacts.ps1" -InputDir "release-assets"
```

## GitHub Actions

- `.github/workflows/release.yml`
  - construye portable + setup en `build_windows`.
  - firma `.exe/.msi` de manera opcional (si hay secretos).
  - publica assets junto con checksums/manifiesto.

Comandos clave ejecutados en `build_windows`:

```powershell
python -m pip install --upgrade pip
pip install -r requirements.txt
pip install pyinstaller
pyinstaller --clean --noconfirm chemuson.spec
$env:CHEMUSON_VERSION = "<version>"
& "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" "packaging\windows\Chemuson.iss"
& "packaging\windows\sign_artifacts.ps1" -InputDir "release-assets"
```

- `.github/workflows/test.yml` (`windows-smoke`)
  - ejecuta pruebas del updater Windows.
  - compila el instalador con Inno en modo smoke.

## Limitaciones de CI (documentadas)

En CI no se ejecuta el instalador sobre el sistema runner para evitar:
- prompts UAC/privilegios,
- mutaciones del host efimero no confiables para regression tests.

En su lugar se simulan pasos criticos:
- seleccion de asset correcto (`installer` vs `portable`),
- construccion de comandos de update Inno,
- compilacion real del setup con Inno Setup.
