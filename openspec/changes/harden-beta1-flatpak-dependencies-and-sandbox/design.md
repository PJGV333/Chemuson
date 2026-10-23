# Diseño — dependencias y sandbox de beta.1

## Baseline y procedencia

La rama local `release/v0.3.0-beta.1-security-prep` parte del SHA exacto `08a90c199f526639453a8ae8b82da957a1458d08`; al crearla, el árbol estaba limpio y coincidía con la rama prep remota. Baseline proporcional registrada en `baseline.md`; pytest completo no se ejecutó por la petición explícita y la política de pruebas vigente.

La inspección estática no encontró imports `PIL`/Pillow en `src/`, `tools`, `packaging` ni tests de runtime; el manual describe erróneamente que Pillow manipula anotaciones. Las imágenes pegadas/persistidas se procesan con Qt (`QImage`, `QPixmap`). `rdkit.Chem.Draw.rdMolDraw2D` importó correctamente bajo un import blocker que rechaza todo `PIL`, y `tests/test_image_clipboard_persistence.py` pasó 9/9 con ese blocker.

Consulta de advisories OSV el 2026-10-23: Pillow 12.0.0 reporta 77 registros, incluidos CVE-2026-42311 y CVE-2026-54058; 12.3.0 todavía reporta CVE-2026-96637, CVE-2026-96643, CVE-2026-96644 y CVE-2026-96778; 13.0.0 no devuelve advisories en la consulta por versión. La metadata de RDKit en PyPI declara Pillow transitivamente. Build Preview #39899351520 resolvió Pillow 13.0.0 para PyInstaller; la inspección de ambos CArchive encontró cinco extensiones nativas `PIL/*` por ejecutable. El Flatpak, que instala RDKit con `--no-deps`, no contiene el paquete top-level `site-packages/PIL` ni metadata Pillow, aunque sí contiene recursos legacy propios de `rdkit/sping/PIL`. Por tanto, se elimina la dependencia directa y Flatpak, y se excluyen/validan los recursos top-level del CArchive; no se afirma que el entorno de build no resuelva Pillow.

## Decisiones

1. Quitar Pillow de `pyproject.toml`, `requirements.txt` y del módulo/wheel específico en Flatpak; añadir `PIL` a exclusiones PyInstaller y comprobar el CArchive Windows/AppImage para impedir recursos transitorios. La instalación normal de RDKit puede resolver Pillow transitivamente en el entorno de build; se audita su versión por separado. No cambiar el comportamiento de imágenes: las rutas relevantes usan Qt/RDKit nativo.
2. La campaña completada `stabilize-ai-molecular-assistant-integration` exigía históricamente conservar Pillow; esta revisión lo supersede con una delta `MODIFIED` para que el contrato vigente refleje el audit actual, sin alterar su baseline/decisión histórica.
3. Actualizar el test de coherencia para prohibir la dependencia directa y el módulo Flatpak; hacer que el validador del ejecutable inspeccione el CArchive y rechace entradas top-level `PIL`/`Pillow`, con cobertura de la exclusión en la spec PyInstaller.
4. Sustituir `--filesystem=home` por `--filesystem=home/.chemuson:create` (autosave, cachés/logs y archivos ChemUSON existentes) y `--filesystem=home/.config/Chemuson:create` (biblioteca de plantillas y fallback de configuración). No cambiar código de portal ni rutas de persistencia.
5. Conservar `--share=network`: Assistant se conecta al endpoint remoto configurado por el usuario y resolución PubChem es opcional/consentida. El permiso sigue siendo amplio a nivel de sandbox y se documenta como riesgo residual; no implica que ChemUSON contacte proveedores automáticamente.
6. Publicar cambios mediante commit/push normal sólo en `release/v0.3.0-beta.1-security-prep`. El push activa el Build Preview del patrón `release/**-prep`; no activa `test.yml` ni la suite completa. No mover `release/v0.3.0-beta.1-prep` ni `main`.
7. Verificar por separado los cuatro grupos de artefactos, sin confundir smoke automatizado con aceptación manual. AppImage y Windows portable ejecutan el smoke 2D→3D disponible; la interacción Qt file chooser en Flatpak requiere retest manual.

## Riesgos y límites

- Las rutas de datos/configuración concedidas contienen estado de ChemUSON y siguen siendo accesibles a la aplicación; esto es deliberado para mantener la compatibilidad actual. El permiso deja de abarcar el resto del home.
- El funcionamiento de QFileDialog con portal y sin acceso general a home no se certifica con la build: el uso real de archivos en la app Flatpak queda como aceptación manual obligatoria.
- El permiso de red permite conectividad general para cualquier código de la aplicación, aunque la funcionalidad normal la use para endpoints configurados y consulta externa explícita. No se modifica la política de privacidad ni se realiza pentest.
- Pillow se elimina como dependencia directa y del Flatpak; la metadata de RDKit puede instalarlo transitivamente en el entorno de build (13.0.0 para el preview auditado). La spec excluye `PIL` y el gate lee el CArchive para verificar cero recursos Pillow en Windows/AppImage. No se prueba el universo de plugins opcionales de RDKit que usan Pillow porque ChemUSON no los llama; un cambio futuro de resolver requiere volver a auditar el paquete.
- Cambiar permisos Flatpak requiere empaquetado exacto y prueba instalada; hasta entonces queda NO-GO para publicar.
