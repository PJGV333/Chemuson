# Diseño — dependencias y sandbox de beta.1

## Baseline y procedencia

La rama local `release/v0.3.0-beta.1-security-prep` parte del SHA exacto `08a90c199f526639453a8ae8b82da957a1458d08`; al crearla, el árbol estaba limpio y coincidía con la rama prep remota. Baseline proporcional registrada en `baseline.md`; pytest completo no se ejecutó por la petición explícita y la política de pruebas vigente.

La inspección estática no encontró imports `PIL`/Pillow en `src/`, `tools`, `packaging` ni tests de runtime; el manual describe erróneamente que Pillow manipula anotaciones. Las imágenes pegadas/persistidas se procesan con Qt (`QImage`, `QPixmap`). `rdkit.Chem.Draw.rdMolDraw2D` importó correctamente bajo un import blocker que rechaza todo `PIL`, y `tests/test_image_clipboard_persistence.py` pasó 9/9 con ese blocker.

Consulta de advisories OSV el 2026-10-23: Pillow 12.0.0 reporta 77 registros, incluidos CVE-2026-42311 y CVE-2026-54058; 12.3.0 todavía reporta CVE-2026-96637, CVE-2026-96643, CVE-2026-96644 y CVE-2026-96778; 13.0.0 no devuelve advisories en la consulta por versión. Como Pillow no es necesario, la decisión es eliminarlo en lugar de fijar un wheel nuevo.

## Decisiones

1. Quitar Pillow de `pyproject.toml`, `requirements.txt` y del módulo/wheel específico en Flatpak. No cambiar el conjunto de comportamiento de imágenes: las rutas relevantes están implementadas con Qt/RDKit nativo.
2. La campaña completada `stabilize-ai-molecular-assistant-integration` exigía históricamente conservar Pillow; esta revisión lo supersede con una delta `MODIFIED` para que el contrato vigente refleje el audit actual, sin alterar su baseline/decisión histórica.
3. Actualizar el test de coherencia de dependencias para que proteja la ausencia de Pillow y extender el contrato existente del manifiesto Flatpak para rechazar Pillow y el permiso total a `home`.
4. Sustituir `--filesystem=home` por `--filesystem=home/.chemuson:create` (autosave, cachés/logs y archivos ChemUSON existentes) y `--filesystem=home/.config/Chemuson:create` (biblioteca de plantillas y fallback de configuración). No cambiar código de portal ni rutas de persistencia.
5. Conservar `--share=network`: Assistant se conecta al endpoint remoto configurado por el usuario y resolución PubChem es opcional/consentida. El permiso sigue siendo amplio a nivel de sandbox y se documenta como riesgo residual; no implica que ChemUSON contacte proveedores automáticamente.
6. Publicar la rama de preparación de forma normal y disparar Build Preview sólo por `workflow_dispatch`, evitando iniciar la suite completa. No mover `release/v0.3.0-beta.1-prep` ni `main`.
7. Verificar por separado los cuatro grupos de artefactos, sin confundir smoke automatizado con aceptación manual. AppImage y Windows portable ejecutan el smoke 2D→3D disponible; la interacción Qt file chooser en Flatpak requiere retest manual.

## Riesgos y límites

- Las rutas de datos/configuración concedidas contienen estado de ChemUSON y siguen siendo accesibles a la aplicación; esto es deliberado para mantener la compatibilidad actual. El permiso deja de abarcar el resto del home.
- El funcionamiento de QFileDialog con portal y sin acceso general a home no se certifica con la build: el uso real de archivos en la app Flatpak queda como aceptación manual obligatoria.
- El permiso de red permite conectividad general para cualquier código de la aplicación, aunque la funcionalidad normal la use para endpoints configurados y consulta externa explícita. No se modifica la política de privacidad ni se realiza pentest.
- Pillow se elimina como dependencia directa; se verifica que Qt y `rdMolDraw2D` funcionan sin el módulo PIL. No se prueba el universo de plugins opcionales de RDKit que usan Pillow porque ChemUSON no los llama.
- Cambiar permisos Flatpak requiere empaquetado exacto y prueba instalada; hasta entonces queda NO-GO para publicar.
