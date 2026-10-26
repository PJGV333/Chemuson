# Tareas

- [x] Leer `AGENTS.md`, política de pruebas, OpenSpec beta.3 y contratos vigentes de release/preview.
- [x] Confirmar `main == origin/main` en `9c171a8267d632e9bd914392680c316c28fe5e87`, registrar status inicial, verificar que beta.1–beta.3 son tags existentes y crear `release/v0.3.0-beta.4-prep` desde main.
- [x] Verificar PNG aprobado: integridad, dimensiones, modo y SHA-256; inspeccionar visualmente el diseño original.
- [x] Capturar baseline focal, tamaños/clasificación de carpetas ignoradas, manifiestos trackeados, dry-run de `git clean -ndX`, clasificación pública beta.3 y disponibilidad de credenciales.
- [ ] Incorporar OpenSpec beta.4, añadir el maestro inalterado y generar variantes raster/ICO/SVG sin vectorización del arte.
- [ ] Integrar icono en QApplication/PyQt6, PyInstaller portable, Inno Setup, AppImage, Flatpak/AppStream, README es/en y manual; validar visibilidad en claro/oscuro.
- [ ] Ajustar únicamente Inno Setup si hace falta para que `CHEMUSON_OUTPUT_DIR` escriba fuera del checkout; conservar paths por defecto y compatibilidad de Actions.
- [ ] Limpiar sólo artefactos generados obsoletos identificados; mantener los manifiestos históricos trackeados de `dist-appimage-diagnostic/`; reportar bytes recuperados y comprobar que no haya ejecutables trackeados.
- [ ] Corregir README es/en, política de versionado y manual para describir la clasificación excepcional de beta.3 sin convertirla en política.
- [ ] Auditar y probar SemVer, selección de canal/preferencias, AppImageUpdate `latest`/`prerelease`, metadatos, Windows, Flatpak, downgrade e integridad; aplicar corrección mínima sólo si se demuestra una regresión.
- [ ] Añadir Jesús David Albino-Espíndola después de los autores actuales, sin asterisco ni cambio de permisos.
- [ ] Actualizar la versión canónica/AppStream a `0.3.0-beta.4` y crear notas públicas bilingües, mantenimiento/logotipo solamente.
- [ ] Ejecutar validaciones focales y `openspec validate prepare-chemuson-0-3-0-beta-4 --strict`; actualizar `validation.md` con comandos y SHA exacto.
- [ ] Commit(s) descriptivos. Publicar rama y abrir/revisar PR si la sesión remota está autenticada y disponible.
- [ ] Ejecutar Build Preview oficial sobre SHA único; verificar preflight, checksums, provenance y los cuatro paquetes. No confundir pruebas locales con Actions.
- [ ] Integrar respetando branch protection únicamente si puede hacerse por el flujo aprobado; dejar el tag sin crear y presentar evidencia para aceptación final del propietario.
