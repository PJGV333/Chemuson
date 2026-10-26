# Propuesta: preparar ChemUSON 0.3.0-beta.4

## Motivación

La campaña incorpora la identidad visual aprobada por el propietario, corrige inconsistencias entre la documentación y la clasificación excepcional de `v0.3.0-beta.3`, reconoce a Jesús David Albino-Espíndola, organiza salidas locales de empaquetado y prepara una beta posterior a beta.3 sin cambiar las capacidades científicas existentes.

La rama `main` estaba sincronizada con `origin/main` en `9c171a8267d632e9bd914392680c316c28fe5e87`. La rama de trabajo es `release/v0.3.0-beta.4-prep`. Los tags anotados beta.1, beta.2 y beta.3 se preservan.

## Alcance

- Incorporar sin rediseño el PNG aprobado `/home/unison-pjgv/Descargas/logo_ChemUSON.png`, verificarlo y producir variantes raster/ICO y un SVG contenedor del PNG original.
- Mostrar el emblema en la ventana y en los paquetes PyInstaller, Inno Setup, AppImage, Flatpak/AppStream, README en español e inglés y manual de usuario. Aclarar que es identidad visual del proyecto y no un escudo institucional oficial.
- Identificar y eliminar sólo salidas locales ignoradas y obsoletas; preservar los manifiestos trackeados de `dist-appimage-diagnostic/`. Ajustar únicamente el destino de Inno Setup para respetar `CHEMUSON_OUTPUT_DIR`; mantener el comportamiento de los workflows oficiales sin dicha variable.
- Actualizar documentación sobre la decisión excepcional del propietario que dejó `v0.3.0-beta.3` con `prerelease=false` y distintivo GitHub **Latest**, sin convertirla en política de las futuras betas.
- Auditar SemVer, selección por canal, preferencias persistidas, Windows, AppImageUpdate, Flatpak, downgrade e integridad. Corregir el updater sólo ante una regresión demostrable.
- Añadir Jesús David Albino-Espíndola a las dos listas documentales de autores, sin asterisco ni cambios de acceso técnico.
- Preparar versión canónica `0.3.0-beta.4`, metadatos AppStream y notas públicas bilingües (español antes de inglés) que destaquen identidad visual y mantenimiento, no capacidades científicas heredadas.
- Corregir, sólo en pruebas, el harness AST de `run_app()` tras el fallo de `test_startup_success_order` en el Build Preview #40371397035: simular `QIcon`, la ruta del recurso y su asignación a la aplicación/ventana; probar además el error de icono ausente sin alterar los contratos de manejo de fallos.
- Validar de forma focalizada y preparar el Build Preview oficial con cuatro paquetes desde un SHA único cuando haya acceso autorizado a GitHub.

## Fuera de alcance y protecciones

- No cambiar Clean2D, ChemName, química, geometría/renderizado molecular, persistencia `.cmsn`, dependencias o funcionalidades científicas.
- No borrar/recrear/mover tags ni Releases históricos, no cambiar preferencias de actualización y no cambiar manualmente la clasificación GitHub de beta.3 o beta.4.
- No incluir `tests/` en ejecutables; no ejecutar `pytest -q` completo.
- No crear el tag `v0.3.0-beta.4` ni una Release. La aceptación manual de los paquetes exactos pertenece al propietario.
- No invitar a Jesús David Albino-Espíndola al repositorio ni modificar miembros/permisos.
- No afirmar que se publicó una rama/PR o se ejecutó Actions sin evidencia autenticada.
