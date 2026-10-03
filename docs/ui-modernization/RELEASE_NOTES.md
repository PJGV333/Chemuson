# Notas de release — modernización de interfaz

**Estado:** borrador de notas para la siguiente publicación; no es una release ni fija versión/tag. La workflow existente `release` genera automáticamente las notas de GitHub (`generate_release_notes: true`); este resumen curado puede servir de descripción de producto al preparar esa publicación.

## Nueva interfaz de ChemUSON

- Barra de aplicación con pestañas de documentos y acciones frecuentes; rail de herramientas con flyouts y panel lateral integrado.
- Temas claro/oscuro, iconografía SVG preparada para HiDPI, onboarding de primera ejecución y exploración de plantillas con selección de un clic.
- Paleta de comandos con `Ctrl+P`. `Ctrl+K` conserva la acción de Clean2D de un paso; `Ctrl+Shift+K` y `Ctrl+Alt+K` conservan las otras acciones Clean2D.
- Acciones de archivos y exportación, tabs/estado de documento y accesos existentes se mantienen disponibles. La modernización no cambia química, Clean2D, nomenclatura ni formato `.cmsn`.

La revisión manual de la interfaz integrada en KDE/Wayland fue aprobada. Evidencia visual: [`after/README.md`](after/README.md). El manual completo, incluida la guía de atajos, está en [`../MANUAL_USUARIO.md`](../MANUAL_USUARIO.md).

## Versión

La fuente actual sigue en `0.3.0-dev`. Dado que esta modernización integrada es un hito visible de producto, `0.4.0-dev` es una propuesta coherente para la siguiente publicación de prueba, sujeta a aprobación y al flujo oficial de release. No se realizó bump, tag, publicación ni cambio de metadatos AppStream en esta fase.
