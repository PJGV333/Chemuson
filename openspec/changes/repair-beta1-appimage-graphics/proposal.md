# Propuesta — reparar el arranque gráfico del AppImage beta.1

## Problema

El AppImage Type 2 construido desde `a38fa8045c858f8ff4b6002950f57e50292aea81` no presenta correctamente la ventana principal en CachyOS/KDE. La evidencia del propietario incluye fallo EGL en Wayland (`EGL not available`, `QOpenGLWidget is not supported`, `QRhiGles2` sin contexto), fallo GLX en XWayland (`qglx_findConfig`, `Could not initialize GLX`, SIGABRT) incluso con `LIBGL_ALWAYS_SOFTWARE=1`, además de errores de Fontconfig al leer configuración del host.

El ejecutable de desarrollo `.venv` sí ha funcionado en el mismo host. Esto apunta a una diferencia de empaquetado/inicialización, pero no determina por sí solo si la causa es colisión de bibliotecas, plugin Qt, dependencia ausente, inicialización del visor 3D o una combinación. Fontconfig se investigará como subsistema separado y sus warnings no se ocultarán.

## Objetivo

Corregir con el menor cambio mantenible el empaquetado o la inicialización gráfica para que el AppImage beta.1 arranque visiblemente en KDE/Wayland, permita usar el editor 2D independientemente de un fallo del visor y conserve la funcionalidad OpenGL/iluminación/renderizado 3D. Mantener soporte X11/XWayland cuando el entorno anfitrión sea compatible.

## Alcance

- Reproducir y preservar el AppImage original y sus logs; no modificarlo.
- Inspeccionar `build_appimage.sh`, `AppRun`, `rthook_qt.py`, `chemuson.spec`, plugins QPA `xcb`/`wayland`, bibliotecas PyInstaller/Qt/GLVND/Mesa/Fontconfig/FreeType y workflows de distribución.
- Comparar la resolución real de bibliotecas en el ejecutable congelado y el entorno de desarrollo, usando `LD_DEBUG`, `QT_DEBUG_PLUGINS`, `/proc/<pid>/maps`, `readelf`, `ldd` y pruebas equivalentes.
- Implementar sólo la causa demostrada y añadir un gate de arranque gráfico real para el AppImage. Conservar los tests offscreen como pruebas complementarias.
- Documentar claramente qué se valida automáticamente y qué aceptación real queda pendiente.

## Fuera de alcance

No modificar controladores, bibliotecas ni configuración del sistema anfitrión; no borrar ni reescribir datos del usuario; no cambiar algoritmos químicos, Clean2D, ChemName o `.cmsn`; no deshabilitar OpenGL/3D ni fijar `offscreen` como solución. No publicar/taggear/fusionar con beta ni tocar `main` o la referencia RC defectuosa.

## Aceptación

1. El AppImage original se reproduce en un entorno gráfico real, con capturas y logs que distinguen Wayland y XWayland, y conservando su SHA.
2. La causa raíz se apoya en rutas/versiones de bibliotecas/plugins cargados, comparación con `.venv`, y pruebas de aislamiento que atribuyen o descartan conflictos GL/EGL/GLX, plugin Qt, inicialización del visor y Fontconfig por separado.
3. La corrección mínima hace aparecer la ventana principal; el canvas 2D permanece utilizable sin requerir que el visor 3D cree antes un contexto.
4. En el paquete corregido se abre/cierra el visor 3D y se representan una molécula sencilla y una compleja, conservando el contexto/iluminación reales. Arranque/cierre no provoca SIGABRT.
5. Una prueba de integración gráfica para el AppImage usa una plataforma Qt real (xcb con Xvfb y Mesa software en CI, no `offscreen`) y verifica el startup/visibilidad de la ventana. Los checks offscreen existentes siguen activos, pero no reemplazan este gate.
6. El host CachyOS se prueba mediante Wayland real; XWayland se prueba por separado. Un entorno automatizado no se presenta como aceptación de ambos desktops.
7. Los cuatro paquetes sólo se reconstruyen por Build Preview existente con provenance del SHA corregido; la RC anterior permanece intacta. No se publica una Release y publicación pública permanece NO-GO.
