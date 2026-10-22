# Diseño — investigación y reparación gráfica del AppImage

## Baseline/evidencia

`baseline.md` fija la rama/SHA de partida, la integridad del AppImage de comparación, el estado de Git, el intérprete de pruebas y la primera reproducción. Los logs/capturas/binarios/inspecciones adicionales se conservan fuera del repo bajo `/tmp/chemuson-appimage-gl-repair-78b0139/`. No se modifica el artefacto original ubicado en `/mnt/HDD_4TB/Descargas_provisional/chemuson-preview-linux-appimage/`.

## Diagnóstico confirmado antes de editar producción

El fallo GL proviene del `libstdc++.so.6` congelado por PyInstaller, no de un plugin Qt ausente ni de GL/EGL/GLX empaquetados. El AppImage carga correctamente `libqwayland.so`/integración EGL o `libqxcb.so`/integración GLX desde su bundle, y resuelve `libGL.so.1`, `libEGL.so.1`, GLX/GLdispatch y Mesa desde CachyOS (`/usr/lib`). Mesa 26.2.4 carga `libgallium-26.2.4-arch3.1.so` y `libSPIRV-Tools.so`; este último requiere `GLIBCXX_3.4.32`. El `libstdc++.so.6` de `_MEIPASS` exporta como máximo `GLIBCXX_3.4.30`; el loader informa `version lookup error ... GLIBCXX_3.4.32 not found`. El fallo provoca la inicialización fallida de EGL en Wayland y la ruta GLX aborta con SIGABRT en XWayland.

Prueba A/B: quitar `LD_LIBRARY_PATH=/home/ccachyavgp/Aplicaciones/Orca_QC` no cambia el fallo. Pre-cargar únicamente `/usr/lib/libstdc++.so.6` hace que Mesa valide `GLIBCXX_3.4.32`, desaparezcan los errores EGL/GLX y la ventana ChemUSON aparezca en XWayland. El mismo bundle con el runtime normal no presenta ventana en Wayland y aborta en XWayland (exit `-6`). El `.venv` usa `libstdc++`, GL/EGL/GLX, Fontconfig y FreeType del host; su contexto OpenGL 4.6 Core/radeonsi funciona. El metadato de los plugins congelados corresponde a Qt 6.11.0; el `.venv` usa Qt 6.11.2. Los plugins cargan de forma completa, por lo que la diferencia de parche Qt no explica el error de símbolo observado.

Fontconfig es una incompatibilidad distinta: PyInstaller congeló Fontconfig 2.13.1 y FreeType 2.11.1, que leen la configuración del host CachyOS Fontconfig 2.18.3 y rechazan `xsi:nil`/aliases nuevos. Pre-cargar sólo las bibliotecas Fontconfig/FreeType del host elimina esos warnings mientras el error GLIBCXX y el SIGABRT XCB permanecen; pre-cargar sólo libstdc++ corrige GL pero deja los warnings. Por tanto, Fontconfig no causa la ventana invisible, pero sí debe resolverse de forma independiente.

## Estrategia de corrección respaldada por evidencia

- Excluir de la colección PyInstaller los ABI de integración del host `libstdc++.so.*`, `libfontconfig.so.*` y `libfreetype.so.*`. El loader usará las versiones del sistema anfitrión; las dependencias químicas/Qt y el runtime de PyInstaller permanecen congelados. No incluir drivers GL/Mesa ni manipular `/etc/fonts`, `LD_PRELOAD`, `LD_LIBRARY_PATH` o preferencias globales.
- Mantener Qt/PyQt, plugins QPA, perfil OpenGL 4.1 Core, renderer, shaders e iluminación intactos. No diferir ni eliminar el dock 3D: la prueba A/B indica que el problema es el runtime ABI que carga Mesa, no la jerarquía del visor.
- Validar el archivo congelado para impedir que estos tres SONAMEs vuelvan al CArchive y añadir un smoke real bajo Xvfb/XCB + Mesa llvmpipe con captura de la ventana, contexto OpenGL 4.1 Core y renderizado de escenas simples/complejas.

## Prueba gráfica distribuible

Añadir un validador que lance el binario del AppImage extraído bajo una sesión gráfica controlada (xcb + Xvfb + Mesa llvmpipe/software, con GL/EGL visibles) y observe que el proceso presenta una ventana y permanece estable; tomar una captura y verificar que no es una imagen vacía/placeholder. No usar `QT_QPA_PLATFORM=offscreen` como prueba gráfica. El launch smoke actual offscreen permanece complementario. Para 3D, una prueba separada debe abrir/cerrar el dock con una molécula sencilla y otra compleja en renderer OpenGL real o en software Mesa, e inspeccionar cierre limpio.

La sesión KDE Wayland de este host está disponible y debe utilizarse para el before/after real. Para XWayland, ejecutar la misma GUI sobre `DISPLAY` con `QT_QPA_PLATFORM=xcb`; no inferir X11 a partir de tests Qt offscreen.

## Límites

Los warnings Fontconfig pueden ser reales pero no necesariamente la causa de la ventana invisible. Un startup Xvfb/llvmpipe es evidencia de integración gráfica automatizada, no de drivers AMD/Wayland reales. La evaluación final debe declarar ambos niveles y no confundirlos.
