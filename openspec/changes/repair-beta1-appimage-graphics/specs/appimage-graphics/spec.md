# APPIMAGE-GRAPHICS-001 — startup gráfico y OpenGL en AppImage

## Purpose

Asegurar que el AppImage distribuido seleccione plugins Qt compatibles con el sistema anfitrión y presente ChemUSON en KDE/Wayland, conservando el editor 2D y el visor OpenGL 3D.

## ADDED Requirements

### Requirement: AppImage presents the main editor on supported desktop sessions
El AppImage MUST mostrar una ventana principal funcional en KDE/Wayland y MUST usar X11/XWayland cuando se solicite o el entorno sea compatible. Su entorno de plugin/biblioteca MUST resolver versiones ABI compatibles de Qt y el stack gráfico de usuario/host. La prueba MUST diferenciar plugin QPA seleccionado, librerías efectivamente cargadas y proceso/ventana visible; no basta con comprobar que el ejecutable sigue vivo.

#### Scenario: Native Wayland desktop startup
- **WHEN** el usuario inicia el AppImage en una sesión KDE Wayland real con Mesa/EGL funcional en el host
- **THEN** ChemUSON muestra su ventana principal sin `EGL not available` ni ausencia de QPA OpenGL que impida presentar el editor
- **AND** el proceso permanece estable y cierra normalmente.

#### Scenario: X11/XWayland startup
- **WHEN** se selecciona Qt `xcb` en un servidor X compatible
- **THEN** el plugin `libqxcb` y su integración GLX/EGL cargan con bibliotecas compatibles del host/bundle
- **AND** la ventana presenta el canvas y el cierre no produce SIGABRT.

### Requirement: Editor 2D does not depend on unused 3D initialization
La ventana principal y el canvas 2D MUST poder presentarse y editar sin requerir que el visor molecular 3D solicitado opcionalmente haya creado un contexto OpenGL. La solución MUST NOT deshabilitar OpenGL/3D, sustituir el renderer, ni rebajar permanentemente el perfil/contexto usado por el visor.

#### Scenario: Open editor before requesting 3D
- **WHEN** el usuario inicia una sesión nueva y no abre el dock 3D
- **THEN** la ventana y el editor 2D son visibles y utilizables
- **AND** el visor 3D conserva su funcionalidad para cuando se abra.

#### Scenario: OpenGL viewer remains available
- **WHEN** el usuario abre/cierra el visor 3D con una molécula simple y otra compleja
- **THEN** ambos modelos se dibujan con el renderer/luz previstos y la ventana sigue respondiendo
- **AND** cerrar el dock o la aplicación no provoca SIGABRT.

### Requirement: AppImage validation includes a real graphical startup gate
El flujo de validación del AppImage MUST mantener los validadores offscreen existentes y añadir un smoke con plataforma Qt gráfica real (xcb sobre Xvfb y Mesa por software/llvmpipe disponible) sobre el ejecutable del AppImage extraído. Debe verificar ventana visible/startup mediante señal/captura específica y límite temporal, y observar cierre limpio. `QT_QPA_PLATFORM=offscreen` MUST NOT usarse para satisfacer este requisito. Este smoke automático MUST quedar separado de retest manual Wayland/XWayland.

#### Scenario: CI validates graphical package startup
- **WHEN** Build Preview construye un AppImage Type 2
- **THEN** el validador ejecuta el binary exacto desde AppDir bajo Xvfb/software Mesa y captura evidencia de interfaz visible
- **AND** cualquier timeout, falta de ventana, error GL fatal o cierre anormal bloquea la carga del artifact.

### Requirement: Fontconfig diagnostics remain independent and evidence-based
El empaquetado MUST NOT silenciar warnings de Fontconfig. Si las bibliotecas bundled y la configuración Fontconfig del host son incompatibles, la validación MUST identificar versiones/rutas/entrada XML y corregir sólo la incompatibilidad demostrada, sin cambiar archivos del host ni datos de usuario.

#### Scenario: Host uses newer Fontconfig configuration syntax
- **WHEN** el AppImage lee `/etc/fonts` con su propia libfontconfig/libfreetype
- **THEN** se demuestra compatibilidad o se resuelve la versión/procedencia del parser
- **AND** no se filtran logs, se quitan warnings a ciegas o se reemplaza la configuración global.

## Compatibility

No cambia versiones de aplicación, química, Clean2D, ChemName, API de persistencia `.cmsn`, IDs de plataforma, nombres de artifacts, driver del host ni updater. El AppImage corregido debe conservar provenance SHA exacta del Build Preview, `publication=false` y no altera la referencia RC anterior.
