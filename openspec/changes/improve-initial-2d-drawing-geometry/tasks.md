# Tareas — geometría inicial 2D

## A. Auditoría, baseline y contratos
- [x] A.1 Crear rama aislada desde `origin/release/v0.3.0-beta.1-prep`, verificar merge `3e0ef606` y SHA base.
- [x] A.2 Leer y trazar eventos de ratón, previews, selección geométrica, comandos, coordenadas MolGraph y renderer.
- [x] A.3 Registrar baseline exacta (`compileall`, colección, pytest, Ruff) y reproducir casos de referencia con coordenadas/capturas.
- [x] A.4 Separar defecto del modelo, defecto visual y orientación libre deliberada; acordar contrato geométrico.
- [x] A.5 Validar OpenSpec estricto antes de modificar producción.

## B. Geometría de colocación
- [x] B.1 Centralizar regla esqueletal sp3 2D de 120° sin cambiar la referencia tetraédrica 109.5° ni el paso global.
- [x] B.2 Aplicar 120° a cadena acíclica libre y sugerencias automáticas locales para sp3; preservar 120° sp2/180° sp y anillos regulares por tamaño.
- [x] B.3 Hacer que una sugerencia por clic corto respete geometría local aun con estructura rotada; conservar snap de arrastre explícito.
- [x] B.4 Corregir resolución del release para que un clic corto sobre átomo existente no se confunda con la arista vecina.
- [x] B.5 Verificar que anillos existentes, aromaticidad, dobles, etiquetas y representación editorial no cambian.

## C. Regresiones y validación
- [x] C.1 Añadir pruebas numéricas de anillos 3–6, transformaciones rígidas, longitudes, ángulos y cierre.
- [x] C.2 Añadir QTests de benceno sustituido, clic corto, arrastre, preview/commit, carbono carbonílico y cadena completa.
- [x] C.3 Cubrir colisiones, `Alt`/`Shift`, ángulo/paso/longitud fijos, aromaticidad, Undo/Redo y persistencia.
- [x] C.4 Renderizar evidencia antes/después con renderer actual a 1×/1.25× y temas disponibles; mantener revisión visual manual separada.
- [x] C.5 Ejecutar familias de tests de dibujo, Visual Lab/proyecciones, arquitectura, persistencia, compileall, Ruff scoped, OpenSpec estricto y `git diff --check` con límites.
- [x] C.6 Completar tabla de métricas, fallos basales/nuevos y pasos de prueba en CachyOS; entregar commits locales sin push/merge/PR/tag (revisión estética manual pendiente).
