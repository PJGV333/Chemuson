# Diseño de la revisión documental

## Principio rector

**El manual solo afirma lo que se ha visto en la aplicación de `0.3.0-beta.4`.** Cuando una función existe pero está incompleta, es experimental o depende de componentes externos, se dice; cuando no se pudo verificar, no aparece. Este principio manda sobre la comodidad de la redacción.

## Decisiones

### 1. Verificación doble, código y ejecución

Los menús, atajos y paletas se comparan entre dos fuentes independientes: el constructor de la interfaz (`main_window_ui_builder.py`, `actions/*.py`, `toolbar.py`, `tool_rail.py`) y el volcado de la ventana real ejecutada sin pantalla. Cuando ambas coinciden, la afirmación entra en el manual. Si difieren, se documenta la diferencia (hallazgo C6, `Ver → Mostrar barra de herramientas`).

### 2. Un capítulo = una explicación completa

Cada función se describe **una sola vez** con su estructura didáctica fija: qué hace, dónde está, cómo se usa, ejemplo sencillo y limitación. El resto del manual enlaza con ancla interna. Se eliminan los capítulos espejo (B1–B6 del baseline).

### 3. Separar conceptos que el manual anterior mezclaba

- **Manipulación geométrica 2D** (Editar → Rotar, `Alt`+arrastre).
- **Rotación de partes** (ramas con `Ctrl+Alt+…`, fragmento con pivote).
- **Perspectiva tridimensional sobre el papel** (Rotación 3D precisa, proyección *trackball*).
- **Modelo tridimensional real** (visor OpenGL y CompChem).

Se dedican a ello los capítulos 6, 15 y la tabla comparativa 6.5, porque la confusión entre estos cuatro niveles era la causa de la mayor parte del texto erróneo anterior.

### 4. Ilustraciones

Solo capturas de widgets de la aplicación ejecutada. El lienzo se rellena por la API interna del editor —los mismos objetos que produce la interfaz al dibujar— y se recorta con `QGraphicsView.fitInView` sobre una región del papel. No hay composiciones, flechas añadidas ni texto sobreimpreso. El script queda versionado (`docs/manual-figures/capturar_figuras.py`) y `docs/manual-figures/README.md` documenta procedencia, entorno y regeneración.

Renumeración por familias (`01` ventana, `03`–`09` menús, `10`–`17` paletas, `20`–`26` lienzo, `30`–`34` contextual, `40`–`47` paneles), independiente de los números de figura del texto, que se asignan por orden de aparición.

### 5. Estructura de capítulos

28 capítulos → 22. Se gana claridad absorbiendo: §16 (barras flotantes) → cap. 7; §17 (TLC/gel) → cap. 9; §18 (paneles) → cap. 12 reducido a navegación; §22 (plantillas) → cap. 13 único; §7.4 (plantillas del menú) → enlace al cap. 13.

### 6. Hallazgos de comportamiento

C1–C10 del `baseline.md` **se documentan y no se arreglan**: arreglarlos exigiría tocar subsistemas críticos (`gui`, `clean2d`) fuera del alcance y sin OpenSpec de comportamiento. El manual describe el comportamiento real, incluida la duplicación de `Limpiar 2D` y las etiquetas en inglés de los resultados, porque ocultarlas volvería a inducir error.

## Archivos afectados

- `docs/MANUAL_USUARIO.md` (reescritura y renumeración)
- `docs/manual-figures/` (37 PNG + `capturar_figuras.py` + `README.md`)
- `openspec/changes/revise-user-manual-beta-4/`

No se tocan `src/`, `tests/`, `packaging/`, `.github/` ni `architecture/modules.yml`.

## Validación

Comprobaciones de rutas de imagen, anclas internas, correlación imagen↔leyenda, consecutividad de números de figura, ausencia de figuras huérfanas en disco y `git diff --check` (política de pruebas para documentación: sin pytest ni compileall).
