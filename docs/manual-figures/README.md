# Ilustraciones del manual de usuario

Carpeta de trabajo de las figuras de [`docs/MANUAL_USUARIO.md`](../MANUAL_USUARIO.md).

## Procedencia

**Todas las imágenes de esta carpeta son capturas auténticas de ChemUSON 0.3.0-beta.4.**
No hay ilustraciones dibujadas a mano, montajes ni capturas de otras versiones.

- **Código capturado:** el del tag [`v0.3.0-beta.4`](https://github.com/PJGV333/Chemuson/tree/v0.3.0-beta.4) (el código de `main` es idéntico en `src/`).
- **Cómo se obtuvieron:** se ejecutó la propia aplicación en modo sin pantalla (`QT_QPA_PLATFORM=offscreen`) y se capturaron los widgets reales —la ventana, los menús, las paletas del *rail*, el lienzo y los paneles— con `QWidget.grab()`. El lienzo se encuadra con `QGraphicsView.fitInView` y se recorta la región de interés; **no se dibuja nada a mano sobre la imagen**.
- **Contenido del lienzo:** las estructuras de ejemplo se crearon con la API interna del propio editor (`add_atom`, `add_bond`, `add_arrow_item`, `add_bracket_item`, `_insert_orbital_item`, etc.), es decir con los mismos objetos que produce la interfaz al dibujar.
- **Fecha de captura:** 28/10/2026.
- **Entorno:** PyQt6 sobre Python 3.14, Linux x86-64, sin pantalla (plataforma `offscreen`).

## Cómo regenerarlas

Desde la raíz del repositorio, con un entorno que tenga PyQt6:

```bash
QT_QPA_PLATFORM=offscreen python docs/manual-figures/capturar_figuras.py \
  ventana menu flyout contenido contextual paneles
```

Para no sobrescribir las capturas existentes, establezca `CHEMUSON_OUT` en un directorio temporal vacío antes de ejecutar el comando. Para usar otro checkout del código de beta.4, indique su ruta en `CHEMUSON_REPO`. Sin estas variables se usa el repositorio y la carpeta actuales.

El script muestra `OK nombre.png anchoxalto` por captura y devuelve un código distinto de cero ante errores registrados. Al pedir todos los grupos comprueba además que se produzcan todos los PNG del inventario versionado (actualmente 37). Compare las imágenes antes de reemplazarlas: las diferencias de fuentes, plataforma o Qt pueden alterar las capturas.

Grupos disponibles: `ventana` (ventana completa), `menu` (menús de la barra), `flyout` (paletas del rail), `contenido` (objetos del lienzo), `contextual` (menú del clic derecho) y `paneles` (vistas del panel lateral). Los nombres de grupos no se interpretan como filtros de figuras. Para ejecutar una figura registrada individualmente use `python docs/manual-figures/capturar_figuras.py fig flechas` (este registro genera los siete archivos 20–26). `fig` sin nombre ejecuta todos los registros individuales; nombres desconocidos causan error.

Las figuras 17 y 25 se corrigieron en la microcampaña: el ciclo de la figura 17 es no aromático; la figura 25 usa posiciones locales de carril para tres manchas de ejemplo y cinco bandas reales. Los valores son ilustrativos, no experimentales.

## Índice de archivos

| Archivo | Figura del manual | Contenido |
|---|---|---|
| `01-ventana-tema-claro.png` | Figura 1 | Ventana principal en tema claro: lienzo, barra de herramientas y panel lateral. |
| `02-ventana-tema-oscuro.png` | Figura 2 | Ventana principal en tema oscuro. |
| `47-paleta-comandos.png` | Figura 3 | Paleta de comandos y búsqueda rápida (`Ctrl+P`). |
| `05-menu-editar.png` | Figura 4 | Menú Editar y acceso a las transformaciones. |
| `08-menu-ver-paneles.png` | Figura 5 | Submenú Ver → Paneles para mostrar u ocultar vistas. |
| `03-menu-estructura.png` | Figura 6 | Menú Estructura y acciones disponibles en beta.4. |
| `09-menu-estructura-plantillas.png` | Figura 7 | Submenú Estructura → Plantillas. |
| `04-menu-ayuda.png` | Figura 8 | Menú Ayuda tal como aparece en beta.4. |
| `25-canvas-ordenes-enlace.png` | Figura 9 | Ejemplos de enlaces sencillos, dobles y triples. |
| `26-canvas-estilos-enlace.png` | Figura 10 | Ejemplos de estilos gráficos de enlaces. |
| `11-flyout-atomos.png` | Figura 11 | Paleta de Átomos del rail. |
| `06-menu-editar-rotar.png` | Figura 12 | Submenú Editar → Rotar. |
| `07-menu-rotar-fragmento-pivote.png` | Figura 13 | Opciones de rotación de fragmentos con pivote. |
| `10-flyout-seleccion.png` | Figura 14 | Paleta de Selección y herramienta de rotación 3D precisa. |
| `13-flyout-corchetes.png` | Figura 15 | Paleta de Corchetes. |
| `12-flyout-flechas.png` | Figura 16 | Paleta de Flechas y sus variantes. |
| `20-canvas-flechas-corchetes.png` | Figura 17 | Anotaciones de flechas, líneas y corchetes sobre un ciclohexano. |
| `14-flyout-simbolos.png` | Figura 18 | Paleta de Símbolos químicos. |
| `22-canvas-simbolos-cargas.png` | Figura 19 | Ejemplos de cargas y anotaciones químicas. |
| `16-flyout-energia.png` | Figura 20 | Paleta de Diagramas de energía. |
| `23-canvas-diagramas-energia.png` | Figura 21 | Ejemplos de diagramas energéticos y subniveles. |
| `17-flyout-orbitales.png` | Figura 22 | Paleta de Orbitales. |
| `24-canvas-orbitales.png` | Figura 23 | Ejemplos de glifos de orbitales. |
| `15-flyout-placas.png` | Figura 24 | Paleta de Placas TLC y geles. |
| `21-canvas-placas.png` | Figura 25 | Placa TLC con Rf 0,25/0,55/0,80 y gel ilustrativo de cinco bandas. |
| `30-contextual-atomo.png` | Figura 26 | Menú contextual sobre un átomo. |
| `31-contextual-enlace.png` | Figura 27 | Menú contextual sobre un enlace. |
| `34-contextual-seleccion.png` | Figura 28 | Menú contextual con varios objetos seleccionados. |
| `32-contextual-lienzo-vacio.png` | Figura 29 | Menú contextual sobre el lienzo vacío. |
| `33-contextual-analisis-submenu.png` | Figura 30 | Submenú Análisis del menú contextual. |
| `40-panel-inspector.png` | Figura 31 | Panel Inspector con datos del objeto seleccionado. |
| `43-panel-plantillas.png` | Figura 32 | Panel Plantillas y biblioteca de estructuras. |
| `42-panel-validacion.png` | Figura 33 | Panel Validación de estructuras. |
| `41-panel-propiedades.png` | Figura 34 | Panel Propiedades químicas. |
| `45-panel-espectroscopia.png` | Figura 35 | Panel Espectroscopía y resultados estimados. |
| `46-panel-compchem.png` | Figura 36 | Panel CompChem y sus opciones de cálculo. |
| `44-panel-apariencia.png` | Figura 37 | Panel Apariencia y opciones visuales. |

## Convención de nombres

`NN-descripcion.png`, donde `NN` agrupa las figuras por tipo:

- `01`–`02` ventana principal (claro/oscuro)
- `03`–`09` menús de la barra superior
- `10`–`17` paletas del *rail*
- `20`–`26` objetos del lienzo
- `30`–`34` menú contextual
- `40`–`47` panel lateral y paleta de comandos

Los números de figura del manual son independientes y se asignan en el orden en que aparecen las imágenes en el texto.
