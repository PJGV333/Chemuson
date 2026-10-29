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
| `01-ventana-tema-claro.png` | Figura 1 | Ventana principal de ChemUSON 0.3.0-beta.4 en tema claro, con un ciclo |
| `02-ventana-tema-oscuro.png` | Figura 2 | La misma ventana en tema oscuro. El tema se cambia con el botón de la |
| `47-paleta-comandos.png` | Figura 3 | Paleta de comandos abierta con `Ctrl+P`. Escriba una palabra (por ejem |
| `05-menu-editar.png` | Figura 4 | Menú |
| `08-menu-ver-paneles.png` | Figura 5 |  |
| `03-menu-estructura.png` | Figura 6 | Menú |
| `09-menu-estructura-plantillas.png` | Figura 7 |  |
| `04-menu-ayuda.png` | Figura 8 | Menú |
| `25-canvas-ordenes-enlace.png` | Figura 9 | Órdenes de enlace creados con la herramienta |
| `26-canvas-estilos-enlace.png` | Figura 10 | Estilos de la paleta de |
| `11-flyout-atomos.png` | Figura 11 | Paleta |
| `06-menu-editar-rotar.png` | Figura 12 |  |
| `07-menu-rotar-fragmento-pivote.png` | Figura 13 |  |
| `10-flyout-seleccion.png` | Figura 14 | Paleta |
| `13-flyout-corchetes.png` | Figura 15 | Paleta |
| `12-flyout-flechas.png` | Figura 16 | Paleta |
| `20-canvas-flechas-corchetes.png` | Figura 17 | Ejemplo de anotaciones sobre un ciclohexano: corchete |
| `14-flyout-simbolos.png` | Figura 18 | Paleta |
| `22-canvas-simbolos-cargas.png` | Figura 19 | Ejemplo con la herramienta |
| `16-flyout-energia.png` | Figura 20 | Paleta |
| `23-canvas-diagramas-energia.png` | Figura 21 | Diagramas creados con la herramienta |
| `17-flyout-orbitales.png` | Figura 22 | Paleta |
| `24-canvas-orbitales.png` | Figura 23 | Doce glifos de la paleta de |
| `15-flyout-placas.png` | Figura 24 | Paleta |
| `21-canvas-placas.png` | Figura 25 | Placa TLC (izquierda) y gel de electroforesis (derecha) insertados con |
| `30-contextual-atomo.png` | Figura 26 | Clic derecho sobre un átomo: |
| `31-contextual-enlace.png` | Figura 27 | Clic derecho sobre un enlace: además de las acciones básicas ofrece |
| `34-contextual-seleccion.png` | Figura 28 | Con varios objetos seleccionados se añade |
| `32-contextual-lienzo-vacio.png` | Figura 29 | Sin objeto bajo el puntero, el menú se reduce a |
| `33-contextual-analisis-submenu.png` | Figura 30 |  |
| `40-panel-inspector.png` | Figura 31 | Vista |
| `43-panel-plantillas.png` | Figura 32 | Vista |
| `42-panel-validacion.png` | Figura 33 | Vista |
| `41-panel-propiedades.png` | Figura 34 | Vista |
| `45-panel-espectroscopia.png` | Figura 35 | Vista |
| `46-panel-compchem.png` | Figura 36 | Vista |
| `44-panel-apariencia.png` | Figura 37 | Vista |

## Convención de nombres

`NN-descripcion.png`, donde `NN` agrupa las figuras por tipo:

- `01`–`02` ventana principal (claro/oscuro)
- `03`–`09` menús de la barra superior
- `10`–`17` paletas del *rail*
- `20`–`26` objetos del lienzo
- `30`–`34` menú contextual
- `40`–`47` panel lateral y paleta de comandos

Los números de figura del manual son independientes y se asignan en el orden en que aparecen las imágenes en el texto.
