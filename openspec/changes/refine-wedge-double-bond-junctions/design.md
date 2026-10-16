# Diseño — unión editorial de wedges, dobles y etiquetas

## Estado inicial y evidencia

Rama `feat/molecular-depiction-excellence`, HEAD inicial `bd597dffc15688615aa7ceb3a519c920e499a5e6`, árbol limpio y sincronizada con `origin/feat/molecular-depiction-excellence`. La galería existente (`tools/visual/render_bond_junction_gallery.py`) se ejecutó sin cambios y guardó paneles baseline en `/tmp/chemuson-depiction-refinement-baseline/before-existing/`.

## Hipótesis de causa raíz por validar

- **Hash + doble:** `canvas_structure._set_bond_item_join_context()` sólo envía vecinos a `BondItem` para `WEDGE`; `HASHED` dibuja marcas transversales aisladas en `BondItem.update_positions()`. El contexto no describe sigma y pi por separado. Confirmar con un junction controlado antes de tocar geometría; no unir la cuña a la línea pi si ésta es un trazo secundario intencionalmente retranqueado.
- **Etiquetas + doble:** `canvas_text._label_open_direction()` suma vectores vecinos sin ponderar orden/grosor, mientras `_configure_bond_rendering()` decide orientación sigma/pi separadamente. `AtomItem._center_label()` centra por el advance/bounding box de texto, no por la silueta visible. Medir si el defecto proviene de side-bearing/altura de glyph, dirección abierta o interferencia específica de pi; no asumir un desplazamiento global.
- **Clearance:** el cálculo de label-shrink fija `1.5` para cuña sólida y usa `hash_stroke/2 + 1.0` para hash, con perfiles de trimming distintos. Medir el gap del trazo realmente pintado (incluidas caps y segmentos hash), no sólo `QGraphicsPathItem.shape()` ni el centro de cada línea.
- **Wedge + bold/Haworth:** el contexto sólido sí incluye vecinos `PLAIN`/`BOLD`, y el miter/fork de extremo ancho ya adapta esquinas. Probar el contorno final y sus intersecciones en pentágono y hexágono antes de ajustar blend, extensión o miter; preservar los casos simples y ramas no cíclicas.

## Método

1. Congelar el comportamiento actual con un caso reproducible por defecto claro/oscuro y métricas geométricas antes del cambio.
2. Añadir primero tests rojos para la causa medida; cubrir no sólo paths centrales, también trazos de pen/caps y la línea pi secundaria.
3. Aplicar el cambio más localizado al cálculo/renderer propietario, sin mover átomos y sin acoplar datos químicos nuevos.
4. Verificar controles: doble normal, wedges sin etiquetas, heteroátomos usuales, anillo sencillo, wedges previos, gallery histórica y Haworth α/β.
5. Guardar una sola comparación versionada de dimensiones idénticas; capturas intermedias a `/tmp`.

## Riesgos

Un miter que parezca cerrado puede ocultar el sigma; contar pi como borde de unión puede producir una cuña visualmente incorrecta. El centrado por contorno puede cambiar etiquetas fuera del caso objetivo. Una clearance menor puede solapar hashes con glifos a otras escalas. La extensión de wedge en una esquina bold puede formar un bulto o saliente. Cualquier regresión de lectura estéreo o de controles bloquea la aceptación.
