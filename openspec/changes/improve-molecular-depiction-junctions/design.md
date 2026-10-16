# Diseño — uniones geométricas de enlaces

## Contexto y baseline observado

La rama experimental parte de `origin/release/v0.3.0-beta.1-prep` en `9a70ee1c908a39e700098bcc0a54d464a606e24c`. El lienzo usa `BondItem.update_positions()` para producir trazos Qt vectoriales. `canvas_structure._bond_neighbor_vectors()` describe cada vecino con dirección, anchura y un punto de extremo; `_set_bond_item_join_context()` inyecta ese contexto sólo en cuñas sólidas. Las cuñas se construyen como un polígono relleno, con adaptadores de punta/base y un miter final acotado. Los enlaces rayados se dibujan como trazos transversales independientes, con posiciones y anchuras lineales.

La galería exploratoria base renderiza casos de 30°, 60°, 90°, 120°, 150° y casi colineales; simple/bold/doble/triple; wedges solid/hash y combinaciones; bifurcaciones; anillos aromáticos/fusionados; N/O/S y cargas; varios grosores y longitudes. Se observa un pequeño saliente transversal en algunas uniones hash→simple/bold donde el trazo rayado terminal coincide con el átomo compartido. Es un defecto visual confirmado para investigar, no autorización para eliminar o desplazar arbitrariamente la raya terminal.

El doble enlace se dibuja con componentes distintos: línea sigma central y línea pi desplazada, normalmente recortada hacia el interior; el doble aromático `bold` también puede tener sigma gruesa y pi fina. El contexto actual de uniones no expresa ese perfil completo. Esto es una hipótesis de insuficiencia de representación, no un defecto visual confirmado: antes de cambiarla hay que demostrar qué contorno debe conectar con la cuña y cuál debe quedar retranqueado.

## Decisiones

1. La galería y la medición preceden a los cambios de producción. El material exploratorio se guarda en `/tmp`; el repositorio contiene sólo un generador reproducible y, si resulta útil, una lámina final única.
2. Todo cálculo nuevo vive dentro de `chemuson.gui` y recibe coordenadas/patrones de trazos inmutables; no lee ni escribe el grafo molecular. El adaptador de `BondItem` conserva su API actual y no reubica átomos.
3. Si se crea `bond_junctions.py`, representa por separado sigma, pi, contornos de cuñas y segmentos hash con geometría explícita; no colapsa doble enlace a una sola línea de grosor compuesto. Las intersecciones se acotan con límites geométricos derivados de las anchuras reales, evitando nuevos factores de compensación arbitrarios.
4. Los trazos `bold` siguen siendo trazos gruesos no estereoquímicos; sólo los estilos `WEDGE`/`HASHED` transmiten estereoquímica. La orientación de la cuña depende de los extremos existentes y nunca se invierte como efecto secundario.
5. Ningún cambio al estilo base/presets de documento. `ChemUSON Publication` queda condicionado al cierre visual de las fases anteriores y revisión de configuración/preferencias.

## Invariantes

- Snapshot exacto de IDs, elementos, coordenadas, conexiones, orden, aromaticidad, carga y metadatos estéreo antes/después de refrescar items.
- Renderizado determinista para iguales datos, estilo, escala, tema y destino.
- Cuñas válidas: orden de esquinas estable, polígono no autointersecado, base finita y relación anchura/longitud acotada.
- Ningún extremo invade el `shape()` de etiqueta/carga/H; uniones locales no crean salientes mayores que el grosor renderizado involucrado.
- Zoom y exportación no dependen del tamaño de captura ni de imágenes rasterizadas.

## Riesgos y límites

- Cambiar el área de unión de una cuña puede ocultar o aparentar invertir estereoquímica; cualquier duda bloquea el parche.
- Los anillos y dobles enlaces poseen líneas pi intencionalmente retranqueadas; conectar todas las líneas al átomo sería incorrecto.
- PyQt puede diferir entre render a pantalla, QImage y PDF; se prueban rutas nativas disponibles y se deja la evaluación impresa/manual al propietario.
- La apariencia actual puede ser correcta en algunos casos: no se homogeneiza la silueta de todos los encuentros.
