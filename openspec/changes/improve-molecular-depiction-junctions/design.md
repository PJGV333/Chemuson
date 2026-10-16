# Diseño — uniones geométricas de enlaces

## Contexto y baseline observado

La rama experimental parte de `origin/release/v0.3.0-beta.1-prep` en `9a70ee1c908a39e700098bcc0a54d464a606e24c`. El lienzo usa `BondItem.update_positions()` para producir trazos Qt vectoriales. `canvas_structure._bond_neighbor_vectors()` describe cada vecino con dirección, anchura y un punto de extremo; `_set_bond_item_join_context()` inyecta ese contexto sólo en cuñas sólidas. Las cuñas se construyen como un polígono relleno, con adaptadores de punta/base y un miter final acotado. Los enlaces rayados se dibujan como trazos transversales independientes, con posiciones y anchuras lineales.

La galería exploratoria base renderiza casos de 30°, 60°, 90°, 120°, 150° y casi colineales; simple/bold/doble/triple; wedges solid/hash y combinaciones; bifurcaciones; anillos aromáticos/fusionados; N/O/S y cargas; varios grosores y longitudes. Se observa un pequeño saliente transversal en algunas uniones hash→simple/bold donde el trazo rayado terminal coincide con el átomo compartido. Es un defecto visual confirmado para investigar, no autorización para eliminar o desplazar arbitrariamente la raya terminal.

El doble enlace se dibuja con componentes distintos: línea sigma central y línea pi desplazada, normalmente recortada hacia el interior; el doble aromático `bold` también puede tener sigma gruesa y pi fina. El contexto actual de uniones no expresa ese perfil completo. Esto es una hipótesis de insuficiencia de representación, no un defecto visual confirmado: antes de cambiarla hay que demostrar qué contorno debe conectar con la cuña y cuál debe quedar retranqueado.

## Resultado del diagnóstico y decisión de implementación

La galería reproducible confirmó que algunos hash marks terminales quedaban centrados exactamente sobre el átomo compartido (`t = 1`), dejando una pequeña prolongación transversal en uniones hash→simple/bold. La comparación antes/después muestra que centrar cada marca en uno de `n` intervalos iguales (`t = (i + 0.5) / n`) reparte el margen en ambos extremos; el ancho sigue creciendo linealmente en esos mismos centros. No se cambiaron recuento, orientación, semántica, estilos ni el trazo vecino. A escala nominal, la captura posterior elimina el saliente aislado; la aceptación final sigue pendiente de revisión visual del propietario.

La inspección de simple/doble/bold frente a cuñas no confirmó un defecto de sigma/pi. Aunque el contexto de uniones no transmite cada componente de un doble, el renderer ya dibuja sigma/pi por separado y retranquea intencionalmente la pi. No se conectaron ni recortaron esas líneas. Tampoco se alteraron contornos de cuña sólida: la matriz nueva comprueba su geometría actual como control.

## Decisiones

1. La galería y la medición preceden a los cambios de producción. El material exploratorio permanece en `/tmp`; el repositorio conserva el generador reproducible y una lámina comparativa antes/después, claro/oscuro.
2. `bond_junctions.py` contiene únicamente el cálculo puro de segmentos hash, sin Qt ni acceso al grafo. El adaptador de `BondItem` conserva su API actual y no reubica átomos.
3. No se introduce un motor genérico de perfiles sigma/pi ni recortes de intersección: la hipótesis de que el tratamiento existente causa un defecto visual no quedó confirmada. Los enlaces múltiples conservan su implementación actual de líneas independientes.
4. Los trazos `bold` siguen siendo trazos gruesos no estereoquímicos; sólo los estilos `WEDGE`/`HASHED` transmiten estereoquímica. La orientación de la cuña depende de los extremos existentes y nunca se invierte como efecto secundario.
5. No se modifica el estilo base ni los presets de documento. El gate para `ChemUSON Publication` no se cumple: no están validados todos los perfiles/exportaciones ni las preferencias/documentos, por lo que se aplaza.

## Invariantes

- Snapshot exacto de IDs, elementos, coordenadas, conexiones, orden, aromaticidad, carga y metadatos estéreo antes/después de refrescar items.
- Segmentos hash deterministas, uniformemente espaciados y de anchura monótona; la API existente y la orientación de extremos se conservan.
- Renderizado determinista para iguales datos, estilo, escala, tema y destino.
- Cuñas válidas: orden de esquinas estable, polígono no autointersecado, base finita y relación anchura/longitud acotada.
- Ningún extremo invade el `shape()` de etiqueta/carga/H; uniones locales no crean salientes mayores que el grosor renderizado involucrado.
- Zoom y exportación no dependen del tamaño de captura ni de imágenes rasterizadas.

## Riesgos y límites

- Cambiar el área de unión de una cuña puede ocultar o aparentar invertir estereoquímica; cualquier duda bloquea el parche.
- Los anillos y dobles enlaces poseen líneas pi intencionalmente retranqueadas; conectar todas las líneas al átomo sería incorrecto.
- PyQt puede diferir entre render a pantalla, QImage y PDF; se prueban rutas nativas disponibles y se deja la evaluación impresa/manual al propietario.
- La apariencia actual puede ser correcta en algunos casos: no se homogeneiza la silueta de todos los encuentros.
