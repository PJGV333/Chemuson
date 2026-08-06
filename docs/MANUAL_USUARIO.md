# Manual de usuario de ChemUSON

**Documento:** manual de usuario de la interfaz gráfica  
**Aplicación:** ChemUSON  
**Estado del proyecto:** desarrollo activo (`0.3.0-dev`)  
**Revisión del manual:** 6 de agosto de 2026  
**Código auditado:** rama `main`, commit `66d2b43c8b677f8d05477156ba453ec182e6e51a`

> ChemUSON es un editor molecular 2D orientado a docencia, investigación y comunicación científica. Este manual describe la interfaz realmente implementada en la revisión indicada. Algunas etiquetas aún aparecen en inglés y el menú **Reacción** contiene por ahora únicamente una entrada deshabilitada, **Próximamente**.

---

## Contenido

1. [Instalación y ejecución](#1-instalación-y-ejecución)
2. [Organización de la interfaz](#2-organización-de-la-interfaz)
3. [Flujo básico de trabajo](#3-flujo-básico-de-trabajo)
4. [Menú Archivo](#4-menú-archivo)
5. [Menú Editar](#5-menú-editar)
6. [Menú Ver](#6-menú-ver)
7. [Menú Estructura](#7-menú-estructura)
8. [Menú Reacción](#8-menú-reacción)
9. [Menú Ayuda](#9-menú-ayuda)
10. [Barra principal](#10-barra-principal)
11. [Barra izquierda de dibujo](#11-barra-izquierda-de-dibujo)
12. [Barra derecha de símbolos científicos](#12-barra-derecha-de-símbolos-científicos)
13. [Barra de formato de texto](#13-barra-de-formato-de-texto)
14. [Selección, transformación y navegación](#14-selección-transformación-y-navegación)
15. [Dibujo de átomos, enlaces, anillos y cadenas](#15-dibujo-de-átomos-enlaces-anillos-y-cadenas)
16. [Flechas, corchetes, símbolos, orbitales y diagramas](#16-flechas-corchetes-símbolos-orbitales-y-diagramas)
17. [Placas TLC y geles de electroforesis](#17-placas-tlc-y-geles-de-electroforesis)
18. [Paneles laterales](#18-paneles-laterales)
19. [Menú contextual](#19-menú-contextual)
20. [Análisis químico, nomenclatura y validación](#20-análisis-químico-nomenclatura-y-validación)
21. [Química computacional 3D](#21-química-computacional-3d)
22. [Plantillas](#22-plantillas)
23. [Archivos, formatos y exportación](#23-archivos-formatos-y-exportación)
24. [Guardado automático y recuperación](#24-guardado-automático-y-recuperación)
25. [Preferencias y actualizaciones](#25-preferencias-y-actualizaciones)
26. [Referencia completa de atajos](#26-referencia-completa-de-atajos)
27. [Limitaciones actuales](#27-limitaciones-actuales)
28. [Resolución de problemas](#28-resolución-de-problemas)

---

## 1. Instalación y ejecución

### 1.1 Linux mediante Flatpak

Instalación estable para el usuario actual:

```bash
sudo pacman -S flatpak
flatpak --user remote-add --if-not-exists flathub https://dl.flathub.org/repo/flathub.flatpakrepo
flatpak install --user https://pjgv333.github.io/Chemuson/flatpak/stable/Chemuson-stable.flatpakref
flatpak run io.github.PJGV333.Chemuson
```

Canal beta:

```bash
flatpak install --user https://pjgv333.github.io/Chemuson/flatpak/beta/Chemuson-beta.flatpakref
```

Actualización de una instalación Flatpak:

```bash
flatpak update
```

### 1.2 Linux mediante AppImage

```bash
chmod +x Chemuson-vX.Y.Z-linux-x86_64.AppImage
./Chemuson-vX.Y.Z-linux-x86_64.AppImage
```

### 1.3 Windows

ChemUSON puede distribuirse como ejecutable portable o mediante instalador. Abra el archivo correspondiente desde el Explorador de archivos. El actualizador de una instalación completa puede aplicar un nuevo instalador al cerrar la aplicación.

### 1.4 Entorno de desarrollo

```bash
python3 -m venv chemuson
./chemuson/bin/pip install -r requirements.txt
./chemuson/bin/pip install -e .
./chemuson/bin/chemuson
```

Consultar la versión desde terminal:

```bash
chemuson --version
```

---

## 2. Organización de la interfaz

La ventana principal se divide en estas zonas:

- **Barra de menús:** Archivo, Editar, Ver, Estructura, Reacción y Ayuda.
- **Barra principal superior:** accesos rápidos a archivo, historial, transformaciones, limpieza 2D y SMILES.
- **Barra izquierda de dibujo:** selección, enlaces, cadenas, anillos, elementos, centros de coordinación y rotación 3D precisa.
- **Barra derecha de símbolos:** texto, corchetes, flechas, placas, símbolos químicos, diagramas electrónicos y orbitales.
- **Barra de formato de texto:** fuente, tamaño, estilos, alineación, color y opacidad.
- **Lienzo:** área de dibujo, selección y manipulación.
- **Paneles acoplables:** Plantillas, Inspector, Validación, Propiedades químicas, Espectros, 3D / CompChem y Apariencia.
- **Barra de estado:** informa la herramienta activa y, cuando procede, el nombre químico calculado o `N/D`.

Los paneles se pueden mostrar u ocultar desde **Ver**. La barra de símbolos también dispone de su propia acción de visibilidad.

---

## 3. Flujo básico de trabajo

1. Cree un documento con **Archivo > Nuevo** (`Ctrl+N`).
2. Elija una herramienta en la barra izquierda o derecha.
3. Dibuje en el lienzo mediante clic, arrastre o la secuencia indicada por cada herramienta.
4. Regrese a **Seleccionar** para mover, escalar, rotar o editar objetos.
5. Revise la estructura con **Estructura > Validar valencias** (`Ctrl+Shift+V`).
6. Aplique **Limpiar 2D** si necesita normalizar la geometría.
7. Guarde el documento editable con **Archivo > Guardar** (`Ctrl+S` o `Ctrl+G`).
8. Genere una figura mediante **Archivo > Exportar como**.

---

## 4. Menú Archivo

| Comando | Atajo | Función |
|---|---:|---|
| **Nuevo** | `Ctrl+N` | Crea un documento vacío. |
| **Abrir...** | `Ctrl+O` | Abre un documento o formato químico reconocido. |
| **Guardar** | `Ctrl+S`, `Ctrl+G` | Guarda el documento actual. Si aún no tiene ruta, solicita un nombre y ubicación. |
| **Archivos recientes** | - | Muestra documentos abiertos recientemente. |
| **Centro de recuperación...** | - | Abre la interfaz para revisar autosaves y sesiones recuperables. |
| **Exportar como > PNG...** | - | Exporta una imagen rasterizada. |
| **Exportar como > SVG...** | - | Exporta una imagen vectorial. |
| **Exportar como > PDF...** | - | Exporta un documento gráfico PDF. |
| **Exportar como > CML...** | - | Exporta información química en Chemical Markup Language. |
| **Salir** | `Ctrl+Q` | Cierra ChemUSON; puede solicitar guardar cambios. |

### 4.1 Consideraciones de guardado

- El formato nativo editable de ChemUSON usa la extensión `.cmsn`.
- El guardado manual limpia autosaves obsoletos asociados con el documento.
- ChemUSON conserva respaldos rotativos recientes.
- Antes de cerrar un documento modificado, atienda cualquier solicitud de guardado.

---

## 5. Menú Editar

| Comando | Atajo | Función |
|---|---:|---|
| **Deshacer** | `Ctrl+Z` | Revierte la última operación registrada. |
| **Rehacer** | `Ctrl+Y` o equivalente de la plataforma | Reaplica una operación deshecha. |
| **Cortar** | `Ctrl+X` | Copia y elimina la selección. |
| **Copiar** | `Ctrl+C` | Copia la selección al portapapeles. |
| **Pegar** | `Ctrl+V` | Inserta contenido químico, texto o imagen compatible. |
| **Duplicar** | `Ctrl+D` | Duplica la selección con desplazamiento. |
| **Eliminar** | `Delete` | Elimina los objetos seleccionados. |
| **Copiar como > SMILES** | - | Copia la selección como SMILES. |
| **Copiar como > Molfile** | - | Copia la selección como Molfile. |
| **Copiar como > CML** | - | Copia la selección como CML. |
| **Copiar como > InChI** | - | Copia la selección como InChI. |
| **Edit Electronic Diagram...** | - | Edita un diagrama electrónico semántico seleccionado. |
| **Rotar > Girar 90° a la izquierda** | - | Rota la selección 90° en sentido antihorario. |
| **Rotar > Girar 90° a la derecha** | - | Rota la selección 90° en sentido horario. |
| **Rotar > Giro 180° horizontal** | - | Refleja horizontalmente la selección. |
| **Rotar > Giro 180° vertical** | - | Refleja verticalmente la selección. |
| **Rotar > Girar rama -60°** | `Ctrl+Alt+Left` | Rota el lado móvil de un enlace seleccionado. |
| **Rotar > Girar rama +60°** | `Ctrl+Alt+Right` | Rota el lado móvil en sentido contrario. |
| **Rotar > Invertir rama (180°)** | `Ctrl+Alt+I` | Invierte la rama elegida. |
| **Rotar > Autoacomodar rama** | `Ctrl+Alt+A` | Busca una disposición menos congestionada. |
| **Rotar > Fragmento con pivote > Definir átomo pivote** | - | Usa el átomo seleccionado como pivote. |
| **Rotar > Fragmento con pivote > Limpiar átomo pivote** | - | Quita el pivote actual. |
| **Rotar > Fragmento con pivote > Girar fragmento -30°** | - | Rota el fragmento alrededor del pivote. |
| **Rotar > Fragmento con pivote > Girar fragmento +30°** | - | Rota el fragmento en sentido contrario. |
| **Rotar > Fragmento con pivote > Invertir fragmento (180°)** | - | Invierte el fragmento alrededor del pivote. |
| **Redimensionar selección...** | `Ctrl+Alt+S` | Abre el diálogo de escala de la selección. |
| **Grosor... > Aumentar** | `Ctrl+Shift+Up` | Aumenta el grosor de enlaces, flechas o corchetes seleccionados. |
| **Grosor... > Reducir** | `Ctrl+Shift+Down` | Reduce el grosor. |
| **Grosor... > Restablecer** | `Ctrl+Shift+0` | Recupera el grosor predeterminado. |
| **Preferencias...** | - | Abre las preferencias de aplicación. |

> En la revisión documentada, **Pegar** puede aparecer dos veces en el menú por la construcción actual de la interfaz; ambas entradas ejecutan la misma acción.

---

## 6. Menú Ver

### 6.1 Visualización química

- **Mostrar carbonos:** alterna las etiquetas de carbono implícitas.
- **Mostrar hidrógenos:** alterna la representación de hidrógenos implícitos.
- **Aromáticos como círculos:** cambia la representación visual de anillos aromáticos compatibles.
- **Dimensiones del dibujo...:** permite ajustar parámetros globales de dibujo, incluida la longitud de enlace.

### 6.2 Zoom y navegación

| Comando | Atajo | Función |
|---|---:|---|
| **Zoom +** | `Ctrl++` estándar | Aumenta la escala del lienzo. |
| **Zoom -** | `Ctrl+-` estándar | Reduce la escala. |
| **Zoom 100%** | `Ctrl+0` | Restablece el zoom. |
| **Mostrar copiar/pegar/zoom en barra superior** | - | Añade o retira esos botones auxiliares de la barra principal. |
| **Reglas** | - | Muestra u oculta reglas horizontal y vertical. |
| **Cuadrícula** | - | Muestra u oculta la cuadrícula del lienzo. |

La rueda del ratón también controla el zoom. `Shift+Rueda` desplaza horizontalmente.

### 6.3 Numeración

El submenú **Numeración** contiene:

- **Mostrar numeración:** activa o desactiva las etiquetas numéricas.
- **Numerar átomos:** enumera átomos.
- **Numerar estructuras:** enumera componentes moleculares.
- **Numerar ambos:** combina los dos modos.
- **Recalcular numeración:** vuelve a generar el orden actual.
- **Incluir numeración en exportación:** decide si las etiquetas se incorporan a PNG, SVG, PDF u otra salida gráfica.

La numeración es visual y no altera la identidad química del modelo.

### 6.4 Tamaño del lienzo

Presets disponibles:

- Carta vertical.
- Carta horizontal.
- A4 vertical.
- A4 horizontal.
- A3 vertical.
- A3 horizontal.
- **Personalizado...** para dimensiones definidas por el usuario.

### 6.5 Paneles y barras

Desde **Ver** se alterna la visibilidad de:

- Barra de simbolismos químicos.
- Plantillas.
- Inspector.
- Validación.
- Propiedades químicas.
- Espectros.
- 3D / CompChem.
- Apariencia.

---

## 7. Menú Estructura

### 7.1 Limpieza y conformación 2D

| Comando | Atajo | Uso recomendado |
|---|---:|---|
| **Limpiar 2D** | - | Ejecuta la ruta principal de limpieza geométrica. |
| **Limpiar 2D (1 paso)** | `Ctrl+K` | Aplica una limpieza completa en una sola acción. |
| **Limpiar 2D para publicación** | `Ctrl+Shift+K` | Prioriza una presentación lista para figura o manuscrito. |
| **Proponer conformero 2D** | `Ctrl+Alt+K` | Genera una propuesta alternativa sin asumir que debe sustituir silenciosamente la geometría existente. |

La operación se aplica a la selección cuando corresponde; si no hay selección, puede actuar sobre la estructura disponible según el contexto.

### 7.2 Validación de valencias

| Comando | Atajo | Función |
|---|---:|---|
| **Validar valencias** | `Ctrl+Shift+V` | Ejecuta el diagnóstico y actualiza el panel Validación. |
| **Siguiente error de valencia** | `F8` | Selecciona el siguiente problema. |
| **Error de valencia anterior** | `Shift+F8` | Regresa al problema anterior. |

Las correcciones deterministas disponibles se aplican desde el panel **Validación** y quedan registradas en deshacer/rehacer.

### 7.3 Polímeros y grupos R

- **Definir repetición de polímero...:** asigna la etiqueta de repetición a una selección compatible.
- **Definir sustituyentes R...:** configura la información asociada con grupos R.

### 7.4 Plantillas

- **Plantillas:** abre la biblioteca organizada por categorías.
- **Guardar selección como plantilla...:** incorpora la selección actual a la biblioteca del usuario.
- **Importar biblioteca...:** importa una biblioteca de plantillas.
- **Exportar biblioteca...:** guarda la biblioteca para respaldo o intercambio.

### 7.5 Conversión y notaciones

- **Nombre a estructura...:** interpreta un nombre mediante las fuentes configuradas y genera una estructura cuando es posible.
- **Importar SMILES...:** construye una estructura a partir de SMILES.
- **Exportar SMILES...:** genera SMILES para la estructura compatible.

### 7.6 Análisis

El submenú **Análisis** contiene:

- **Nombre (SMILES):** obtiene una identificación basada en la notación disponible.
- **Fórmula química.**
- **Masa exacta.**
- **Peso molecular.**
- **m/z.**
- **Análisis elemental...**
- **Todo:** ejecuta el conjunto de análisis disponibles.

---

## 8. Menú Reacción

En esta revisión, el menú contiene únicamente **Próximamente**, deshabilitado. Las flechas de reacción, equilibrio, retrosíntesis y mecanismo sí pueden dibujarse como objetos gráficos, pero aún no constituyen un motor completo de construcción, balance o análisis de reacciones.

---

## 9. Menú Ayuda

- **Guía rápida...:** abre la introducción breve integrada.
- **Buscar actualizaciones...:** consulta el canal configurado.
- **Acerca de Chemuson...:** muestra versión e información del proyecto.

---

## 10. Barra principal

Botones visibles de forma predeterminada:

1. **Nuevo**.
2. **Abrir**.
3. **Guardar**.
4. **Deshacer**.
5. **Rehacer**.
6. **Girar 90° a la izquierda**.
7. **Girar 90° a la derecha**.
8. **Reflejar horizontalmente**.
9. **Reflejar verticalmente**.
10. **Limpiar 2D**.
11. **Dibujar desde SMILES...**.

Al activar **Ver > Mostrar copiar/pegar/zoom en barra superior**, se agregan:

- Copiar.
- Pegar.
- Zoom +.
- Zoom -.

---

## 11. Barra izquierda de dibujo

### 11.1 Selección

El botón desplegable ofrece:

- **Seleccionar:** selección rectangular y manipulación convencional.
- **Selección libre:** lazo para contornos no rectangulares.

### 11.2 Enlaces

La paleta incluye:

| Tipo | Comportamiento |
|---|---|
| **Enlace sencillo (incremental)** | Al hacer clic sobre un enlace existente, incrementa el orden; después del máximo soportado vuelve al inicial. |
| **Enlace grueso** | Línea de mayor peso visual. |
| **Enlace doble** | Fija orden 2. |
| **Enlace triple** | Fija orden 3. |
| **Enlace aromático** | Marca el enlace como aromático. |
| **Cuña sólida (wedge)** | Estereoquímica hacia el observador. |
| **Cuña rayada (hashed wedge)** | Estereoquímica alejándose del observador. |
| **Ondulado (wavy)** | Estereoquímica indefinida. |
| **Flexible** | Enlace curvo; `Shift` permite la variante tipo S durante el ajuste. |
| **Intermolecular** | Interacción no covalente o asociación gráfica. |
| **Coordinativo** | Enlace de coordinación. |

### 11.3 Cadena

**Cadena** permite crear una secuencia de enlaces desde un átomo existente o desde un espacio vacío. El número de segmentos se ajusta durante el gesto; también existen atajos numéricos sobre un átomo.

### 11.4 Anillos

La paleta contiene:

- Benceno aromático.
- Anillos no aromáticos de 3 a 12 miembros.
- **Tamaño personalizado...**, dentro del intervalo 3-12.

Un anillo puede insertarse libremente, fusionarse sobre un enlace o unirse a un átomo. La orientación se previsualiza antes de finalizar.

### 11.5 Elementos

Accesos directos:

`C`, `N`, `O`, `S`, `P`, `F`, `Cl`, `Br`, `I`, `H`.

**Tabla periódica...** abre el selector completo. Con la herramienta de elemento:

- Clic en espacio vacío: crea un átomo.
- Clic en un átomo: cambia su elemento.
- Clic sobre un hidrógeno implícito visible: sustituye ese hidrógeno cuando la geometría lo permite.

### 11.6 Centro de coordinación

**Centro de coordinación (esfera)** crea o convierte un átomo en centro coordinativo con representación esférica. Si el elemento activo es carbono, la creación libre usa `Pd` como valor inicial de seguridad. El estilo de la esfera se modifica desde el menú contextual.

### 11.7 Rotación 3D precisa

Selecciona una estructura o un punto de referencia y abre el control de rotación exacta. Para rotación interactiva temporal también puede usarse `Alt` con la herramienta de selección.

---

## 12. Barra derecha de símbolos científicos

### 12.1 Texto

Crea una anotación editable. Tras hacer clic en el lienzo se activa inmediatamente el cursor de texto. Salga de la edición haciendo clic fuera o cambiando de herramienta.

### 12.2 Corchetes y marcos

Paleta disponible:

- Corchetes `[]`.
- Corchete izquierdo `[`.
- Corchete derecho `]`.
- Esquinas.
- Llaves `{}`.
- Llave izquierda `{`.
- Llave derecha `}`.
- Marco rectangular.
- Marco con esquinas redondeadas.
- Paréntesis `()`.

Para crear un objeto, arrastre desde una esquina hasta la opuesta. Después puede moverlo, escalarlo, rotarlo y cambiar su grosor.

### 12.3 Flechas y líneas

Paleta completa:

- Línea.
- Línea discontinua.
- Flecha directa.
- Flecha directa abierta.
- Flecha directa discontinua.
- Flecha retro.
- Flecha retro abierta.
- Flecha retro discontinua.
- Flecha doble.
- Flecha doble abierta.
- Flecha doble discontinua.
- Equilibrio.
- Equilibrio discontinuo.
- Flecha de retrosíntesis.
- Flecha curva de dos electrones.
- Flecha curva tipo fishhook de un electrón.

Las flechas rectas se crean con dos clics: inicio y final. Las curvas usan tres etapas: inicio, final y ajuste de curvatura. Los extremos se ajustan automáticamente a átomos, enlaces, radicales o pares solitarios cercanos. Mantenga `Alt` para desactivar temporalmente ese ajuste.

### 12.4 Placas

- **Placa TLC**.
- **Gel de electroforesis**.

Consulte la sección específica para sus controles contextuales.

### 12.5 Símbolos químicos

- Carga positiva.
- Carga negativa.
- Carga alterna `+/-`.
- Signo más.
- Signo menos.
- Electrón desapareado.
- Par solitario.
- Ancla ondulada.
- Radical catión.
- Radical anión.
- Carga parcial `δ+`.
- Carga parcial `δ-`.

Las cargas positivas y negativas incrementan o disminuyen la carga formal del átomo. La carga alterna recorre `0 -> +1 -> -1 -> 0`. Los electrones y símbolos se anclan al átomo elegido o al átomo cercano cuando procede.

### 12.6 Diagramas de energía

Presets de cajas:

- Subnivel `s`.
- Nivel personalizable.
- Subnivel `p`.
- Subnivel `d`.
- Subnivel `f`.
- Híbrido `sp`.
- Híbrido `sp2`.
- Híbrido `sp3`.

El menú también ofrece constructores de:

- **Atomic diagram...**
- **Diatomic MO diagram...**
- **Ligand field diagram...**

Y bibliotecas de presets atómicos, orbitales moleculares diatómicos y campo de ligandos.

Desde el menú contextual de un diagrama se puede modificar, según el tipo:

- Etiqueta y lado de etiqueta.
- Número de cajas.
- Ocupación electrónica.
- Color de contorno, relleno, etiqueta y flechas `↑`/`↓`.
- Visibilidad del fondo, contorno o base.
- Restauración del estilo.

### 12.7 Orbitales

La paleta vectorial incluye familias `s`, `p`, `d`, `d(z²)`, varios orbitales `f`, híbridos `sp`, `sp²` y `sp³`, orbitales enlazantes `σ` y `π`, y toroides, con variantes visuales de contorno, sombreado o sólido según el glifo. El objeto insertado puede moverse, redimensionarse, rotarse y estilizarse; algunos orbitales permiten editar individualmente lóbulos, colores y opacidades desde el menú contextual.

---

## 13. Barra de formato de texto

| Control | Intervalo/atajo | Función |
|---|---:|---|
| **Familia tipográfica** | Lista de fuentes del sistema | Cambia la fuente del texto seleccionado o la predeterminada. |
| **Tamaño** | 6-144 pt | Ajusta el tamaño. |
| **Negrita** | `Ctrl+B` | Activa/desactiva negrita. |
| **Cursiva** | `Ctrl+I` | Activa/desactiva cursiva. |
| **Subrayado** | `Ctrl+U` | Activa/desactiva subrayado. |
| **Subíndice** | `Ctrl+=` | Activa subíndice y desactiva superíndice. |
| **Superíndice** | `Ctrl+Shift+=` o `Ctrl++` | Activa superíndice y desactiva subíndice. |
| **Alinear a la izquierda** | - | Alinea el bloque de texto. |
| **Centrar** | - | Centra el bloque. |
| **Justificar** | - | Justifica el bloque. |
| **Color de texto** | - | Abre el selector de color. |
| **Opacidad** | 0-100 %, pasos de 5 % | Con selección afecta solo los elementos seleccionados; sin selección se aplica al lienzo/contenido global según el contexto. |

El menú de la herramienta **Texto** también expone acciones de fuente, tamaño, estilos, subíndice, superíndice, alineación y modo de color de etiquetas por elemento o negro.

---

## 14. Selección, transformación y navegación

### 14.1 Selección básica

- Clic en un objeto: selección única.
- Arrastre en vacío con **Seleccionar**: marco rectangular.
- Herramienta **Selección libre**: lazo.
- `Ctrl+Arrastre`: fuerza selección libre con la herramienta convencional.
- `Shift+Clic` o `Ctrl+Clic`: agrega, retira o alterna objetos de la selección.
- `Ctrl+A` o `Ctrl+E`: selecciona todos los objetos visibles.
- Clic derecho: selecciona el objeto bajo el cursor y abre su menú contextual.

### 14.2 Controles visuales de selección

La caja de selección incluye controles para:

- **Mover:** arrastre el control central o el propio contenido.
- **Escalar:** arrastre el control de escala.
- **Rotar:** arrastre el control de rotación.

También puede abrir una escala numérica con `Ctrl+Alt+S`.

### 14.3 Movimiento preciso

- Flechas de dirección: desplazan 1 unidad/píxel de escena.
- `Shift+Flecha`: desplaza 10 unidades.

Se admiten átomos, textos, flechas, corchetes, diagramas, orbitales e imágenes seleccionados.

### 14.4 Panorámica y zoom

- Rueda: zoom.
- `Shift+Rueda`: desplazamiento horizontal.
- Botón central + arrastre: panorámica.
- Mantener `Space` + arrastrar con botón izquierdo: panorámica.
- `Ctrl+0`: zoom 100 %.

### 14.5 Rotación 3D

- `Alt+Arrastre` sobre una estructura con la herramienta de selección: rotación 3D interactiva.
- `Alt+0`: restablece la perspectiva de rotación 3D.
- **Rotación 3D precisa:** solicita parámetros exactos.
- `Alt+Clic` sobre un átomo, cuando no inicia una rotación, recorre la orientación manual de su etiqueta.

### 14.6 Cancelación y edición

- `Esc`: cancela la inserción pendiente de una plantilla.
- Durante el ajuste de un enlace flexible: clic izquierdo confirma; clic derecho cancela.
- `Delete` o `Backspace`: elimina selección; si no existe selección, intenta eliminar el objeto bajo el cursor.
- Al editar texto, `Delete` y `Backspace` actúan sobre el texto, no sobre el objeto completo.

---

## 15. Dibujo de átomos, enlaces, anillos y cadenas

### 15.1 Atajos numéricos contextuales

Con el cursor sobre un átomo:

- `1` a `9`: crea una cadena con la longitud indicada.
- `Shift+3` a `Shift+9`: crea un anillo del tamaño indicado unido al átomo.

Con el cursor sobre un enlace:

- `1` a `9`: fija o cambia el orden de enlace admitido por el modelo.
- `Shift+3` a `Shift+9`: crea un anillo del tamaño indicado fusionado sobre el enlace.

### 15.2 Cambio rápido de etiqueta atómica

Seleccione exactamente un átomo, sin enlaces seleccionados, y escriba una letra. ChemUSON abre el editor de etiqueta inicializado con ese carácter. Este flujo admite símbolos de elemento y etiquetas personalizadas según la validación interna.

### 15.3 Dobles enlaces

ChemUSON orienta automáticamente la línea `π`, especialmente en anillos. Para un único doble enlace seleccionado:

- `Ctrl+Alt+Rueda`: alterna manualmente el lado de la línea `π`.

La operación participa en deshacer/rehacer.

### 15.4 Enlaces y geometría

- Inicie sobre un átomo para extender la estructura.
- Inicie en vacío para crear un nuevo fragmento.
- Finalice cerca de otro átomo para conectar cuando el ajuste lo reconoce.
- Un clic sobre un enlace existente modifica su orden o estilo, según la herramienta activa.
- Los enlaces de interacción y coordinación se conservan como estilos distintos de un enlace covalente normal.

---

## 16. Flechas, corchetes, símbolos, orbitales y diagramas

### 16.1 Restricción y ajuste de flechas

Los extremos se aproximan automáticamente a objetivos químicos. Use `Alt` para colocación libre. Durante una flecha curva, mueva el puntero en la tercera etapa para ajustar la curvatura; `Shift` invierte el signo de la curvatura.

### 16.2 Objetos gráficos transformables

Flechas, corchetes, textos, diagramas de energía, diagramas semánticos, orbitales e imágenes pueden:

- seleccionarse junto con estructuras;
- moverse por arrastre o teclado;
- escalarse y rotarse mediante controles;
- copiarse, cortarse, duplicarse o eliminarse;
- participar en deshacer/rehacer.

### 16.3 Ancla ondulada

La herramienta **Ancla ondulada** necesita un átomo de referencia. Se utiliza para representar un punto de unión indeterminado o flexible asociado a la estructura.

---

## 17. Placas TLC y geles de electroforesis

### 17.1 Placa TLC

Inserte una placa desde la paleta de placas. Las manchas permanecen dentro de su carril y entre la línea base y el frente de disolvente.

Menú contextual de una mancha:

- **Cambiar color...**
- **Opacidad:** 20, 40, 60, 80 o 100 %.
- **Tamaño:** pequeña, normal o grande.
- **Mostrar Rf / Ocultar Rf.**
- **Set Rf...:** introduce un valor entre 0 y 1.
- **Eliminar mancha.**

El movimiento vertical actualiza el valor `Rf`; el arrastre queda registrado para deshacer.

### 17.2 Gel de electroforesis

El objeto de gel admite carriles y bandas. La escala vertical puede representar distancia normalizada, `log(Masa/kDa)` o masa en kDa. Las bandas se mantienen dentro del intervalo de corrida. Use el clic derecho sobre bandas o el objeto para acceder a sus controles de estilo y edición disponibles.

---

## 18. Paneles laterales

### 18.1 Plantillas

Botones:

- **Nueva categoría.**
- **Importar.**
- **Exportar.**

Doble clic o activación sobre una plantilla inicia su inserción. Menú contextual:

- Plantilla: Insertar, Renombrar plantilla, Eliminar plantilla.
- Categoría: Renombrar categoría, Eliminar categoría.

### 18.2 Inspector

Muestra información de la selección:

- Átomo: elemento, ID, carga y coordenadas.
- Enlace: orden, estilo y aromaticidad.
- Texto: fuente, tamaño y estado sub/superíndice.
- Diagrama de energía: preset, etiqueta, cajas y ocupación.
- Selección múltiple: conteo de átomos, enlaces y textos.

### 18.3 Propiedades químicas

Presenta propiedades calculadas del documento o selección compatible. Puede incluir:

- fórmula;
- masas;
- análisis elemental;
- problemas de valencia;
- `logP`;
- TPSA;
- donadores y aceptores de hidrógeno;
- enlaces rotables;
- alertas de Lipinski.

**Copiar propiedades** coloca la tabla en el portapapeles como TSV.

### 18.4 Validación

Columnas del reporte:

- objeto afectado;
- severidad;
- código;
- desglose de valencia;
- mensaje;
- sugerencia.

Botones:

- **Validar**.
- **Anterior**.
- **Siguiente**.
- **Copiar reporte**.
- **Exportar reporte...** (`TSV` o `JSON`).
- Selector de corrección.
- **Aplicar corrección**.

Seleccionar una fila enfoca el átomo afectado en el lienzo.

### 18.5 Espectros

Pestañas:

- `¹H NMR`.
- `¹³C NMR`.
- `MS`.

Botones:

- **Copiar tabla**.
- **Exportar tabla...** (`TSV` o texto).

Las predicciones muestran fuente y confianza. Seleccionar una fila de NMR puede seleccionar el átomo asociado en el lienzo.

> Las predicciones son heurísticas y no sustituyen mediciones experimentales.

### 18.6 3D / CompChem

Se detalla en la sección 21.

### 18.7 Apariencia

Control actual:

- **Bordes de enlace:** redondeados o rectos.

El cambio se aplica al estilo de dibujo activo.

---

## 19. Menú contextual

El clic derecho adapta las opciones al tipo de selección.

### 19.1 Acciones generales

- Cortar.
- Copiar.
- Pegar.
- Eliminar.
- Redimensionar selección.

### 19.2 Enlaces, flechas y corchetes

- Incrementar, disminuir o restablecer grosor.
- Para líneas/flechas rectas: igualar longitudes y establecer longitud numérica.
- Para enlaces: color personalizado y restablecer color.
- Para un enlace apto: reordenar una de las dos ramas, girar `+/-60°`, invertir 180° o autoacomodar.

### 19.3 Átomos y cargas

- Elegir átomo de unión para una etiqueta compleja.
- Carga formal: incrementar, disminuir, establecer o restablecer a cero.
- **No implícitos:** impide completar automáticamente hidrógenos/atributos implícitos para el átomo.

### 19.4 Coordinación

- Convertir un enlace normal en coordinativo o viceversa.
- Activar o desactivar la esfera de coordinación.
- Hacer la esfera transparente o visible.
- Cambiar/restablecer color.
- Mostrar/quitar relleno.
- Cambiar tamaño.
- Distribuir ligandos alrededor del centro cuando existen al menos dos ligandos compatibles.

### 19.5 Diagramas y orbitales

Según el objeto aparecen controles de etiqueta, ocupación, número de cajas, títulos, niveles, carriles, colores, opacidad, lóbulos y restauración de estilo. **Edit Electronic Diagram...** abre el constructor correspondiente para diagramas semánticos.

---

## 20. Análisis químico, nomenclatura y validación

### 20.1 IUPAC-lite

ChemUSON incluye un motor de nomenclatura de cobertura parcial. Atiende cadenas, algunas insaturaciones, sustituyentes comunes, anillos simples, varios heteroaromáticos y determinados sistemas fusionados, además de algunos casos de coordinación, carbohidratos y esteroides predefinidos.

Cuando una estructura excede la cobertura o la interpretación no es confiable, la interfaz debe degradar a `N/D` en lugar de bloquearse.

### 20.2 Nombre a estructura

**Estructura > Nombre a estructura...** puede usar fuentes offline y, cuando está habilitado y disponible, una consulta a PubChem. La operación usa caché y trabajo en segundo plano para evitar bloquear la interfaz.

### 20.3 Descriptores y RDKit aislado

Las tareas que dependen de RDKit se ejecutan preferentemente en un proceso aislado. Esto protege la interfaz frente a fallos nativos o cálculos lentos. La disponibilidad de determinados descriptores depende del backend instalado.

### 20.4 Validación

La validación reporta la valencia observada como suma de órdenes de enlace, hidrógenos asignados e hidrógenos implícitos. No toda advertencia tiene una corrección automática; el botón solo se habilita cuando existe una acción definida.

---

## 21. Química computacional 3D

El panel **3D / CompChem** ofrece:

### 21.1 Backend

- **RDKit aislado**.
- **Open Babel**, si el ejecutable `obabel` está disponible.

### 21.2 Campo de fuerza

- UFF.
- MMFF94.
- MMFF94s.

### 21.3 Botones

| Botón | Función |
|---|---|
| **Generar conformero 3D** | Construye coordenadas tridimensionales. |
| **Optimizar** | Minimiza la geometría con el backend y campo de fuerza elegidos. |
| **Reset/regenerar** | Descarta el estado calculado y vuelve a generar. |
| **Proyectar a 2D** | Sustituye, con confirmación, la disposición 2D por una proyección del conformero; la acción es deshacible. |
| **Exportar XYZ** | Guarda coordenadas cartesianas. |
| **Exportar input > ORCA** | Genera un archivo de entrada para ORCA. |
| **Exportar input > Gaussian** | Genera un archivo de entrada para Gaussian. |
| **Exportar input > NWChem** | Genera un archivo de entrada para NWChem. |

La tabla de resultados muestra paso, energía, convergencia y mensaje. La energía es dependiente del método y no debe compararse indiscriminadamente entre backends o campos de fuerza.

---

## 22. Plantillas

### 22.1 Inserción

1. Abra **Ver > Plantillas**.
2. Expanda una categoría.
3. Active o haga doble clic en una plantilla.
4. Mueva la previsualización al lienzo.
5. Haga clic en espacio vacío o sobre un átomo para insertarla/conectarla.
6. Pulse `Esc` para cancelar.

### 22.2 Crear una plantilla

1. Seleccione la estructura o fragmento.
2. Use **Estructura > Guardar selección como plantilla...**.
3. Asigne nombre y categoría.

### 22.3 Administración

Las bibliotecas pueden importarse y exportarse. El menú contextual permite renombrar o eliminar plantillas y categorías del usuario. Las plantillas base distribuidas con la aplicación pueden estar protegidas o regenerarse según la versión.

---

## 23. Archivos, formatos y exportación

### 23.1 Documento editable

- `.cmsn`: formato nativo de ChemUSON.

### 23.2 Intercambio químico

- SMILES: importación, exportación y copia.
- Molfile: apertura/intercambio y copia según la ruta usada.
- CML: exportación y copia.
- InChI: copia desde la selección compatible.

### 23.3 Exportación gráfica

- PNG: raster.
- SVG: vectorial.
- PDF: documento gráfico.

La numeración solo se incluye si está activada **Incluir numeración en exportación**.

### 23.4 Datos auxiliares

- Validación: TSV o JSON.
- Espectros: TSV o TXT.
- 3D: XYZ.
- Química computacional: entradas ORCA, Gaussian y NWChem.

### 23.5 Portapapeles

ChemUSON puede copiar objetos de dibujo y representaciones químicas. Al pegar, prioriza los datos que reconoce; también admite texto e imágenes compatibles en las rutas implementadas. Las imágenes pegadas se ajustan inicialmente para no exceder una fracción razonable del papel y luego pueden transformarse como cualquier anotación.

---

## 24. Guardado automático y recuperación

ChemUSON realiza autosaves:

- aproximadamente cada dos minutos;
- también después de unos segundos de inactividad cuando existen cambios.

Al iniciar, si hay sesiones pendientes, se ofrece recuperar o descartar cada una. El **Centro de recuperación** permite revisar estos elementos de forma explícita.

Rutas habituales:

```text
~/.chemuson/autosave/
~/.chemuson/autosave/old/
~/.chemuson/crash_logs/crash_YYYYMMDD_HHMMSS.txt
```

- `autosave/`: sesiones pendientes.
- `autosave/old/`: autosaves archivados después de recuperar o descartar.
- `crash_logs/`: registros de fallos inesperados.

El autosave complementa el guardado manual; no sustituye la creación de un archivo `.cmsn` en una ubicación conocida.

---

## 25. Preferencias y actualizaciones

Abra **Editar > Preferencias...**.

### 25.1 Nomenclatura

Preferencias implementadas en el proyecto:

- **Nombre avanzado (fase 4/6)**.
- **Usar RDKit aislado**.

### 25.2 Actualizaciones

Parámetros persistentes:

- Actualizaciones activadas/desactivadas.
- Canal `stable` o `beta`.
- Modo `notify` o `silent`.
- Intervalo de comprobación en horas.
- Fecha/hora de la última consulta.

Chequeo manual: **Ayuda > Buscar actualizaciones...**.

En Flatpak, la actualización dentro de la aplicación se delega al remoto y se aplica mediante `flatpak update`. En instalaciones Windows, el actualizador puede descargar y ejecutar el instalador adecuado al cerrar. Las descargas se verifican mediante metadatos de integridad; si falla el reemplazo, existe una ruta básica de restauración.

---

## 26. Referencia completa de atajos

### 26.1 Atajos globales y de menú

| Atajo | Acción |
|---|---|
| `Ctrl+N` | Nuevo documento. |
| `Ctrl+O` | Abrir. |
| `Ctrl+S` | Guardar. |
| `Ctrl+G` | Guardar, alternativa explícita. |
| `Ctrl+Q` | Salir. |
| `Ctrl+Z` | Deshacer. |
| `Ctrl+Y` o equivalente estándar | Rehacer. |
| `Ctrl+X` | Cortar. |
| `Ctrl+C` | Copiar. |
| `Ctrl+V` | Pegar. |
| `Ctrl+D` | Duplicar selección. |
| `Delete` / `Backspace` | Eliminar selección u objeto bajo el cursor. |
| `Ctrl+A` / `Ctrl+E` | Seleccionar todo en el lienzo. |
| `Ctrl++` | Zoom +; también es superíndice cuando el foco/estado de texto usa esa acción. |
| `Ctrl+-` | Zoom -. |
| `Ctrl+0` | Zoom 100 %. |
| `Ctrl+K` | Limpiar 2D (1 paso). |
| `Ctrl+Shift+K` | Limpiar 2D para publicación. |
| `Ctrl+Alt+K` | Proponer conformero 2D. |
| `Ctrl+Shift+V` | Validar valencias. |
| `F8` | Siguiente error de valencia. |
| `Shift+F8` | Error anterior. |
| `Ctrl+Alt+S` | Redimensionar selección. |
| `Ctrl+Shift+Up` | Aumentar grosor. |
| `Ctrl+Shift+Down` | Reducir grosor. |
| `Ctrl+Shift+0` | Restablecer grosor. |
| `Ctrl+Alt+Left` | Girar rama -60°. |
| `Ctrl+Alt+Right` | Girar rama +60°. |
| `Ctrl+Alt+I` | Invertir rama. |
| `Ctrl+Alt+A` | Autoacomodar rama. |

### 26.2 Formato de texto

| Atajo | Acción |
|---|---|
| `Ctrl+B` | Negrita. |
| `Ctrl+I` | Cursiva. |
| `Ctrl+U` | Subrayado. |
| `Ctrl+=` | Subíndice. |
| `Ctrl+Shift+=` | Superíndice. |
| `Ctrl++` | Superíndice, alternativa. |

### 26.3 Lienzo y gestos contextuales

| Gesto/atajo | Acción |
|---|---|
| Rueda | Zoom. |
| `Shift+Rueda` | Desplazamiento horizontal. |
| Botón central + arrastre | Panorámica. |
| `Space` + arrastre izquierdo | Panorámica. |
| Flechas | Mover selección 1 unidad. |
| `Shift+Flechas` | Mover selección 10 unidades. |
| `Alt+Arrastre` | Rotación 3D interactiva en selección compatible. |
| `Alt+0` | Restablecer perspectiva 3D. |
| `Alt+Clic` en átomo | Recorrer anclaje/orientación de etiqueta cuando no inicia rotación. |
| `Ctrl+Alt+Rueda` | Invertir lado de la línea π en un único doble enlace seleccionado. |
| `Ctrl+Arrastre` en vacío | Selección libre. |
| `Shift/Ctrl+Clic` | Alternar objeto en la selección. |
| `Esc` | Cancelar inserción pendiente de plantilla. |
| Letra con un átomo seleccionado | Abrir editor de etiqueta atómica. |
| `1`-`9` sobre átomo | Crear cadena de longitud indicada. |
| `Shift+3`-`Shift+9` sobre átomo | Crear anillo unido al átomo. |
| `1`-`9` sobre enlace | Cambiar orden de enlace. |
| `Shift+3`-`Shift+9` sobre enlace | Crear anillo fusionado al enlace. |
| `Alt` al colocar flecha | Desactivar snapping del extremo. |
| `Shift` al ajustar flecha curva | Invertir curvatura. |
| Clic derecho durante enlace flexible | Cancelar ajuste. |
| Clic izquierdo durante enlace flexible | Confirmar ajuste. |

> Los atajos creados mediante `QKeySequence.StandardKey` se adaptan a la plataforma. Los atajos escritos expresamente como `Ctrl+...` son los definidos literalmente por la aplicación.

---

## 27. Limitaciones actuales

- Proyecto en desarrollo activo; menús, formatos y resultados pueden cambiar.
- El menú **Reacción** no implementa aún construcción o balance químico de reacciones.
- IUPAC-lite no cubre toda la nomenclatura IUPAC ni todos los sistemas fusionados, heterociclos o casos multifuncionales.
- No existe soporte actual para CDXML, MRV o PDB.
- Algunas propiedades requieren RDKit; la optimización alternativa requiere `obabel`.
- Las predicciones de NMR y MS son estimaciones heurísticas.
- Determinadas operaciones usan subconjuntos de elementos o aproximaciones.
- ChemUSON no sustituye software regulatorio validado ni análisis experimental.
- Algunas etiquetas de interfaz permanecen en inglés, entre ellas partes de los constructores de diagramas electrónicos.

---

## 28. Resolución de problemas

### La estructura no se nombra

1. Ejecute **Validar valencias**.
2. Corrija errores deterministas desde el panel.
3. Pruebe **Limpiar 2D**.
4. Verifique que la estructura esté dentro del alcance de IUPAC-lite.
5. Active RDKit aislado en preferencias si la ruta avanzada lo requiere.

### Clean 2D produce una disposición no deseada

- Deshaga con `Ctrl+Z`.
- Pruebe **Proponer conformero 2D** para evaluar una alternativa.
- Limpie solo la selección relevante.
- Reacomode una rama concreta con el menú contextual o `Ctrl+Alt+A`.

### No aparece un backend 3D funcional

- Confirme que RDKit está instalado en el entorno de ejecución.
- Para Open Babel, compruebe:

```bash
obabel -V
```

- Revise el mensaje del panel **3D / CompChem**.

### Recuperación tras un cierre inesperado

1. Reinicie ChemUSON.
2. Use el diálogo inicial o **Archivo > Centro de recuperación...**.
3. Revise `~/.chemuson/crash_logs/` si necesita diagnosticar el fallo.

### Una figura exportada no incluye numeración

Active **Ver > Numeración > Incluir numeración en exportación** antes de exportar.

### La línea π de un doble enlace está en el lado incorrecto

Seleccione únicamente ese doble enlace y use `Ctrl+Alt+Rueda`.

### Una flecha se pega al objetivo equivocado

Mantenga `Alt` mientras coloca el extremo para desactivar el ajuste automático.

---

## Nota de mantenimiento

Este manual debe actualizarse cuando se incorporen nuevas acciones, atajos o herramientas. Las fuentes de verdad principales para una revisión futura son:

- `src/chemuson/gui/actions/`
- `src/chemuson/gui/main_window_ui_builder.py`
- `src/chemuson/gui/toolbar.py`
- `src/chemuson/gui/text_toolbar.py`
- `src/chemuson/gui/canvas/`
- `src/chemuson/gui/docks.py`
- `README.md`

La organización general se inspiró en convenciones de manuales de editores químicos consolidados, pero todas las funciones descritas aquí se verificaron contra la implementación de ChemUSON indicada al inicio.
