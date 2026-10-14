# Delta de especificación — pulido del panel lateral

## ADDED Requirements

### Requirement: La presentación científica no altera los datos

La tabla de Propiedades químicas SHALL representar subíndices de fórmula y superíndices de isótopos/cargas/aductos mediante texto Unicode legible, aplicando reglas únicamente a esos campos. Las celdas de magnitudes, porcentajes, concentraciones y desplazamientos químicos SHALL conservar el texto original. La selección/copia SHALL devolver TSV legible y químicamente equivalente sin HTML ni caracteres de marcado de presentación.

#### Scenario: Fórmula, NMR y MS
- **WHEN** se muestran `C66H75Cl2N9O24`, `C2H6O`, `^1H NMR est.`, `^13C NMR est.`, `M+` y `[M+H]+` en campos aplicables
- **THEN** se ven `C₆₆H₇₅Cl₂N₉O₂₄`, `C₂H₆O`, `¹H NMR est.`, `¹³C NMR est.`, `M⁺` y `[M+H]⁺`
- **AND** masas, logP, porcentajes, intensidades y ppm no se transforman.

#### Scenario: Copiar y seleccionar
- **WHEN** el usuario selecciona/copía texto en tema claro u oscuro
- **THEN** texto/alineación siguen siendo legibles y copiables
- **AND** el TSV conserva los datos fuente originales, sin HTML ni glifos agregados de presentación.

### Requirement: Lienzo y SidePanel se dividen con un splitter ajustable

La shell SHALL mantener intacto el rail y alojar lienzo y SidePanel en un divisor horizontal Qt. El panel SHALL iniciar con ancho cercano a 340 px, permitir al menos 300 px y crecer al arrastrar el separador a la izquierda; a la derecha SHALL reducirse hasta el mínimo. El lienzo SHALL ocupar el espacio restante y mantener área positiva en tamaños de ventana compactos y grandes. Ocultar el panel SHALL devolver todo su espacio al lienzo y mostrarlo SHALL recuperar un ancho razonable y la página seleccionada.

#### Scenario: Redimensionar en ambos sentidos
- **WHEN** se arrastra el handle del splitter hacia la izquierda y hacia la derecha
- **THEN** el panel cambia de ancho en el sentido esperado, el rail permanece fijo y el canvas usa el espacio restante.

#### Scenario: Ocultar y restaurar
- **WHEN** se oculta y vuelve a abrir el SidePanel desde Ver
- **THEN** el lienzo recupera el espacio completo y el panel vuelve con ancho utilizable, misma instancia y página, sin scroll ni pérdida de controles.

### Requirement: El SidePanel inicia oculto, con navegación explícita

Cada nueva sesión SHALL iniciar con SidePanel oculto independientemente del último valor persistido de visibilidad. SHALL conservar `ui/side_panel/active_tab` y las demás preferencias. Las siete acciones de página del menú Ver SHALL abrir el panel directamente en la página solicitada; no se eliminará ni sustituirá ninguna acción.

#### Scenario: Nueva sesión
- **GIVEN** cualquier valor previo de visibilidad y una página válida guardada
- **WHEN** se construye una ventana nueva
- **THEN** el SidePanel está oculto y la página guardada permanece seleccionada, sin reescribir otra preferencia.

#### Scenario: Navegación Ver
- **WHEN** se activa cada una de las siete acciones de página desde Ver
- **THEN** SidePanel se muestra en la página correspondiente y la acción de visibilidad refleja el estado real.

### Requirement: Onboarding no depende de la geometría de un panel oculto

En primera ejecución el onboarding SHALL conservar sus tres pasos y no SHALL mostrar el SidePanel automáticamente. El tercer paso SHALL resaltar un acceso visible al menú Ver y explicar cómo abrir el panel; SHALL NOT usar la geometría de un target oculto.

#### Scenario: Primera ejecución con panel oculto
- **WHEN** aparece el onboarding y se avanza al paso lateral
- **THEN** el panel continúa oculto, el foco corresponde al acceso visible al menú, la ventana y el rail no se desplazan y no aparece una superposición anclada a geometría inválida.

### Requirement: Las pestañas conservan etiqueta completa accesible

Cada botón de pestaña SHALL conservar su texto completo y SHALL exponer dicho texto como tooltip y nombre accesible aunque el viewport lo recorte. El comportamiento de navegación/overflow SHALL permanecer igual.

#### Scenario: Pestaña recortada
- **WHEN** se posa el cursor sobre `Prop...ades` o `Vali...ción`
- **THEN** el tooltip muestra `Propiedades` o `Validación` completas, sin cambiar el orden ni el acceso por overflow.
