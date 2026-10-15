# Spec Delta — guardado, cierre de pestañas y acceso a menús

## ADDED Requirements

### Requirement: Guardar usa un solo flujo para menú y atajo

Ctrl+S y Archivo → Guardar SHALL invocar la misma `QAction` de guardar y el mismo `FileController` sobre el documento activo exclusivamente. Ctrl+G SHALL conservar su comportamiento preexistente. El QAction SHALL seguir disponible en modo compacto aunque el `QMenuBar` esté oculto; no SHALL existir un handler/atajo paralelo.

#### Scenario: Guardar documento existente desde teclado o menú
- **WHEN** el documento activo está dirty y tiene ruta nativa `.cmsn`
- **AND** el usuario invoca Ctrl+S o Archivo → Guardar
- **THEN** se actualiza ese archivo sin un diálogo de ruta, se serializa el contenido actual y el documento queda clean
- **AND** ningún documento de otra pestaña se guarda o marca clean.

#### Scenario: Guardar documento nuevo
- **WHEN** el documento activo no tiene una ruta nativa
- **THEN** Ctrl+S o Archivo → Guardar solicita una ruta nativa y sólo tras guardar con éxito fija la ruta y limpia dirty.

#### Scenario: Cancelación o error al guardar
- **WHEN** el usuario cancela el diálogo o la escritura CMSN falla
- **THEN** el documento permanece abierto y dirty, y el contenido/ruta previa no se pierde ni se marca guardado.

#### Scenario: Exportación MOL
- **WHEN** el usuario exporta MOL
- **THEN** sigue usando la acción de exportación de estructura existente, sin reemplazar el guardado nativo CMSN ni alterar ruta/dirty/autosave.

### Requirement: Cerrar pestañas conserva la decisión y el documento objetivo

El botón X SHALL estar visible y accionable para pestañas existentes y SHALL delegar al flujo de cierre/confirmación existente. Una pestaña dirty SHALL cerrarse por Save únicamente después de que el canvas objetivo quede clean; Discard requiere elección explícita; Cancel o Save cancelado/fallido SHALL conservar el documento íntegro y abierto. El resultado SHALL ser correcto también al cerrar una pestaña que no está activa. Una pestaña limpia SHALL poder cerrarse sin prompt de descarte.

#### Scenario: Guardar antes de cerrar una pestaña no activa
- **WHEN** Save se elige para una pestaña dirty no activa y su `.cmsn` se escribe correctamente
- **THEN** sólo esa pestaña se guarda y se cierra, y se restaura la pestaña que estaba activa.

#### Scenario: Descartar o cancelar
- **WHEN** el usuario elige explícitamente Discard
- **THEN** sólo se cierra la pestaña objetivo
- **WHEN** elige Cancel, o el diálogo/guardado se cancela o falla
- **THEN** la pestaña objetivo continúa abierta con datos y dirty intactos.

### Requirement: La barra de menús clásica es opcional y persistente

Ver SHALL contener una opción checkable titulada exactamente **Mostrar barra de menús**. Activarla SHALL mostrar permanentemente Archivo, Editar, Ver, Estructura, Reacción y Ayuda; desactivarla SHALL restaurar el modo compacto con acceso por hamburguesa/Alt. Ambas presentaciones SHALL reutilizar los mismos `QMenu` y `QAction`, atajos y callbacks. La preferencia SHALL persistir por el sistema de configuración UI existente. La barra tradicional es el **default de la beta** (`fix-beta-tab-close-and-toolbar-visibility`): el default anterior (compacto) se migra una sola vez mediante el marcador de versión de defaults y toda elección posterior del usuario se respeta.

#### Scenario: Cambiar modo de menús
- **WHEN** el usuario alterna la opción en Ver
- **THEN** la visibilidad de la barra cambia y los mismos menús/acciones funcionan, sin duplicar handlers ni atajos.

#### Scenario: Reiniciar aplicación
- **WHEN** se reinicia tras elegir una de las dos opciones
- **THEN** la preferencia elegida vuelve a aplicarse.

#### Scenario: Temas y tamaños
- **WHEN** se usa la barra compacta o clásica en tema claro/oscuro y en ventanas de tamaños diferentes
- **THEN** los seis menús y sus acciones siguen accesibles, y el rail de dibujo, Command Palette y estructura de documentos permanecen sin cambios.

#### Scenario: Default de la beta sobre configuración previa
- **GIVEN** una configuración escrita cuando el default era compacto (`ui/show_menu_bar=false`) y sin marcador de defaults
- **WHEN** se abre la ventana
- **THEN** la barra de menús aparece visible una vez (defaults nuevos) y la siguiente elección explícita del usuario persiste.