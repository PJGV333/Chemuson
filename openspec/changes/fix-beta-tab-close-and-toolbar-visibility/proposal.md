# Propuesta — cierre visible, barra de menús tradicional y barra de herramientas compacta

## Problema

La aceptación visual de la rama `fix/beta-save-shortcut-menu-access` encontró tres
defectos restantes en el shell de la beta:

1. Algunas pestañas se creaban con el botón Cerrar sin icono hasta el siguiente
   cambio de tema.
2. `main_toolbar` continuaba ocultándose incondicionalmente aunque existía una
   opción distinta para acciones auxiliares.
3. Tras la primera corrección, el propietario verificó que la interpretación de
   las barras no era la pedida: se mostró una barra de acciones con nueve
   botones y el menú tradicional seguía oculto por defecto. Además, al pasar una
   pestaña de limpia a modificada, el indicador dirty desplaza la X hacia el
   borde e incluso fuera del rect de la pestaña.

## Alcance

- Inicializar el icono de cierre en cada pestaña al construirla; garantizar
  contraste, legibilidad, área de interacción cómoda, tooltip y accesibilidad,
  manteniendo separado el punto dirty y reutilizando el flujo de cierre existente.
- **A1**: mostrar permanentemente, de forma compacta, la barra de menús
  tradicional con los seis menús existentes (`Archivo | Editar | Ver | Estructura
  | Reacción | Ayuda`), visible por defecto y ocultable desde
  `Ver → Mostrar barra de menús`, conservando el acceso hamburguesa/Alt cuando
  esté oculta. Sin duplicar `QMenu`, `QAction`, atajos ni controladores, y sin
  una segunda fila excesivamente alta.
- **A2**: reducir la barra de herramientas superior a tres botones pequeños,
  uniformes y theme-aware (Nuevo, Abrir, Guardar) usando las mismas QActions.
  Deshacer/Rehacer quedan representados únicamente por los botones de la app bar;
  rotaciones, volteos, Clean2D y SMILES siguen en sus rutas actuales (menús,
  atajos, rail y paleta de comandos) y no se eliminan del programa.
  `Ver → Mostrar barra de herramientas` controla exclusivamente esa barra.
- **A3**: la geometría de la pestaña no cambia al pasar de limpia a modificada:
  el espacio del indicador dirty queda reservado de forma permanente y sólo varía
  su pintura. La X permanece visible, alineada, dentro del rect de su pestaña y
  sin solaparse con el botón `+` ni con otras pestañas; el cierre sigue delegando
  al flujo Save/Discard/Cancel existente.
- Cubrir temas claro/oscuro, ventanas pequeñas, DPI elevado (métricas lógicas),
  documentos nuevos y reordenados, guardado/Undo/Redo y persistencia con tests
  Qt acotados.

## Fuera de alcance

No cambiar Ctrl+S, handlers/controladores, persistencia CMSN, `CanvasTabManager`,
flujos de cierre, Undo/Redo, rail de dibujo, panel lateral, iconos de rail,
química, Clean2D, ChemName, empaquetado, canales, workflow de release, beta,
`main` o tags. No añadir dependencias ni duplicar QActions de comandos o toolbars.
No se elimina ninguna acción del programa: sólo se retiran de la barra superior.
