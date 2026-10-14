# Propuesta — pulido final del panel lateral beta

## Por qué

La aceptación final de `v0.3.0-beta.1` identificó cuatro defectos de presentación/uso en el panel lateral: notación científica sin jerarquía tipográfica, ancho rígido, panel visible al iniciar y nombres de pestaña inaccesibles cuando se recortan.

## Alcance

- Dar formato Unicode de subíndices/superíndices a las celdas científicas aplicables de Propiedades químicas, sin cambiar los valores originales ni el TSV.
- Hacer redimensionable el divisor horizontal entre lienzo y panel, conservando el rail, un ancho mínimo y el espacio completo del lienzo al ocultarlo.
- Iniciar el panel oculto en cada ventana/sesión, restaurar la pestaña activa y mantener los siete accesos explícitos desde Ver.
- Adaptar el tercer paso del onboarding para resaltar el acceso visible a Ver, no un panel oculto.
- Añadir tooltips con el nombre completo de cada pestaña.
- Actualizar el contrato del panel, la nota beta y la matriz manual con estados explícitos.

## Fuera de alcance

No cambiar química, cálculos de propiedades, ChemName, Clean2D, ChemIO, RDKit, Assistant, serialización, exportación ni capacidades nuevas. No publicar tags/releases ni integrar automáticamente a la rama de preparación.
