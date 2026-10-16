# Propuesta — cierre conjunto de cabecera beta y Clean2D con cuñas

## Motivo

Las ramas independientes `fix/beta-save-shortcut-menu-access` y `clean2d/wedge-bond-campaign` se consolidan aquí sobre `release/v0.3.0-beta.1-prep` sin integrarlas a beta. La prueba simultánea confirmó que el shell clásico y los nuevos tests de Clean2D coexisten, pero la ejecución real volvió a reproducir `Estructura 2D compleja preservada; no se redibujó globalmente` con `C[C@H](F)CO`, estructura estereoquímica documentada en la matriz beta. La firma cíclica corregida no llega a evaluar candidatos: la clasificación recibe dos bloques `stereo_center` por un único enlace wedge (`BondStereo` en ambos extremos) y activa `multiple_stereo_centers`.

La inspección confirma además que la cabecera seguía mostrando ChemUSON, el matraz no oficial y la versión dentro de la AppBar, con otra fila de acciones encima de la barra de documentos. Eso contradice la aceptación final del propietario.

## Alcance

- Mantener visible por defecto el QMenuBar tradicional, usando los seis menús existentes.
- Reducir la zona superior a dos franjas compactas: menú y documentos/acciones. Retirar de la zona de pestañas el nombre grande, el matraz y la etiqueta de versión; conservar metadatos, título de ventana y Acerca de.
- Presentar Nuevo, Abrir, Guardar y Pegar en la franja de documentos/acciones mediante las QAction originales; ocultar el botón `+` redundante en la AppBar real. Undo/Redo siguen únicamente a la derecha. No duplicar el elemento Pegar en Editar.
- Preservar Ctrl+S, cierre Save/Discard/Cancel, X y punto dirty, temas, panel lateral oculto al inicio y redimensionable.
- Diagnosticar una molécula con estereoquímica especificada de referencia, centros con wedge/hash, moléculas desconectadas y selección total/parcial. Registrar selección efectiva, complejidad/política, candidatos/fuentes/rechazos, candidato elegido, operación, coordenadas y resultado.
- Corregir solamente la representación duplicada de un centro de cuña en el conteo estructural de centros Clean2D. El cambio de política es indispensable: un enlace dirigido `UP/DOWN` aporta un centro en su extremo de origen, no uno por cada extremo; desduplicar por ID conserva el conteo de centros distintos. No rebajar los gates de estéreo ni forzar candidatos.
- Reforzar la prueba de cuñas con referencia previa y comparación R/S independiente/round-trip de dos enantiómeros cuando RDKit esté disponible. Incluir un control con cuña dibujada que no constituye un centro R/S asignable.
- El controlador GUI verifica R/S visible en el worker RDKit aislado antes de aplicar candidatos sobre cuñas/hash; rechaza el fallback cacheado si cambia R/S aunque la firma cíclica se conserve. No altera ranking ni fuerza un candidato.

## Fuera de alcance

No cambiar ChemName, ChemIO de producción, persistencia `.cmsn`, conectividad, identidad molecular, reglas de cuñas, renderizado químico, ranking de candidatos ni políticas de preservación distintas del conteo defectuoso demostrado; el gate GUI R/S sólo rechaza candidatos que la asignación independiente demuestra que son inseguros. No integrar a `release/v0.3.0-beta.1-prep` o `main`, crear tags/releases ni publicar paquetes oficiales. No declarar aceptación beta; falta prueba manual del propietario.
