# Proposal — estabilización química y del visor 3D (fases 1–4)

## Propósito

Cerrar defectos de fidelidad geométrica y química detectados durante la consolidación del editor 3D, y refinar el viewport existente sin cambiar su arquitectura. La geometría documental continúa separada de `MolGraph`; los hidrógenos creados por RDKit siguen siendo derivados efímeros.

## Alcance

- Aplicar a los hidrógenos derivados las mismas transformaciones rígidas que reciben sus fragmentos al trasladar, rotar, ajustar longitudes/ángulos/dihedros y optimizar.
- Llevar la estereoquímica explícita y soportada del grafo al worker RDKit y verificarla después de generar y relajar conformeros. Si la representación no define inequívocamente una configuración, no inventarla ni afirmar preservación.
- Preservar la semántica aromática y órdenes de enlace en las rutas RDKit, escena y render; no inferir aromaticidad indiscriminadamente.
- Compactar los controles del dock y mejorar legibilidad/contraste del viewport con cambios acotados; la aceptación visual requiere display real y puede quedar pendiente si el entorno bloquea el lanzamiento.
- Investigar por separado los SIGSEGV de `pytest` y del cierre Qt/SIP con evidencia reproducible disponible; no ocultarlos ni atribuir causa común sin prueba.

## Fuera de alcance

- Reescribir `MolGraph`, la arquitectura GUI, el sistema de persistencia o el formato `.cmsn`.
- Persistir átomos de H generados, iniciar fase 5, añadir dependencias runtime, cambiar Clean2D/ChemName, ni resolver estereoquímica avanzada no representada inequívocamente.
- Eludir restricciones de lanzamiento gráfico, instalar/modificar el sistema o usar privilegios elevados.
- Tratar el SIGSEGV conocido de la suite completa como resuelto a partir de pruebas focalizadas u offscreen.

## Módulos

- M05 `geometry3d`: transformaciones y metadatos efímeros.
- M04 `chemio`: payload y reconstrucción RDKit aislada.
- M08 `gui/viewport3d`: malla, render, dock y pruebas de ciclo de vida.
- OpenSpec y tests focalizados correspondientes. No se prevén nuevos paquetes ni dependencias.

## Criterio de éxito

Los H derivados conservan geometría respecto a su fragmento; la conformación respeta estereoquímica explícita inequívoca o rechaza con diagnóstico estable; aromaticidad y enlaces múltiples se interpretan/renderizan consistentemente; el visor ofrece espacio de trabajo útil; y los dos SIGSEGV quedan corregidos mediante regresión probada o documentados separadamente con reproducción mínima y limitaciones claras.