# Proposal — corrección de docks, edición 3D y conformeros complejos

## Propósito

Corregir tres defectos de uso observados en ChemUSON-3D sin reescribir su arquitectura: legibilidad del `MolecularViewportDock` al flotar, coherencia entre herramienta/modo de edición y gestos reales de edición, y diagnósticos robustos cuando la generación 3D de una estructura compleja no termina en una escena visible.

## Hitos independientes

- **A — Dock y estado de herramienta:** tema claro/oscuro aislado al dock 3D, reinyección al acoplar/flotar, controles contrastados, escena/contexto OpenGL y selección de herramienta preservados; ninguna herramienta aparece activa si el modo de edición está apagado.
- **B — Edición directa:** selección visual de átomos, H derivados y enlaces; construcción/cambio de elemento y orden de enlace con validación química; movimiento 3D de átomos y fragmentos con control de profundidad; handles guiados de longitud, ángulo y dihedro; preview cancelable y una transacción undoable por gesto.
- **C — Generación compleja:** instrumentación por etapas, fallos tipados, intentos ETKDG acotados y distintos, conservación de geometría válida aunque MMFF/UFF no puedan minimizarla, cancelación/timeout y publicación solo para el documento vigente. Verificar tetrandrina real cuando haya documento disponible y usar una referencia molecular exacta en regresión sin confundirla con el dibujo original del usuario.

## Invariantes

`MolGraph` y los IDs documentales siguen siendo autoridad; H derivados son metadata de escena; no se confunde una geometría 2D con una 3D; la configuración estereoquímica explícita se preserva y la no especificada no se inventa; el cálculo RDKit continúa aislado fuera del hilo UI; snapshots de pestaña/firma impiden aplicar resultados obsoletos; no cambian persistencia `.cmsn`, dependencias, paleta KDE ni otros docks.

## Evidencia inicial

La ruta actual solo escribe el fallo de generación en `CompChemDock`; el viewport puede quedar vacío con su estado previo. La generación aislada de una referencia tetrandrina con 46 átomos pesados y 52 enlaces sí devolvió 46 coordenadas en ~0.33 s en el baseline. No se encontró un documento `.cmsn` de la tetrandrina dibujada por el usuario; existe una prueba Clean2D con una estructura “tetrandrine-like”, que no es sustituto del caso real. Por ello la causa concreta del documento original permanece abierta.

## Fuera de alcance

No ejecutar la GUI durante esta campaña, no elevar privilegios, no modificar configuración del sistema, no añadir dependencias, no iniciar cambios Clean2D/ChemName ni persistencia, no copiar código Avogadro, no afirmar validación visual con pruebas offscreen y no repetir la suite completa por el SIGSEGV Qt/SIP conocido.

## Criterio de éxito

Cada hito queda en un commit local independiente en `feature/chemuson-3d`. Las regresiones reproducen comportamiento y conexión de señales, no solo presencia de widgets. Las limitaciones de aceptación gráfica y del archivo original quedan explícitas; no se publica ni integra el branch.