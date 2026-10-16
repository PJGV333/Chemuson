# Propuesta — robustez de depiction y dibujo molecular

## Por qué

La primera iteración de uniones corrigió únicamente las marcas terminales de cuña rayada. La revisión manual posterior detectó otras áreas que requieren diagnóstico independiente: Clean2D puede dejar un benceno dibujado a mano visualmente abierto, algunas cuñas quedan demasiado separadas de las etiquetas o no se integran bien en bifurcaciones, el dibujo se siente limitado por el borde visible del canvas y los anillos tipo Haworth necesitan una experiencia más natural.

## Qué cambia

- Reproducir primero Clean2D sobre anillos bencénicos construidos manualmente frente a plantillas; aislar candidatos, causas de rechazo y cualquier deformación antes de cambiar código.
- Ajustar cuñas sólidas/rayadas y sus uniones bifurcadas sólo tras capturar evidencia visual comparable, con pruebas de clearance y orientación.
- Permitir construir enlaces/estructuras junto al borde visible y ampliar el área navegable sin recortar ni limitar coordenadas.
- Investigar el flujo explícito de anillos/Haworth y mejorar únicamente defectos de geometría o ergonomía reproducibles.
- Mantener una galería Qt reutilizable y una sola lámina before/after versionada; mantener capturas exploratorias en `/tmp`.

## Invariantes y fuera de alcance

No cambiar identidad/conectividad del grafo, IDs, órdenes/estilos de enlaces, cargas, aromaticidad, metadatos estereoquímicos, persistencia, import/export, ni el significado de `bold`. Clean2D sólo puede aceptar coordenadas seguras y undoables; ante un candidato inseguro debe conservar el estado inicial. No alterar mixins, orden de despacho de Qt o interacciones fuera de los casos documentados. Sin dependencias nuevas, cambios de ChemName/ChemIO, preset Publication, merges a `main`/`beta-prep`, tags ni releases.

## Aceptación

Cada corrección requiere reproducción base, prueba de regresión y evidencia before/after a escala comparable. La cobertura automatizada no sustituye la revisión visual manual del propietario.
