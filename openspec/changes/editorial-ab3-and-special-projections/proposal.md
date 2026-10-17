# Propuesta — Preset editorial AB3 y proyecciones especiales

## Por qué

El propietario aprobó los valores experimentales de depiction derivados del corpus químicamente válido A1–B3, pero aún no deben convertirse en defaults de producción. El laboratorio necesita cargarlos como perfil nombrado y reproducible; B4 debe conservarse como fixture de estrés sin participar en calibración editorial. Además, la plantilla Haworth vigente produce cruces/colisiones y mezcla reglas de proyección con la construcción de grafos en un módulo GUI general.

## Alcance

- Añadir el preset experimental `ChemUSON Editorial AB3` al Molecular Depiction Visual Lab, con carga desde UI y exportación reproducible, sin alterar `CHEMDOODLE_LIKE` ni los defaults de `DrawingStyle`.
- Marcar B4 y los demás casos fuera del corpus de selección como no elegibles para calibración; conservar B4 en render y reportes con etiqueta explícita de fixture geométrico no editorial.
- Establecer constructores puros y parámetros tipados para Haworth, Fischer y silla, desacoplados del renderer Qt y de `DrawingStyle` global.
- Reconstruir Haworth con geometría de anillo, aristas frontales y sustituyentes explícitos, preservando la API heredada y sin añadir metadatos estereoquímicos al grafo.

## Fuera de alcance

No promover estilos a defaults globales, cambiar Clean2D, el editor general, persistencia ni formato `.cmsn`; no inferir estereoquímica de moléculas; no implementar UI completa para proyecciones ni expandir el corpus editorial aprobado.

## Módulos probables

`tools/visual/depiction_lab_app/` y sus pruebas; `src/chemuson/gui/templates.py`; nuevo `src/chemuson/gui/projections/` perteneciente a M08; `architecture/modules.yml`; OpenSpec y evidencia del laboratorio.
