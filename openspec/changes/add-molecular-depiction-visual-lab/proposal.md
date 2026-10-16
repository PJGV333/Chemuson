# Propuesta — ChemUSON Molecular Depiction Visual Lab

## Por qué

ChemUSON ya dispone de un renderer Qt de enlaces, átomos y etiquetas, además de galerías deterministas para cuñas y uniones. Se necesita una herramienta de desarrollo aislada para comparar ese renderer sin modificar la geometría química, explorar ajustes con feedback inmediato, mantener perfiles reproducibles y entregar evidencia para revisión visual del propietario.

## Qué cambia

- Añadir `tools/visual/depiction_lab/` como aplicación Qt independiente que presenta el `ChemusonCanvas`, `AtomItem` y `BondItem` de producción en dos paneles sincronizados.
- Añadir un modelo JSON validado de parámetros de representación, perfil Current de solo lectura y perfil Working editable. Los valores Current reproducen exactamente `CHEMDOODLE_LIKE`; toda extensión experimental tiene semántica opcional y el valor `None` conserva la ruta actual.
- Añadir una galería con al menos doce casos estables, comparación, zoom/orientación, temas, diagnósticos de geometría, carga/guardado de perfiles y exportación reproducible.
- Añadir pruebas cortas, aisladas y offscreen para el modelo, renderer real, interfaz, fixtures, invariantes químicos y exportación.
- Documentar las limitaciones de calibración y un camino futuro de promoción del perfil aprobado, sin modificar las preferencias ni el estilo del editor principal.

## Invariantes y fuera de alcance

No se cambian Clean2D, ChemName, ChemIO, persistencia, conectividad, IDs, coordenadas almacenadas, órdenes de enlace, aromaticidad, cargas ni metadatos estéreo. La orientación de cuña y el grafo se conservan al cambiar parámetros. No se aplican perfiles al editor de producción, no se añade dependencia externa y no se reorganiza la GUI. Las escenas se construyen desde fixtures deterministas; no se ejecuta Clean2D.

## Aceptación

La herramienta debe arrancar con un comando documentado; Current y Working deben renderizar los mismos fixtures mediante el motor de producción; los controles habilitados deben modificar propiedades efectivas, y los opcionales no calibrables se presentan como no disponibles con su motivo. Perfiles, exportaciones y resumen deben permitir reproducir cada resultado. La aceptación visual final no forma parte de esta propuesta.
