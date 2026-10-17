# Propuesta — Estabilizar CI y perfeccionar Molecular Depiction Visual Lab

## Por qué

La ejecución GitHub Actions `38996089270` del SHA `87e298c` agotó el límite de 300 s en shards 3 y 5, aunque los seis shards restantes completaron. Los logs muestran progreso continuo y trabajo Qt/Clean2D costoso concentrado en archivos específicos; el plan balanceó únicamente node IDs. Además, el centrado diferido de `ChemusonCanvas` puede recentrar los dos canvas aislados del Visual Lab sobre una hoja de trabajo que el laboratorio no usa. La proyección Haworth necesita coordenadas editoriales por vértice, y las geometrías Fischer/Haworth/silla deben retener sus convenciones estereoquímicas fuera de `MolGraph` cuando éste no puede almacenarlas.

## Alcance

- Balancear los ocho shards de manera determinista por tiempos históricos por archivo, conservando la colección exacta, límites finitos y terminación del grupo de procesos. Instrumentar y publicar duraciones por archivo; mantener 300 s hasta que datos nuevos justifiquen otro valor.
- Centrar el Visual Lab mediante una subclase/adaptador local al laboratorio, con centrado diferido y tras resize; no monkey-patchear ni cambiar el canvas normal.
- Refinar geometría canónica de piranosa/furanosa, con reglas direccionales de sustituyentes, aspecto/perspectiva/rotación validados y política explícita de aristas frontales; beta-D-glucopiranosa es referencia.
- Añadir semántica tipada e inmutable para Haworth α/β y D/L, Fischer horizontal/vertical y silla axial/ecuatoral-arriba/abajo. Mantenerla en un contenedor de proyección en memoria y hacer explícito que la conversión a `MolGraph`/`.cmsn` es sólo gráfica y no persiste esos descriptores.
- Añadir pruebas reales de viewport y renderer, invariantes de orientación, reportes/evidencia anterior/nueva y verificación OpenSpec.

## Fuera de alcance

No alterar los valores aprobados de `ChemUSON Editorial AB3`, su condición opt-in, el estado especial de B4 ni la exclusión de A4 de la calibración automática. No cambiar Clean2D, reglas químicas generales, `MolGraph`, la serialización `.cmsn`, ni la geometría/centrado predeterminados del editor de producción. No hacer merge.
