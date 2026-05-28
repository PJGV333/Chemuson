# Arquitectura 3D y Química Computacional

ChemUSON mantiene Python/PyQt6 como aplicación principal. El canvas 2D, herramientas, docks, empaquetado y controladores siguen siendo la base de la app. La capa 3D se añade como módulo especializado para no acoplar coordenadas 3D al modelo 2D principal.

## Módulos

`chemuson.geometry3d.model` define los contratos puros:
- `CoordinateSet3D`: coordenadas asociadas a IDs de átomo de `MolGraph`.
- `OptimizationSettings`: force field, iteraciones, timeout y pasos por actualización.
- `OptimizationFrame`: frame intermedio para una UI o visor futuro.
- `OptimizationResult`: resultado final con energía, convergencia y mensaje.
- `ForceField`: `UFF`, `MMFF94`, `MMFF94s`, `GAFF`, `Ghemical`.
- `SceneMolecule3D`: escena desacoplada para render 3D, con átomos/enlaces listos para visor y soporte de átomos generados.

`MolGraph` sigue siendo el modelo 2D principal y editable. Los hidrógenos o átomos generados por RDKit/Open Babel pueden existir en `SceneMolecule3D` sin mutar el documento.

`chemuson.geometry3d.service` conserva la API existente:
- `conformer_3d_for_graph`
- `conformer_3d_for_graph_async`
- `project_conformer_to_2d`

El servicio usa `CoordinateSet3D` internamente y mantiene compatibilidad con `Conformer3DResult` para no romper llamadas existentes.

## Backends

### RDKit

`chemuson.geometry3d.rdkit_backend` usa el worker aislado en `chemuson.chemio._rdkit_worker`. No importa RDKit en el proceso principal.

Soporta:
- generación de conformación vía `graph_conformer3d`;
- optimización vía `graph_optimize3d`;
- force fields `UFF`, `MMFF94`, `MMFF94s`.

La optimización acepta un `CoordinateSet3D` inicial. Si se envían coordenadas, el worker crea un conformer RDKit con esas posiciones antes de minimizar; si no se envían, usa ETKDG como fallback.

`steps_per_update` ejecuta la minimización en bloques y devuelve frames con:
- `step`
- `energy`
- `converged`
- `positions`

Esto prepara una optimización interactiva estilo Avogadro sin bloquear la UI. La UI debe consumir `rdkit_backend.optimize_iter(...)` desde un worker/hilo, no desde el hilo principal.

### Open Babel

`chemuson.geometry3d.obabel_backend` es opcional y se activa si existe el ejecutable `obabel` en `PATH`.

Soporta de forma declarativa:
- `GAFF`
- `Ghemical`
- `MMFF94`
- `MMFF94s`
- `UFF`

Open Babel usa MOL como formato principal para preservar conectividad, orden de enlace, carga y coordenadas. XYZ queda como fallback si MOL falla. El backend valida que el número de átomos de salida coincida con el input e incluye `stdout`/`stderr` en metadata cuando están disponibles.

Si no está instalado, devuelve `OptimizationResult` con `ok == False` y un mensaje claro. No rompe la app ni la suite.

## Caché

`chemuson.geometry3d.cache.cache_key_for_3d(graph, settings, backend)` combina:
- hash químico base (`chemical_graph_hash`);
- backend;
- force field;
- `max_iters`;
- `steps_per_update`;
- semilla.

Esto evita mezclar resultados de UFF/MMFF/Open Babel o corridas con distinta configuración.

## Exportación

`chemuson.geometry3d.export_xyz` exporta `MolGraph + CoordinateSet3D` como XYZ.

`chemuson.compchem.exporters` prepara entradas para:
- ORCA
- Gaussian
- NWChem

Los exportadores aceptan molécula, coordenadas 3D, carga, multiplicidad, método, base, memoria, núcleos y tipo de cálculo.

## Visor Rust/wgpu Futuro

La frontera recomendada para Rust es `CoordinateSet3D` más metadatos de render:
- posiciones;
- radios/colores derivados;
- enlaces;
- selección/picking;
- mediciones de distancia, ángulo y diedro.

Python debe seguir orquestando documentos, edición 2D, generación de coordenadas y exportación. Rust/wgpu debe entrar como visor/renderizador 3D embebible o proceso especializado, no como reemplazo de la GUI principal.
