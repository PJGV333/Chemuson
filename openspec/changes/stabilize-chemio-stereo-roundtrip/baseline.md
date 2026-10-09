# Baseline — ChemIO stereo round-trip

Captura previa a cualquier edición de código/documentación, rama `fix/chemio-stereo-roundtrip` creada desde `release/v0.3.0-beta.1-prep`.

## Git

- Rama de origen antes de aislar: `chemname/iupac-robustness`, HEAD `054f6c7fc378975289c15e850a950a257942773b`, upstream `origin/chemname/iupac-robustness`, divergencia `0 0`.
- Beta local y remota coincidían en `6aeef19028ffd047f8a39bc8f6063ea0b57210bf`; divergencia `0 0`. Coincide con la referencia suministrada.
- Árbol limpio (`git status --short` sin salida). La rama ChemName no tenía commits locales pendientes de publicación. La rama nueva parte exactamente del SHA beta; no tiene upstream configurado aún.
- `git status --short`: sin salida.

## Comandos baseline obligatorios / limitados

- `timeout 120s python -m compileall src tests tools packaging`: exit 0.
- `timeout 300s python -m pytest --collect-only -q`: exit 0; **2120 tests collected in 1.02s**.
- Suite `python -m pytest -q`: **NO EJECUTADA**, por límite explícito de esta campaña contra suite monolítica/bloqueos Qt; no se declara aprobada.
- `timeout 120s ruff check src tests tools packaging --select F401,F811,F821,E722,E741`: exit 1 por un F401 preexistente fuera de alcance: `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3` (`math`). No se modificó.
- Entorno: Python 3.14.7, pytest 9.1.1, RDKit 2026.03.6, PyQt6 local.

## Reproducción independiente de fallos

Procesos pytest separados, cada uno con timeout externo 120 s y `QT_QPA_PLATFORM=offscreen`:

- `tests/test_smiles_stereo_import.py::test_chiral_smiles_import_creates_wedge_or_hash`: **FAIL** en `assert any(atom.stereo_cip ...)`; la cuña sí existe.
- `tests/test_smiles_stereo_import.py::test_amino_acid_chiral_smiles_import_creates_wedge_or_hash`: **FAIL** en la misma aserción; la cuña/hash sí existe.

Comparación RDKit independiente: `C[C@H](O)F` es R y su opuesto `C[C@@H](O)F` es S; `N[C@@H](C)C(=O)O` es S y `N[C@H](C)C(=O)O` es R. En la ruta aislada actual SMILES→MolBlock→ChemIO graph, el grafo contiene wedge/hash pero no `stereo_cip`; grafo→worker→SMILES produce `CC(O)F` y `CC(N)C(=O)O`, que RDKit ve sin configuración asignada. No se usa esa salida defectuosa como referencia esperada.

## Trazado de rutas

- Ruta A (SMILES→ChemIO graph): el worker RDKit conoce los tags y genera cuñas MOL. El parser fallback importa estilo/paridad en `Bond`, pero el grafo no recibe descriptores CIP/E/Z. En el deduplicador MOL se ordena el par de índices antes de volver a guardar el enlace, aunque el sentido de la cuña corresponde al orden CTAB original; `C[C@H](O)F` evidencia el caso invertido.
- Ruta B (graph→ChemIO→SMILES): `rdkit_io.molgraph_to_rdkit_with_map` ignora `stereo_cip` y estilo/paridad wedge/hash. `_graph_request_payload` envía estilo, pero `_rdkit_worker._build_mol_from_graph_payload` no lo aplica ni lee `stereo_cip`; por ello ambos enantiómeros colapsan a la misma salida aquiral. El escritor fallback SMILES tampoco escribe stereo.
- Ruta C (SMILES→ChemIO→MOL/SDF→ChemIO): el parser local reconoce MOL stereo codes `1/6/4`, pero normaliza el orden endpoint; el MOL fallback escribe siempre code `0`, perdiendo wedge/hash. Para `F/C=C/F` y `F/C=C\\F`, el bloque original contiene coordenadas 2D pero no un campo E/Z que el parser interno copie a `stereo_ez`; el worker de exportación además no recibe coordenadas, así que el SMILES de salida pierde E/Z. El formato/backends se deben probar por código/descriptor, no por aspecto.
- Ruta D (representación 2D): `smiles_depict_candidates` llama `Chem.WedgeMolBonds`; importa en grafo wedges solo para centros especificados. Tetrandrina no debe heredar cuñas: RDKit detecta 0 centros especificados, 2 potenciales `?`, y 0 enlaces E/Z explícitos en su SMILES histórico.

Evidencia de dirección MOL: bloque para `C[C@H](O)F` escribe el enlace `2  1 ... stereo=1` (centro en átomo 2); el parser lo guarda como `(1,2, WEDGE)`, invirtiendo el extremo de referencia. El aminoácido escribe `2  3 ... stereo=6` y su par ya está ordenado, por eso la dirección parece conservarse.

## Alcance de baseline

No se modificó código, reglas ChemName, Clean2D, GUI, Persistencia ni packaging. Los 2 fallos P0 son reales y reproducibles; el criterio tetrandrina requiere corregir el test, no el motor.
