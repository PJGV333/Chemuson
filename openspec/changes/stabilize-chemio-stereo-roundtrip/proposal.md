# Propuesta — ChemIO stereo round-trip

## Por qué

La importación de SMILES quirales dibuja una cuña o un enlace discontinuo, pero los dos casos P0 observados pierden su configuración química al exportar/reimportar. El grafo interno guarda estilos de enlace dirigidos, aunque no siempre conserva la orientación de sus extremos; el exportador MOL fallback escribe paridad cero y las rutas RDKit/worker omiten las anotaciones tetraédricas del grafo. La salida puede convertirse silenciosamente en una estructura aquiral o en un estereoisómero distinto.

## Alcance

- Preservar centros tetraédricos especificados y E/Z explícitos en las rutas ChemIO SMILES y MOL/SDF soportadas.
- Mantener conectividad, identidad de átomos/enlaces, órdenes, cargas formales e isótopos.
- Corregir las expectativas de tetrandrina si su SMILES carece de estereoquímica especificada.
- Añadir una matriz de regresión ChemIO con comparación química independiente de RDKit y evidencia OpenSpec.

## Fuera de alcance

No se modifican Clean2D, ChemName, geometría/heurísticas de dibujo, GUI general, persistencia `.cmsn`, CompChem, packaging ni workflows; no se inventa estereoquímica ni se añade una dependencia runtime. No se declaran resueltos otros fallos Qt, Clean2D o CompChem.

## Módulos probables

`src/chemuson/chemio/rdkit_io.py`, `src/chemuson/chemio/rdkit_safe.py`, `src/chemuson/chemio/_rdkit_worker.py` y tests ChemIO/SMILES.
