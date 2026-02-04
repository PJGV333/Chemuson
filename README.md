# Chemuson

Editor de dibujo molecular 2D libre y open source.

## Stack Técnico

- **GUI**: PyQt6
- **Química**: RDKit
- **Lenguaje**: Python 3.10+

## Instalación y Ambiente

He creado un ambiente virtual llamado `chemuson`. Si necesitas recrearlo o instalarlo en otro lugar:

```bash
python3 -m venv chemuson
./chemuson/bin/pip install -r requirements.txt
```

## Ejecución

Para iniciar la aplicación usando el ambiente virtual:

```bash
./chemuson/bin/python src/main.py
```

## IUPAC-Lite (soporte actual)

**Soportado**
- Lineales: alcanos acíclicos, halógenos, alquilos lineales.
- Lineales con una insaturación: `-ene` o `-yne`.
- Lineales con un grupo funcional: `-ol` o `-amine` (uno solo).
- Cicloalcanos simples (3–8) con sustituyentes simples (halógenos, alquilos lineales).
- Benceno con sustituyentes simples (halógenos, alquilos lineales, hydroxy, nitro).
- Heteroaromáticos monocíclicos simples (furan, thiophene, pyrrole, pyridine).
- Heteroaromáticos con 2 heteroátomos: diazinas (pyridazine/pyrimidine/pyrazine) y imidazole/oxazole/thiazole.
- Naftaleno y monosustituidos simples (halógenos, alquilos lineales).

**No soportado aún**
- Anillos fusionados o puenteados.
- Heterociclos aromáticos/no aromáticos.
- Múltiples insaturaciones.
- Múltiples grupos funcionales (dioles, diaminas, etc.).
- Sustituyentes ramificados o con insaturación.
- Heteroaromáticos con 3+ heteroátomos.
- Naftaleno polisustituido o fusionados mayores.

**Ejemplos (entrada → salida)**
- `C1CCCCC1` → `cyclohexane`
- `CC1CCCCC1` → `methylcyclohexane`
- `ClC1CCCCC1` → `chlorocyclohexane`
- `C1=CC=CC=C1` → `benzene`
- `Clc1ccccc1` → `chlorobenzene`
- `Cc1ccccc1` → `methylbenzene`
- `n1ccccc1` → `pyridine`
- `o1cccc1` → `furan`
- `s1cccc1` → `thiophene`
- `c1ccc2cccc3` → `naphthalene`
- `Br-(CH2)12-Cl` → `1-bromo-12-chlorododecane`
