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
- Lineales con insaturaciones múltiples: `-diene`, `-triene`, `-diyne` y `en-yne`.
- Lineales con un grupo funcional: `-ol`, `-amine`, `-al`, `-one`, `-oic acid`, `-nitrile` (uno solo).
- Cicloalcanos simples (3–8) con sustituyentes simples (halógenos, alquilos lineales).
- Benceno con sustituyentes simples (halógenos, alquilos lineales, hydroxy, nitro).
- Sustituyentes simples en cadena/anillo: alkoxy (C1–C3), nitro, amino.
- Heteroaromáticos monocíclicos simples (furan, thiophene, pyrrole, pyridine).
- Heteroaromáticos con 2 heteroátomos: diazinas (pyridazine/pyrimidine/pyrazine) y imidazole/oxazole/thiazole.
- Heteroaromáticos con 3 N: triazine y triazole.
- Naftaleno con mono/di-sustituidos simples (halógenos, alquilos lineales).
- Sustituyentes alquilo ramificados simples (isopropyl, sec-butyl, tert-butyl, neopentyl).
- Anillos como sustituyentes: phenyl, benzyl, cyclohexyl.
- Aromáticos fusionados: anthracene, phenanthrene (mono-sustituidos) y pyrene (mono/di-sustituidos).
- Hetero-fusionados: quinoline, isoquinoline, indole (mono-sustituidos).

**No soportado aún**
- Anillos fusionados complejos (fuera de naphthalene/anthracene/phenanthrene).
- Hetero-fusionados fuera de quinoline/isoquinoline/indole.
- Heterociclos no aromáticos.
- Heteroaromáticos fuera del subset listado.
- Múltiples insaturaciones.
- Múltiples grupos funcionales (dioles, diaminas, etc.).
- Sustituyentes ramificados fuera del set soportado o con insaturación.
- Heteroaromáticos con 3+ heteroátomos.
- Fusionados mayores (tres+ anillos aromáticos) con di-sustitución.
- Numeración alternativa para fused systems fuera de pyrene.

**Pyrene (numeración)**
- Por defecto: IUPAC 2004.
- Alternativa: CAS. Usa `NameOptions(fused_numbering_scheme="cas")`.

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
- `c1ccc2ccccc2n1` → `quinoline`
- `c1ccc2ccccc2[nH]1` → `indole`
