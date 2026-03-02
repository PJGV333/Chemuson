# RDKit en Chemuson (modo aislado)

Chemuson usa RDKit de forma opcional para import/export, validación y estereoquímica.

## Recomendación de instalación

Usar `conda`/`micromamba` con canal `conda-forge` para reducir problemas binarios:

```bash
micromamba create -n chemuson-rdkit -f packaging/conda/environment.yml
micromamba activate chemuson-rdkit
```

## Modo aislado para estereoquímica

El motor de nombres soporta extracción estereo en subproceso para evitar que un fallo nativo de RDKit tumbe la app:

- API: `NameOptions(rdkit_isolated=True)` (default actual).
- Implementación: `chemuson.chemio.rdkit_safe` + worker `chemuson.chemio._rdkit_worker`.
- Si el worker falla por código de salida, señal o timeout, el motor degrada a fallback sin crashear.

## Smoke test recomendado

```bash
pytest -q tests/test_rdkit_smoke.py
```

Si RDKit no está disponible o falla al importar en subproceso, el test se marca como `skip`.
