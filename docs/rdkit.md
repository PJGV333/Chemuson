# RDKit en Chemuson (modo aislado)

Chemuson usa RDKit de forma opcional para import/export, validación y estereoquímica.

Por defecto, Chemuson evita importar RDKit en el proceso principal. Las extensiones nativas de RDKit pueden abortar el intérprete si fueron compiladas para otra versión de Python o si provienen de artefactos empaquetados incompatibles. La ruta recomendada es ejecutar RDKit en el worker aislado de `chemuson.chemio.rdkit_safe`.

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
- Cobertura aislada:
  - `stereo_descriptors_for_chain(...)` para `R/S` y `E/Z`.
  - `advanced_stereo_descriptors_for_chain(...)` para `M/P`, `R_a/S_a`, `endo/exo`, `si/re` (best-effort).
- Si el worker falla por código de salida, señal o timeout, el motor degrada a fallback sin crashear.

## Importación directa de RDKit

La importación directa en `chemuson.chemio.rdkit_io` está deshabilitada por defecto. Para entornos controlados donde RDKit está instalado y validado en el mismo intérprete, se puede habilitar explícitamente:

```bash
CHEMUSON_ENABLE_DIRECT_RDKIT=1 pytest -q
```

Sin esa variable, `_rdkit_available()` devuelve `False` y las rutas de UI/exportación usan workers aislados o fallbacks internos cuando corresponde.

## Pruebas

`pytest` está configurado para buscar pruebas solo en `tests/` y no recorrer artefactos de distribución como `dist/`, `dist-appimage/` o `dist-flatpak/`. Esto evita importar paquetes empaquetados incompatibles durante la suite local.

En GUI, estas rutas se controlan desde:

- `Editar -> Preferencias -> RDKit -> Nombre avanzado (fase 4/6)`
- `Editar -> Preferencias -> RDKit -> Usar RDKit aislado`

## Smoke test recomendado

```bash
pytest -q tests/test_rdkit_smoke.py
```

Si RDKit no está disponible o falla al importar en subproceso, el test se marca como `skip`.
