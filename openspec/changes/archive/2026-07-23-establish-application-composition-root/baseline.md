# Baseline Report

- **Fecha**: 2026-07-23
- **Rama base**: `main`
- **HEAD base**: `7ea76b2983742bff862e4fea155b6b13b6ba1564`
- **Rama de trabajo**: `architecture/phase6-establish-application-composition-root`
- **Entorno**: Python 3.12; imagen Work sin `libEGL.so.1`

## Comandos y resultados

### `git status --short`

Sin salida antes de crear la rama y los artefactos: árbol limpio.

### `python -m compileall -q src tests tools packaging`

Exit code 0.

### `pytest --collect-only -q`

Exit code 2. Pytest alcanzó **1025 pruebas recolectadas**, pero interrumpió la
colección con **56 errores**. Todos los errores observados son
`ImportError: libEGL.so.1: cannot open shared object file` al importar PyQt6 o
módulos GUI que lo cargan transitivamente.

### `pytest -q`

Exit code 2 por los mismos 56 errores de colección relacionados con
`libEGL.so.1`. No se ejecutaron tests después de la interrupción de colección.

### `ruff check src tests tools packaging --select F401,F811,F821,E722,E741`

Exit code 1 por el único diagnóstico histórico:

```text
F401 `math` imported but unused
tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3:8
```

### `openspec validate --all --strict`

Exit code 0: **18/18 especificaciones válidas**.

## Confirmación

No se modificaron snapshots, baselines ni excepciones para ocultar resultados.
La ausencia de EGL es una limitación conocida de esta imagen y no un fallo
introducido por este cambio. La implementación deberá añadir pruebas de
bootstrap que puedan ejecutarse sin importar Qt.
