# Validación — refinamiento de wedges y dobles enlaces

## Estado

OpenSpec y baseline inicial registrados en `bd597dffc15688615aa7ceb3a519c920e499a5e6`. Baseline completo monolítico: `pytest -q` agotó 480 s aproximadamente al 34%; se conserva ese resultado, no se declara verde. `compileall` pasó, colección 2275, Ruff selectivo conserva sólo el F401 histórico `math` documentado en `baseline.md`.

La galería original antes de cambios está en `/tmp/chemuson-depiction-refinement-baseline/before-existing/`; la matriz dedicada de 15 casos antes del cambio está en `/tmp/chemuson-depiction-refinement-baseline/focused-before/` (hashes en `baseline.md`). La matriz confirma visualmente la desconexión terminal hash↔sigma, el bulto/spike en ambos rings wedge↔bold y los estados de clearance. La medición óptica del glyph O simple da un desplazamiento <0.2 unidades en las tres orientaciones; aún se validará en tests y no se hará corrección global sin un fallo reproducible. Pruebas de regresión, fixes, comparación versionada y shards finales están pendientes.

## Entrega

Pendiente. Sin cambios de producción aún. Se permanecerá en `feat/molecular-depiction-excellence`, sin merge/tag/release; la revisión visual del propietario se informará por separado de las pruebas automatizadas.
