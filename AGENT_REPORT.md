# AGENT_REPORT — higiene del repositorio (en curso)

## Alcance y punto de partida

Rama `maintenance/repository-hygiene-closure`, creada desde `origin/main` en
`1db4f63b52af79247745b3a8a220fb728348218c`. No se integra a `main`, no se
reescribe historia y no se cambia comportamiento del producto. El checkpoint
`afe1245` guarda OpenSpec, baseline, memoria de campañas y política antes de
retirar artefactos.

El OpenSpec activo es
`openspec/changes/repository-hygiene-closure-2026-10-03/`. La baseline
preliminar registra 1816 tests recolectados, 1760 passed/55 skipped/1 fallo
CompChem ya conocido, 269 tests de arquitectura y un F401 histórico Clean2D.
Las medidas completas están en `baseline.md` del OpenSpec.

## Decisiones y cambios en curso

- `docs/history/CAMPAIGNS.md` consolida campañas de arquitectura, UI y Clean2D,
  decisiones, resultados, experimentos descartados y ramas únicas protegidas.
- `docs/history/REPOSITORY_POLICY.md` fija el proceso de baseline, auditoría,
  retención y ciclo de ramas.
- Se retiran únicamente outputs sin consumidor verificado: PostScript/PNG de
  depuración, parches sueltos ya históricos, reporte orbital generado, código
  de demostración del spike, capturas intermedias y generadores one-shot de UI.
- Se retienen OpenSpecs/contratos Markdown, evidencia final UI, capturas
  KDE/Wayland aprobadas, checks JSON, fixtures actuales, Clean2D productivo,
  la referencia normativa `pyqt6-spike/theme.py` y ramas no-ancestro.
- Se actualizaron comentarios/docstrings de tema para no enlazar prototipos
  retirados. No se modifica implementación Clean2D ni código/test químico.
- No se añadieron reglas `.gitignore`: los outputs retirados no se regeneran
  por tests/CI ni requieren una regla amplia.

## Verificación

**Pendiente al redactar esta actualización:** tests y validadores después de la
poda, medición final, inventario post-prune y push normal de la rama. No se
reporta ningún gate como aprobado hasta ejecutarlo. Los fallos de baseline no
se cambian ni se silencian.

Reporte final y lista exacta de archivos/ramas se completarán en
`docs/history/REPOSITORY_CLEANUP_2026-10-03.md` al cerrar las verificaciones.
