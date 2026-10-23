# Tareas — aminas ordinarias con NH₂/H₂N/HN horizontales

## A. Diagnóstico y baseline
- [x] A.1 Rama `fix/primary-amine-label-orientation-beta1` desde
      `origin/stabilize/v0.3.0-beta.1` @ `33ed5223`; worktree limpio; `main` y ramas de
      release intactas.
- [x] A.2 Reproducir el defecto con script diagnóstico offscreen (etilamina, anilina,
      dietilamina, 1,1-diamino-3-fenilpropan-2-ona, acetamida, N,N′-dietilurea, amonio)
      y capturas `before/`.
- [x] A.3 Registrar baseline (status, compileall, pytest completo, Ruff) en `baseline.md`.

## B. Corrección
- [x] B.1 Reutilizar `_nitrogen_geometry_context`; restringir candidatos verticales a
      contextos distintos de `amine` en `_refresh_atom_label` /
      `_choose_implicit_hydrogen_position`.
- [x] B.2 Añadir refresco de N adyacentes a los extremos de un enlace en
      `add_bond_item`/`update_bond_item`/`remove_bond_item`
      (`_nitrogen_context_dependent_refresh_ids`).
- [x] B.3 No modificar layout de glifos, anclaje óptico, recorte, geometría del grafo,
      persistencia ni Clean2D.

## C. Regresiones y evidencia
- [x] C.1 Nuevo `tests/test_amine_label_orientation.py`: horizontalidad con rects reales,
      no solape con enlaces, reproducción exacta del dicetona diamino, controles
      amida/urea/amonio, edición+undo/redo, arrastre+undo, temas claro/oscuro.
- [x] C.2 Actualizar el contrato de amina secundaria de
      `tests/test_amide_urea_2d_continuity.py` (`above` → `left`) documentado en
      `design.md`.
- [x] C.3 Capturas comparativas claro/oscuro en `evidence/`.

## D. Verificación y entrega
- [x] D.1 Tests dirigidos de etiquetas/H/enlaces en verde; `compileall`; Ruff focal;
      `openspec validate --strict`; `git diff --check`.
- [x] D.2 Clasificar los 63 fallos del baseline monolítico como ruido de entorno
      (pasan aislados antes y después del cambio); ver `validation.md`.
- [x] D.3 Commits y push normal a la rama de fix; sin merge a `stabilize/v0.3.0-beta.1`
      (revisión visual pendiente del propietario); sin paquetes nuevos.
