# Validación — CI, viewport y proyecciones

## Implementación

- Rama: `feat/molecular-depiction-visual-lab`.
- Baseline renderer: `87e298c46535e1dffb0477f36412b3dedc87bd98`.
- SHA de implementación reproducido en evidencia: `0f501cc70388b429d34c7bf7132ad227849a8aad`.
- Se preservaron el timeout externo de **300 s**, los defaults del editor y el preset Editorial AB3 aprobado; no se editó `.cmsn`, `MolGraph` ni `BondStereo`.

## CI local por shards

Colección exacta: **2.358 node IDs**, SHA `a55f94c1035aed5641188abeb7a986a9694e2d4311ba15c1a89dfaff92f08d1f`. El perfil de duraciones tiene SHA `fb12b2cbb3e519dbc27270d6e849d5ac5881f7dd92eafda92cf79c24af1fe830`. LPT por archivo produjo distribución `[275,231,262,316,442,282,264,286]` y cargas previstas entre 145,427–145,618 s.

| Shard | Asignados/timed/JUnit | Resultado pytest | Tiempo |
|---:|---:|---|---:|
| 0 | 275/275/275 | 271 passed, 4 skipped | 51,428 s |
| 1 | 231/231/231 | 231 passed | 71,263 s |
| 2 | 262/262/262 | 259 passed, 3 skipped | 76,676 s |
| 3 | 316/316/316 | 316 passed | 78,770 s |
| 4 | 442/442/442 | 442 passed | 108,427 s |
| 5 | 282/282/282 | 278 passed, 3 skipped, **1 failed** | 136,318 s |
| 6 | 264/264/264 | 260 passed, 4 skipped | 66,889 s |
| 7 | 286/286/286 | 280 passed, 6 skipped | 71,028 s |
| **Total** | **2.358/2.358/2.358** | **2.337 passed, 20 skipped, 1 failed** | todos < 300 s |

Se ejecutó el verificador fail-closed. Confirmó asignación única, cobertura JUnit y timings exactos de los ocho shards; salió **1** únicamente porque el estado de la matriz es failure y el shard 5 terminó con exit 1. No se omitió ni se marcó como exitoso ningún nodo.

### Falla basal reproducida

El único fallo local fue `tests/test_molecular_assistant_reference_resolution.py::test_ai_success_with_equivalent_nonidentical_smiles_is_verified`, que agota el `QSignalSpy.wait(8000)` tras el orden UI `test_beta_compact_shell_and_tabs.py` → `test_beta_save_shortcut_menu_access.py`. El mismo prefijo de 69 pruebas anteriores y el nodo fallido, ejecutado desde un worktree limpio del SHA basal `87e298c`, reprodujo **69 passed, 1 failed en 79,38 s**. El nodo aislado pasa **1/1 en 0,75 s**; la familia 3D/Molecular Assistant pasa **30 passed, 1 skipped**. Esos tests y `MolecularAssistantController` no cambiaron en esta propuesta. Se conserva y reporta como una falla preexistente de lifecycle/orden Qt; no se aumentó su timeout, no se editaron tests ni se cambió el sharder para ocultarla. Evidencia y solicitud de tratamiento separado: `AGENT_REPORT.md`.

La baseline completa anterior también registró `tests/test_compchem3d_dock.py::test_compchem_controller_generates_async_with_fake_backend` como único fallo de su `pytest -q` monolítico. No se repitió el monolito: el baseline demoró 1.310,93 s, por encima del límite local de diez minutos.

## Pruebas y verificaciones

- Focales de proyección, geometría real, viewport/export, preset, arquitectura y contrato CI: **118 passed**.
- `pytest -q tests/architecture`: **283 passed**.
- `python -m compileall src tests tools packaging`: PASS.
- Ruff focal de todos los archivos Python modificados: PASS.
- Ruff global: conserva sólo el F401 basal ajeno `math` en `tests/test_clean2d_para_disubstituted_aromatic_layout_v1.py:3`.
- `openspec validate stabilize-visual-lab-ci-and-projections --strict --no-interactive`: PASS.
- `git diff --check`: PASS (el log completo de Actions se versiona como gzip determinista para conservar bytes originales sin whitespace de diff).
- La actualización de `tests/ci/expected_pytest_nodeids.txt` fue explícita y revisada; el manifiesto cubre exactamente 2.358 IDs.

## Evidencia de renderer

- Comparación Haworth C3, igual renderer/tema/tamaño 1200×820, baseline/refined: `evidence/haworth/before-after-C3.png` (SHA-256 `806eb83f2e18ed7dc6a5d4599e78e358bc27c2deea41821b984b3d117a497a5b`). Los bundles de reproducción y sus hashes están en `evidence/haworth/{baseline,refined}/`.
- Ventana Visual Lab real, Current/Working centrados y sincronizados: `evidence/viewport/depiction-lab-C3-synced.png` (SHA-256 `1994ec67608aaef7ac4587dde7b37a49a8f5a834a29da8f52b1c2674d2630706`). Probe: 12 observaciones, cerrada correctamente.
- Haworth piranosa/furanosa, Fischer y silla A/B con renderer real y metadatos semánticos: `evidence/projections/projection-geometry-gallery.png` y `.md` (PNG SHA-256 `681a35cba206250b9135d3a1341998efea14a7c10f26baa5d7a3b1deb5a95642`).
- Esta evidencia verifica geometría y representación, **no constituye aceptación estética del propietario**. C3 permanece fuera de AB3; A4/B4 siguen excluidos y B4 continúa siendo fixture no editorial.

## Estado remoto

Se publicó mediante push normal la rama `feat/molecular-depiction-visual-lab`. El SHA local/remoto confirmado es `9e7e92a5fd5e31168632b79226ae3f9ac4b26da8`; el worktree quedó limpio. No se hizo merge.

GitHub Actions run [`39002887955`](https://github.com/PJGV333/Chemuson/actions/runs/39002887955), sobre ese SHA, **falló**. La planificación, los smoke de Windows/Flatpak y siete de ocho shards pasaron. Shard 5 falló en el mismo nodo basal `tests/test_molecular_assistant_reference_resolution.py::test_ai_success_with_equivalent_nonidentical_smiles_is_verified`: `QSignalSpy.wait(8000)` venció; el shard registró 278 passed, 3 skipped, 1 failed de 282 IDs, todos ejecutados y cronometrados. El verificador final rechazó correctamente la matriz por ese exit 1; no faltó cobertura. Resultado concuerda con la reproducción local limpia sobre el SHA basal `87e298c` y no se corrigió fuera de alcance.

Por tanto, la ejecución remota no está verde, la falla Qt de lifecycle/orden permanece abierta y no se hará merge. Esta evidencia visual tampoco implica aceptación estética del propietario.
