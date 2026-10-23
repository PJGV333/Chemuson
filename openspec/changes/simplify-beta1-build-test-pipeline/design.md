# Diseño — triggers, salidas y contratos

## Decisiones

1. **`test.yml` solo `workflow_dispatch`.** Los jobs (plan, 8 shards acotados,
   summary fail-closed, windows-smoke, flatpak-smoke) no cambian: la capacidad
   de auditoría completa permanece intacta y sigue verificando la colección
   contra `tests/ci/expected_pytest_nodeids.txt`. No hay ramas protegidas ni
   status checks requeridos que referencien estos jobs (verificado por API),
   por lo que quitar los triggers no deja checks colgados.
2. **La suite se conserva íntegra en `tests/`.** `tests/test_test_workflow_contract.py`
   no aserta sobre los triggers, solo sobre la estructura de jobs y el
   fail-closed del sharder: sigue vigente sin modificación.
3. **`CHEMUSON_OUTPUT_DIR` como raíz de salidas.** Por defecto se preserva el
   comportamiento actual (`dist/`, `dist-appimage/`, `dist-flatpak/` en el
   repositorio). Si se define, cada script cuelga sus subdirectorios de esa
   raíz; los argumentos posicionales y variables específicas tienen prioridad.
   `build_linux_all.sh` además dirige PyInstaller con `--distpath`/`--workpath`
   dentro de la salida elegida.
4. **`.gitignore`**: `dist-appimage-*/`, `dist-preview*/`, `dist-flatpak-preview/`
   y `pyinstaller-work/`. Los manifiestos de procedencia/checksums ya trackeados
   como evidencia (diagnostic build) permanecen trackeados; no se reescriben.
