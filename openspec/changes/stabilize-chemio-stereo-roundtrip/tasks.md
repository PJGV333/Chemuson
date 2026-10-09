# Tareas — ChemIO stereo round-trip

## 1. Alcance, Git y baseline
- [x] 1.1 Verificar ramas/upstream reales, árbol limpio y ausencia de commits pendientes; crear `fix/chemio-stereo-roundtrip` desde el SHA beta indicado sin transportar ChemName.
- [x] 1.2 Leer OpenSpec activa y registrar baseline limitado; no ejecutar suite monolítica.
- [x] 1.3 Reproducir cada test P0 en proceso aislado y confirmar el estado estereoquímico de tetrandrina con RDKit.
- [x] 1.4 Trazar las rutas SMILES→grafo, grafo→SMILES, SMILES→MOL→grafo y grafo→MOL; determinar el contrato visual de cuñas.

## 2. Mecanismo ChemIO
- [x] 2.1 Mantener dirección y paridad de enlaces wedge/hash al leer/escribir MOL V2000.
- [x] 2.2 Reconstruir centros tetraédricos explícitos en las conversiones RDKit directa y worker aislado sin depender del texto literal `@`/`@@`.
- [x] 2.3 Preservar E/Z explícito cuando se represente y fallar explícitamente ante degradación no representable.
- [x] 2.4 Mantener conectividad, orden/identidad de átomos/enlaces, órdenes, carga e isótopos.

## 3. Regresión química independiente
- [x] 3.1 Añadir matriz para los dos P0, enantiómeros opuestos, formas SMILES con distinto orden de vecinos, aminoácidos y moléculas multicentro soportadas.
- [x] 3.2 Añadir controles aquirales, centros potencialmente quirales no asignados, y tetrandrina no estereoespecificada; no generar cuñas artificiales.
- [x] 3.3 Añadir E/Z representativo y round-trips SMILES/MOL/SDF con comparación independiente de RDKit, fórmula, conectividad, enlaces, cargas/isótopos y estereoquímica.
- [x] 3.4 Corregir la expectativa histórica de tetrandrina con evidencia RDKit; mantener el caso vancomicina explícito como control positivo.

## 4. Validación y entrega
- [x] 4.1 Ejecutar red/green en procesos aislados y grupos focales, cada comando con timeout externo ≤300 s (máximo absoluto 600 s); no repetir un bloqueo sin diagnóstico cambiado.
- [x] 4.2 Ejecutar compileall, colecta acotada, Ruff focal, OpenSpec estricto disponible y `git diff --check`; documentar full-suite no ejecutada y la única F401 baseline de Ruff global.
- [x] 4.3 Registrar archivos, límites de backend, regresiones/resultados y estado Git/push en `validation.md`.
- [x] 4.4 Crear commits pequeños; sin merge/PR/tag ni publicación.
- [ ] 4.5 Hacer push normal sólo a `fix/chemio-stereo-roundtrip`; pendiente porque GitHub no dispone de credenciales en este entorno.
