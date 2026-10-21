# Propuesta — estabilización final de ChemUSON 0.3.0-beta.1

## Por qué

La candidata parte de `release/v0.3.0-beta.1-prep` en `61f7647b683db3eca19ea06e5ec82ab2897d640a`. El workflow `test` de ese SHA ejecutó la colección completa en ocho shards, pero fallaron los shards 0 y 6. La evidencia disponible señala una demora en la finalización de un QThread de CompChem y un doble de prueba obsoleto que no implementa `progress_ready`; un wrapper de subprocess reproduce el segundo fallo. El teardown/cierre Qt debe verificarse separadamente, sin confundir fallos de harness con crashes del paquete.

Además, antes de decidir GO/NO-GO hay que volver a auditar los cuatro formatos existentes y la matriz manual de aceptación, conservar la procedencia/hash de artefactos y reconciliar la documentación histórica del AppImage. La campaña está congelada en funcionalidades.

## Qué cambia

- Abrir esta campaña OpenSpec en la rama local `stabilize/v0.3.0-beta.1`, creada en el mismo directorio desde el SHA beta esperado.
- Reproducir y clasificar los tres nodos fallidos del run Actions `39679607676`, conservar logs/JUnit/timings/plan y corregir sólo la regresión demostrada.
- Mantener los contratos de señales de cada worker simulado sincronizados con `CompChem3DWorker`; no cambiar producción por inferencia ni aumentar timeouts de forma general.
- Ejecutar los nodos focalizados y la matriz existente de ocho shards con `tools/ci_pytest_shards.py`; documentar los resultados reales aunque no sean verdes.
- Auditar los workflows/scripts de preview y release, AppImage Type 2, Flatpak, Windows portable/Inno, recursos, checksums, procedencia y documentación; reutilizar el pipeline actual.
- Actualizar la matriz manual existente y notas/documentación sólo con evidencia exacta; separar automatización de aceptación visual/química del propietario.
- Presentar dictamen GO/NO-GO. Se permite sólo commit y push normal a `stabilize/v0.3.0-beta.1` con autorización explícita del owner; sin merge, tag, Release, PR, worktree ni canal público.

## Fuera de alcance

No agregar funcionalidades, rediseñar la GUI, cambiar Clean2D, modelos o algoritmos químicos, modificar serialización `.cmsn`, cambiar permisos/canales del pipeline ni ocultar pruebas. No crear paquetes/dependencias. No cambiar producción GUI salvo que una regresión de ciclo de vida quede reproducida y una corrección mínima se justifique en esta especificación antes de editar.

## Aceptación

1. Los tres node IDs indicados tienen causa raíz clasificada con referencia a evidencia del run exacto; se distinguen doble obsoleto, comportamiento de producción y condiciones del entorno.
2. Los dos casos de cierre ejecutan sus aserciones de cierre tardío y señal pendiente; no fallan por un contrato ficticio de worker. Ninguna falla se convierte en skip/xfail ni queda fuera del manifiesto.
3. La colecta coincide exactamente con el manifiesto y los ocho shards se ejecutan con JUnit/log/timing/proceso conservados. Cero fallos, errores, timeout, señales fatales o crashes es el criterio de CI verde; skips preexistentes se enumeran y justifican individualmente.
4. Los cuatro paquetes, cuando exista una ejecución autorizada para el nuevo SHA, se verifican por separado con checksum/provenance; las comprobaciones automatizadas no se presentan como aceptación manual.
5. La matriz/manual, notas y guías reflejan sólo resultados trazables. Windows nativo, instalación, OpenGL visual u otros entornos no disponibles figuran NOT TESTED.
6. No se publica ni se afirma estabilidad sin evidencia. La decisión final identifica bloqueos P0/P1, advertencias, artefactos y pruebas pendientes.
