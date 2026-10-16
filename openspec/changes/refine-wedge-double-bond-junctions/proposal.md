# Propuesta — refinamiento geométrico de wedges y dobles enlaces

## Por qué

La galería y la revisión manual posteriores a las primeras mejoras detectaron cuatro problemas visuales que requieren reproducción específica: cuña hash junto a líneas sigma/pi, etiqueta ópticamente descentrada en un doble enlace, clearances inconsistentes entre etiquetas y cuñas, y encuentros abultados de cuña sólida con aristas bold de proyecciones Haworth.

## Qué cambia

- Extender la galería Qt con ejemplos explícitos de hash↔doble, O junto a dobles en ejes horizontal/vertical/diagonal, etiquetas en cuña sólida/hash y pentágonos/hexágonos Haworth con aristas bold y wedges sólidos.
- Medir la geometría y capturar baseline antes de modificar rendering.
- Corregir solamente la causa raíz reproducida en el renderer de BondItem, la colocación óptica de labels o el contexto de unión; mantener cambios acotados y conservar controles de enlaces sencillos y anillos.
- Añadir pruebas de geometría, clearance, render determinista, identidad del grafo y regresiones de depiction ya aceptadas.
- Versionar evidencia before/after comparable en temas claro y oscuro y conservar artefactos exploratorios en `/tmp`.

## Invariantes y fuera de alcance

No cambiar átomos/coordenadas, conectividad, IDs, órdenes, estilos, aromaticidad, cargas, metadatos estereoquímicos, semántica de `bold`, persistencia, undo/redo, Clean2D, ChemName, import/export ni el modelo químico. Sin nuevas dependencias, refactor amplio, merge a `main`/`release/v0.3.0-beta.1-prep`, tags o releases.

## Aceptación

Los casos antes/después deben usar idénticos datos y escala. Las pruebas deben demostrar que mejora la integración hash↔doble y wedge↔bold, centra ópticamente etiquetas sin invadir trazos y equilibra clearances; los controles comunes no deben degradarse. La revisión visual automatizada no sustituye la del propietario.
