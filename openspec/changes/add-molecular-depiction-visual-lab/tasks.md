# Tareas — Molecular Depiction Visual Lab

## 1. Diagnóstico y OpenSpec
- [x] 1.1 Hacer fetch de la rama experimental, verificar SHA remoto, árbol limpio y crear rama aislada desde ese SHA.
- [x] 1.2 Leer catálogo de módulos, contratos relevantes, cambios de depiction, renderer y galerías existentes.
- [x] 1.3 Ejecutar y registrar baseline exacto (con timeout externo) antes de modificar código/documentos.
- [x] 1.4 Definir controles que usan parámetros existentes y separar overrides opt-in de unsupported.

## 2. Núcleo renderer compartido
- [ ] 2.1 Añadir modelo tipado, serialización JSON estricta, rangos y migración versionada.
- [ ] 2.2 Crear adaptador que inicializa dos `ChemusonCanvas` desde copias del mismo grafo y aplica estilos sólo al canvas Working.
- [ ] 2.3 Implementar los dos primeros casos (cuña sólida + etiqueta; hash + doble) y probar que sólo el renderer calibrado cambia.
- [ ] 2.4 Parametrizar los overrides opt-in con defaults `None` de comportamiento equivalente y pruebas de paridad.
- [ ] 2.5 Construir el corpus estable de mínimo 12 casos y snapshots de preservación.

## 3. Interfaz y diagnóstico
- [ ] 3.1 Implementar ventana Qt de dos paneles con selección, zoom/ángulo sincronizados, temas y ocultación de controles.
- [ ] 3.2 Añadir sliders/spinboxes reales, opciones discretas, categorías, tooltips y reset individual/global.
- [ ] 3.3 Añadir overlays de geometría y leyenda, aislados de exportación.
- [ ] 3.4 Implementar carga/guardado Current/Working, comparar perfiles y almacenamiento propio del laboratorio.

## 4. Reproducción y validación
- [ ] 4.1 Implementar exportación PNG actual/calibrado/comparativa/contact sheet, perfil y resumen Markdown SHA/versionado.
- [ ] 4.2 Añadir README y comando de arranque comprobado; probar apertura/cierre offscreen aislado con timeout.
- [ ] 4.3 Añadir pruebas focalizadas de parámetros, UI, renderer, invariantes y exportación.
- [ ] 4.4 Ejecutar tests focales y shards CI delimitados; reconciliar contra baseline y fallos conocidos sin ocultarlos.
- [ ] 4.5 Validar OpenSpec y documentar evidencia/limitaciones en `validation.md` y `AGENT_REPORT.md`.

## 5. Entrega de rama
- [ ] 5.1 Hacer commits pequeños y coherentes en `feat/molecular-depiction-visual-lab`.
- [ ] 5.2 Si autenticación disponible, push únicamente a esta rama; verificar SHA remoto y árbol limpio.
- [ ] 5.3 Declarar `READY FOR OWNER VISUAL REVIEW` sólo si todos los criterios objetivos pasan; nunca aprobar el preset estético.
