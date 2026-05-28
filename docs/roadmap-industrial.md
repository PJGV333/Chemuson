# Hoja de ruta industrial para ChemUSON

Este documento consolida una estrategia técnica para llevar ChemUSON hacia capacidades comparables con suites como ChemDraw, MarvinSketch o ChemDoodle, aprovechando la arquitectura actual basada en PyQt6, `MolGraph` y RDKit.

## Objetivo

Evolucionar ChemUSON desde un editor 2D académico robusto hacia una plataforma profesional de dibujo, validación e inteligencia química, con interoperabilidad industrial y herramientas analíticas avanzadas.

## 1) Mejoras críticas sobre capacidades existentes

### 1.1 Rotación pseudo-3D → conformaciones 3D reales

**Situación actual**
- `project_3d_rotation` en `geom.py` proyecta desde un plano 2D con restricción de tilt para preservar legibilidad.

**Meta**
- Migrar a un pipeline híbrido:
  1. Generación de coordenadas 3D reales con RDKit (`EmbedMolecule` + minimización MMFF/UFF).
  2. Control trackball sobre matriz de rotación 3D.
  3. Proyección 3D→2D con depth cueing (ancho/opacidad según Z).

**Entregables técnicos**
- `chemuson/geometry3d/` con servicio asíncrono de conformaciones.
- Cache por hash químico para evitar recomputar conformaciones.
- Modo de interacción `Alt + Arrastre` sin límite rígido de ±60°.

---

### 1.2 Limpieza 2D (Ctrl+K) con heurísticas profesionales

**Situación actual**
- Limpieza 2D funcional orientada a normalización geométrica base.

**Meta**
- Introducir layout dirigido por fuerzas con restricciones químicas:
  - Longitudes objetivo según tipo de enlace.
  - Penalización angular para `sp2`/aromáticos (~120°).
  - Simetría en cadenas y anillos.
  - Resolución de colisiones etiqueta-enlace.

**Entregables técnicos**
- Nuevo motor `clean2d_v2` con parámetros configurables.
- Pruebas visuales de regresión (snapshots SVG/PNG).
- Opción de “limpieza rápida” y “limpieza de publicación”.

---

### 1.3 Validación de valencia con diagnóstico interactivo

**Situación actual**
- `validate` marca errores visuales, pero con retroalimentación limitada.

**Meta**
- Añadir capa explicativa por átomo/enlace:
  - Tooltips de hover (explicación matemática del error).
  - Panel de inspección con sugerencias de corrección (cargas, orden de enlace, protonación).

**Entregables técnicos**
- Modelo de `ValidationIssue` con código, severidad, mensaje y hint.
- Navegación “siguiente/anterior error”.
- API reutilizable para UI y exportes.

---

### 1.4 Nomenclatura bidireccional (Name→Structure)

**Situación actual**
- IUPAC-lite con degradación segura (`N/D`) y fallback RDKit aislado.

**Meta**
- Agregar conversión inversa desde texto:
  - Nombre común/sistemático → estructura.
  - Integración con OPSIN y/o APIs públicas (PubChem).
  - Construcción automática de `MolGraph` + limpieza inicial.

**Entregables técnicos**
- Cuadro de búsqueda “Nombre a estructura”.
- Conectores desacoplados (offline/online).
- Registro de procedencia y nivel de confianza.

## 2) Funciones nuevas para competir en productividad

### 2.1 Dock de inteligencia química en tiempo real

- Fórmula, masa exacta, peso molecular, análisis elemental.
- `logP`, TPSA, HBD/HBA, enlaces rotables, alertas de regla de Lipinski.
- Recalculo incremental para no bloquear UI.

### 2.2 Superátomos y grupos abreviados

- Etiquetas como `Me`, `Et`, `Ph`, `Boc`, `Ts`, `Ac`.
- Vista colapsada en canvas + representación expandida en `MolGraph`.
- Preservación correcta en validación, exportación y cálculo de propiedades.

### 2.3 Mecanismos de reacción (arrow pushing)

- Flechas curvadas de par electrónico y anzuelo radicalario.
- Snapping inteligente a átomos, enlaces y centros de pares solitarios.
- Capas/objetos editables independientes del esqueleto molecular.

### 2.4 Predicción espectral integrada (NMR/MS)

- Estimación base de desplazamientos de ^1H y ^13C con fingerprints de entorno.
- Vista de espectro con interacción cruzada (click pico ↔ highlight átomo).
- Arquitectura de plugins para futuros predictores más avanzados.

### 2.5 Soporte polímeros y Markush

- Corchetes de repetición (`n`, `x`) y definición de grupos `R`.
- Restricciones lógicas para listas de sustituyentes válidos.
- Exportes compatibles con formatos de intercambio químico.

### 2.6 Interoperabilidad de formatos

- Estructurales: `.mol`, `.sdf`, `.pdb`, `.cml`.
- Suites: `.cdxml`, `.mrv`.
- Publicación: SVG/PDF/EPS mejorados (resolución y perfiles de color).

## 3) Plan de ejecución por fases

### Fase 1 (0–3 meses): Fundaciones
- `clean2d_v2` + `ValidationIssue` interactivo.
- Dock de propiedades básicas en tiempo real.
- Refactor de servicios asíncronos para tareas químicas pesadas.

### Fase 2 (3–6 meses): Productividad avanzada
- Superátomos y abreviaturas.
- Arrow pushing.
- Conversión Name→Structure con conector OPSIN/PubChem.

### Fase 3 (6–12 meses): Diferenciadores
- Pipeline 3D real con proyección robusta.
- Predicción espectral base.
- Polímeros/Markush y expansión mayor de import/export.

## 4) Arquitectura recomendada

- **Núcleo químico desacoplado**: servicios puros para layout, validación, descriptores y 3D.
- **Orquestación asíncrona**: workers/subprocesos para RDKit pesado y APIs externas.
- **UI reactiva**: señales Qt granulares para actualizaciones incrementales.
- **Trazabilidad**: telemetría local opcional de desempeño (tiempos de layout, validación, render).

## 5) Métricas de éxito

- Tiempo medio de limpieza 2D por estructura.
- Errores de valencia detectados/corregidos por sesión.
- Latencia de actualización del panel de propiedades.
- Porcentaje de import/export round-trip sin pérdida semántica.
- Tasa de éxito Name→Structure por fuente (offline/online).

## 6) Riesgos y mitigación

- **Dependencia de servicios externos**: modo offline y cache local.
- **Bloqueo de interfaz**: aislamiento en subprocesos + límites de tiempo.
- **Compatibilidad entre formatos**: estrategia de mapeo por prioridad semántica.
- **Complejidad de UI**: feature flags y despliegue gradual.

## 7) Principios de implementación

1. No degradar estabilidad del editor 2D existente.
2. Mantener degradación segura ante casos no soportados.
3. Medir rendimiento antes/después en cada entrega.
4. Diseñar APIs internas orientadas a pruebas automatizadas.
5. Priorizar valor de usuario por iteración (shipping incremental).

## 8) Ruta de cierre basada en el estado actual

Esta ruta parte del estado observado en mayo de 2026: `clean2d_v2` ya existe, está integrado y fue validado funcionalmente en uso manual, por lo que no se prioriza como siguiente bloque salvo regresiones puntuales. La prioridad pasa a cerrar MVPs ya iniciados y convertirlos en capacidades robustas.

### Bloque A: estabilización técnica antes de ampliar superficie

**Objetivo**
- Asegurar que la suite completa pueda ejecutarse sin fallos nativos ni contaminación de rutas de build.

**Trabajo**
- Aislar `dist-flatpak/`, `dist/` y artefactos empaquetados del `PYTHONPATH`/descubrimiento de pruebas.
- Documentar entorno recomendado de pruebas con Python soportado y RDKit compatible.
- Añadir smoke test explícito del worker RDKit aislado, sin importar RDKit en el proceso principal.

**Criterio de terminado**
- `pytest -q` ejecuta sin `segmentation fault`.
- Las pruebas RDKit fallan de forma controlada o se saltan cuando RDKit no está disponible.
- No se importa RDKit directamente desde builds empaquetados durante pruebas locales.

### Bloque B: validación química interactiva completa

**Objetivo**
- Convertir `ValidationIssue` en una herramienta de corrección, no solo resaltado.

**Trabajo**
- Añadir panel/dock de validación con lista de errores, severidad, mensaje, fórmula del cálculo y sugerencia.
- Permitir click en issue para seleccionar y centrar el átomo/enlace afectado.
- Añadir acciones seguras de corrección donde aplique: ajustar carga formal, reducir orden de enlace, limpiar H explícitos/asignados.
- Exponer `ValidationIssue.as_dict()` para exportes/reportes.

**Criterio de terminado**
- `Ctrl+Shift+V` abre/actualiza el panel además de resaltar.
- `F8`/`Shift+F8` sincronizan selección, panel y tooltip.
- Hay pruebas unitarias del modelo y pruebas GUI mínimas del flujo de navegación.

### Bloque C: Name->Structure profesionalizado

**Objetivo**
- Pasar de resolución básica a un flujo trazable offline/online.

**Trabajo**
- Mantener conector offline actual para nombres comunes.
- Añadir conector OPSIN cuando esté disponible localmente o vía configuración opcional.
- Mantener PubChem como conector online con cache, procedencia y confianza.
- Añadir diálogo con fuente, confianza, SMILES resultante y opción de cancelar antes de insertar.
- Ejecutar resolución en worker para no bloquear UI.

**Criterio de terminado**
- La UI no se congela durante consultas online.
- Cada inserción conserva fuente, nombre resuelto, confianza y SMILES en estado/documento cuando sea posible.
- Hay tests con conectores falsos y cache; PubChem real no es obligatorio en CI.

### Bloque D: superátomos y abreviaturas robustas

**Objetivo**
- Hacer que grupos como `Me`, `Et`, `Ph`, `Boc`, `Ts`, `Ac` sean objetos químicos coherentes en dibujo, cálculo, validación y exportación.

**Trabajo**
- Formalizar un modelo interno de abreviatura con etiqueta colapsada, grafo expandido y átomo de anclaje.
- Unificar la expansión usada por RDKit con validación, propiedades y exportación.
- Añadir edición segura: colapsar/expandir, cambiar etiqueta, preservar anclaje.
- Definir límites claros para abreviaturas no soportadas.

**Criterio de terminado**
- Propiedades, espectros MVP y exportación SMILES/Molfile usan la expansión correcta.
- La vista colapsada no rompe selección, undo/redo ni persistencia `.cmsn`.
- Tests cubren abreviaturas principales y casos no soportados.

### Bloque E: 3D real como modo de interacción sólido

**Objetivo**
- Completar el salto de pseudo-3D a conformaciones 3D reales sin perder legibilidad 2D.

**Trabajo**
- Mantener generación RDKit aislada con cache por hash químico.
- Hacer la generación completamente asíncrona desde UI.
- Aplicar depth cueing visible en render: opacidad/ancho/color por Z donde no perjudique publicación.
- Añadir controles de reset, regenerar conformación y modo fallback explícito.
- Evitar bloquear o degradar estructuras grandes con timeout y mensajes claros.

**Criterio de terminado**
- `Alt + Arrastre` usa conformación real cuando existe y fallback pseudo-3D cuando falla.
- La UI informa si está usando RDKit, cache o fallback.
- Hay pruebas de proyección, cache, fallback y no bloqueo.

### Bloque F: espectros MVP útil y extensible

**Objetivo**
- Convertir la predicción heurística actual en una herramienta claramente marcada, navegable y extensible.

**Trabajo**
- Mostrar confianza y fuente por predicción.
- Mejorar reglas heurísticas para carbonilos, aromáticos, heteroátomos y equivalencias simples.
- Añadir interfaz de plugin documentada para predictores externos.
- Permitir exportar tabla de picos.

**Criterio de terminado**
- El dock indica explícitamente que la predicción es estimada.
- Click pico <-> átomo funciona para ^1H y ^13C.
- Tests cubren familias simples: alcano, alcohol, carbonilo, aromático.

### Bloque G: polímeros y Markush semánticos

**Objetivo**
- Pasar de anotaciones detectables a entidades semánticas exportables.

**Trabajo**
- Modelar corchetes de repetición como parte persistente del documento químico.
- Definir grupos R con listas de sustituyentes permitidos.
- Añadir panel de edición de R-groups y repetición (`n`, `x`, rangos simples).
- Preparar exportación semántica cuando el formato lo soporte y fallback gráfico cuando no.

**Criterio de terminado**
- `.cmsn` conserva polímeros/Markush sin pérdida.
- Exportes no compatibles degradan de forma visible y documentada.
- Tests cubren persistencia, edición y resumen semántico.

### Bloque H: interoperabilidad industrial

**Estado mayo 2026**
- Primer corte implementado con CML import/export/copy y matriz de compatibilidad en `docs/interoperability.md`.
- Quedan pendientes CDXML/MRV/PDB y metadatos avanzados de publicación.

**Objetivo**
- Ampliar import/export sin perder información química ni visual esencial.

**Trabajo**
- Consolidar round-trip `.mol`/`.sdf` como base confiable.
- Añadir CML como primer formato semántico XML por menor complejidad que CDXML/MRV.
- Evaluar CDXML/MRV por fases: primero exportación parcial, luego importación.
- Mejorar exportación SVG/PDF con metadatos, perfiles y opciones de publicación.

**Criterio de terminado**
- Matriz de compatibilidad documentada por formato y entidad química.
- Tests de round-trip para los formatos soportados.
- Exportes parciales declaran explícitamente qué información se pierde.

### Orden recomendado inmediato

1. Bloque A: estabilización de pruebas/RDKit.
2. Bloque B: panel de validación interactiva.
3. Bloque C: Name->Structure asíncrono y trazable.
4. Bloque D: superátomos coherentes en cálculo/exportación.
5. Bloque E: 3D real completo en UI.
6. Bloque F: espectros mejorados.
7. Bloque G: polímeros/Markush semánticos.
8. Bloque H: formatos industriales.

### Siguiente implementación sugerida

El siguiente bloque práctico debe ser **Bloque A**, porque desbloquea verificación confiable antes de tocar más funcionalidad. Después, el mayor valor de usuario está en **Bloque B**, ya que aprovecha `ValidationIssue` ya implementado y lo convierte en una herramienta visible de diagnóstico/corrección.
