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
