# Design — estabilización química y del visor 3D

## Invariantes

1. `MolGraph` sigue siendo autoridad de conectividad, propiedades y estereoquímica documental. `CoordinateSet3D` conserva coordenadas por IDs documentales estables.
2. Los H/átomos derivados RDKit son metadata transitoria de escena y nunca se agregan al grafo ni a `.cmsn`.
3. Cada transformación geométrica conserva la relación padre-enlace de H derivados y aplica la misma transformación rígida a todo fragmento movido. Una edición deformante de un único átomo mantiene al H unido a su padre, pero no simula una reorientación química global no especificada.
4. El worker recibe un payload explícito. RDKit debe asignar/interpretar los descriptores soportados y la configuración observada se compara con la solicitada después del embedding y de la minimización. Configuración ambigua/no soportada produce diagnóstico; nunca se marca como preservada por conservar sólo texto en `MolGraph`.
5. Aromaticidad explícita y Kekulé son representaciones admitidas y diferenciadas. Sanitización RDKit puede percibir aromaticidad en un Kekulé válido; ChemUSON no promueve anillos conjugados arbitrarios. La escena retiene orden/aromaticidad para render consistente.
6. Cambios GUI son acotados a los widgets/malla existentes. No se crean recursos ni arquitectura de render nueva sin métricas; una verificación offscreen no certifica contraste, encuadre ni layout visual.
7. El SIGSEGV pytest y el SIGSEGV Qt/SIP son fallos distintos hasta que una reproducción diferencial demuestre vínculo. No usar `os._exit`, excepciones silenciosas, destrucción indiscriminada o sleeps.

## Ruta de investigación del crash

La evidencia existente incluye un core de la demo en el que la UI muere durante `sip_api_visit_wrappers` → `release_QApplication` al finalizar el script. El hilo de UI está en la destrucción del wrapper; otros hilos observados pertenecen a XCB/Qt y Mesa. Una nota previa reporta una reproducción con dock oculto, que debilita al renderer como condición necesaria pero no identifica causa. El trabajo separará: salida normal/orden de referencias Python, cierre de ventana sin trabajo, trabajo Qt activo, renderer/contexto y teardown de aplicación. No se vuelve a lanzar GUI si requiere aprobación adicional; no se diagnostica estabilidad sólo con offscreen.

## Cambio mínimo previsto

- Extraer/añadir un helper puro que transforme coordenadas de H metadata junto a la transformación aplicada a los IDs parentales; reutilizarlo desde edición y Live MM, sin incluir IDs derivados en comandos/persistencia.
- Completar el payload RDKit y la importación de texto para no perder etiquetas CIP/EZ explícitas: el worker puede devolver descriptores por índice de átomo/enlace y el adaptador los aplica a IDs documentales estables. Para R/S absolutos inequívocos se asignan tags y se valida CIP. Para E/Z se fijan los átomos de referencia ordenados a ambos lados del doble enlace; si falta representación suficiente, devolver limitación explícita. Validar tras embedding y relajación.
- Propagar `is_aromatic` y orden a `SceneBond3D`; malla de orden doble/triple/aromático se representa sin cambiar el picking/documento 2D.
- Reorganizar controles en secciones compactas y ajustar iluminación/paleta mediante contrato de malla/shader existente; conservar Live MM accesible.

## Verificaciones

- Tests first RED/GREEN de rotación derivada y H en etano/etanol/butano; undo/redo, cancelación y posterior Live MM.
- Tests RDKit R, S, E, Z y controles sin especificación, incluyendo pos-embedding y pos-optimización; casos ambiguos rechazados.
- Tests de benceno aromático explícito/Kekulé, ciclohexano, eteno y etino; las estructuras siguen distintas químicamente cuando corresponde y órdenes quedan en escena/malla.
- Tests puros de multiplicidad de enlace y wiring del dock. Revisión gráfica manual pendiente cuando no haya autorización.
- Reproducción de shutdown no gráfica, análisis del core disponible y comparación de ciclo de vida con/sin trabajo activo. No repetir suite completa conocida; correr shards seguros.

## No cambios

Sin cambios de formato `.cmsn`, sistema de documentos, Clean2D/ChemName, jerarquía de mixins o dependencias runtime. Si la investigación exige una modificación transversal de ownership, dependencia nueva o movimiento de paquetes, detenerse y pedir decisión de diseño conforme a `AGENTS.md`.