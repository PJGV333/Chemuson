## ADDED Requirements

### Requirement: La optimización en vivo emite progreso geométrico real

El control `Optimización en vivo` SHALL ejecutar minimización MM incremental fuera del hilo GUI y SHALL mostrar frames recibidos mientras el proceso sigue trabajando. El sistema SHALL NOT presentar como progreso una lista de frames almacenada y emitida únicamente después de terminar la optimización.

#### Scenario: Frame llega antes del término
- **GIVEN** una molécula parametrizable, coordenadas iniciales y Live MM activo
- **WHEN** el worker inicia bloques de minimización
- **THEN** el dock aplica secuencialmente frames de coordenadas/energía antes del evento terminal del worker
- **AND** cada frame tiene secuencia creciente y coordenadas finitas identificadas por `atom_id`

#### Scenario: Cambio de coordenadas en una operación real
- **GIVEN** butano con torsión inicial no mínima
- **WHEN** se activa Live MM y se registra el flujo de frames
- **THEN** se observan al menos dos snapshots geométricos distintos antes del resultado terminal
- **AND** se registra energía en kcal/mol para cada frame disponible

### Requirement: Los cálculos se ejecutan con RDKit aislado y force field explícito

La minimización SHALL usar el worker RDKit existente, preferir MMFF94 cuando parametrice toda la molécula, permitir MMFF94s y ofrecer UFF como alternativa visible. La interfaz SHALL mostrar el force field efectivamente utilizado, energía en kcal/mol, iteración y estado. Si falta parametrización, SHALL mostrar `Sin parámetros`/error y SHALL NOT simular frames, borrar la geometría ni modificar conectividad.

#### Scenario: MMFF94 aplicable
- **GIVEN** un grafo para el que MMFF94 proporciona parámetros
- **WHEN** comienza Live MM
- **THEN** se muestran frames y el estado identifica MMFF94 con energía en kcal/mol

#### Scenario: MMFF94 no parametriza
- **GIVEN** una estructura no soportada por MMFF94
- **WHEN** comienza Live MM bajo política automática con UFF disponible
- **THEN** UFF se selecciona una sola vez y se informa explícitamente
- **AND** si tampoco hay parámetros UFF, no se emiten frames y las coordenadas previas se conservan

#### Scenario: RDKit no disponible o falla
- **GIVEN** el worker no arranca, expira o devuelve error
- **WHEN** el controlador recibe el fallo
- **THEN** el dock lo presenta con diagnóstico estable
- **AND** la geometría committed, el canvas 2D y el undo stack siguen válidos

### Requirement: El planificador limita trabajo y coalesce ediciones

La interfaz SHALL usar debounce/coalescencia para no iniciar una minimización completa por cada evento de ratón. SHALL mantener como máximo un proceso de Live MM y un último snapshot pendiente por documento. SHALL permitir interacción con la molécula mientras el worker progresa. Los frames del worker SHALL actualizar un preview en UI sin bloquear input ni render.

#### Scenario: Ráfaga de eventos de arrastre
- **GIVEN** el usuario arrastra un átomo generando múltiples posiciones en menos de 100 ms
- **WHEN** el planificador procesa la ráfaga
- **THEN** sólo el snapshot más reciente se programa después del debounce
- **AND** no se acumulan procesos o colas ilimitadas

#### Scenario: Worker y UI concurrentes
- **GIVEN** un worker de minimización activo
- **WHEN** el usuario navega, edita o cambia un valor numérico
- **THEN** los eventos de GUI siguen respondiendo sin esperar el cálculo
- **AND** el nuevo estado se conserva como fuente de la siguiente revisión

### Requirement: La minimización sujeta el átomo bajo control del usuario

Durante un arrastre, el worker SHALL fijar mediante restricciones los IDs que el usuario está manipulando, utilizando el mecanismo `AddFixedPoint` del force field RDKit. El preview SHALL mantener esos átomos en las posiciones enviadas por el usuario mientras otros átomos se relajan. Al terminar/cancelar el gesto, la restricción temporal se quita en el siguiente snapshot.

#### Scenario: Arrastre de átomo con Live MM
- **GIVEN** un átomo bajo el cursor y el resto de la molécula relajándose
- **WHEN** llegan frames durante el arrastre
- **THEN** el átomo fijado no se desplaza por la optimización y los vecinos sí pueden cambiar
- **AND** al soltar se permite relajar la molécula completa

#### Scenario: Pinned atom ausente del snapshot
- **GIVEN** una revisión donde el ID fijado ya no pertenece al grafo
- **WHEN** se construye la petición del worker
- **THEN** la restricción se omite de manera controlada y el snapshot obsoleto se descarta

### Requirement: Revisiones impiden que resultados viejos sobrescriban ediciones nuevas

Cada trabajo SHALL llevar `document_token`, `job_id`, `graph_revision` y `coordinate_revision`. El hilo UI SHALL validar estos campos antes de aplicar cada frame. Una edición topológica o geométrica SHALL invalidar inmediatamente el trabajo anterior y SHALL programar sólo el snapshot más reciente. Los frames de una pestaña inactiva SHALL actualizar únicamente ese documento si su revisión sigue válida; los frames de una pestaña cerrada SHALL descartarse.

#### Scenario: Edición concurrente invalida trabajo
- **GIVEN** un trabajo activo en revisión N
- **WHEN** el usuario edita y crea revisión N+1
- **THEN** los frames de N no cambian ni el preview ni la persistencia de N+1
- **AND** se cancela/reinicia el worker con el snapshot más reciente

#### Scenario: Cambiar de pestaña
- **GIVEN** dos documentos con trabajos/revisiones independientes
- **WHEN** el usuario cambia la pestaña activa
- **THEN** cada resultado sólo se aplica a su documento token y revisión coincidentes
- **AND** el renderer muestra exclusivamente el estado de la pestaña activa

### Requirement: Cancelar y cerrar detiene el proceso y cierra recursos

El usuario SHALL poder detener y aceptar el último frame válido, o cancelar y restaurar el snapshot committed previo. Al cerrar/cambiar documento, el supervisor SHALL cancelar el hijo, esperar un límite acotado, terminarlo forzosamente si no responde, cerrar pipes y suprimir emisiones tardías. No SHALL bloquear el hilo GUI durante espera indefinida.

#### Scenario: Parar y aceptar
- **GIVEN** varios frames válidos ya aplicados
- **WHEN** el usuario selecciona `Parar/aceptar`
- **THEN** el proceso se detiene y la geometría preview final se guarda mediante una sola acción undoable
- **AND** los frames anteriores no agregan entradas al historial

#### Scenario: Cancelar
- **GIVEN** preview de Live MM aún no aceptado
- **WHEN** el usuario cancela
- **THEN** se restauran las coordenadas committed anteriores sin una entrada undo
- **AND** la química y coordenadas 2D permanecen intactas

#### Scenario: Cerrar ventana con worker activo
- **GIVEN** un worker Live MM en ejecución
- **WHEN** empieza el apagado de la ventana
- **THEN** el subproceso se termina y sus señales tardías se ignoran
- **AND** la liberación del visor y sus recursos OpenGL sigue la ruta existente
