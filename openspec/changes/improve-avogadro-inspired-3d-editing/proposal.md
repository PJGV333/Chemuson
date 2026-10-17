# Proposal — edición molecular 3D directa inspirada en Avogadro

## Propósito

Convertir el visor 3D existente en un espacio de edición molecular directa, no solo en una vista de conformeros y controles numéricos. Recuperar contratos de interacción útiles de Avogadro 1.2 mediante una implementación original en Python/PyQt6, OpenGL 4.1 y RDKit aislado, integrada con el editor 2D de ChemUSON.

## Alcance por hitos

- **A — Sincronización 2D/3D:** detectar cambios químicos efectivos de `MolGraph`, invalidar geometrías y trabajos incompatibles, aplicar debounce y generar automáticamente la escena vigente; preservar conformeros ante cambios exclusivos del layout 2D.
- **B — Hidrógenos editables:** permitir picking, resaltado y movimiento de H derivados mediante identidad visual tipada, separada de IDs documentales. Al confirmar, promover el H seleccionado a átomo `H` explícito con enlace, coordenada y undo/redo en una sola transacción química. Cancelar conserva el estado documental.
- **C — Edición directa:** ofrecer selección, construcción, manipulación de átomos/fragmentos y enlaces, medidas/indicadores geométricos, ayuda contextual y navegación diferenciada; mantener controles numéricos como ajuste preciso.
- **D — Lectura visual:** estabilizar iluminación en espacio de cámara y añadir representaciones básicas compatibles con picking: alambre, bolas y varillas y espacio lleno aproximado.
- **E — Live MM integrado:** usar el worker RDKit aislado y streaming existente durante gestos/construcción, fijar temporalmente el conjunto manipulado, descartar resultados obsoletos y confirmar solo una transacción undoable por acción.
- Crear matriz comparativa basada en fuente Avogadro 1.2.0 y hoja de ruta separada para superficies científicas avanzadas.

## Invariantes

`MolGraph` es la autoridad química compartida; coordenadas 2D y 3D permanecen separadas; las identidades derivadas no se convierten implícitamente en IDs documentales; persistencia, estereoquímica, aromaticidad y valencias conservan sus contratos existentes; cada pestaña y cada resultado asíncrono permanecen aislados. No copiar, traducir ni adaptar literalmente código GPL de Avogadro.

## Fuera de alcance

No añadir dependencias ni paquetes, no cambiar el formato `.cmsn`, no modificar Clean2D/ChemName, no ejecutar binarios Avogadro, no construir superficies electrónicas/energéticas, no elevar privilegios y no declarar aceptación visual o cierre Qt/SIP a partir de pruebas offscreen. No repetir sin causa nueva la suite monolítica cuyo baseline conocido termina en SIGSEGV.

## Criterio de éxito

Cada hito queda cubierto por regresiones que ejercitan comportamiento y wiring de UI, no solo una API. Los cambios de conectividad actualizan el viewport sin reabrirlo; un H derivado puede seleccionarse y promoverse con undo/redo sin duplicación; las herramientas directas resuelven gestos de usuario; iluminación y picking concuerdan con el estilo; Live MM no bloquea ni aplica frames obsoletos. La validación visual manual y el SIGSEGV Qt/SIP se informan por separado.