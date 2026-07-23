## Context

`ChemusonWindow.__init__` coordina actualmente cinco grupos de trabajo:

1. estado, preferencias y controllers;
2. pestañas y canvas inicial;
3. docks y trabajos en segundo plano;
4. toolbars, menús y conexiones;
5. status bar, tema y chequeo diferido de actualizaciones.

Las factorías de acciones y `MainWindowUiBuilder` ya separan parte del detalle,
pero el orden completo permanece embebido en la clase de 2603 líneas.

## Goals / Non-Goals

**Goals:**

- Hacer visible y revisable la composición del shell en `gui/shell/`.
- Reducir el constructor de `ChemusonWindow` a inicialización base y una única
  delegación.
- Preservar orden, atributos, parentage Qt, conexiones y efectos observables.
- Mantener `MainWindowUiBuilder` como colaborador existente.
- Crear una costura estable para extraer regiones posteriores sin ampliar esta
  fase.

**Non-Goals:**

- Rediseñar la interfaz o corregir recortes visuales.
- Cambiar nombres públicos o handlers de `ChemusonWindow`.
- Mover operaciones de archivo, Clean2D, nomenclatura, compchem, templates o
  validación.
- Dividir los monolitos del canvas o `items.py`.
- Introducir DI, service locator, mixins nuevos o un framework de componentes.

## Decisions

### M08 owns a shell package

`src/chemuson/gui/shell/` pertenecerá a M08. Su API será interna y recibirá la
instancia existente de `ChemusonWindow`; no creará una segunda subclase ni
cambiará la ruta pública de la ventana.

### Extract composition, not behavior

El ensamblador moverá literalmente las asignaciones y conexiones del
constructor. Los métodos invocados y callbacks seguirán definidos en
`ChemusonWindow` o en colaboradores ya existentes. Ningún handler de química,
documento o canvas se trasladará en esta fase.

### Preserve ordering

El orden del constructor es un contrato: acciones antes de widgets que las
consumen; canvas antes de docks dependientes; menús antes de sincronización;
señales antes del estado inicial; tema antes del chequeo diferido. La revisión
del diff y pruebas de caracterización impedirán reordenamientos accidentales.

### One delegation seam

Después de `QMainWindow.__init__`, `ChemusonWindow.__init__` creará o invocará
un único ensamblador del shell. No conservará una copia parcial del bloque ni
añadirá aliases de compatibilidad.

### Existing builders remain authoritative

`MainWindowUiBuilder` seguirá construyendo menús y barra superior; las
factorías en `gui/actions/` seguirán creando acciones. El nuevo paquete sólo
coordina regiones y conexiones, evitando duplicar responsabilidades.

## Risks / Trade-offs

- [Orden alterado] → comparación AST y pruebas GUI existentes conservan
  atributos y secuencia.
- [Acoplamiento circular shell↔main_window] → el shell no importa
  `ChemusonWindow`; trabaja contra el objeto recibido.
- [God object trasladado] → sólo se extrae ensamblaje; los handlers no migran y
  las siguientes extracciones deberán tener OpenSpecs propios.
- [Cambio visual] → se conservan constructores, parents, áreas Qt y orden.

## Migration Plan

1. Capturar baseline e inventariar el constructor.
2. Añadir caracterización arquitectónica del seam.
3. Crear `gui/shell/` y trasladar literalmente la composición.
4. Reducir `ChemusonWindow.__init__` a delegación.
5. Actualizar M08 y documentación.
6. Ejecutar pruebas focalizadas, arquitectura, OpenSpec y suite disponible.

Rollback: revertir el cambio devuelve el bloque al constructor; no existe
migración de datos ni configuración.
