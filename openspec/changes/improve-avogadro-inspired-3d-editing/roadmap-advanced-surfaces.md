# Hoja de ruta — superficies moleculares avanzadas

Este documento separa capacidades científicas avanzadas de las herramientas de edición de moléculas. Ninguna superficie se simula mediante un color arbitrario ni por reutilizar esferas vdW.

| Capacidad futura | Contrato requerido | Datos/backend necesarios | Representación/validación pendiente |
|---|---|---|---|
| Superficie vdW/solvente | Definir radios, densidad de muestreo, probe solvente, cavidades y normalización | Malla de superficie a partir de geometría y radios; considerar algoritmo de solvent-accessible/surface | Control de resolución, orientación de normales, picking/selección de superficie, referencias visuales |
| Densidad electrónica | Identificar orbitales/función de onda y unidad espacial; documentar aproximaciones | Backend de química cuántica capaz de producir grid volumétrico con metadatos | Isovalores, mapas de color con escala, convergencia y validación frente a datos del backend |
| Potencial electrostático | Definir malla, unidades, carga/estado electrónico y convención de signo | Cargas/función de onda y solver electrostático con trazabilidad | Mapa continuo, leyenda numérica, casos de referencia y no confundir con colores CPK |
| Superficie de energía potencial | Declarar coordenadas colectivas, rango, algoritmo y campo/estado químico | Cálculos repetidos de backend explícito; probablemente molécula/solvente/temperatura | Malla/contours etiquetados, incertidumbre, convergencia y comparación numérica |

Antes de implementar cualquiera: especificar unidades, procedencia, precisión/costo, cancelación asíncrona, almacenamiento/serialización, derechos de dependencias, pruebas numéricas de referencia y contrato de picking. No forma parte de los hitos A–E.
