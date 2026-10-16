# Propuesta — uniones editoriales de enlaces moleculares

## Por qué

La representación 2D ya conserva la geometría y la semántica química del grafo, pero el acoplamiento de trazos en encuentros con cuñas puede producir extremos visualmente abruptos. La lógica de `BondItem` incluye adaptaciones específicas para la cuña sólida y tramos de cuña rayada, mientras que el contexto de uniones del canvas reduce cada enlace vecino a un vector y un ancho escalar. Hace falta una campaña visual reproducible que distinga uniones correctas de defectos reales y permita corregir solamente los encuentros demostrablemente defectuosos.

## Qué cambia

- Añadir una galería determinista de moléculas sintéticas y estructuras cíclicas con uniones simples, dobles, triples, bold, cuñas sólidas/rayadas, bifurcaciones y proximidad a etiquetas/cargas.
- Medir y documentar los casos base en temas claro/oscuro antes de ajustar geometría.
- Añadir un módulo geométrico interno únicamente si el diagnóstico demuestra que el cálculo actual necesita representar explícitamente perfiles de trazos (incluidas las líneas sigma/pi de dobles enlaces).
- Realizar ajustes estrechos de los contornos de cuña y terminaciones rayadas sólo si una comparación visual y pruebas geométricas confirman una mejora sin regresiones.
- Añadir pruebas de renderizado/invariantes que comprueben determinismo, seguridad de etiquetas y preservación del grafo químico.

## Fuera de alcance

No se mueven átomos ni se cambia Clean2D; no se modifica la química, estereoquímica, conectividad, órdenes, cargas, aromaticidad, persistencia, undo/redo, import/export o interacción del editor. No se cambia la semántica wedge/hash ni el significado independiente de `bold`. No se reescribe `BondItem`. No se crea el preset `ChemUSON Publication` salvo que las mejoras geométricas anteriores queden validadas en los cuatro estilos de contexto y exportación; esa condición no se presume.

## Criterio de aceptación

Una mejora geométrica debe tener un ejemplo base reproducible, una comparación antes/después a igual escala en ambos temas, pruebas de límites y conectividad, y ausencia de regresiones en renderizado, anillos aromáticos/Haworth, Clean2D y arquitectura. El pase automatizado no constituye aceptación visual del propietario.
