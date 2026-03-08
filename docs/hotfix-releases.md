# Hotfix Releases

Fecha: 2026-03-08

## Objetivo

Permitir publicar correcciones pequenas con la menor friccion posible para:

- testers del canal `beta`
- usuarios finales del canal `stable`

sin compilar localmente ni preparar artifacts a mano.

## Limitacion importante

El updater de Chemuson **no distribuye commits sueltos**. Solo ve **releases publicadas** con una **version mas nueva** que la instalada.

Eso significa:

- un commit en Git no llega por si solo a los usuarios;
- un asset corregido con la misma version tampoco se ofrecera como update;
- cada hotfix necesita una version nueva y artifacts nuevos.

## Flujo recomendado

### 1. Hotfix para testers

1. Haces commit/push del arreglo a la rama que quieras publicar.
2. En GitHub Actions ejecutas manualmente la workflow `release`.
3. Usas por ejemplo:
   - `version = 0.2.2-beta.1`
   - `channel = beta`
4. La workflow:
   - sincroniza `src/chemuson/_version.py` dentro de CI,
   - compila Windows/Linux/Flatpak,
   - publica una **prerelease** `v0.2.2-beta.1`.
5. Los testers con canal `beta` podran recibirla desde `Ayuda -> Buscar actualizaciones...`.

### 2. Promocion a estable

Cuando el hotfix ya fue validado por testers:

1. Ejecutas otra vez la workflow `release`.
2. Usas por ejemplo:
   - `version = 0.2.2`
   - `channel = stable`
3. Se publica la release estable y los usuarios del canal `stable` la veran.

## Como lo usan los testers

Los testers deben tener el canal `beta` en preferencias de actualizacion.

Comportamiento esperado:

- `stable` solo ve releases estables;
- `beta` ve tanto prereleases beta como releases estables mas nuevas.

## Que automatiza ahora CI

La workflow `release` ya no depende de que la version del repositorio coincida manualmente con el nombre del release.

Durante la ejecucion:

- toma la version indicada en la workflow o en el tag,
- actualiza `_version.py`,
- construye artifacts con esa misma version,
- publica la release/tag correspondiente.

## Convencion sugerida

- Beta para testers: `X.Y.Z-beta.N`
- Hotfix estable: `X.Y.Z`

Ejemplo:

- testers: `0.2.2-beta.1`
- testers con ajuste extra: `0.2.2-beta.2`
- salida estable final: `0.2.2`

## Recomendacion operativa

Si quieres mover cambios pequenos muy rapido:

- usa `beta` como canal de entrega inmediata,
- valida con testers,
- promociona el mismo hotfix a `stable` cuando quede confirmado.
