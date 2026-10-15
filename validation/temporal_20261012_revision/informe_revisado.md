# Cierre operativo de la prueba temporal

Fecha: 12 de octubre de 2026.

## Conclusión

Se conserva el paso de transporte de **0.015625 ms**. Su comparación con **0.0078125 ms** satisface equivalencia práctica en los tres estados examinados bajo **±2.5 puntos porcentuales para S(100 s)** y **±3% relativo para RMST(100 s)**.

El margen de supervivencia se revisó después de conocer los resultados. Por tanto, esta conclusión no equivale a aprobar el criterio prefijado original de ±2 pp: con ese margen, el estado B sigue inconcluso. No se modificaron simulaciones ni intervalos para obtener la aceptación.

## Diseño y resultados

Tres estados congelados de dos réplicas mecánicas; 8192 trazadores por estado y paso; tres pasos temporales (0.03125, 0.015625 y 0.0078125 ms); observación hasta 100 s. Las comparaciones comparten posiciones iniciales e incrementos brownianos anidados dentro de cada estado. Se completaron 1152 bloques y 73728 trayectorias, sin errores numéricos.

| Estado | Diferencia S, pp [IC 97.5%] | Diferencia RMST, % [IC 97.5%] | Criterio original | Criterio revisado |
|---|---:|---:|---|---|
| A: R1, f₀ = 0, estado inicial | 0.024 [−1.398, 1.447] | 0.338 [−0.629, 1.317] | Cumple ambos | Cumple ambos |
| B: R1, f₀ = 2.2 nN, 10 s | 0.806 [−0.782, 2.392] | 0.318 [−0.682, 1.339] | S inconcluso; RMST cumple | Cumple ambos |
| C: R9, f₀ = 8.8 nN, 15 s | −0.049 [−1.424, 1.326] | 0.149 [−1.067, 1.385] | Cumple ambos | Cumple ambos |

Las diferencias son paso actual menos paso más fino; la diferencia relativa de RMST usa el valor fino como denominador. Se exige que **todo el intervalo** quede dentro del margen. La comparación secundaria entre 0.03125 y 0.015625 ms también cumple los márgenes originales en los tres estados; no demuestra convergencia monótona ni autoriza por sí sola a cambiar el paso de producción.

## Alcance

- Los IC son condicionales a cada estado. Son marginales del 97.5%, con cobertura conjunta nominal del 95% para los dos indicadores de una comparación, no para todas las geometrías a la vez. El bootstrap de RMST tiene cobertura aproximada.
- No se ha demostrado equivalencia uniforme de toda S(t), un límite matemático cuando el paso tiende a cero, ni estabilidad para todos los estados de producción.
- La prueba usa 50 capas de 100 nm, radio del trazador 0.25 µm y difusión 0.4 µm²/s. No valida automáticamente otros radios, parámetros o transporte durante contracción.
- El umbral de 2.5 pp es operativo, no biológico: alrededor de S = 25%, supone una tolerancia relativa del 10%. No hay una calibración científica independiente del resultado que determine ese número.
- No resuelve las dudas sobre número de plaquetas, estereología, parámetros mecánicos ni certificación espacial de conectividad.

## Cambios documentales

Se actualizaron el resumen y la sección de verificación del minipaper, la discusión y la cohorte S8 del suplemento. La tabla ampliada conserva la diferencia entre el criterio original y el revisado.

Los [resultados originales](../temporal_20261009/results.json) y el
[informe original](../temporal_20261009/informe.md) se conservan por separado de
[results_revised_margin.json](results_revised_margin.json). No se repitió mecánica
ni transporte para revisar el margen.

## Figuras y documentos

- [Supervivencia](01_supervivencia.png).
- [Equivalencia con margen revisado y referencia al original](02_equivalencia.png).
- [Suplemento compilado](../../paper/validation_supplement.pdf).
- [Protocolo original](../temporal_20261009/protocolo.md).
