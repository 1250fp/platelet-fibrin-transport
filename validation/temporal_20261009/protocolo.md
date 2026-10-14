# Prueba temporal confirmatoria del transporte

## Pregunta y alcance

¿El paso utilizado en el minipaper, 0.015625 ms, produce S(100 s) y RMST(100 s) prácticamente equivalentes a 0.0078125 ms en tres geometrías representativas ya guardadas?

Es una verificación numérica del transporte congelado. No calibra el tiempo mecánico, la densidad de plaquetas, la representación 3D–2D ni la validez biológica de las capas.

## Diseño fijado antes de ejecutar

| Estado | Motivo |
|---|---|
| rep1_f0_0_tm0 | Estado inicial, antes de la evolución mecánica. |
| rep1_f0_2.2_tm10 | Misma geometría que en el piloto temporal anterior. |
| rep9_f0_8.8_tm15 | Estado mecánico a 15 s, fuerza mayor y otra réplica. |

Son tres estados geométricos procedentes de **dos réplicas mecánicas**, no tres réplicas independientes. No se agrupan los trazadores de diferentes estados para estrechar sus intervalos.

- Tres pasos: 0.03125, 0.015625 y 0.0078125 ms.
- 8192 posiciones válidas por geometría; las mismas en los tres pasos. 73728 trayectorias en total.
- Nuevas semillas para las posiciones y el ruido; no se reutilizan resultados exploratorios como si fueran confirmatorios.
- Incrementos brownianos anidados: el paso doble suma dos incrementos del fino; el cuádruple suma cuatro.
- Geometrías, radio 0.25 µm, D=0.4 µm²/s, 50 capas, H=100 nm y frontera lateral intactos.
- Horizonte de transporte fijo: 100 s. Ocho procesos externos, sin paralelismo interno.
- Ninguna evolución mecánica nueva ni modificación del número de plaquetas.

## Criterios

$$\Delta S=S_{0.015625}(100)-S_{0.0078125}(100).$$

$$\delta_R=\frac{\mathrm{RMST}_{0.015625}(100)-\mathrm{RMST}_{0.0078125}(100)}{\mathrm{RMST}_{0.0078125}(100)}.$$

Los intervalos completos deben quedar dentro de ±2 puntos porcentuales y ±3%, respectivamente. Se conservan los IC marginales del 97.5% del piloto previo: cotas binomiales exactas conservadoras para la discordancia y 20000 remuestreos pareados para RMST. Cobertura conjunta nominal del 95% para los dos indicadores de cada comparación; no se afirma cobertura simultánea de todos los IC entre geometrías.

El paso actual solo se considerará apoyado en este conjunto si ambos criterios primarios pasan en cada una de las tres geometrías. El contraste 0.03125–0.015625 ms es diagnóstico. No se exige que un paso más grueso pase para aceptar el paso intermedio. Las curvas S(t) serán descriptivas, no una prueba de equivalencia uniforme.

## Por qué no solo 1024

La discordancia observada en el piloto terminado fue 59/256. Manteniendo el método conservador de IC y suponiendo diferencia real nula, la semianchura orientativa es 5.08 pp para 1024, 2.52 pp para 4096 y 1.77 pp para 8192. Por eso se fija 8192; esto no asegura aprobación si existe sesgo o mayor discordancia. Es una planificación aproximada de precisión, no un cálculo que garantice potencia estadística.

## Decisión y límites

- **Cumple:** respaldar 0.015625 ms solo para las geometrías, radio, difusión, fronteras y horizonte comprobados.
- **Inconcluso:** distinguir falta de precisión de una diferencia relevante; no aumentar automáticamente muestra ni márgenes.
- **Fuera del margen:** revisar el integrador/refinar el paso; no actualizar resultados de producción automáticamente.
- Un fallo numérico detiene el análisis y nunca cuenta como censura.
- Las incertidumbres de estereología y conectividad permanecen pendientes por separado.

## Trazabilidad

`protocol.json` fija posiciones, semillas, parámetros, criterios y hashes. `controls.json` recoge replays bit a bit del ejecutor anterior y controles del acoplamiento browniano. `status.json` informa solo del estado operativo. Los resultados se analizan una sola vez al completar todos los bloques.

Fuentes metodológicas: [Lakens, 2017](https://doi.org/10.1177/1948550617697177), sobre equivalencia y márgenes prefijados; [Giles, métodos Monte Carlo multinivel](https://people.maths.ox.ac.uk/gilesm/mlmc.html), sobre acoplar resoluciones mediante incrementos brownianos compartidos. Aquí se utiliza ese acoplamiento para comparar pasos; no se implementa un estimador multinivel.
