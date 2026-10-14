# Comparación numérica: seis y nueve réplicas, y extensión a 15 s

Esta comparación separa dos cambios: ampliar la cohorte y avanzar el estado mecánico.

## Cohorte en tiempos comunes

| Tiempo mecánico (s) | f0 (nN) | S seis (%) | S nueve (%) | RMST seis (s) | RMST nueve (s) |
|---:|---:|---:|---:|---:|---:|
| 0 | Todas | 27.86 | 24.22 | 59.63 | 55.97 |
| 5 | 0 | 29.30 | 25.48 | 60.18 | 56.23 |
| 5 | 2.2 | 21.61 | 21.35 | 53.10 | 51.63 |
| 5 | 4.4 | 19.66 | 17.53 | 52.47 | 50.04 |
| 5 | 8.8 | 20.05 | 18.01 | 50.82 | 48.74 |
| 10 | 0 | 30.14 | 27.08 | 59.88 | 56.91 |
| 10 | 2.2 | 21.16 | 19.92 | 53.07 | 50.56 |
| 10 | 4.4 | 20.05 | 18.45 | 50.72 | 48.82 |
| 10 | 8.8 | 16.73 | 15.54 | 47.43 | 45.57 |

## Evolución mecánica con las mismas nueve réplicas

| f0 (nN) | S tm10 (%) | S tm15 (%) | RMST tm10 (s) | RMST tm15 (s) |
|---:|---:|---:|---:|---:|
| 0 | 27.08 | 27.34 | 56.91 | 56.47 |
| 2.2 | 19.92 | 18.53 | 50.56 | 49.69 |
| 4.4 | 18.45 | 16.97 | 48.82 | 47.46 |
| 8.8 | 15.54 | 14.32 | 45.57 | 43.97 |

Las medias dan el mismo peso a cada réplica mecánica. El transporte está congelado en cada estado; su horizonte sigue siendo 100 s.

Los intervalos por contraste son bootstrap percentil de índices pareados, condicionales a cada red y al muestreo observado. No miden la incertidumbre entre coágulos ni el error numérico.

Las nuevas réplicas no resuelven por sí mismas la equivalencia temporal inconclusa ni la conectividad no resuelta. Este archivo contiene comparaciones descriptivas, no una nueva declaración de validación.
