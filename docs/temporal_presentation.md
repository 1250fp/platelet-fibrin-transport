# Convergencia temporal: guion para cinco minutos

## Mensaje central

El paso de transporte de **0.015625 ms** satisface equivalencia práctica frente a **0.0078125 ms** en los tres estados examinados, **bajo el criterio revisado** de ±2.5 puntos porcentuales para supervivencia y ±3% para RMST. El margen de supervivencia se revisó después del análisis: con el margen original de ±2 pp, un estado permanece inconcluso.

## 0:00-0:45 | Qué estamos verificando

«La pregunta de esta prueba es numérica: si reduzco el paso temporal, ¿cambia de manera relevante el escape de los trazadores? El paso temporal no es la duración de la contracción ni el tiempo total de observación. Es el intervalo con que el algoritmo avanza el movimiento browniano. Aquí las redes permanecen congeladas y observamos el escape durante cien segundos. Buscamos separar un cambio asociado al método numérico de un cambio que podría atribuirse a la geometría del coágulo».

## 0:45-1:35 | Diseño de la comparación

Mostrar el gráfico de supervivencia, sin interpretar todavía diferencias entre estados.

«Reutilizamos tres estados geométricos: uno inicial de la réplica uno; otro de esa misma réplica a diez segundos, con fuerza basal de 2.2 nanonewtons; y un tercero de la réplica nueve a quince segundos, con fuerza basal de 8.8 nanonewtons. Son tres estados de dos réplicas, no tres réplicas independientes. Probamos tres pasos: 0.03125, 0.015625 y 0.0078125 milisegundos. En cada estado utilizamos 8192 trazadores por paso. No repetimos la evolución mecánica».

## 1:35-2:15 | Por qué la comparación es pareada

«Dentro de cada estado, las tres versiones parten de las mismas posiciones. Además, el desplazamiento browniano de un paso grueso se construye sumando incrementos del fino. Así comparamos las resoluciones usando un mismo camino browniano subyacente, aunque las trayectorias reflejadas pueden separarse después de tocar obstáculos. La comparación principal es entre el paso usado en el trabajo, 0.015625 milisegundos, y uno dos veces menor. El paso mayor sirve como comparación adicional».

## 2:15-3:00 | Qué medimos

«La supervivencia S(t) es la fracción que todavía no escapó por ninguna cara axial. El tiempo medio restringido, o RMST, es el promedio del tiempo observado hasta escapar, limitado a cien segundos. Los trazadores que no escapan aportan cien segundos: no los excluimos del promedio. Las curvas de supervivencia de los tres pasos están próximas, pero mirar curvas parecidas no basta. Necesitamos estimar las diferencias y su incertidumbre».

## 3:00-4:10 | Resultado y lectura de las barras

Mostrar el gráfico de equivalencia.

«Cada punto es la diferencia entre el paso actual y el más fino; cada barra es su intervalo de confianza marginal del 97.5%. La zona verde muestra las tolerancias aceptadas: ±2.5 puntos porcentuales para supervivencia y ±3% para RMST. En todos los estados las barras quedan dentro de esas zonas. El caso limitante es el estado a diez segundos: su diferencia de supervivencia es de 0.806 puntos y su intervalo va de menos 0.782 a más 2.392 puntos. Para RMST, los tres intervalos quedan ampliamente dentro del tres por ciento».

«Es importante aclarar que el margen de supervivencia original era ±2 puntos, indicado por las líneas discontinuas. Con ese criterio el estado intermedio era inconcluso. La aceptación con ±2.5 es una revisión posterior del criterio operativo; no modificamos los datos ni los intervalos».

## 4:10-5:00 | Conclusión y límites

«Concluimos que el paso utilizado presenta estabilidad temporal práctica, dentro del criterio revisado, para estos estados y estas condiciones. No significa que el error sea cero ni que hayamos probado el límite matemático de paso nulo. Tampoco certifica automáticamente todos los estados de producción, otros radios de trazador o redes que cambien durante el transporte. Por último, esta prueba no valida el número de plaquetas ni la conversión de tres a dos dimensiones: esas son decisiones físicas del modelo. Por ahora conservamos el paso actual y reportamos de manera explícita el margen bajo el que lo aceptamos».

## Ecuaciones de apoyo (no hace falta explicar todas en cinco minutos)

Para N trazadores con tiempo de escape τ y horizonte T = 100 s:

$$
\widehat S(t)=\frac{1}{N}\sum_{i=1}^{N}\mathbf 1\{\tau_i>t\}.
$$

$$
\widehat{\operatorname{RMST}}(T)=\frac{1}{N}\sum_{i=1}^{N}\min(\tau_i,T)
=\int_0^T\widehat S(t)\,dt.
$$

La expresión del promedio se aplica aquí porque el único cierre administrativo es el horizonte común T; no se usa la media solo entre quienes escaparon.

$$
\Delta S_{\rm pp}=100\,[\widehat S_{\Delta t}(100)-\widehat S_{\Delta t/2}(100)].
$$

$$
\delta_{\rm RMST}(\%)=100\,
\frac{\widehat{\operatorname{RMST}}_{\Delta t}(100)-\widehat{\operatorname{RMST}}_{\Delta t/2}(100)}
{\widehat{\operatorname{RMST}}_{\Delta t/2}(100)}.
$$

El criterio revisado exige contener el **intervalo completo**, no solo el punto:

$$
\operatorname{IC}(\Delta S_{\rm pp})\subset(-2.5,2.5),\qquad
\operatorname{IC}(\delta_{\rm RMST})\subset(-3,3).
$$

## Datos que conviene tener a mano

| Estado | ΔS e IC (pp) | ΔRMST e IC relativo (%) |
|---|---:|---:|
| A: R1, inicial | 0.024 [-1.398, 1.447] | 0.338 [-0.629, 1.317] |
| B: R1, 2.2 nN, 10 s | 0.806 [-0.782, 2.392] | 0.318 [-0.682, 1.339] |
| C: R9, 8.8 nN, 15 s | -0.049 [-1.424, 1.326] | 0.149 [-1.067, 1.385] |

- 50 capas; H = 100 nm; espesor total = 5 µm.
- Radio del trazador = 0.25 µm; D = 0.4 µm²/s.
- Reflexión en sólidos y laterales; absorción del centro en z = 0 y z = L.
- Liberación central condicionada a espacio válido, no a conectividad.
- 73728 trayectorias = 3 estados × 3 pasos × 8192; no son 73728 redes independientes.
- IC de supervivencia: cotas binomiales exactas conservadoras para pares discordantes.
- IC de RMST: 20000 remuestreos de índices pareados; cobertura aproximada.
- Cobertura conjunta nominal del 95% para los dos indicadores de una comparación; no se afirma cobertura simultánea de todos los estados.
- Ninguna de estas bandas demuestra equivalencia uniforme de toda la curva S(t).

## Preguntas previsibles

**¿Por qué 2.5 pp?** Es una tolerancia operativa revisada después de conocer los resultados, no un umbral biológico ni una recomendación numérica universal. Se conserva la conclusión original con 2 pp. A una supervivencia del 25%, 2.5 pp representa un 10% relativo.

**¿Por qué una barra contiene cero?** La estimación admite diferencias de ambos signos. Para equivalencia importa que todo el intervalo entre en el margen, no demostrar una diferencia respecto a cero.

**¿Las barras representan la variabilidad entre coágulos?** No: la incertidumbre es condicional a cada geometría. Los tres estados proceden de dos réplicas mecánicas.

**¿Por qué no basta con mirar curvas superpuestas?** La escala visual puede ocultar diferencias y no cuantifica su incertidumbre.

**¿Se puede usar un paso todavía mayor?** No se cambió el paso de producción. El contraste grueso-intermedio es favorable en estos estados, pero el objetivo principal era comprobar el paso actual contra una referencia más fina.

## Fuentes y archivos

- [Lakens (2017), Equivalence Tests: A Practical Primer](https://doi.org/10.1177/1948550617697177): explica la equivalencia y la elección explícita de márgenes; no prescribe nuestros valores ni convierte una revisión posterior en una prueba prefijada.
- [Giles, métodos Monte Carlo multinivel](https://people.maths.ox.ac.uk/gilesm/mlmc.html): referencia para el acoplamiento de resoluciones mediante incrementos brownianos compartidos. Aquí no se implementa un estimador multinivel.
- [Gráfico de supervivencia](../validation/temporal_20261012_revision/01_supervivencia.png).
- [Gráfico de equivalencia](../validation/temporal_20261012_revision/02_equivalencia.png).
- [Suplemento compilado](../paper/validation_supplement.pdf).
- [Informe del criterio revisado](../validation/temporal_20261012_revision/informe_revisado.md).
- El [protocolo original del 9 de octubre](../validation/temporal_20261009/protocolo.md) y sus resultados permanecen intactos en `validation/temporal_20261009`. La revisión se guarda por separado.
