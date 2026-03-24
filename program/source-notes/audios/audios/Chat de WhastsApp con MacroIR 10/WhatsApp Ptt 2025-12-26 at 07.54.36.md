# Transcripción

- Archivo: `WhatsApp Ptt 2025-12-26 at 07.54.36.ogg`
- Duración: 2:52
- Motor: `whisper.cpp` (`whisper-cli`)
- Modelo: `ggml-large-v3-turbo-q5_0`
- Fecha: 2026-02-10

---

Bueno, esta mañana otra vez me puse a pensar en el gráfico este de la relación entre la actualización de la distribución de probabilidad con una sola medida en el caso de los dos algoritmos, el MR y el IR, porque me di cuenta que en realidad yo podía hacer una aproximación a este, basada en simular la distribución usada por el algoritmo, es decir, P minúscula por P mayúscula, es decir, no, P min por P, es decir, la media del vector el estado de probabilidad por la probabilidad de transición, o sea, la diagonal del vector de probabilidad por la probabilidad de transición. O sea, ahí me quedaría una especie de matriz que te dice la probabilidad de cada combinación. Y claro, yo en realidad podría aplicarlo a eso para cada combinación de números de canales abiertos. y bueno, en generar las muestras, así que sé yo. Pero igual, de todos modos, termino recreando el microscopic recursive, el algoritmo microscopic recursive, así que decidí abandonar eso. La intuición es que vos te quedas así, o sea, vos cuando quieras hacer el microscopic recursive, claro, va a estar sobre, va a ser aplicado el microscopic recursive, no sobre el estudio, sobre el metaestado. O sea que tengo que expandir el metaestado de aplicar el microscopic recurso y sobre el metaestado y no lo puedo... no puedo hacer el truco que hago con las instrucciones normales de marginalizar porque la marginalización de las multivariadas, de las multinomiales es una pesadilla. Lo cual tendría que preguntar, a ver...
