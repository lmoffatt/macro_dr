# Transcripción

- Archivo: `WhatsApp Ptt 2026-01-15 at 16.01.41.ogg`
- Duración: 2:42
- Motor: `whisper.cpp` (`whisper-cli`)
- Modelo: `ggml-large-v3-turbo-q5_0`
- Fecha: 2026-02-10

---

Bueno, las variables relevantes, como dijimos, son fundamentalmente dos. En realidad, a ver, serían para... o tres, tres o cuatro, a ver, para serían. Una es el residuo, por cada medición, residuo normalizado. El otro es el gradiente de la log likelihood. Y después necesito también, para calcular la ficha de information matrix, necesito el gradiente de la media, el gradiente de la varianza y la varianza por medición. Entonces, claro, esas serían todas las variables que necesito por cada medición. Y a partir de eso yo puedo calcular todos los parámetros que necesito. digamos que calculo para una serie de, como se dice, de muestras de simulaciones. Ahora bien, si yo ahora lo reverso, que es tratar de samplear el espacio de parámetros para una simulación fija, claro, ahí yo no me queda otra que, Sí, tener un... ¿Cómo se llama? Un Monte Carlo Markov Chain. O sea, yo podría samplear el espacio, pero no serían samples de los parámetros. Sería una exploración del espacio de parámetros y ver la log likelihood por ahí. No sería exactamente un sample. el espacio de parámetros. Para eso tengo que hacer la Montecarlo Markov Chain. Porque ahí me queda un poco la duda. Si yo me quiero limitar a esto nada más, sí tengo que quedarme ahí, me parece. Bueno.
