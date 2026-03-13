# Transcripción

- Archivo: `WhatsApp Ptt 2025-12-26 at 08.05.27.ogg`
- Duración: 2:16
- Motor: `whisper.cpp` (`whisper-cli`)
- Modelo: `ggml-large-v3-turbo-q5_0`
- Fecha: 2026-02-10

---

Cambios en macro IR en cuanto a la interfaz. Entonces la idea es que ahora tengo que habilitar la posibilidad de repeticiones. Es decir, que en lugar de simular una vez, simulo n veces. O sea, si yo tengo un experimento, hago una simulación, en lugar de hacer una simulación única, hago n simulaciones, cada una con su propio seed value. O sea, ahí una de las cuestiones es, claro, sí, Cada una genera un seed value independiente. Entonces, hay una cuestión más que estaba pensando, que estaría buena, y es para tener una especie de cadena de reproducibilidad que, en realidad, cuando yo haga cualquier tipo de medición no trivial, como por ejemplo el likelihood, que tenga que referirme a un file siempre. O sea, que yo tomo el likelihood y tomo el nombre de un file. ¿Cuál es la ventaja de eso? La ventaja de eso es que luego me queda una trazabilidad de dónde vienen los datos, que es independiente de tener exactamente lo que se llama el protocolo, el script. O sea, la idea sería que me puedo independizar del script y yo solamente con las data frames pueda reproducir qué es lo que hago. Es una especie de reproducibilidad a partir de data frames. Sería una especie de principio. Yo creo que eso estaría bueno. A ver, podría tratar de pensarlo.
