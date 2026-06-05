**Archivo:** WhatsApp Ptt 2026-04-17 at 10.27.06.mp3
**Duración:** ~09:24
**Hablante:** Luciano

[00:00] Bueno, Macro IR, digamos, cierre, capítulo final, no sé cómo llamarlo.

[00:07] Estoy acá en la reserva, en el parque de Andes, mirando una aviota.

[00:16] Pasar. Estoy en el muelle circular.

[00:22] Ahí pasa un avión.

[00:28] Bueno, vine a la facultad y estaba dudando si iré al parque, la memoria del parque, la reserva para hablar al celular o no.

[00:41] Y bueno, al final dije que sí, porque tengo que pensar cómo hacer el cierre del paper.

[00:50] Entonces, nada, pensándolo al principio pensé, bueno, tengo que empezar a filosofar acerca de qué trataba este paper.

[00:58] Bueno, el paper trataba de, digamos, es la continuación, si se quiere, de mi paper original de Bury's Biblical Journal.

[01:09] Y el que saqué en el Communication Biology.

[01:16] En el sentido de que en el Microsoft Biology mostré que el algoritmo de este macro IR sirve para obtener datos y ahora lo que quiero es estudiar su validez y además presentarlo bien, presentar su deducción.

[01:34] qué significa este algoritmo de meta estado no sé cómo llamarlo que tengo que hacer una de las cosas que tengo que hacer es pulir mi nomenclatura

[01:46] o sea tengo que presentar mi nomenclatura de digamos del estado contorno no me acuerdo como se llamaba el estado

[01:59] el "boundary state", que es el estado que tiene que ver con los límites del intervalo de mención.

[02:14] Está representado el estado en los dos estados frontera, donde empieza y termina el intervalo de mención,

[02:27] que ahí sí es un punto, no es una combinación de todos los puntos posibles, sino que vos tomas

[02:36] dos estados puntuales, inicial y final, y dados esos dos estados puntuales, bueno, vos consideras

[02:44] todas las trayectorias posibles, no tratás de obtener información acerca de cuáles

[02:59] de todas las trayectorias posibles es la que ocurre, sino que vos te interesa solamente

[03:05] el estado inicial y el final.

[03:07] Porque claro, el estado inicial y el estado final, especialmente el estado final, te va

[03:12] va a determinar después cómo continúa la trayectoria, por la propiedad marcoviana,

[03:22] de que en un estado puntual está definido todo el estado del sistema, o sea, que no

[03:28] tenés memoria, que tu pasado está codificado en tu presente, o no sé si codificado, pero

[03:40] futuro en realidad está codificado tu presente, para predecir el futuro, vos necesitas solamente

[03:48] el estado presente, eso tenes una buena predicción.

[03:54] Y bueno, planteo eso como una estrategia de análisis de sistemas marcovianos, digamos

[04:09] continuos y en los cuales tus mediciones no son mediciones instantáneas sino que

[04:17] representan el estado promedio en el intervalo de tiempo

[04:22] una medición promedio en el intervalo

[04:32] Ahí tenemos la salvedad de que, bueno, no estamos suponiendo un tipo de integración exacta, ¿no?

[04:42] Desde un punto hasta el otro punto, cuando en realidad las mediciones usan filtros y bueno,

[04:50] y para entender los filtros y todo hace falta otro tipo de análisis, que eso va a ser otro manuscrito,

[04:58] que es muy importante, es otro manuscrito, pero bueno, porque se va a permitir estudiar con precisión,

[05:08] o por lo menos sí, precisión algorítmica o lo que ocurre en alta frecuencia.

[05:21] Bien, pero bueno, en principio presento el algoritmo, presento el tema de que, bueno,

[05:29] el lado de estado inicial en estado final, digamos que tengo como una especie de meta estado

[05:34] que lo puedo mapear, o sea, es equivalente a un modelo con más estados

[05:42] y entonces a ese modelo puedo aplicarle el algoritmo macro R clásico

[05:49] y entonces obtengo el update y con eso obtengo las fórmulas del posterior de lo que sería

[06:07] el estado boundary. Pero después lo proyecto hacia el final, no me interesa la parte inicial,

[06:18] o sea, la puedo calcular igual para una cuestión después, si realmente me interesa

[06:25] de alguna manera reconstruir la trayectoria real del sistema, pero si solamente me interesa obtener

[06:34] la likelihood para poder después estimar los parámetros, bueno, no son los necesarios,

[06:39] ahora proyectas el estado final y al final te queda, digamos, un algoritmo muy parecido a macro R,

[06:46] pero con algunas diferencias de términos.

[06:52] Y bueno, y entonces, en realidad vos podés plantear dos algoritmos,

[06:59] uno que simplemente que sea el estado inicial, o sea que vos corregís tu estado inicial

[07:06] y después aplicás la transición marcoviana al estado inicial del intervalo,

[07:14] O el estado medio, si querés también intervalo, también podría ser.

[07:17] Pero bueno, en principio no, el estado inicial.

[07:19] Y bueno, y entonces tendríamos tres algoritmos, macro R clásico, macro MR y macro IR.

[07:32] Y bueno, y entonces los comparo con...

[07:37] Entonces la idea es ahora, sí, mostrar en qué es mejor este algoritmo.

[07:43] O sea que entonces tengo que ahora desarrollar una teoría para comparar el logaritmo.

[07:51] La teoría entonces se basa en básicamente dos propiedades.

[07:59] Una es que el score se hace, o sea la esperanza del score tiene que ser cero.

[08:05] Yo lo que te daría es si hay un bias.

[08:09] Y después la otra es que, digamos, dos estimaciones de la covarianza del gradiente y el gesiano deberían coincidir.

[08:26] y si no coinciden, digamos vos tenés un factor que te indica una inflación de varianza

[08:34] que estaría indicando de alguna manera una especie de pseudo-réplica, por llamarlo de alguna manera,

[08:40] o sea, una especie de sobreestimación de la likelihood.

[08:53] Entonces, bueno, yo calculo las dos estimaciones, veo que entonces ahora,

[09:00] y bueno, ahora yo tomo un sistema mínimo para poder entenderlo bien,

[09:05] que es una vida todo encerrado nada más,

[09:08] cosa de tener bien fija, tener una sola constante de tiempo

[09:13] y relacionar esa constante de tiempo con la constante de integración.

[09:19] y bueno y yo
