**Archivo:** WhatsApp Ptt 2026-04-21 at 17.22.23.mp3
**Duración:** ~14:37
**Hablante:** Luciano

[00:00] En la 9 de julio, bueno, ayer creo que fue, pasé uno de los días más intensos con Macro IR que recuerde,

[00:12] porque claro, me puse a ver una pregunta que era,

[00:17] ¿cómo dependía del largo del intervalo?

[00:21] Sí, el largo del intervalo de integración, es decir, el número de mediciones que tomaba por intervalo de medición.

[00:30] cómo aumentar ese número me aumentaba la resolución, la sensibilidad en los parámetros,

[00:39] es decir, cuánto podía reducir la incerteza en los parámetros.

[00:44] Y bueno, para mi sorpresa absoluta, el error de los parámetros corregido por la covarianza del score me daba constante.

[01:05] es decir, básicamente, con un error importante, porque es una medida bastante ruidosa, pero

[01:14] bueno, se veía que vos tomaras una medida, una medición cada un tau o 100, el error

[01:25] es prácticamente el mismo.

[01:27] O sea, con lo cual, digamos, lo que es la variable esta que yo estaba estimando,

[01:36] que es la, ¿cómo se llama?

[01:39] La distorsión de la información, la matriz de distorsión de la información

[01:46] que te indica cómo se infla la covarianza con, digamos,

[01:51] cómo se infla la covarianza.

[01:55] en realidad sale de dividir una constante con la matriz de información de Fisher.

[02:04] Es decir, como que hay un mecanismo de compensación perfecto,

[02:09] pero que no pude formular matemáticamente.

[02:13] Es como que, claro, digamos, al principio yo decía, no puede ser.

[02:17] Si yo mido más veces, tengo que tener más información.

[02:20] Pero claro, hay casos donde esto no ocurre.

[02:24] Por ejemplo, si vos tenés un proceso de poisson y vos me dicen en un intervalo de tiempo,

[02:32] que yo, cien cuentas, por ejemplo, en una hora,

[02:34] vos a ver si tenés en los primeros diez minutos, tuviste diez cuentas, ocho, cinco, lo que sea,

[02:40] no te cambia nada.

[02:42] Lo único que importa es cuántas cuentas tuviste al final.

[02:45] Es decir, que lo que te vale es el acumulado.

[02:50] como eso se distribuye a lo que es estimar el parámetro no te afecta.

[02:56] Aunque sí te afecta es si vos comparás ese modelo con otros modelos más complicados.

[03:01] O sea que en este caso habría lo mismo, ¿no?

[03:04] Porque tenemos un modelo que es un poco más complejo que un proceso de poisson,

[03:09] lo que tenemos es un canal que se abre y se cierra,

[03:12] pero bueno, hay solamente una constante de tiempo,

[03:15] o sea, hay dos rates, uno de la abertura y otro de clausura,

[03:17] pero te usan constante de tiempo y bueno, tenés cuatro parámetros más que dos ya se sacan por el

[03:26] período inicial que es la línea de base y el ruido y después antes el número de canales y la

[03:32] conductancia por canal y claro digamos en realidad no necesitas para obtener seis puntos necesitas

[03:40] siete puntos, ¿no? O sea, más puntos serían como innecesarios.

[03:46] No sé. Claro, traté después de reducir eso e incluso tomando, digamos, fracciones de tau,

[03:59] es decir, y tampoco me daba igual constante. Pero claro, si vos tomás fracciones de tau,

[04:04] vos decís, bueno, mi primera medición es después de 5 TAUS, claro, te va a dar prácticamente 99%,

[04:13] que se yo, no sé, pero claro, como vos reducís el ruido también, ese 99% va a tener una decisión bastante grande,

[04:22] Entonces, terminás teniendo de alguna manera la misma información.

[04:28] Es un tanto, digamos, contraintuitivo y, digamos, no sé, muy extraño como resultado.

[04:41] Pero, no sé, parecía muy robusto.

[04:45] La verdad que yo quedé medio atónito con eso, porque digo, digamos, seguro que se me escapa algo muy obvio y quedo como un boludo que no me di cuenta de eso.

[04:57] Pero bueno, no sé.

[05:01] Supongo que lo pondré como una opción.

[05:06] Honestidad intelectual ante todo.

[05:09] No, no estaría mal decirlo.

[05:15] Bueno, creo que esto sí tiene que ser un poco el mensaje del paper, porque, a ver, qué sé yo,

[05:23] tenemos que ser un paper más o menos humano, ¿no? Entonces mostrar que uno no entiende eso,

[05:28] me parece que está bueno y llama un poco la atención.

[05:32] Entonces lo voy a poner bastante adelante, este tema.

[05:42] Y entonces, a ver, digamos, porque la onda acá es un poco quedarse con lo esencial.

[05:48] Entonces la figura 1 sí está perfecta, o sea, porque nosotros estamos mostrando el algoritmo de IR.

[05:57] Claro, después, lo que pasa es que ahora el algoritmo de IR y el MR,

[06:05] Si al final de cuentas, vos lo que te importa es, la variable importante es la resolución de los parámetros.

[06:17] Si la resolución es más o menos la misma, alpiste.

[06:20] Pero igual creo que no es exactamente la misma.

[06:23] Porque lo que yo veo, digamos, y acá es una cosa importante, es que hay una diferencia en la matriz de distorsión.

[06:30] La matriz de distorsión es más importante en el caso del MR que en el IR.

[06:39] Pero bueno, claro, digamos...

[06:43] O sea, igual queda un poco la duda.

[06:50] Lo que queremos probar es, bueno, usemos el IR que es mejor que el MR.

[06:56] ¿Eso queda completamente probado?

[07:00] Y no del todo, porque bueno, si vos al final tenés que corregir, bueno...

[07:05] ...qué es lo que vale, ¿no?

[07:10] Es la corrección. Igual se ve que es marginalmente mejor siempre.

[07:16] Yo creo que digamos, no es que son idénticos, o sea...

[07:24] El caso, por ejemplo, para el IRT es más difícil, ¿no? Ahí ya...

[07:30] Ya, digamos, las diferencias son más pequeñas y la verdad que no es tan fuerte.

[07:38] El caso, digamos, lo más fuerte sí es quizás el gráfico de la correlación cruzada, ¿no?

[07:49] O sea, de cómo una medición depende de la subsiguiente, depende de la anterior, ¿no?

[07:57] eso ahí sí claramente se nota que es mejor el IR que el MR.

[08:06] Entonces quizás, yo creo que eso lo voy a tener que incluir,

[08:11] quizás esa podría ser la figura 2 en vez de la 3.

[08:16] O sea, después me queda el caso del MR.

[08:25] El MR, por ejemplo, es bastante bueno.

[08:28] No hay tanta historia para la cosa.

[08:32] Tenés como una caída, pero es un error bastante pequeño para este tipo de casos.

[08:44] Hay una cuestión que es clara.

[08:53] Si nos ponemos a medir cuál es el mejor de los algoritmos en todas circunstancias, creo que el mejor siempre es el IR, ¿no?

[09:00] No habría situación donde el otro algoritmo sea mejor.

[09:06] Eso estaría bueno, ves, medir, digamos, cuál es el algoritmo que es mejor.

[09:12] O sea, si se demuestra, claro, eso está bueno, ves.

[09:20] Claro, está bien. Eso me parece que es un resultado honesto.

[09:26] Ahora se me ocurrió una especie de cuadrícula donde vos decís, bueno, hay alguna circunstancia donde el IR no sea mejor que los demás, por ejemplo.

[09:39] Eso no es mala idea, ¿no?

[09:50] Después la otra es, bueno, cuando el algoritmo deja de ser confiable, ¿no?

[09:57] O sea, a ver, y ahí tenemos dos cosas, ¿no?

[10:01] O sea, en realidad, digamos, el problema es que ninguno es confiable por debajo de tau

[10:07] por este asunto de la cross-correlación temporal, ¿no?

[10:12] O sea, todos necesitan corrección.

[10:15] Eso es así y no hay vuelta, ¿no?

[10:20] pero claro

[10:25] digamos, cuál necesita más corrección

[10:29] digamos, siempre creo que el IR

[10:31] es el que menos corrección necesita

[10:33] igual estaría bien ver

[10:35] si hay otros casos

[10:37] yo creo que, digamos, ese

[10:38] por ejemplo es un buen

[10:40] un buen caso, porque bueno, vos podés decir

[10:43] bueno, cuál es mejor, entonces ya está

[10:44] siempre es mejor ese

[10:46] yo lo que necesito es un mensaje claro

[10:49] Pepe tiene que dar un mensaje claro

[10:51] en estas circunstancias

[10:54] usemos siempre IR

[10:56] por ejemplo

[10:56] pero bueno, si quieren usar IR

[10:59] úsenlo en estas circunstancias

[11:02] si quieren usar IR

[11:03] usen nuestras circunstancias

[11:04] es decir, uno podría decir

[11:07] para cada algoritmo en qué circunstancias

[11:10] se pueden usar

[11:11] entre comillas

[11:12] con qué correcciones

[11:14] y cuál es el recomendado

[11:17] y si uno podría recomendar

[11:19] siempre IR y chau, pero habría que factorizar el costo computacionales.

[11:26] Y eso no sé si es tan simple, bueno habría que, eso es algo que tengo que hacer,

[11:33] bien, sí, tomar el costo computacional de cada uno.

[11:42] si eso eso tengo que hacerlo

[11:44] y se lo voy a pedir a claude que me lo haga eso es un punto importante pues claro tenemos eso

[11:52] claramente vos tienes un costo computacional y claro digamos o sea vos vos podrías sacarlo

[11:57] con el tiempo computacional pero bueno no no lo tengo optimizado nr y mnr digamos no

[12:07] o no los tengo súper optimizados entonces es un poco dudoso podría tratar también decirle a

[12:16] cloud que me los optimice o algo así o sea claramente lo que es un nr digamos vos podés

[12:32] optimizarlo mucho porque no tenés que calcular matrices

[12:38] para calcular

[12:41] la corriente todavía una ya no es infinito más rápido que esté bien hecho

[12:49] todos los recursos y vos tienen un costo más importante pero claro digamos si

[12:56] vos al final lo que está haciendo es promediando entonces ese costo reduce

[13:02] el tema de que promediar

[13:04] vos no perdés información

[13:06] es un poco

[13:08] el mensaje fuerte

[13:10] del paper

[13:11] y claro, la pregunta acá

[13:13] más importante es

[13:15] que yo me tengo que hacer ahora

[13:17] me basta con

[13:21] digamos, para este paper

[13:24] con un sistema

[13:25] de dos estados

[13:27] los reviewers los van a comprar

[13:30] o sea, ¿por qué no hago

[13:31] por lo menos un sistema de tres estados, o cuatro estados, o qué sé yo, o los 48 estados que tengo.

[13:41] Yo, digamos, algo en mí un poco se resiste a avanzar en la complejidad de esto,

[13:50] que me quiero quedar en esta simplicidad, ¿no?

[13:53] Como dándole, digamos, ganas de que, bueno, a ver ahora qué pasa cuando vos metes dos estados, ¿no?

[14:00] O sea, como forzar una segunda temporada del paper, algo así, como que estoy guardando cartas para la segunda temporada.

[14:09] Entonces la pregunta es, ¿el reviewer aceptará eso? Bueno, me puede pedir, dame, no, dame un adelanto de la segunda temporada o no te acepto, digamos.

[14:19] O sea, yo creo que eso es un...

[14:23] Creo que es una estrategia razonable,

[14:26] dado que no estoy tan urgido de tiempo.

[14:30] Que estoy urgido de tiempo, digamos, ya lo tengo que publicar.

[14:32] Y además, bueno, sí, ya lo hago con esto solo y ya cierra.
