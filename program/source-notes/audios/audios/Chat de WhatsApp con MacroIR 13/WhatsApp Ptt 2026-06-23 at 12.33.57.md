**Archivo:** WhatsApp Ptt 2026-06-23 at 12.33.57.mp3
**Duración:** ~32:40
**Hablante:** Luciano

[00:00] Entonces, primero decir que los procesos marcovianos son muy importantes para modelar procesos biológicos de cualquier naturaleza, son bastante universales.

[00:22] El segundo punto es decir que los procesos marcovianos son, los tan independents, ¿no?

[00:30] Se supone que son continuos, instantáneos, pero que las mediciones nunca son instantáneas.

[00:41] Siempre vos tenés un tiempo de medición.

[00:50] pero que eso, digamos, generalmente no se tiene en cuenta a la hora de modelar estos procesos,

[00:59] o sea que una versión, digamos, naíf del uso de estos modelos implica un observable instantáneo,

[01:15] representante del modelo marcoviano.

[01:18] Después me especializo en canales, digo que ese observable sería la corriente.

[01:24] Y que el algoritmo que se ha presentado para modelar esto es macro R,

[01:41] Pero que ya, digamos, en el momento en que lo presenté, decía que tenía ese problema, que no servía para corrientes promediadas.

[01:57] Entonces, ese problema se solucionó con Macro IR, que ya lo usamos en un trabajo anterior,

[02:09] donde se encontró además una asimetría en la activación de los canales pulinérgicos.

[02:21] Pero bueno, en este caso lo que quiero hacer es, digamos, en ese trabajo no se puso énfasis en el algoritmo,

[02:31] sino en los resultados, y ahora quiero poner énfasis en estudiar este algoritmo.

[02:38] En cuanto a por qué es necesario hacer, digamos, por qué un algoritmo más naif,

[02:53] en el cual vos simplemente supongas que, bueno, yo tengo un intervalo, bueno,

[02:59] tomo la... lo modelo digamos la conductancia como si fuera un... lo llaman el "open channel

[03:14] noise" o sea con una corriente promedio y una varianza de esa corriente promedio y tomo

[03:20] eso, y lo tomo como si fuera instantáneo, frente a que eso no funciona, voy a mostrar,

[03:29] y que es necesario realmente hacer, digamos, tomar un modelo en el cual el Estado incluye

[03:47] el estado inicial y el estado final durante el intervalo y que haciendo ese moderado,

[03:57] o sea, haciendo el boundary state, lo llamo yo el estado boundary, se puede generar un

[04:10] algoritmo que tiene la misma complejidad algorítmica de macro R pero que usa la

[04:19] información durante el intervalo de medición para obtener una mayor precisión en el estado

[04:34] final del algoritmo y con eso tenés una mejor, una más precisa determinación de la likelihood.

[04:46] Con lo cual acá entre paréntesis tengo que mostrar que la likelihood es mejor.

[04:51] Y después lo estaba evitando eso.

[04:56] Bien, entonces este trabajo lo que va a hacer es presentar macro IR y mostrar los regímenes en los cuales funciona y dónde deja de funcionar y por qué deja de funcionar.

[05:12] o sea digamos el problema genérico de esto

[05:21] y

[05:24] claro o sea en macro r digamos ya es una

[05:31] una aproximación en cuanto a que una distribución discreta de estados

[05:41] discretos lo tomo como un continuo.

[05:46] O sea, el espacio de estados, los cuales vos podés tener distintas combinaciones

[06:01] de número de canales en distintos estados, que es un espacio discreto,

[06:06] o sea vos tenés la combinatoria de isquieta, la llevo a una combinatoria continua donde vos podés tener 2,5 canales en un estado y 0,5 en otro

[06:21] y represento el probability state, o sea el estado de probabilidades con una distribución

[06:37] continua en lugar de una distribución discreta, o sea con una gaussiana multivariada en lugar

[06:45] de una multinomial, digamos. Y entonces esa aproximación debería ser válida, bueno,

[06:55] para n mayores que ponerle 10, 30, no sé, una cosa así, habría que ver bien cuál es

[07:02] el criterio. Y después la otra simplificación está en que la likelihood para un, digamos,

[07:15] para un intervalo solo.

[07:19] También yo la aproximo con una gaussiana,

[07:23] cuando en realidad la verdadera probabilidad es bastante más compleja.

[07:27] Es una distribución gamma,

[07:38] no, como se llama una distribución de telegráfica estocástica,

[07:44] una cosa así, que vendría a ser como la proyección de una gama multinomial,

[07:50] pero que además está en un simplex y que vos la proyectás sobre un vector de conductancias.

[08:04] Es una distribución bastante compleja, no para nada normal,

[08:12] en caso de que vos estés en intervalos muy pequeños,

[08:18] en intervalos bastante más pequeños que la constante de tiempo que estás considerando.

[08:28] Entonces tenemos como dos regiones, que sería la multinomial para número de canales

[08:37] y la telegráfica estocástica para intervalos cortos.

[08:47] Y una distribución normal en todo sentido,

[08:51] normal para la distribución de canales y normal para la ECTICU en el centro

[08:57] y ahí es donde vive Macro IR.

[08:59] Entonces la idea es un poco explorar cómo es la transición entre el régimen donde la aproximación macro-R y R funciona, donde ya no funciona tan bien.

[09:14] Y entonces, ¿cómo yo defino que funciona bien o no funciona bien?

[09:20] y entonces para eso uso tres estimadores

[09:25] el primero es simplemente la macro IR

[09:30] supone un producto de distribuciones normales

[09:35] entonces cada distribución normal tiene un residuo

[09:38] entonces puedo poner los residuos

[09:40] o sea, puedo ver la esperanza del residuo

[09:47] que debería ser cero, por ejemplo, eso es una cosa

[09:49] y después es la varianza del residuo.

[09:52] El residuo está estandarizado.

[09:57] Ese es el primer estimador.

[09:59] El segundo estimador es el score,

[10:04] es decir, la derivada de la log-tyclicum,

[10:10] estimada en dos puntos, en realidad.

[10:15] En el primer punto, lo que hago es simulo un experimento con parámetros conocidos,

[10:26] simulo muchas veces y entonces tomo el promedio de la ecriculación evaluada en ese punto.

[10:35] Y eso sería el score, y ese debería ser cero.

[10:45] Ahora, para poder entender cuánto se aparta, porque la likelihood puede ser muy grande,

[10:51] muy chica, no sé, yo eso lo proyecto, digamos, usando la matriz de covarianza de los parámetros,

[11:09] obtenida, bueno, como vamos a ver después. La proyecto y me da, digamos, cómo sería

[11:15] el bias, ¿no? Se lo transforma en un bias, en un error en los parámetros, y veo, digamos,

[11:26] eso. Ese es el segundo estimado. Y el tercer estimado tiene que ver con la propiedad de

[11:35] que la esperanza del score sea cero, bueno la tercera es que la varianza del score es igual al

[11:45] negativo del Gessiano. Y entonces yo el Gessiano lo mío directamente, o sea tomando la derivada del

[11:53] score, o sea tomo el score que lo calculo analíticamente, tengo un sistema de derivadas

[12:00] automáticas en el código, entonces con una derivada numérica de... yo tengo el Gessiano y a ese lo

[12:10] comparo con la varianza del score. La varianza del score se la toma... ah y ahí está el tema. La

[12:23] varianza del score está en realidad bien, el problema es con el Gessiano, porque a veces en

[12:27] algunos algoritmos me da que puede darse

[12:33] entonces los algoritmos no son perfectos es decir si los algoritmos

[12:42] son imperfectos o sea vos no te va a dar que el score sea cero, te va a dar un bias

[12:50] y entonces lo lógico sería medir y entonces lo que puede ocurrir es que el

[12:57] gesiano te dé indeterminado

[13:01] entonces en esos casos conviene evaluar el gesiano en el óptimo en el punto

[13:07] donde la esperanza de del escorce a cero realmente no sé que es el punto

[13:14] donde vas a ir a parar, digamos, si vos usás esta aproximación de la Lanklin-Hood,

[13:22] donde vas a ir a parar para los parámetros que vos elegiste originalmente,

[13:32] esos van a estar distorsionados por una cantidad y eso se determina, digamos,

[13:35] midiendo el punto donde el score se hace cero.

[13:43] optimizas el score y además que la likelihood es máxima, pero bueno, esa es otra historia.

[13:50] O sea, el punto de máxima likelihood estimation, en ese punto yo calculo el gesiano y lo comparo

[14:00] con el mismo punto de valorado, la varianza del score. Bueno, son dos matrices y entonces

[14:12] ahora para cuantificar cómo se distorsiona la varianza en cada uno de los parámetros.

[14:21] Lo que hago es una corrección que se usa en geodesia, no sé qué,

[14:35] que es una visión sándwich que mantiene simetría.

[14:41] Entonces me queda una, digamos, calculo una matriz de dispersión de la información.

[14:48] Y esa, digamos, matriz de dispersión de la información junto con el vector de bias,

[14:55] junto con los dos escalares de, ¿cómo se llama?

[15:05] El residuo estándar y la disovarianza, esos cuatro indicadores van a indicar qué tan bueno es el algoritmo,

[15:17] o sea, qué tanto distorsiona la realidad.

[15:20] ¿Y distorsiona en qué sentido?

[15:22] Bueno, distorsiona en el sentido de que tengo que corregir el parámetro obtenido

[15:27] y por otro lado que tengo que corregir la varianza obtenida.

[15:35] esas son las dos cosas

[15:36] bueno, entonces

[15:39] hice el estudio

[15:42] el estudio lo hago

[15:44] tomé, digamos

[15:49] un monedlo minimalista

[15:52] solamente con dos estados

[15:54] para que quede bien claro

[15:57] porque después la idea es como que

[15:59] digamos, la esperanza es que

[16:02] las conclusiones se puedan

[16:04] extrapolar a más estados

[16:07] porque bueno, yo como voy a tener

[16:09] digamos la

[16:09] ¿cómo se dice?

[16:12] voy a tener por

[16:14] digamos la idea es como que

[16:19] los patrones que yo encuentre los pueda

[16:21] extrapolar a más estados, eso quiero decir

[16:23] bien

[16:24] entonces es de dos estados

[16:27] abierto y cerrado

[16:30] y tengo dos

[16:31] un K, una K, dos y una K, dos

[16:32] y entonces yo lo que estudio, digamos, este parámetro

[16:41] son básicamente tres parámetros, el número de canales

[16:44] que eso me va a permitir ver la transición entre la gaussian y la multinomial

[16:54] y el tiempo de integración

[17:01] desde un tiempo de integración igual al tiempo de relajación del sistema

[17:07] o parecido, en realidad no sino igual a uno sobre la constante cinética

[17:15] a, digamos, y hacerlos cada vez más reducidos, hasta 7, entre 17, entre dos décadas nomás estudié

[17:26] y entonces se ve que a tiempos más cortos estaríamos más en el régimen, a tiempos cortos y a poco ruido, que sería el tercer factor de ruido,

[17:46] estaríamos en el régimen donde la aproximación gaussiana a la likelihood no es buena,

[17:56] sino que tendría que ser una aproximación, digamos,

[17:59] sino que estaría en un régimen más complejo que es el ruido telegráfico o stocástico.

[18:08] Entonces, nada, y la idea es estudiar todos los algoritmos que yo presento, a ver cómo se comportan en estas regiones.

[18:25] Lo primero que se ve es que todos los algoritmos, o sea, primero que los no recursivos andan muy mal,

[18:34] digamos en cuanto a

[18:39] más que nada en cuanto

[18:41] no en cuanto al bias

[18:42] sino en cuanto a la

[18:44] inflación de la covarianza

[18:47] básicamente

[18:50] la covarianza se infla tantas veces

[18:57] como N tenés

[19:00] tenés tantas pseudo-réplicas

[19:02] sería como parecido al problema de hacer do-réplicas. Si vos tomas 100 medidas por cada tau,

[19:11] es como que te incrementa 100 veces la varianza. La varianza deducida por el algoritmo es 100 veces

[19:27] más grande de 100 veces más chicas que las gradas

[19:34] o sea en la distorsión es casi proporcional

[19:42] al número de mediciones para un intervalo

[19:48] y en macro R ya no es así, se corrige en gran medida, pero queda un residuo importante,

[19:59] o sea, se te llegan a inflar las cobrianzas en factor 2, 2 o 3, no sé.

[20:09] Mientras que en macro R lo máximo que llega es a 1.3, una cosa así, en estas circunstancias.

[20:18] O sea que es claramente el ganador.

[20:22] Y bueno, lo que se ve después es que, digamos, hasta 100 canales la cosa anda bien en cuanto al número de canales.

[20:33] Con 10 canales se distorsiona bastante.

[20:37] O sea, se nota una distorsión importante.

[20:42] lo cual es esperable porque ya no estás en el régimen donde la distribución normal sea una buena aproximación.

[20:51] Y en cuanto al río telegráfico, también lo mismo.

[20:58] O sea, vos ves que sí, que hay dos factores que lo hacen aumentar.

[21:08] O sea, por un lado tenés este factor de correlación que aumenta de la misma manera.

[21:17] Que eso tiene que ver más que nada con las, supuestamente, que un canal está correlacionado consigo mismo.

[21:29] O sea, es un canal que está abierto y estuvo abierto durante todo el intervalo de medición.

[21:35] tiene una corrección 1

[21:37] y

[21:39] y eso

[21:41] digamos infla artificialmente

[21:44] infla un poco

[21:44] si yo miro dos veces el mismo

[21:48] canal que es abierto

[21:49] tengo una

[21:50] correlación

[21:52] que no es descrita

[21:55] con este modelo

[21:57] continuo

[21:58] eso por un lado

[22:04] y por otro lado lo que tengo sí es el efecto este de que la mala especificación de la función de likelihood hace que sea muy lejana a la normal

[22:23] y eso se ve en la distorsión geométrica, que es la otra variable que yo he visto.

[22:34] esa distorsión compite, vos tenes dos cosas, a intervalos más cortos el ruido se hace

[22:50] más importante, porque cuanto más corto ese intervalo el ruido blanco se hace más

[22:58] grande pero y a la vez el ruido digamos de

[23:04] gating, la probabilidad de que los canales se cierren se hacen más chicas

[23:09] entonces a tiempos muy cortos gana el régimen gaussiano y a tiempos

[23:16] muy largos también regalma el régimen gaussiano porque ya

[23:21] muchos procesos guasonianos es equivalente a lo normal

[23:27] como la suma de muchos procesos independientes,

[23:29] si te hace el rima central del límite,

[23:33] si te hace normal.

[23:35] Y en periodos intermedios vos tenés una región

[23:39] donde estos procesos poasonianos o telegráficos

[23:47] le ganan a la normal y entonces se ve una distorsión

[23:50] en el...

[23:52] en la... o sea, se ve

[23:58] una distorsión

[24:00] de 1.5, no sé,

[24:02] una cosa de 1.2, no me acuerdo.

[24:04] O sea, es una distorsión, no es muy grande,

[24:06] pero está...

[24:08] Y bueno, y todo eso baja

[24:10] con aumentar el ruido instrumental,

[24:12] digamos, todos esos efectos se hacen

[24:14] más pequeños, pero también, digamos,

[24:16] claramente tenés menos información

[24:18] porque está tapada por el ruido instrumental, la información importante ya no aparece.

[24:30] Bien, entonces en definitiva tenemos un algoritmo macro IR que se demuestra que funciona bien,

[24:44] en este caso solamente para dos estados, pero funciona bien para dos estados,

[24:53] los otros algoritmos fallan, digamos, catastróficamente o fallan,

[25:00] de manera que lo podemos poner como un algoritmo para usarlo, para esta región gaussiana, por llamar de alguna manera.

[25:11] Para la región multinomial estaría el algoritmo microscopic recursive, micro IR,

[25:17] que presentaremos oportunamente.

[25:23] Y para el otro no hay nada.

[25:29] no sé, habría que...

[25:31] digamos, es como altamente...

[25:35] bueno, también, digamos, o sea,

[25:37] ambos algoritmos serían

[25:39] muy costosos

[25:40] y...

[25:43] no sé, o sea, no...

[25:49] posiblemente, qué sé yo, pueda ser

[25:51] uno de partículas, no sé, un filtrado

[25:53] de partículas, no sé qué se podría hacer

[25:55] para

[25:55] para

[25:57] intervalos muy breves

[26:00] algo tiene que existir

[26:02] bien

[26:04] y después finalmente

[26:06] me queda decir en la discusión

[26:08] que bueno

[26:10] que

[26:11] que bien, si bien yo lo

[26:14] probo solamente, o sea lo estoy probando

[26:16] para dos estados, pero eso

[26:18] es muy probable que se pueda

[26:20] extrapolar

[26:21] a más estados en el sentido

[26:24] que bueno

[26:24] por ejemplo

[26:26] cinéticas más

[26:29] lo importante acá es

[26:30] digamos vos

[26:32] tener una idea de

[26:33] de cómo es tu cinética

[26:36] comparada con

[26:38] el ruido

[26:39] eso te da esa idea

[26:43] el número de canales

[26:44] y el... sí

[26:45] no, cómo es tu cinética no comparada con el ruido

[26:48] cómo es tu cinética comparada con

[26:50] con el tiempo

[26:52] con el...

[26:54] perdón, no es con el ruido, con el tiempo de integración, o sea vos tenés dos, la cinética comparada con el tiempo de integración

[27:01] y después el ruido también comparado con la cinética, y el número de canales, esas son las tres cosas

[27:09] y eso digamos vos lo podés ver para, vos podés digamos llevar esta reglita que está medida para dos estados nada más

[27:21] podés, digamos, con total impunidad, como una aproximación grosera, tomarla para cualquier tipo de...

[27:32] ¿cómo se dice?

[27:37] de estado de más estados, bueno, un modelo de muchos más estados vos podés...

[27:47] decir bueno, nada, o sea, cuál es el más rápido, el más lento y ver cómo se deforma.

[27:56] O sea, acá lo que queda claro es que los que más se deformarían serían los más lentos,

[28:04] los procesos más lentos, o sea que sí,

[28:09] claramente esos podrían tener un poco de distorsión en su...

[28:18] en la...

[28:20] especialmente si vos los medís con intervalos regulares

[28:26] ahora si vos usás intervalos espaciados logarítmicamente

[28:31] entonces ahí estarías como siempre en un punto óptimo

[28:37] O sea, esa sería también como una ventaja de esa forma de trabajar.

[28:44] Esa podría ser una de las conclusiones también, una de las discusiones.

[28:48] La otra discusión es cómo se puede estimar la corrección a la evidencia.

[28:54] Después otra discusión de que todo este proceso se puede aplicar para modelos más complejos,

[29:02] lo aplicas solamente para el punto óptimo y vos optimizas el modelo cinético, obtenes

[29:11] el punto óptimo y después a partir del punto óptimo simulas y podes ver digamos

[29:25] es esta matriz de distorsión de la información en ese punto por ejemplo y ahí digamos este quizás

[29:34] corregir un poco la estimación de los parámetros o bueno ver ver cuánto te está haciendo la mierda

[29:42] y también este puedes corregir la evidencia o sea uno podría corregir ex facto las evidencias

[29:53] con esto, podría ser.

[30:01] Y bueno, este modelo es como... y la cuestión es por qué nadie lo planteó antes esto y

[30:13] yo creo que es porque en el fondo no... o sea, si bien es mejor la "like" y "judi"

[30:19] todo no hay un aumento sensible en la disminución sensible en la coherencia de los parámetros

[30:31] es bastante parecida con lo cual digamos no sería tan útil en algún sentido o sea que no

[30:43] si vos usás, qué sé yo, usás Bustra, usás otros métodos, vos podrías no usar esto.

[30:51] O sea, digamos, la necesidad de esta corrección surge más que nada de un pensamiento teórico,

[31:01] de un análisis de que, bueno, no, eso no son.

[31:06] y no surgió de que si yo aumentaba el número de canales y trataba de recoger no volvía al mismo punto

[31:22] me daban vallas, especialmente el número de canales y todas esas cosas

[31:28] se discriminaba. Entonces, nada, ahí me doy cuenta del problema teórico y encuentro la

[31:37] solución, ¿no? ¿Cierto? Y bueno, entonces, irá algo que no es mucho más costoso computacionalmente,

[31:49] o sea, se puede decir que se reemplace completamente macro R por macro IR, o sea,

[31:56] se usa este modelo y se puede usar para todo tipo de obtención de la Eclix-Hood de procesos

[32:07] marcovianos. Esto es completamente general. En principio la teoría es aplicable a todos,

[32:20] es bastante universal, en ese sentido está pensado en canales, pero en principio es bastante

[32:33] universal. Y bueno, eso sería más o menos todo el paper.
