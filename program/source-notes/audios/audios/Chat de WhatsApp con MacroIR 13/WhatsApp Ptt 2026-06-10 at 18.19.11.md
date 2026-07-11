**Archivo:** WhatsApp Ptt 2026-06-10 at 18.19.11.mp3
**Duración:** ~13:36
**Hablante:** Luciano

[00:00] Sí, no sé, a esta altura ya estoy un poco con el cerebro hecho puré con esto.

[00:09] Pero bueno, vamos a ver, digamos, cuáles son las preguntas que yo me había hecho,

[00:16] qué es lo que puedo responder, qué es lo que no.

[00:20] Bueno, uno de mis objetivos era mostrar que Macro IR no funciona,

[00:26] o sea que es diferente, que falla por algún lado.

[00:30] Y eso lo conseguí, por lo menos hasta, creo que hasta para mil canales,

[00:36] a ciertas circunstancias.

[00:38] Y eso es importante porque, digamos, si vos mostrás cuánto falla,

[00:45] sabés cuánto no falla.

[00:46] Es decir, podés cuantificar la magnitud del error.

[00:53] Bien, eso es un punto, ¿no?

[01:00] Claro, es el menos malo de los modelos, sería algo así.

[01:07] Entonces, sí, yo creo que la salida tiene que ser algo que...

[01:11] Yo no puedo hacer la matriz de dispersión de la información porque eso es medio teórico y todo.

[01:18] yo pongo los parámetros recuperados y la matriz de correlación de los parámetros

[01:31] recuperados comparado con lo que se mide, ahí me parece que eso tiene un poco de sentido,

[01:41] queda queda eso bien bien verificado

[01:47] bien eso por un lado

[01:54] quizás si te haga hacer un gráfico tipo

[02:03] si no son tal son 6 o por ahí vamos a evitar algunos parámetros

[02:11] que se yo pero pero bueno eso eso es una cuestión

[02:17] que ya lo tengo que hacer y si hacer con a medida que aumenta el número de datos

[02:23] también

[02:26] o sea que el Gauss-Newton es parte fundamental del tema

[02:36] y bueno y me quedan como indicadores tempranos y por ejemplo el R estándar y el R cuadrado estándar

[02:47] como indicadores tempranos me parecen muy piolas

[02:56] y bueno y después bueno tengo

[03:01] ya indicadores más sofisticados como bueno

[03:05] descomponer bueno qué es lo que anda mal, que bueno ahí ya se ve que son las dos cosas

[03:11] que es tanto el...

[03:15] o sea que tengo cambios en el Gaussian... no, el Gaussian

[03:20] Information distortion, no sé cómo se llamaba,

[03:23] Gaussian Fischer Distortion,

[03:26] me parece que anda básicamente bien,

[03:31] salvo para algunas cosas,

[03:33] para, digamos, pocos canales.

[03:36] Después el...

[03:47] Los otros variables, bueno me queda el sample distortion y el correlation, no, el sample

[03:58] distortion, si el correlation, distortion, lo llamaba, cross correlation distortion,

[04:03] eso los tengo que ver de a pares porque puede ser que tengan signos diferentes y eso un poco

[04:16] jodan en el fondo porque en realidad los dos tienen que ser uno para que la lectura sea buena

[04:26] eso también es un tema a tratar

[04:37] y el correlation distortion es bastante caro de simular, es extremadamente caro

[04:45] el otro... a ver...

[04:52] no es tan caro, no es tan caro, no, el que es caro es el que te mide la constante de tiempo

[05:01] ese sí, a ese lo tengo que volver a ver también

[05:07] Bien, porque tengo que ver por qué falla, ¿no?

[05:12] Sí.

[05:12] Eh, sí.

[05:14] Eh...

[05:16] O sea, el petro va a ser medio largo en este sentido, porque, bueno, tengo que analizar

[05:24] todas esas cuestiones.

[05:25] Eh...

[05:27] Sí.

[05:31] claro, si, si me da que hay correlaciones

[05:37] si cuando vos tenés poco ruido

[05:40] si tengo más ruido probablemente esa correlación

[05:43] desaparezca, aunque puede ser que esté

[05:45] pero claro, posiblemente no pese tanto

[05:49] posiblemente sea esa correlación la que

[05:56] haga que no

[05:59] bueno, tiene sentido esa correlación

[06:02] ahora que pienso porque justamente

[06:04] es como

[06:05] claro que no hay

[06:10] no hay realmente

[06:12] información extra, ¿no?

[06:13] porque si vos estás midiendo el mismo canal una y otra vez

[06:16] que está cerrado con él

[06:19] bueno, ahí la correlación es total

[06:21] eso es interesante

[06:26] porque creo que eso igual, no sé

[06:27] parece que no se daba o si se daba con la micro IR, no sé, esa es una pregunta, pero bueno.

[06:36] Pero en definitiva estoy como cerca, yo creo que sí que necesito un par de días para explorar todos

[06:44] los datos y todo, o sea que lo que tengo que hacer básicamente fundamentalmente es, bueno,

[06:50] volver a correr todo esto con más más amplios

[06:56] posiblemente con más ruido también con los otros otros más ruidos de volver

[07:03] cuánto tardaron en correr estas cosas pero pero yo creo que sí con

[07:10] 4.000 8.000 más o menos vamos a estar bien

[07:15] no había tanto error. Me parece que no, me equivoco. Lo que tengo que hacer ahora es implementar el Gaussian

[07:26] Newton. Lo tengo que implementar del todo, cosa de analizar eso.

[07:44] Eso me parece que es fundamental.

[07:46] Y, claro, ahora yo en realidad estoy guardando las samples.

[07:54] O sea que podría usar las mismas samples para todos los casos, para distintos casos.

[08:02] Incluso, la verdad que no estaría mal usar las mismas samples y promediarlas a distintos intervalos de dimensión también.

[08:14] Entonces, sí, ahí sería todo más o menos igual.

[08:20] Es más, hasta podrían ser todos los mismos canales.

[08:24] Entonces, es la misma información que vos la vas subdividiendo.

[08:32] si no es mala idea

[08:36] el tema es que estamos yendo un poco al carajo con todo

[08:43] pero bueno

[08:45] pero bueno ya está

[08:49] ahora lo que sí te voy a hacer es eso

[08:51] implementar el gaussian

[08:52] y bueno

[08:54] y ver un poco las cosas que tengo que guardar

[08:57] porque eso tampoco las puedo recuperar

[08:59] y bueno

[09:01] listo

[09:02] Ah, y la evolución, porque ya no estaba mostrando cómo era la evolución, y eso sí es medio

[09:09] un bardo porque son muchos datos, es pesado pesado, pero no sé, ahí no sé, es interesante

[09:22] interesante porque bueno, vos podías tener... si, que creo que era por ejemplo... si vos

[09:30] no tenés agonista, obviamente no tenés nada de la cinética, si no hay agonista

[09:39] o digamos si no hay agonista y no hay... está cerrado, si está hay agonista bueno vas a tener el ON

[09:46] pero si una vez que lo sacaste el ON no tenés más información inmediatamente

[09:52] Eso es interesante.

[09:56] Claro, y ahí posiblemente, claro, vos no tenés...

[10:04] Si, ahí no...

[10:06] Si tenés información de la no podés tenerla nunca, ¿no?

[10:09] La verdad que no sé hasta qué punto eso tiene sentido analizar.

[10:21] Intratrace

[10:23] Casi salía demasiado largo

[10:27] Puedo mostrar así más

[10:31] Como hice otras veces

[10:34] Era un poco para ver dónde estaba

[10:37] La pregunta era de dónde salía la

[10:41] La extorsión, pero ya lo tengo

[10:44] Tengo respondido que es de la

[10:47] la correlación es bastante

[10:52] tengo eso también

[10:55] lo que llaman el Tauese también

[10:59] lo cual creo que estaríamos más o menos bien con todo

[11:07] creo que voy a tener que tomar un par de días para pensar en otra cosa

[11:13] también puedo hacer

[11:15] una serpiente o lo que sea

[11:17] no sé

[11:18] igual yo lo que

[11:21] si no sé

[11:23] algo así voy a tener que hacer

[11:25] estoy un poco saturado

[11:26] ya no doy más

[11:31] ya un poco me aburrí

[11:42] todo pero bueno

[11:42] Ya estamos.

[11:47] Creo que ya está, cierro el capítulo.

[11:51] Bien.

[11:52] Increíble como si son la mayor cantidad de tiempo se perdió en encontrar pequeños errores

[12:11] que muchos los...

[12:13] Ah, bueno, no siempre le echa la culpa a otros,

[12:16] pero sí, digamos, un poco la IA tuvo que ver,

[12:19] en el sentido que aparecieron polizones

[12:22] que no sabía lo que eran,

[12:24] pero bueno, uno también comete errores.

[12:26] No sé, creo que no es.

[12:28] Más de todo justo echarle la culpa

[12:31] toda a la IA, porque eso...

[12:34] Bueno, no sé, quizás sí tenga que volver a pensar en Luthier,

[12:41] eso un poco me va a dar un poco de energía

[12:44] y bueno, tengo que ver lo de la solnaturalización

[12:48] y especialmente ya si se cae el atabón

[12:51] como parece que se caiga, digo, conseguir un laburo

[12:54] o sea, estoy en el horno

[12:57] ya estoy en menos 1200 dólares

[13:02] la concha de su madre

[13:03] no sé cómo voy a ser

[13:05] en fin

[13:07] adiós

[13:11] Ya de empezar a aplicar todo tipo de cuentas ahí en Amazon.

[13:17] Bueno, vamos a ver.

[13:19] Sí, no sé.

[13:22] La verdad tengo que sacar este paper y después dedicarme a Amazon.

[13:26] Sacarlo cuanto antes y después buscar la bura.

[13:30] No me quedo otra.
