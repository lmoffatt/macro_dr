**Archivo:** WhatsApp Ptt 2026-04-29 at 13.10.35.mp3
**Duración:** ~13:26
**Hablante:** Luciano

[00:00] Bueno, creo que estamos llegando a la luz. Estoy yendo hacia el año de julio, ¿no?

[00:06] Para darle un poco de contexto a esto.

[00:10] Bueno, creo que sí, que estoy ya llegando a la luz.

[00:21] Fue una excelente decisión incluir a Micro IR en el análisis.

[00:30] Yo la verdad que no se podía usar, porque se puede usar hasta 100 o 200 canales sin problema,

[00:43] ya con más. Obviamente eso es porque son solo dos estados, si tienes más estados es

[00:55] un quilombo. Pero para eso, para dos estados, funciona bien. Y bueno, me alcanza para ver

[01:05] los fundamentos, digamos, cómo están. Bueno, entonces, lo primero es que micro IR efectivamente

[01:16] no tiene correlación temporal, es decir, las medidas son completamente independientes

[01:23] en el residuo estándar.

[01:27] Y también no hay distorsión en la...

[01:31] cambios en la información de...

[01:35] no hay cambios en la matriz de distorsión de la información.

[01:42] Entonces, lo cual me indica que los cambios que sí veo en macro R

[01:48] son producto de la aproximación normal del espacio de probabilidades.

[01:56] O sea, que era lo que me decía la guía, pero bueno, yo no le creía al todo.

[02:00] O sea, que digamos, tomando solamente dos momentos, la media y la covarianza

[02:07] no alcanza como para tener una buena, una feliz representación de todo,

[02:11] necesitamos más momentos y bueno, eso no es posible.

[02:17] Entonces, eso cierra un poco ese capítulo en sentido que, bueno, sí se muestra que es una aproximación.

[02:25] Entonces, la segunda cosa es que micro IR funciona con, puede hacer lo que llamo la QDT,

[02:36] que es el cálculo de la matriz exponencial y la conductancia, etc.

[02:47] el tema es que eso

[02:50] lo puedo calcular con

[02:52] autovalores, con eigenvalues

[02:54] o con una expansión

[02:56] de Taylor

[02:57] y bueno

[02:59] yo estoy usando

[03:03] los eigenvalues porque son más

[03:04] rápidos, más precisos, etc

[03:06] pero claro, ya con 100 canales

[03:09] ya los eigenvalues

[03:11] se dan la mierda

[03:12] no sirve, entonces ahí sí tengo que usar

[03:14] la expansión de Taylor y bueno

[03:16] es lo que estoy ahora compilando

[03:19] que ahora me fui a caminar porque está compilando eso

[03:21] entonces con eso podría hasta llegar

[03:23] a 50 canales

[03:27] los resultados andan bien

[03:31] hasta 10 canales, en 20 ya hay problemas

[03:33] en 50 es un desastre

[03:35] entonces tengo que

[03:37] lograr eso, que eso funcione

[03:40] bien, entonces

[03:43] eso lo haría con 100 canales

[03:47] entonces queda

[03:49] realmente solucionado

[03:50] y lo que hice

[03:51] la

[03:53] ¿cómo se llama?

[03:56] la genialidad de ayer

[03:58] fue que

[04:00] claro, que la

[04:02] matriz de la

[04:04] Fisher Information Matrix

[04:06] no la puedo estimar

[04:09] con una aproximación

[04:10] gaussiana, sino que tengo que

[04:12] usar la derivada del gradiente.

[04:16] Entonces simplemente eso, y bueno, funciona y está bien.

[04:20] Y entonces lo que se ve es que, digamos,

[04:23] la matriz de la Gaussian Fisher Information

[04:28] funciona bastante bien para macro R

[04:32] y bastante mal, o sea,

[04:35] dependiendo de cuánto para micro IR.

[04:39] O sea que, en principio, podría usar la matriz esta para corregir la función de likelihood de macro R.

[04:53] O sea que sería un poco la propuesta o el camino a seguir.

[04:58] O sea que eso funcionaría bien.

[05:01] Así que bueno, creo que ya está.

[05:06] ahora sí el paper se cierra realmente. Hay un tema que es que, bueno, micro R con

[05:12] con intervalos largos, intervalos del orden de tau, ahí tiene vallas, no sé por qué

[05:23] mierda, creo que para un canal, así que eso no sé si indica errores en GEMIN o algo así.

[05:30] Vamos a ver ahora qué pasa con la expansión de Taylor, que a ver si se lo soluciona.

[05:36] O si, quizás, justamente, necesitamos más, estamos con el mismo problema, ¿no?

[05:46] Necesitamos más momentos de la distribución de la corriente, ¿no?

[05:57] O sea, porque ahí también están los mismos.

[05:58] Ahora estamos usando solamente dos momentos, la media y la varianza,

[06:02] quizás para intervalos largos del TAU vos necesitás más momentos y todo el edificio

[06:16] de IR, digamos, de alguna manera, no te digo que colapsa, pero bueno, no es que funciona

[06:23] tan bien, en fin, pero bueno, pero está bien, eso tiene sentido.

[06:32] yo pensaba que al contrario

[06:35] que iba a tener problemas

[06:36] con intervalos más cortos

[06:39] porque

[06:39] ahí son de poisson

[06:42] no es

[06:44] normal

[06:46] pero se ve que eso no calienta tanto

[06:49] ah, tengo que ver también

[06:51] con el ruido

[06:52] ah, porque puede ser por eso

[06:55] porque claro, al ser el ruido

[06:57] más elevado

[06:58] comparado con

[07:00] con la conductancia

[07:03] al ser los intervalos más chiquititos

[07:05] de medición, entonces

[07:06] claro, por eso, al ser más caos

[07:09] yo creo que es eso, es espoison

[07:10] no es, es que digamos

[07:12] si vos ya tenés

[07:14] sí, el ruido muy pequeño

[07:17] comparado con

[07:18] el ruido, digamos

[07:20] el ruido instrumental

[07:22] es mucho más chico

[07:24] que el ruido

[07:26] de cinética

[07:30] entonces ahí, claro, no te vale la aproximación de Gaussian,

[07:35] tendría que ser una aproximación de Poisson, y bueno, eso sería otro tema, ¿no es cierto?

[07:40] Está bien, eso es importante decirlo.

[07:43] Bueno, entonces ahora queda el tema del paper,

[07:46] porque el paper ya pasaría a presentar también los microscópicos recursivos,

[07:53] y bueno, sí, va a haber que presentarlos, no queda más remedios,

[08:00] los tomo como los...

[08:03] como se dice, como la...

[08:06] lo que se llama la fuente de verdad, no sé, el gold standard serían.

[08:11] Y bueno, voy a analizar cómo funcionan los dos,

[08:20] o sea, con este tema de Poisson y Gauss

[08:25] para la... ¿cómo se llama?

[08:30] el ruido, la composición de Gauss y Poisson cuando Gauss es grande es como Gauss, entonces

[08:40] ahí no jode, pero claro si vos tenés la otra situación entonces si tenés un errorcito

[08:49] ahí, bien, que no vamos a tratar, pero bueno está bueno señalar.

[08:55] y bueno, entonces

[08:59] se presentan

[09:01] digamos estos seis o siete conceptos

[09:04] todos

[09:05] encadenados

[09:07] para mostrar

[09:09] digamos

[09:10] las limitaciones de este método

[09:14] básicamente es eso

[09:16] o entender

[09:17] en qué contexto funciona

[09:19] esa era un poco la idea

[09:20] y la verdad que yo creo que en eso estoy muy contento

[09:24] que se ve que no funciona ningún contexto, en el sentido de que siempre tenés que hacerla,

[09:29] digamos, para tener una idea verdadera de la likelihood, tenés que hacer la corrección.

[09:38] Pero la corrección es posible, digamos, no es mucho más cara, no es más cara que la optimización,

[09:48] que el cálculo de la evidencia, así que no creo que sea, no es un problema.

[09:54] eso no alula el método para nada

[09:57] así que bueno, más o menos

[09:59] ya queda cerrado

[10:02] este capítulo de Macro IR

[10:04] creo que voy a poder publicar este paper

[10:09] es muy importante, es un paper bueno

[10:11] la verdad que sí, estoy muy contento

[10:13] porque vi un problema que no había visto

[10:18] creo que más o menos lo puedo justificar

[10:20] y bueno

[10:24] ya queda

[10:26] algo sólido

[10:29] algo monolítico

[10:31] como a mi me gusta

[10:32] un paper que no

[10:34] que no quede

[10:36] que sea

[10:38] fundacional

[10:40] que tenga buenos cimientos

[10:43] esa es un poco

[10:44] mi idea y creo que está bien

[10:46] que eso lo cumple

[10:48] y bueno

[10:49] y ahora lo voy a ver si lo escribo ya esta semana, lo que me quedaría hacer es ver un poco cómo terminar esto,

[10:58] o sea, está bien, coincidan macro R y micro R en el rango de casos, y bueno, listo, puede ser que pueda conseguir,

[11:16] o sea, yo puedo correr una vez, ponerle un par de veces

[11:20] un micro IR para 1000 canales

[11:24] sería interesante ver si

[11:32] si poder ahí recuperar más información

[11:37] que el macro IR, ¿no?

[11:38] te da una idea, por ejemplo, de cuánta información se pierde

[11:43] así que bueno eso sería todo más o menos en este cuasi cierre de este paper que me ha contado bastante sangre

[11:57] pero bueno creo que es un paper que me ha dejado muy contento

[12:01] ahora sí siento que casi no tiene huecos

[12:08] no

[12:10] un trabajo sólido así que bueno veremos ahora bueno vamos a apuntar a

[12:18] a ilay veremos qué onda ahí supongo que puede andar y este y bueno eso es como sigo

[12:29] pero bueno primero tengo que terminar esto después ver cuál es el siguiente

[12:36] espacio pero pero si en principio bueno incorporar micro r realmente cierra cierra el tema mucho más

[12:46] que por ejemplo meter más estados que también sería interesante manejar y todo pero luego fíjate que

[12:53] ahí digamos y abrir más preguntas no es que cerradas acá estamos contestando preguntas abiertas

[13:02] por eso. Con tres estados no abriría, abro nuevas preguntas, pero no tengo posibles

[13:13] respuestas desde acá. Está bien, no hace falta hacerlo, me parece que eso sí. Está

[13:19] bien. Eso es una cosa para decir en el CAI, pero digamos, ¿por qué me limito a este

[13:25] análisis?
