**Archivo:** WhatsApp Ptt 2026-06-22 at 11.44.24.mp3
**Duración:** ~09:22
**Hablante:** Luciano

[00:01] Bueno, estoy de vuelta con macro IR atascado.

[00:05] Ya no sé, casi revivo un tema más para no concluir, que es ver en estado estacionario.

[00:16] Porque, claro, ¿qué es lo que vi?

[00:19] Lo que vi es que no hay diferencia casi entre macro IR y macro IR en cuanto a la covarianza de los parámetros.

[00:29] no se ganaría en disminución de error de los parámetros

[00:34] solamente en no tener un bias

[00:41] y tener una mejor estimación del error de los parámetros

[00:46] en forma directa

[00:49] con menor error

[00:52] esa sería la gran ventaja del macro IR

[00:57] pero no aparentemente no se demuestra una ventaja en cuanto al error de los parámetros

[01:06] pero podría ser que se llegue más rápido al mínimo

[01:10] porque al tener menos vallas, partiendo del valor de parámetros que usas en la simulación

[01:21] llegas al valor óptimo en menos pasos

[01:24] si eso implica también que tenés una mejor estructura de la matriz de covarianza

[01:32] del gesiano y todo como para llegar al mínimo más rápido no lo sé

[01:37] es posible o no, no lo sé

[01:39] quizás sí, quizás no

[01:42] en realidad, sí

[01:45] uno querría suponer que podría ser que sea mejor

[01:48] pero bueno, también uno podría haber supuesto que

[01:52] y vamos a tener menor error y no tenemos menor error en la cobrianza de los parámetros.

[01:58] Bueno, entonces eso me llevó a pensar, bueno, pero tengo que mostrar que es mejor en algún sentido,

[02:04] y bueno, digo, por ahí, si lo pongo en un régimen donde realmente, que sea puramente estocástico,

[02:13] que sería el régimen estacionario, donde vos no tenés más cambios de agonistas,

[02:23] solamente tenés el ruido estocástico y la recuperación del ruido estocástico en sí misma,

[02:28] entonces puede ser que ahí sí alguna variante de macro IR brille más que las otras.

[02:39] Y además, bueno, está el tema de caracterizar bien en esos regímenes

[02:46] donde no tenés otra herramienta más, qué es lo que podés obtener.

[02:53] O sea, yo creo que sí, lo de estacionario es importante porque hace la caracterización de macro IR

[03:02] y dejarla para otro paper

[03:05] digamos, es un poco

[03:07] puede ser un poco salame slicing

[03:09] además yo ya lo hice

[03:12] en el otro trabajo

[03:14] o sea, me la podrían pedir

[03:16] los reviewers

[03:17] entonces lo que podría hacer es enviarla

[03:20] y bueno, después

[03:21] pedir

[03:22] ponerla si los reviewers me la piden

[03:25] podría ser

[03:26] la verdad que no es mucho trabajo

[03:29] hacerla, o sea, serían unos días más

[03:32] de corrida y yo puedo seguir escribiendo el paper en entre mientras. O sea que eso no

[03:37] es un problema, o sea, lo voy a hacer. Igual ahora estoy corriendo con un poquito más de

[03:42] ruido, ¿no? Bueno, en definitiva, ¿qué es lo que estoy pensando respecto de la estructura

[03:48] del paper? O sea, yo tenía una idea que era ver en qué regímenes era mejor un algoritmo

[03:54] que el otro y lo que estoy encontrando es que todos los otros algoritmos son malos en

[04:00] todos los regímenes, o sea, no encuentro un solo régimen donde sean buenos los otros

[04:06] algoritmos, vos tenés distorsión por correlación para, digamos, taus muy cortos, en macro R,

[04:20] y después para taus largos también te da muy mal la parte de, también, sí, te da

[04:29] homogénicamente mal

[04:30] o sea que macro R

[04:33] no funcionaría nunca

[04:34] macro MR no soluciona

[04:37] ninguno de los problemas de macro R

[04:39] y es más, hasta los empeoran

[04:41] es una cosa loca

[04:42] con un algoritmo que no

[04:44] sirve, o sea

[04:46] claramente macro IR sí soluciona

[04:49] parcialmente

[04:51] no del todo, pero sí

[04:52] lo mejora bastante, mejora mucho

[04:54] eso no cabe duda

[04:56] pasa de una distorsión

[04:58] de 2 a 1,2 o 1,4

[05:01] con el EI

[05:02] y hay distorsiones inexistentes

[05:05] para ciertos casos

[05:07] ah, tendría que haber el bias

[05:12] a ver si es mejor el bias

[05:13] aunque sea el macroNR

[05:15] macroNR

[05:16] claro, en cuanto a

[05:19] a distorsión de la información

[05:21] ya son valores estratosféricos

[05:24] o sea, mientras el 1 era valor 2

[05:26] factor 2

[05:27] acá estamos en factores cientos

[05:30] o sea

[05:31] no es ni comparable

[05:33] tendría que ponerlo en escala logarítmica

[05:36] porque sea comparable

[05:36] ¿qué más?

[05:41] lo cual me hace pensar

[05:45] si vale la pena ponerlo o no

[05:47] porque también el otro problema con macro NR

[05:49] es que

[05:50] muchas veces ni siquiera

[05:52] llega al

[05:55] al máximo

[05:56] o sea, el algoritmo de Newton no siempre llega, entonces es medio ahorita.

[06:04] Entonces, bueno, en definitiva la pregunta es, yo tengo dos maneras de presentar el paper,

[06:09] una es simplemente caracterizar macro IR y decir que los otros son un desastre

[06:14] y no ahondar demasiado en ello.

[06:20] Claro, podría poner lo que tengo acá, lo que obtuve en material suplementario,

[06:27] y en el principal pongo macro IR.

[06:31] Y lo que me quedaría como sustancia del paper es justificar por qué falla en las condiciones en que falla,

[06:44] o sea, en qué condiciones no falla nada.

[06:48] lo que llega a ver es que siempre falla un poco o sea que sería lo esperable

[06:53] o sea vos y acá está un poco está un poco digamos la paradoja del planteo que me hice

[07:03] por el planteo que me hice es nada tengo que encontrar donde falla y después de

[07:09] y falla y falla entonces no sirve no si sirve que justamente la idea es

[07:18] claro encontrar en qué régimen se haya y tener alguna explicación de por qué fallaría

[07:24] yo creo que esa sería un poco la idea entonces las fallas son bueno por un lado la normalidad

[07:36] o sea, porque yo tengo, estoy aproximando una multinomial con una normal,

[07:45] o sea, tengo una aproximación que falla por ahí,

[07:48] y por el otro lado tengo la nueva aproximación que es,

[07:51] que si bien yo estimo, claro, también, digamos, estoy aproximando, perdón,

[07:59] una poasón con una normal, ¿no?

[08:01] porque yo tengo una especie de convulación, una poasón medio extraña de cuánto, digamos,

[08:10] porque yo tengo, sí, la distribución de la conductancia tiene una probabilidad finita,

[08:25] grande de que la conductancia sea constante y que tenga una varianza cero y yo estoy considerando

[08:34] una distribución normal o sea cuando tendría que ser una distribución de ésta de pozo no

[08:40] ya estoy aproximando una de pozo con una normal entonces ahí es la otra el otro fallo no sería

[08:50] una multinomial con una poasón y una multi poasón

[08:58] serían las dos aproximaciones que fallan, entonces tendría que ver

[09:04] esas dos como fallan y esos límites y a ver si eso explica la falla, yo creo que más o menos eso sería un poco la

[09:19] La centralidad del paper, ¿no?
