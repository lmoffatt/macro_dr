**Archivo:** WhatsApp Ptt 2026-05-31 at 20.07.53.mp3
**Duración:** ~07:07
**Hablante:** Luciano

[00:00] Bueno, hoy es domingo, estoy después del domingo familiar, me puse un poquito a trabajar con la computadora

[00:10] porque, bueno, estaba bastante atascado con un tema que es que me aparecían, digamos, para la information

[00:20] distortion matrix valores que se iban a la mierda, o sea, intervalos muy grandes

[00:29] y no entendía lo que pasaba

[00:31] y bueno, lo que pasa es que

[00:32] el fin, la

[00:34] Fischer Information Matrix

[00:36] se hace

[00:38] o singular

[00:41] o incluso peor

[00:43] todavía

[00:44] no positiva

[00:48] o sea, tenés algún autovalor

[00:50] hasta negativo

[00:51] y eso puede ser porque como yo estoy

[00:54] directamente

[00:56] midiendo al GCA

[00:58] Algesia no puede ser no positivo

[01:01] O sea, puede ser mixto con autovalores positivos y negativos

[01:05] Nada lo impide

[01:07] Y eso de hecho pasa con macro R y todo eso

[01:10] Entonces, no sé, en algún momento se me prendió la lamparita

[01:15] Ah, se me prendió, ya entiendo cómo fue

[01:17] Entonces, claro, me puse a ver, bueno, cómo es que yo estoy

[01:21] Porque el objetivo de todo esto es tener una corrección

[01:27] para la evidencia.

[01:30] Y claro, en un momento,

[01:32] es alto que pasa

[01:34] cuando vos tenés

[01:36] una

[01:40] una

[01:43] una ficha de información matrix

[01:46] que sea estructuralmente

[01:48] digamos

[01:50] no se dice

[01:52] que no

[01:55] no invertible, que tiene algún

[01:59] autovalor que sea cero

[02:00] o negativo incluso

[02:02] entonces bueno

[02:04] el determinante de eso

[02:06] te da cero

[02:08] la corrección es infinita

[02:10] pero claro, en realidad

[02:12] si vas a hacer la cuenta bien

[02:14] tenés que incluir el prior

[02:15] yo estaba obviando el tema

[02:18] valleseano de tener el posterior

[02:21] likelihood

[02:21] claro, si vos incluís el prior

[02:24] entonces ya está, no hay problema

[02:27] porque vos sumás el prior

[02:29] te va a dar un autovalor

[02:31] positivo

[02:32] que va a, digamos, a lo que pasa

[02:35] es que

[02:36] el likelihood no te aporta

[02:39] nada, el prior te queda igual

[02:41] para esa variable que no

[02:43] ganás ninguna información y bueno

[02:45] está todo bien, no es

[02:47] que se te va a infinito

[02:49] la evidencia, no, vos tenés

[02:51] un prior

[02:52] que justamente delimita eso

[02:54] entonces ahí me cayó la lamparita

[02:56] me he prendido la lamparita que todo el análisis

[02:58] que tengo que hacer

[03:01] lo tengo que hacer sobre los posterios

[03:03] o sea, lo que tengo que

[03:05] definir una posterior

[03:06] information distortion

[03:09] posterior

[03:10] todo, digamos

[03:12] o sea, tengo dos versiones, una que es

[03:14] like click y otra que es posterior

[03:16] lo que llamaba yo information distortion

[03:19] matrix se va a llamar

[03:20] likelihood information distortion

[03:24] y posterior information distortion, etc.

[03:28] Y cuando yo tengo,

[03:33] entonces el análisis verdadero sería sobre la posterior,

[03:38] pero puedo también hacer un análisis sobre la likelihood,

[03:41] pero solo lo hago cuando la matriz no es singular,

[03:48] cuando la matriz de la Fisherman's Formation Matrix no es singular.

[03:52] Cuando es singular, listo, no hago ningún análisis y chao.

[03:56] Y bueno, y ahí, digamos, Claude me convenció de hacer, digamos,

[04:02] de separar en distintas categorías, que es que, bueno,

[04:07] vos podés ser que todos tus bootstrap sean positivos,

[04:12] sean positivos definidos, o que algunos no lo sean,

[04:18] que ninguno lo sea

[04:20] que vos tengas una

[04:22] una matriz que sea safe

[04:24] una que sea aceptable y otras

[04:26] categorías que puso

[04:28] así que bueno, la cuestión es que

[04:30] la teoría se me duplicó

[04:32] porque ahora tengo

[04:33] dos versiones, una de la

[04:35] y otra de los posteriores

[04:38] y ya me queda algo

[04:40] que anda

[04:41] y además defino

[04:43] una nueva

[04:45] categoría

[04:47] que es la de information gain o algo así

[04:51] que es cuánto se reduce el prior con los datos

[04:57] y este sí depende del número de réplicas

[05:00] o sea, ese no, hasta ahora estaba haciendo cosas que

[05:03] digamos, este

[05:05] que el número de réplicas lo único que aumentaba

[05:09] era la precisión de lo que yo estaba calculando

[05:12] en este caso no, me cambié el valor

[05:14] pero bueno, digamos, eso está bien

[05:17] y lo puedo presentar como un resultado

[05:20] o sea que voy a tener que hacer

[05:21] para distintos números de réplicas

[05:26] así que bueno, eso es todo un avance

[05:29] yo estaba bastante atascado con esto

[05:32] porque con las matrices singulares

[05:35] lo que hacía era proyectaba

[05:37] sobre los autovalores que existían

[05:42] y todas esas cosas son

[05:44] digamos medio como

[05:45] digamos este

[05:47] para la que te criaste

[05:49] no es claro cuál es la respuesta correcta

[05:52] en este caso directamente obvio

[05:53] todo lo que sea proyectar

[05:55] no proyecto nada, simplemente si la matriz

[05:57] es singular, no se hace

[06:00] el análisis láctico y se pasa directamente

[06:02] al análisis del posterior

[06:03] que el posterior siempre existe

[06:05] y todo, digamos, yo todas las interpretaciones

[06:08] las puedo sacar porque puedo ver por ejemplo

[06:09] si no hay ninguna información acerca de

[06:11] el número de análisis y el ruido

[06:13] y la conductancia, bueno, simplemente

[06:16] va a quedar que no hay reducción

[06:19] en el prior de esas variables

[06:21] está todo bárbaro, no hay ningún problema

[06:24] así que es bien, creo que con esto ya

[06:28] liberé el obstáculo

[06:31] que me estaba atascando con el paper

[06:33] y creo que ahora va a poder andar

[06:35] digamos como por un tubo, vamos a ver

[06:37] ahora estoy haciendo que Claude me escriba

[06:41] todas estas cosas, o sea, iba a ser un

[06:43] una reescritura de código importante

[06:45] y todo lo que hice hasta ahora

[06:47] no me sirve, se lo tengo que rehacer

[06:49] pero bueno, ya todo está en marcha

[06:51] y probablemente no necesite 16.000

[06:54] samples que por ahí con 4.000

[06:55] me alcanza, no sé, vamos a ver

[06:57] qué pasa con eso, pero

[06:58] la cuestión es que estoy contento

[07:01] que creo que esto

[07:02] ya lo logré, digamos,

[07:04] desatascar.
