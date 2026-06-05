**Archivo:** WhatsApp Ptt 2026-04-01 at 10.31.07.mp3
**Duración:** ~03:55
**Hablante:** Luciano

[00:00] Bueno, no estuve dando muchas vueltas a la arquitectura de Indexed y no, realmente todas las soluciones que había propuesto antes estaban de alguna manera u otra mal.

[00:15] Entonces, a ver, vamos a ver si escribo ahora la arquitectura que sí sería.

[00:20] bueno, primero

[00:23] se mantiene

[00:26] digamos la invariante de que

[00:30] los types de expression de T

[00:32] devuelven un T, entonces lo que voy a hacer es voy a especializar

[00:36] un type de expression de indexed T

[00:38] de manera que a esa type de expression

[00:42] de indexed T tenga dos métodos

[00:46] adicionales, uno es

[00:47] este my indexed my index space que vuelve el index space del indexed que está dentro

[00:58] y otra que es el run by coordinate que corre en la función por coordenadas o sea que te da digamos

[01:13] el valor de una coordenada del index t.

[01:17] Bien.

[01:18] Entonces la idea es especializar todos los types

[01:25] de expression t para el

[01:29] index t.

[01:30] Bien, eso por un lado.

[01:34] Y por el otro lado, entonces, a ver,

[01:38] donde se decide, ah, y después lo que voy a tener es

[01:42] el objeto

[01:44] TypedIndexed

[01:46] TypedIndexed

[01:49] de T

[01:50] que lo que va a hacer es tener

[01:52] contener adentro un

[01:54] puntero a un

[01:56] TypedExpression de un

[01:58] IndexedT

[01:59] entonces ahí yo tengo

[02:02] un IndexedT

[02:03] y ahí me puede dar

[02:05] ahí yo especializo el

[02:08] run

[02:08] con Environment y Coordinate

[02:11] y lo que hago es acceder al runByCoordinate del objeto que tengo adentro.

[02:20] Y también el getIndexed, va a ser el getMyIndexedSpace.

[02:30] Bien, eso es con respecto a la arquitectura de los TypedExpressions.

[02:37] Después, la decisión de cómo creo yo los Typed Indexed es en el momento en que compilo los argumentos.

[02:50] Entonces el Compile Argument lo que va a hacer es, si es T, lo manda como un Typed Expression T,

[02:59] Y si es un indexed T, lo manda como un typeed expression.

[03:06] Es typeed indexed T y manda ese objeto adentro.

[03:11] Entonces, ahí está.

[03:13] Y después lo que tengo que hacer es especializar los distintos de typeed.

[03:17] Digamos, voy a tener que tener dos funciones.

[03:21] Bueno, de typeed function, typeed function evaluation.

[03:29] va a estar el t y el indexed t.

[03:32] Entonces cuando yo hago el compile function evaluation,

[03:39] no sé cómo se llama, ahí yo o lo mando a un t o a un indexed t.

[03:46] Bueno, eso más o menos sería la forma.

[03:51] Vamos a ver si esto un poco lo puede interpretar el códex o no.
