**Archivo:** WhatsApp Ptt 2026-03-31 at 16.30.27.mp3
**Duración:** ~08:18
**Hablante:** Luciano

[00:00] Bueno, acá en macro IR ahora estoy con un problema que no es macro IR, es de, si queremos, de Luthier, que tiene que ver con una idea muy interesante que es la de hacer, de indexar, digamos, o sea, transformar una función en una especie de producto tensorial de funciones.

[00:25] es decir, vos por ejemplo tenés una función

[00:28] yo tomo una simulación

[00:29] y tengo un parámetro que puede ser los modelos

[00:31] otro parámetro que puede ser

[00:33] los parámetros

[00:37] cada vez hay parámetros para su modelo

[00:40] o por ejemplo

[00:42] los experimentos

[00:46] y entonces

[00:48] yo tengo un

[00:51] entonces en lugar de hacer

[00:53] una simulación para cada combinación

[00:55] yo le meto en la simulación

[00:57] un vector

[00:59] de modelos, de experimentos

[01:02] de parámetros y entonces

[01:03] combino todo

[01:05] y la idea de eso

[01:07] es usar

[01:08] un tipo

[01:11] diseñado especialmente para eso

[01:13] que lo llamo tipos indexados

[01:15] los cuales tienen ejes

[01:17] que indican

[01:20] ejes por los cuales pueden variar

[01:22] por ejemplo, el eje modelos, el eje parámetros, el eje experimentos.

[01:26] Y entonces vos, cuando tengas una función para argumentos indexados,

[01:37] te sale, digamos, el producto vectorial o tensorial de todos los ejes de los distintos argumentos indexados que hay.

[01:54] Entonces estoy tratando ahora de armar, digamos, la forma de que esto funcione, ¿no?

[02:02] en el código que haga todo este indexado solo, mágicamente.

[02:07] Al principio yo pensaba que simplemente poniendo un vector donde va un elemento alcanzaría,

[02:14] pero bueno, claro, ahí yo tengo que tener para cada vector,

[02:17] tengo que tener un nombre asociado al eje.

[02:20] Con lo cual, finalmente decidí que me conviene hacer una función

[02:28] que construya ejes con un identificador del eje y una lista de valores.

[02:36] Y entonces ahí se diferencia.

[02:40] Entonces los ejes te sirven tanto como eje como el elemento.

[02:47] O sea, si vos aplicas una función sobre el eje, es como que te queda el eje sobre ese elemento.

[02:53] Es decir, como que se te va para arriba la función.

[02:58] automáticamente. Y bueno, esto lo tengo bastante implementado, y me tengo implementado a nivel de

[03:14] las clases, pero me falta implementarlo a nivel del DSM, lo que sería la interfaz.

[03:23] y ahí en la interfaz en la compilación el tema es donde descubro que es un eje o un elemento y

[03:31] este y bueno el lugar el chico que tiene que tomar la decisión es compile argument y compile

[03:41] argument esté teniendo el tipo del elemento ve si lo que está compilando es corresponde al tipo

[03:50] del elemento o si corresponde a índices a ejes del elemento. Entonces ahí es donde puede

[03:58] ser uno u otro y luego dependiendo de si alguno de los argumentos es un eje o no hago

[04:09] dos tipos de Function Evaluation que sería Function Evaluation o Function Indexed Function

[04:17] evaluation o algo así. Así que esa es más o menos la forma de hacerlo. Y bueno, lo único que me quedaría

[04:29] a decidir es si le agrego un método para diferenciar. O sea, cómo diferencio si estoy compilando hacia un

[04:41] un argumento o una lista de argumentos, o sea, un eje de argumentos.

[04:48] Tengo que ver cómo...

[04:51] O sea, lo más fácil es agregar un método que sea isAxis o algo así,

[05:00] o isIndexed y listo.

[05:04] Y entonces ahí aplico una u otra y listo.

[05:11] Creo que lo tengo bastante liquidado el asunto.

[05:14] Y está bien.

[05:18] Y bueno, y ahí ya casi estaría tentado a hacer índices dependientes.

[05:26] Que lo traté de hacer y no pude.

[05:28] Porque ahora me doy cuenta que, por ejemplo, para lo que sean los modelos...

[05:36] Ah, acá por ahí puedo hacer con el algoritmo y chao.

[05:38] Sí, o con Lightly.

[05:41] y listo

[05:43] pero igual

[05:45] si likelihood igual

[05:46] si no sé

[05:49] yo creo que me convendría

[05:51] tener un

[05:52] una forma de generar argumentos

[05:55] dependientes

[05:56] en realidad la idea es

[05:59] simplemente que vos

[06:01] tengas un eje que sea un índice

[06:03] un índice puro

[06:06] y después bueno

[06:07] vos derivas de ese

[06:08] O sea, yo tendría que poder generar ejes, ¿no?

[06:19] O sea, formas de generar ejes más o menos arbitrariamente.

[06:24] Pero claro, a ver...

[06:26] Claro, no, lo que es es si el eje, si el tamaño del eje depende de otro eje.

[06:35] ahí es donde la cosa se pone peliaguda, así el size depende del eje, entonces ahí sí yo necesito la dependencia, ¿no es cierto?

[06:51] Entonces ahí es un problema y bueno, tengo el problema de cómo eso lo manejo, ¿no?

[07:02] O sea, si hago vector de vector, esa es la mejor manera de hacer vector de vector, la verdad.

[07:08] Me parece que medio como que no hay mucha historia más.

[07:18] Lo que pasa es que yo no lo quería hacer.

[07:20] La ventaja de hacerlo, como se dice, templado, que vos podés hacer vector de vector de vector, sin drama, recursivamente.

[07:39] Mientras que cómo hacerlo con cosas que sean dinámicas, bueno, tengo que hacerlo con árboles, ¿no?

[07:48] Es la manera de hacerlo, con trees.

[07:50] Trees podés hacer cualquier cosa, creo.

[07:54] Sí, la idea de los trees es fundamental, ¿no?

[08:00] Trees en computación es el concepto de datos más importante que tienen, creo.

[08:10] No sé cómo es eso, si se separó o no.

[08:16] Ahora voy a averiguar.
