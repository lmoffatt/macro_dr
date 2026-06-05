**Archivo:** WhatsApp Ptt 2026-04-28 at 14.27.43.mp3
**Duración:** ~04:44
**Hablante:** Luciano

[00:00] desde hace varios días estoy ahora estoy en la reserva caminando pero bueno desde varios días

[00:05] estoy implementando micro IR con la ayuda de Claude. Claude es muy cabezadura y hace lo que

[00:16] quiere y bueno yo tengo que andar corrigiéndolo para que haga el micro IR que yo quiero y bueno

[00:26] Ahora el tema que pasó es que da cualquier cosa la matriz de dispersión de información, de distorsión de información, o sea, me da, digamos, con la covarianza corregida, como si la covarianza fuera más chica con intervalos más largos, lo cual ya es demasiado, no tiene ningún sentido.

[00:49] o sea que algún error, bueno el error era que no estaba bien calculando bien la información de Fischer

[00:57] ni la varianza, no sé, o sea, hacía lo que se le cantaba, o sea, tomándolo como una instrucción normal

[01:04] y ahora digamos yo pensando cómo hacer el cálculo del GCA, no?

[01:11] Te doy cuenta que es un quilombo porque vos lo que tenés es una normal pesada por un prior

[01:18] Entonces el Gessiano de eso es más complicado, no es solamente el Gessiano anormal multiplicado por el Prior,

[01:28] porque tenés un factor del score multiplicado por la derivada primera del prior

[01:44] y también tenés un factor que sería la normal multiplicada por la derivada segunda del prior

[01:54] que eso no sé qué carajo es

[01:58] es un quilombo, tengo que sacar la derivada segunda del Gessiano del Prior me pegó un tiro en las pelotas, no sé cómo hacerlo.

[02:08] Entonces estoy pensando en directamente calcular el Gessiano como la derivada numérica del score y chao,

[02:17] por eso me cago de risa.

[02:20] digamos me cuesta un poco pero lo infinito, entonces lo podría ser.

[02:27] Pero claro, digamos el tema es...

[02:29] Bueno, pero igual la pregunta es si justamente yo digamos lo que sospecho es que mi aproximación del

[02:40] el heciano con el heciano-gaussiano, quizás no esté bien, digamos, que sea insuficiente,

[02:52] que tenga algún factor más, y bueno, por eso eso explicaría, digamos, este fenómeno

[03:00] de la covarianza plana para macro IR.

[03:08] bueno estoy con eso te voy a ver cómo lo hago claro quizás lo que podría hacer si es hacer

[03:15] lo mismo para para ambos y chao porque ahí ya no tengo más drama así bueno no es tan dramático

[03:24] sea es un factor 6 este si me reduce en 6 el número de réplicas que puedo obtener pero bueno

[03:33] más o menos igual puedo reducir otras cosas, o sea por ejemplo puedo trabajar con menos,

[03:41] o sea no llegar a 100 muestras sino un poquito menos y con eso más o menos

[03:46] llego a lo mismo, o sea en vez de hacer a partir del intervalo de tau de 0,01 empiezo con 0,05

[03:56] ponerle y ahí ya tengo un factor 5

[03:58] que lo compenso con

[04:00] con lo otro y listo

[04:03] y más o menos tendría

[04:04] unos tiempos más o menos parecidos

[04:07] podría ser que yo

[04:09] sí, un factor 10

[04:10] y ya con eso estoy seguro

[04:12] sí, podría ser algo así

[04:15] y después voy a hacer la otra parte

[04:17] con más tiempo

[04:19] pero me parece que

[04:21] que voy a ir por ese lado

[04:24] por el tema de calcular la herida

[04:26] numérica de gradiente como el hessiano y chao, ahí ya no hay ningún problema, me tiene que dar lo mismo.

[04:35] Entonces claro, por lo menos lo hago como un test case para probar un par de casos a ver si me da lo mismo.
