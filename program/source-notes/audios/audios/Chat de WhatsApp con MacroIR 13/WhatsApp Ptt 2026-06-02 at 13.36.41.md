**Archivo:** WhatsApp Ptt 2026-06-02 at 13.36.41.mp3
**Duración:** ~07:56
**Hablante:** Luciano

[00:00] Bueno, estoy un poco en crisis con Macro IR porque me puse a tratar de ver por qué tenía una gran variabilidad para Macro IR,

[00:16] con muchas samples, 16.000 réplicas, pero con 10.000 canales y 0.01 de intervalo,

[00:29] que me daba una gran variabilidad en la Information Dispersion Matrix.

[00:40] Y bueno, eso un poco me sacó de eje porque pasé de tener un poco de dispersión a mucha dispersión, entre 4.000 y 6.000 y 16.000 samples.

[00:52] No entiendo bien por qué.

[00:54] La guía me decía que puede ser una cola larga de la distribución de cocientes.

[01:09] porque es una distribución de cociente de varianza, ¿no?, en la Information Distortion Matrix.

[01:17] Y bueno, te toca o no te toca, y bueno, eso alteró, y realmente es como que

[01:29] alguna de estas matrices de la Fischer Information Matrix

[01:39] es como que tenía una gran dispersión.

[01:45] Entonces la pregunta es si hay una falla catastrófica

[01:50] o si es una tendencia.

[01:58] Y bueno, y entonces me puse a ver que bueno que por ahí, sumándole el prior, yo podía lograr que una matriz que por ahí era singular se vuelva determinada y chau.

[02:13] pero bueno el problema es si yo tengo en el promedio la matriz esta de Fischer

[02:20] tiene autodolores negativos lo cual es posible para un Gessiano no para una

[02:28] matriz de varianzas, porque las varianzas todas son positivas, pero sí para una

[02:34] matriz de covarianzas, quiero decir, las varianzas en la diagonal es positiva

[02:40] pero sí puede ser para una matriz, un hessiano puede ser tranquilamente indeterminado, si está mal especificado el modelo.

[02:51] Entonces, nada, ahí me quedé como atrancado con un tema que no, como un auto que no puede salir, viste, sin agujero.

[03:00] Y entonces, bueno, esta mañana me puse a pensar que podía ser que estas matrices indeterminadas,

[03:08] que especialmente son prevalentes en macro R, no IR, cuando tenés un tau igual a 1,

[03:16] interval tau igual a 1,

[03:24] Puede ser que sea producto del error con un factor de corrección que lo llamo el trust coefficient o algo así,

[03:38] que lo que hace es evitar que las prioridades sean negativas.

[03:43] Y estuve viendo este tema y bueno, entre otras cosas, cuando me pongo a ver

[03:53] cuánto tiene cuánto vale ese valor ese es el trasko efficient para 10.000

[04:00] canales veo que para 10.000 canales jamás pisa jamás este entra digamos se

[04:06] corrige ese factor o sea que es es prácticamente uno es un poquito menos

[04:10] que uno por diseño pero es básicamente uno mientras que bueno para

[04:14] número de canales menores y hay un efecto que aumentar y vamos a aumentar

[04:20] la diversión de intervalos. O sea que puede ser, o sea, todavía queda posible que para

[04:26] que la matriz está negativa de macro R, sí sea por eso, pero no la de macro IR. La de

[04:34] macro IR no sé qué pasa. Hay algo raro ahí. Y yo no sé, yo tengo miedo que sea, no sé,

[04:44] algún error de paralelismo que escriban una thread sobre la otra o no sé qué, uno de esos

[04:53] errores raros que es difícil de saber. Entonces estoy un poco ahí medio como que no sé bien

[05:05] qué hacer, quiero como tener un poco más de confianza en mis datos, entenderlos y creo que

[05:11] esto de presentar los bootstraps directamente no sirve, que tengo que ver los datos individuales

[05:18] y quizá trabajarlos en R directamente, porque puede ser que el bootstrap esté metiendo quilombos.

[05:26] Así que lo que voy a hacer probablemente es ahora salvar los datos de hacer menos

[05:33] o sea, menos corridas haré yo 1000 o 500 nomás y verlas esas directamente en R, a ver si puedo.

[05:44] El tema es que claro, si con eso tenés mucho error en cuanto a la varianza y todas esas cosas,

[05:51] pero bueno, no sé. Supongo que haré algo así, vamos para ver si detecto, o sea,

[05:59] En realidad las corridas individuales serían para detectar outliers.

[06:03] O sea, puntos, sí.

[06:05] Las corridas individuales son para detectar outliers y ver si correlacionan con algo.

[06:10] O sea, sí, exactamente.

[06:13] Outliers es la Fisher más que nada, que es la que es más sospechosa en estos casos.

[06:24] Porque me parece que sí, la Fisher.

[06:25] Y además que tengo que ver los autovalores, los tengo que mirar en todos casos

[06:32] y posiblemente tendría que tener también los autovectores salvados.

[06:36] No solamente los autovalores, sino los autovectores también.

[06:41] De las matrices importantes.

[06:42] Las matrices importantes son claramente la numerical fischer y después las otras.

[06:51] la Gaussian Fisher, la covarianza de la Lock likelihood, y después tengo la misma, la covarianza de la Lock likelihood,

[07:06] pero de la suma y la suma de la covarianza, las covarianzas individuales suma.

[07:18] Después la que tengo, claro, son las, como llaman, lo que llaman los lags.

[07:26] En sonidos, mire, la verdad, que lleva mucho tiempo de código, pero nino, mire, habría que verlo.

[07:36] Pero bueno, pero en principio sí, creo que querría hacer es analizar un poco, ver si hay variabilidad,

[07:44] tiene una variabilidad, digamos, de tipo outliers y los datos.

[07:48] Eso me parece que sería la manera determinada que se haya incluido.

[07:52] Esa es una cosa rara.
