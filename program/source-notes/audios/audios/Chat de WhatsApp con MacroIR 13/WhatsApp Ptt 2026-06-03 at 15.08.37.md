**Archivo:** WhatsApp Ptt 2026-06-03 at 15.08.37.mp3
**Duración:** ~08:13
**Hablante:** Luciano

[00:00] Bueno, muy contento estaba yo con mi plan, ya que me sentía la computadora y me aburrí.

[00:09] A ver, digamos, porque no, le pedí un resumen a Claude, una recapitulación de todo lo que había emprendido en los últimos días.

[00:24] y la verdad es que sí que necesito un poco repensar las ideas que estuve trabajando,

[00:34] que no llegaron a concretarse porque bueno, un nuevo análisis ataba o obturaba una idea y sugería una nueva y así.

[00:46] Bueno, entonces, ¿cuál es el motivo mío de esas obras?

[00:54] bueno, que como estaba diciendo, que tengo una gran variabilidad en la matriz de distorsión de la información.

[01:08] Para algunos casos, que es una cosa medio difícil de manejar.

[01:14] Bueno, eso me llevó a pensar que bueno, que tengo problemas con esa matriz

[01:19] y que una manera de solucionarlo es si esa matriz es singular es sumarle un prior.

[01:26] Entonces, después me di cuenta que el prior tiene un cierto poder de fuego,

[01:34] pero sobrepasado ese poder de fuego, que es bastante bajo para un prior no muy informativo.

[01:43] Entonces, necesito otra cosa y ahí caí en la solución.

[01:49] que es usar, digamos, medir el Gessiano,

[01:54] que no está garantizado que sea definido, positivo definido,

[02:01] sino que es lo que es, medirlo en el máximo.

[02:05] Entonces, si está lo mencionado en el máximo,

[02:06] por definición de máximo tiene que ser positivo definido,

[02:10] si no, no es máximo.

[02:14] y bueno, si lo encontraba en máximo, bueno, ya no es ninguna likelihood de nada, digamos,

[02:20] que es un máximo.

[02:21] Bueno, y eso solucional también el problema de, bueno, decir que es una likelihood.

[02:27] O sea, que eso, eso es, está bien, es fundamental hacerlo.

[02:32] Bien.

[02:33] Aquí hay dos problemas más que son.

[02:41] Bueno, primero, para llegar al máximo, entonces no usaría yo el gesiano numérico,

[02:46] sino el gesiano gaussiano, o sea, el gesiano por aproximación, la aproximación gaussiana,

[02:52] que ya está garantizado que sea semi positivo definido, pero podría ser semi, no positivo,

[02:58] o sea, podría ser singular.

[02:59] Podría ser que alguno de los parámetros no aporte un sorongo,

[03:05] y entonces no lo puedas usar para obtener el máximo.

[03:11] Ahí sí necesitas sí o sí usar el prior.

[03:14] Entonces uno podría estar tentado a usar el prior de movida, listo.

[03:20] Y entonces no tenés que andar diferenciando entre los casos que están determinados y los que no.

[03:29] Si uso el prior tengo el problema de decir, bueno, cuánto prior uso.

[03:33] también puedes usar un prior genérico de una década y "Chao" alrededor del punto verdadero

[03:41] dice "Chao Pinela" y se acabó la historia y bueno, vos lo ves como

[03:47] como disminución del prior en cuanto a la sensibilidad

[03:57] Podés igual identificar sistemas inidentificables como que no se reduce el volumen del pliego,

[04:06] en la dirección que no es identificable.

[04:09] Bueno, es un poco, digamos, sí, o sea, claro, sí, me parece que está bien, o sea,

[04:22] además bueno como estoy en un enfoque bayesiano está bien usar prior

[04:28] no sé me parece que que cerraría así

[04:34] y bueno y después me queda el tema de que no solucione porque me aparece

[04:43] variabilidad y bueno ese problema lo soluciono como ya lo dije analizando

[04:49] las réplicas individuales.

[04:51] Entonces, ¿qué es lo que está pasando?

[04:56] Que bueno, yo tengo todo un montón de análisis

[05:01] de la matriz de dispersión de la información, etcétera, etcétera,

[05:05] que están basados en la likelihood,

[05:09] los tengo que basar también ahora en la posterior likelihood,

[05:11] y además tengo que agregar toda una sección

[05:15] de la optimización de la posterior likelihood.

[05:19] no sería MAP. Me parece que acá hay dos puntos que podrían ser, que podría ser yo hago un MLE,

[05:25] que es maximum likelihood estimation, o un maximum posterior, un MAP.

[05:33] Hay razones para hacer las dos cosas, o sea, la ventaja del MAP es que bueno, sirve en el caso

[05:41] en que vos tengas titular y clic en parámetros que no estén determinados.

[05:50] El otro tema que me queda por solucionar es si yo hago, si optimizo

[06:01] hago una optimización general o hago optimizaciones para cada réplica en particular.

[06:11] Yo creo que la ventaja de hacer las optimizaciones de las distintas réplicas es que ahí puedo estimar directamente cuál es la variabilidad de los parámetros,

[06:22] digamos, sin tener que usar, digamos, deducciones estadísticas, ¿no? La miro directamente.

[06:31] y eso me parece que está bueno comparar eso con, me da como una especie de test,

[06:41] comparar eso contra la variabilidad de los parámetros óptimos

[06:51] con la matriz, la deducida matriz de covarianzas de los parámetros.

[07:03] O sea, tengo dos formas de obtener la matriz de covarianzas de los parámetros.

[07:07] O sea, yo creo que eso, digamos, por el costo de un 2x, o sea, que hago dos optimizaciones,

[07:18] una optimización de todos los samples y otra de cada una en particular.

[07:25] A partir de eso, digamos, puedo obtener, digamos, dos...

[07:34] Puedo, sí, obtener eso, una estimación directa de la matriz de covarianza

[07:41] de los parámetros deducidos.

[07:48] Yo creo que eso está bien, o sea, lo tengo que hacer.

[07:52] Entonces, básicamente, bueno, tendría que todo este audio,

[07:59] generar la documentación que tenga que generar, y bueno, y generar el código.

[08:06] Así que bueno, más o menos creo que con esto lo tendría todo hecho.
