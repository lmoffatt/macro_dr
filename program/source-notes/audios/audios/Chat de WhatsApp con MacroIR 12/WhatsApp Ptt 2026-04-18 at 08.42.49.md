**Archivo:** WhatsApp Ptt 2026-04-18 at 08.42.49.mp3
**Duración:** ~17:11
**Hablante:** Luciano

[00:00] A ver, repaso lo que tengo que hacer.

[00:02] Bueno, repaso uno, en realidad tendría que agarrar con lo que yo ya tengo,

[00:12] y observar los determinantes, todas esas cosas, para ver si funcionan.

[00:17] O sea, ver todos los parámetros que me di, que todos aporten algo,

[00:26] o que por lo menos sean funcionales,

[00:29] y ver digamos cuál es reducir, cosa de poder hacer un análisis mucho más global, más

[00:41] amplio de condiciones. Ese es el punto 1 y tratar de ver, por ahí podría hasta hacerlo

[00:54] de cero la figura 2, ya van tantas veces que Dios me libre, pero sí, cosa que, bueno,

[01:05] en principio la figura 2 lo que tiene que tener es, como ya he dicho yo, el bias, la

[01:16] distorsión y la sensibilidad. Esas tres cosas son parámetros específicos, creo, sí, todas.

[01:33] El Baías es pictorial, nosotros son matriciales, pero todos puedo dar una versión escalar,

[01:45] o sea puedo tomar el determinante y esto si es importante tengo que tomar lo que sería

[01:54] la norma del bias, cosa de tener un único parámetro.

[01:58] Bien, o sea que ya tengo tres números, tres números y tengo el cubo este entre condiciones

[02:11] digamos ahí yo podría más o menos mostrar todo

[02:19] muy bien, eso está

[02:20] digamos la...

[02:23] después me queda la evolución

[02:25] de la evolución

[02:27] me quedan

[02:30] digamos serían cuatro

[02:34] a ver

[02:35] serían cinco etapas

[02:39] si se querés, sí, que serían, digamos, la inicial, que no tiene que darte nada, la del pulso y la del pospulso.

[02:50] Ahí, digamos, tenés parámetros que están presentes siempre

[03:03] y parámetros que están presentes a veces. El que está presente siempre es el ruido

[03:14] y la baseline, después todos los otros empiezan desde el ON y el ON solamente está en el

[03:23] durante el pulso y en el post pulso solo me quedan los otros, el off y etcétera.

[03:32] Quizás puede ser interesante

[03:49] mostrar las matrices no es fácil

[03:53] y

[03:56] confuso o sea yo ya

[04:02] si es confuso pero

[04:19] claro

[04:21] o sea

[04:29] la verdad que la evolución

[04:32] no sé si vale la pena algo

[04:34] o sea, sí queda lindo

[04:36] en la evolución

[04:37] el

[04:38] el lag de cross-correlación

[04:43] porque ahí te permite

[04:44] ver bien cómo decaen

[04:46] unos y otros

[04:48] y ese gráfico si es bastante diagnóstico y permite comparar los algoritmos de una manera

[05:01] muy clara, o sea en ese sentido está bien, yo me parece que es muy visual, o sea me parece

[05:13] que va para el paper, no sé. Y así como lo planteé, me parecería que estaría bien.

[05:23] Bien, o sea que entonces me quedarían las relaciones de, a ver, vamos pensando para

[05:42] el paper que es lo que va a la bahia, bueno esos tres, mis tres guerreros y la cross correlación

[05:53] temporal un poco para explicar el mecanismo.

[06:09] y basta y ahí

[06:11] y ahí yo me quedaría me parece porque...

[06:19] ah!

[06:20] si bueno puedo poner el coso corregido y no corregido

[06:23] no solo la covarianza

[06:26] no corregida y corregida obvio

[06:29] si, el terminante, son todas cosas

[06:37] Está bien, o sea que está bien.

[06:40] Ahora, el otro es la figura 3, es decir, que la cross correlación temporal sería una figura más.

[06:48] Y esa figura hasta podría ser, no sé, me parece que estaría bien.

[07:07] y a ver no sé

[07:12] qué más

[07:17] me parece que ya está y ya puedo dejar la figura 2 nomás, figura 2 está terminada

[07:33] figura 3 sería la figura 2 no sé si sería eso pero bueno para los datos para

[07:40] generar la figura 2 datos 2

[07:43] y después claro bueno otras figuras ya es este de la misma pero

[07:52] con estacionario

[08:02] me parece que estoy bien

[08:09] me parece que estoy bien

[08:16] bueno tendría que probar eso también pero bueno

[08:21] tengo que hacer un test del simulador

[08:31] de ver los intervalos de simulación

[08:39] a ver si afectan en algo una forma de ver si afectan la likelihood

[08:47] bueno la likelihood es un tema

[08:55] Sí, tendría que ver por algoritmo, ¿no?

[08:58] O sea, por algoritmo, ¿cuánto es la diferencia?

[09:01] Yo creo que eso es un tema uno, ¿no?

[09:03] Porque, a ver, yo lo que quiero ver es, yo estoy comparando por algoritmo.

[09:08] Entonces, claro, sería...

[09:12] Claro, y tomo el IR como máximo.

[09:20] Entonces sería salto con este y con este y con este.

[09:25] Y de ahí ver un poco el resto de las cosas, está bien eso.

[09:37] Claro, pasa que hay cosas absolutas.

[09:40] Por ejemplo, el likelihood...

[09:46] el "Leglish" depende de todo, depende de...

[09:53] ahí si depende del ruido, depende de...

[09:59] varía muchísimo con todos los otros parámetros

[10:05] entonces ahí sí la única comparación que vale es intraparámetros

[10:09] el valor absoluto es medio difícil de interpretar

[10:15] que te dice 200 o 1500 que se yo que significa

[10:20] y no lo que vale es la diferencia

[10:25] y

[10:28] pero bueno lo podes si, o sea que tienes que hacer

[10:34] digamos ahí tendría que ser con rap y con rap y con...

[10:43] pero lo mejor ahí es hacer un Delta con el likelihood

[10:50] después está la pregunta si el Delta likelihood sirve para algo

[11:04] la verdad

[11:06] no sé

[11:09] me parece que no es muy útil

[11:13] O sea, yo lo que está bien, este análisis así global está bien, es ir bien a los cosas

[11:30] que son diagnósticos, entonces son esos tres más un cuarto que te daría un poco idea

[11:37] de mecanismo que es la correlación temporal que sería el otro grado. Y después haría

[11:47] el mismo análisis con el modelo que se llama estacionario y luego me quedaría probar con

[12:07] modelos un poco más complejos como había dicho pero pero yo no sé ahí

[12:13] si vale la pena hacerlo lo haría digamos lo tiro pues se hace un pedo

[12:20] pero bueno ahí

[12:23] ahí tenés más una mayor complejidad

[12:29] porque tenés como dos variables más

[12:37] bueno acá tendríamos una variable más que es la relación on off

[12:43] pero no sé si eso realmente importa mucho podría meterlo

[12:55] para ver en realidad es como que te desmiente como si tuvieras menos canales

[13:01] pero básicamente

[13:06] Sí.

[13:08] Sería como la P de equilibrio por coso o algo así.

[13:13] Sí, eso podría servir para eso de ahí, para el "vayas".

[13:21] Sí, ahí podría servir eso.

[13:24] Y ahí es un poco importante entender el...

[13:36] el bias, claro que no es la constante sino el número de aperturas, pero sí, ahí podría

[13:47] meter el tema este, de la eficacia sería, o sea, mayor o menor eficacia. Claro, lo podés

[13:59] poner con concentración, ¿no? O sea, en realidad... Podríamos manejar con concentración.

[14:10] No, la eficacia acá no hay, eficacia pues no. Eso sería en otro tipo de modelo.

[14:19] si

[14:25] e

[14:34] yo creo que está bien así no me

[14:37] me meta mucho en eso

[14:41] yo lo dejaría como para un reviewer

[14:45] me pregunte algo

[14:46] Está bien, así.

[14:52] Pero podría ser, la verdad.

[14:57] Podría ser un experimento específico para eso, atacar el tema de la violencia.

[15:08] Bueno en realidad todos estarían un poco afectados.

[15:13] El valor es lo que te da error, el valor que te da esta mal, de todos modos un poco el

[15:29] mensaje es que tenes que hacer la corrección, si el valor es diferente tenes que hacer una

[15:35] corrección, punto. Calculas y listo y bueno, qué sé yo. Después, después, te va a poner

[15:46] una cota en el error de alguna manera, qué sé yo, supongo que es la mitad de la corrección,

[15:50] está mal, no sé, algo así. Y listo.

[16:14] si tengo que ver eso ahí que pasó pero bueno

[16:23] si tengo que ver qué pasó con los con los de con el nr

[16:34] a ver qué pasó después la B corta creo que la voy a volar

[16:40] no tiene mucho sentido

[16:47] me queda eso ahí

[16:51] y bueno

[16:57] estamos

[17:02] no hay mucho más que

[17:07] Descero.
