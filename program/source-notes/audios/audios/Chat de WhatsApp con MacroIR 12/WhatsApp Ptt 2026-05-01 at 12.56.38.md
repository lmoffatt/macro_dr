**Archivo:** WhatsApp Ptt 2026-05-01 at 12.56.38.mp3
**Duración:** ~15:57
**Hablante:** Luciano

[00:00] Ya llegando a 9 de julio, estoy cruzando en nuestro viejo estacionamiento.

[00:06] Bueno, no sé, tuve que salir a caminar, ya estoy tratando de surcir los finales detalles de micro IR.

[00:19] O sea, había una versión que hizo Claude.

[00:25] ahora esto es una versión

[00:27] toda diseñada por mí

[00:29] porque no me gusta

[00:30] la arquitectura de cloud

[00:32] quiero hacer algo que yo mismo

[00:34] entienda

[00:36] o que siga mi estilo

[00:39] mis prejuicios

[00:43] y bueno entonces

[00:45] vuelvo a decir

[00:48] que la gran

[00:49] contribución

[00:54] sino la gran idea es cómo unir, cómo generar la matriz de transición microscópica a partir

[01:06] de la matriz de transición individual. Que con una simple operación de unir dos matrices,

[01:16] o sea, donde vos las multiplicás, digamos, elemento a elemento, pero después unís dos elementos que son indistinguibles,

[01:26] básicamente por la cantidad de estados, digamos, por la distribución de estados.

[01:36] Entonces, bueno, así uno puede generar recursivamente la matriz de transición microscópica de muchos, la matriz ensamble de transición.

[01:54] y con una serie de cosas como por ejemplo el hacer el square, no sé cómo se llama,

[02:05] que vos, por ejemplo, para generar 8 canales, primero generás 1 más 1,

[02:10] después haces, con eso tenés 2, después haces 2 veces 2, tenés 4, y haces 2 veces 4 y tenés 8,

[02:16] o sea, con eso minimizás el número de operaciones.

[02:19] Bueno, eso es para generar lo que sería la QDT.

[02:24] Bien, o sea, que ahora lo que estoy viendo es ya dentro del ciclo de recursión de microarras,

[02:33] digamos, como voy generando los diagnósticos.

[02:36] Entonces, yo tengo ahí un montón de variables que puedo tomar,

[02:47] que básicamente tienen que ver con la prioridad posterior y dentro del intervalo,

[02:57] la distribución de prioridades, que es con el P-mini y el P-cov,

[03:05] que en este caso sí lo generaría a partir de justamente el P mayúscula,

[03:13] que sería el P-micro, no sé cómo llamarlo, P-micro,

[03:17] Entonces la idea es que el Pmicro lo transformo en Pmini-Quepeco.

[03:24] Voy a tener un Pmicro Square que sería Pmicro Inicial-Final.

[03:36] permite el intervalo

[03:39] o no sé

[03:40] voy a llamarlo

[03:41] que es lo que yo no

[03:43] no logro

[03:45] no instancio

[03:49] en el caso del macroscópico

[03:51] eso es una cosa virtual

[03:54] que nunca la instancio

[03:55] acá la instancio la tengo que instanciar

[03:57] también voy a tener que instanciar

[03:59] bueno, tengo instanciada

[04:01] la

[04:01] la matriz de

[04:05] se llama de conductancias

[04:08] entonces

[04:09] lo que voy a hacer es

[04:12] esa matriz de conductancias

[04:15] transformarla en

[04:17] una matriz de

[04:20] likelihood individuales

[04:22] de likelihood por estado

[04:23] y después con eso

[04:25] sacar

[04:26] el posterior

[04:28] o sea, posterior likelihood

[04:31] y la suma de eso

[04:34] sería el posterior pemin posterior y después de con... después puedo sacar la... está bien ahí puedo... todas las variables las voy a crear, así tengo nombres y...

[05:04] y está todo más o menos como se dice el abeleado.

[05:09] Podría ponerles un nombre ahí nomás, pero eso es una buena pregunta, ¿no?

[05:15] Si podría dar la pena una clase para eso.

[05:19] Y bueno, entonces ya tengo todo, digamos, porque ese es el proceso original, ¿no?

[05:26] O sea, claro, lo que no tengo es, a ver, tendría que ver si yo puedo generar

[05:32] el IGS o no sé cómo se llama que es la dirección de update. Yo creo que sí la puedo generar

[05:50] a partir supongo yo de la media...

[05:58] Claro, tengo las dos, la GS0 y la GS.

[06:05] Las puedo generar a partir de la PEMIN posterior y PEMIN prior.

[06:21] y el otro cojo. O sea ahí lo que puedo hacer es ver digamos, puedo proyectar cada partecita,

[06:39] todo el mecanismo de Macro R, de Calman, yo lo proyecto, o sea, comparo digamos, la

[06:57] aplicación digamos sería

[07:01] la aplicación ingenua de

[07:11] el algoritmo de Kalman en el estado digamos macroscopizado del estado

[07:19] microscópico comparado con el coso verdadero.

[07:26] Yo tengo un functor que me manda de

[07:32] microscópico a macroscópico y entonces yo lo que puedo hacer es comparar ambos

[07:37] caminos y ver un poco cómo divergen.

[07:41] Eso es fácil. Yo tengo por ejemplo el calpe mililisial, entonces puedo comparar

[07:47] el camino macroscópico con el microscópico a ver cómo proyectan.

[07:54] Claro, eso está bueno, es fácil de hacer y es diagnóstico de alguna manera.

[08:07] yo creo que esto un poco podría superar lo que sería este paper solo porque ya se va

[08:17] muy a la mierda con los detalles yo creo que eso sería un estudio aparte

[08:26] de ver

[08:29] cómo

[08:32] cómo lo macaroscópico altera

[08:36] digamos, perturba las cosas

[08:40] es una herramienta poderosa esto

[08:43] la pregunta es si está dentro del

[08:48] paper o no

[08:49] y ahí lo que importa es

[08:52] si se puede

[08:54] yo puedo eliminarlo

[08:55] que el paper siga teniendo fuerza.

[08:58] ¿Cuánto yo estoy agregando

[09:04] al paper

[09:06] y cuánto estoy obteniendo?

[09:08] Es el retorno

[09:10] al investment.

[09:11] Yo creo que

[09:16] digamos,

[09:17] desde el punto de vista de la carga

[09:19] conceptual

[09:21] para un lector

[09:23] el Skysnes y el Kurtosis está bien. Bueno, puede ser la distancia de KL también podrían

[09:35] dar, tenerlos como indicadores. Pero esto ya es irse muy a la mierda, de ver los detalles

[09:44] Especialmente el GS, el GS sería un poco irme un poquito a la mierda.

[09:54] Pero bueno, digamos, el detalle no está mal, o sea, yo creo que lo tengo que implementar

[10:01] porque lo implemento ahora y ya está, ya me queda como un paper subsiguiente que sería

[10:09] simplemente bueno analizar esa parte, o sea, cómo digamos la...

[10:17] como es el mecanismo con el cual se pierde la...

[10:28] ¿cómo se llama esto?

[10:31] la ergodicidad sería, no sé cómo se llama, pero bueno el hecho de que tengo una memoria,

[10:38] Sí, sería la...

[10:40] ¿Cómo se llama?

[10:41] La propiedad de Markov.

[10:43] De alguna manera se pierde la propiedad de Markov

[10:45] en el sentido que

[10:47] tiene memoria.

[10:49] O sea,

[10:50] el macro R tiene memoria.

[10:54] No es como micro R que todo está

[10:55] en la memoria

[10:58] del Estado.

[11:00] No.

[11:01] Vos hay una par de graper de eso.

[11:04] Y eso lo pagas

[11:05] en correlación temporal

[11:10] bueno, eso hay que ser directamente abiertos

[11:16] me parece que no tiene sentido ocultarlo

[11:20] se paga, se paga

[11:23] me parece que está bien

[11:29] en el fondo lo que queremos ver

[11:32] más que nada

[11:33] es el tema de la resolución de parámetros, porque nuestra idea no es estudiar estadística

[11:41] acá, sino de verlo como método para recuperar parámetros. Entonces la pregunta es, ¿cuánto

[11:52] pierdo en la covarianza? Claro, pero acá en este caso si estaría solamente analizando

[11:59] parámetros en el siguiente paper, estoy analizando evidencia.

[12:04] Ya es otro tema más fuerte.

[12:09] Bien, entonces está bien.

[12:14] Yo creo que el límite lo ponemos ahí en el MKL y el GS si quedará para un próximo

[12:25] Voy a dejarlo como una especie de superdiagnosis, vamos a ver.

[12:39] ¿Qué más?

[12:41] Claro, porque es el tema ese.

[12:44] El tema es, el microscópico es un diagnóstico del macroscópico, básicamente.

[12:51] Básicamente la idea es esa.

[12:54] ¿Por qué tengo la ilusión, entre comillas, verdadera?

[13:02] Que no es la verdadera porque, bueno, estamos haciendo esa aproximación de la ilusión normal,

[13:12] o sea, de la varianza, ¿no?

[13:15] O sea, para hacerlo realmente en serio no sé cómo sería.

[13:21] Eso quizás es lo que no se, tendría que ver qué es lo que quise hacer con el microscópico

[13:30] estocástico que...

[13:32] Sí, eso es lo que tengo que ver, tengo que ver qué carajo quise hacer con el microscópico

[13:36] estocástico, porque también sería otra cosa más para ver.

[13:43] Bien, entonces ahí cierra más o menos todo.

[13:50] Y ahora concentrarme en terminar esta mierda de una vez y ya enviarlo.

[14:06] Esto que he implementado, corrido, ya lo tengo que correr hoy mismo, esto tiene que

[14:13] corriendo y ya vio de generar las figuras y enviarlo lunes o martes o sea tiene que por

[14:25] lo menos estar las figuras ya seguro listas y ya digamos el paper básicamente es bien

[14:37] iniquitado.

[14:38] O sea, sí, hasta uno podría decir, bueno, excluyo el microscópico recursivo, podría

[14:45] ser.

[14:46] Me quedo hasta ahí porque ya estoy, digamos, seguro de que esto está bien, digamos, un

[14:53] método está bien.

[14:56] Pero me parece que está bien, mejor le prometo, entonces por esto tengo que entrar en eLife

[15:01] seguro.

[15:02] y el tema más que nada es de que no tiene datos experimentales pero bueno ya fue

[15:09] eso no es el tema digamos

[15:15] acá la idea es

[15:18] es otra, digamos, surfesear, ¿no? Sería sacar a la superficie la metodología para testear

[15:34] estos algoritmos de la forma central eléctrica de una manera robusta y más o menos clara,

[15:43] con una metodología que permite avanzar.

[15:45] Eso es un tema muy importante, ¿ves?

[15:48] O sea, la ciencia buena es la ciencia que avanza.

[15:53] Ah, ahí hay una idea que se me ocurrió cuando estaba...
