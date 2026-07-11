**Archivo:** WhatsApp Ptt 2026-06-09 at 14.14.42.mp3
**Duración:** ~05:22
**Hablante:** Luciano

[00:01] Bueno, siguen las aventuras de Macro IR.

[00:05] Estuve bastante dándole, dándole, porque tratando de ver dónde está el supuesto bug

[00:16] que me genera esa variabilidad enorme en la matriz de Fischer calculada numéricamente.

[00:23] Bueno, al final, usando la guía, logré varios comandos que me permiten detectar

[00:32] cuándo es que hay un salto en la derivada segunda, y bueno, logré encontrar que el salto

[00:40] está dado por la función Truss Coefficient.

[00:47] Y bueno, yo lo que esta mañana llegué a mi oficina y sin abrir la IA nada,

[00:55] me puse a grabatear en el pizarrón, como debía ser la fórmula,

[01:03] O sea, que bueno, que sería, en realidad es el máximo entre el cociente, entre la distancia a recorrer y la probabilidad para cada estado y uno.

[01:27] es el mínimo entre 1 y eso

[01:32] entonces te da x para menores que 1

[01:36] y 1 para mayores que 1

[01:39] pero claro, el tema después me di cuenta

[01:44] es que yo tengo que ver para todas las direcciones

[01:47] ver cuál es el mínimo, entonces tengo lo que llaman el soft mínimo

[01:51] y al final la ecuación

[01:55] porque yo traté de hacer una ecuación que sea diferenciable y todo,

[01:59] y bueno, lo que me quedó es una cosa muy muy bonita,

[02:03] muy parecida a una función de partición, no sé cómo llamarla,

[02:09] que es lo que llaman el log sum exp,

[02:14] es decir, el logaritmo de la suma de los exponenciales,

[02:19] y entonces vos tenés un factor K,

[02:22] que si el K es negativo te da el mínimo,

[02:24] si K es positivo te da el máximo y es 1/K por la sumatoria de la exponencial de K por la variable que tenés.

[02:39] Y lo que voy a hacer es que una de las variables sea 1.

[02:43] O sea, me quedan todos los ratios de distancia dividido, o sea, desplazamiento dividido,

[02:53] probabilidad actual, que obviamente no tiene que ser mayor que 1,

[03:01] o digamos 1 menos p, o sea, dependiendo de si el desplazamiento es positivo o negativo.

[03:09] Y a eso, bueno, lo exponencio, lo multiplico por K, lo exponencio, lo sumo,

[03:16] y además le sumo uno que sea a la K.

[03:19] Y después todo eso lo divido, lo saco al logaritmo y lo divido por K.

[03:27] Y entonces eso me da el valor mínimo.

[03:32] Y es una hermosa ecuación.

[03:35] ahora lo que tengo que hacer es implementarla con derivadas y mi derivada la verdad que está muy linda

[03:41] ahora lo que me queda ver es porque también hay un trast

[03:47] positive

[03:50] definite semi definite

[03:55] que es positive

[03:58] symmetric definite

[04:00] que tengo que ver

[04:03] si también tengo que hacer una corrección

[04:07] porque

[04:08] me obligó a hacer la corrección

[04:10] pero yo la verdad que no estoy muy seguro

[04:12] de que haga falta

[04:14] creo que decía que hacía falta

[04:16] por una cuestión numérica

[04:17] porque ahí digamos

[04:20] supuestamente eso no funciona mal

[04:22] pero bueno, quiero hacerlo de una cosa

[04:24] más linda, más prolija con primeros principios

[04:26] como logré hacerlo

[04:28] para la media

[04:30] así que bueno, salí un poco a caminar

[04:32] porque ya estaba un poco saturado el cerebro de tanta información

[04:38] y tanta cosa que no había definido cómo manejar.

[04:45] Y bueno, ahora que ya lo tengo definido, no sé,

[04:47] daré una vuelta por acá por el parque y me iré a trabajar y a definir esto.

[04:53] Estoy un poco cansado, la verdad,

[04:57] Todo esto de compilar todo este raíz para encontrar dónde era el punto que había problemas, me llevó un montón de tiempo.

[05:09] La verdad que esto se lo podría haber pensado de una, pero bueno, me siento más cómodo teniendo algo en el cual yo pueda identificar dónde es que hay un salto en la derivada.
