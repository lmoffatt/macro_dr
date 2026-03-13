Aquí tienes la organización y transcripción de los archivos de audio proporcionados.

### 1. Organización de los Archivos

Basado en la nomenclatura de los archivos (secuencia numérica `WA0004` a `WA0017` del mismo día 2025-10-31) y la progresión lógica del contenido (desde detalles técnicos específicos sobre distancias KL y priors, pasando por el estado actual del código, hasta la decisión estratégica de cerrar el proyecto), el orden cronológico y lógico es el siguiente:

1. **PTT-20251031-WA0004.opus**: Planteo inicial sobre la distancia KL y la simplificación de cálculos con muchos datos.
2. **PTT-20251031-WA0005.opus**: Hipótesis sobre *priors* anchos y la posibilidad de un tercer *paper*.
3. **PTT-20251031-WA0008.opus**: Resumen de los audios anteriores (menciona "hice como tres audios"), consolidación del hallazgo teórico sobre *KL distance* y definición del tercer punto del trabajo.
4. **PTT-20251031-WA0009.opus**: Diagnóstico del estado actual ("monstruoso") del software MacroIR y plan de acción inmediato (tests).
5. **PTT-20251031-WA0010.opus**: Estrategia de estabilización ("clavos") y la intención futura de "cerrar el boliche".
6. **PTT-20251031-WA0017.opus**: Reflexión final sobre el cierre del proyecto (analogía con cerrar una empresa), limpieza del código para usuarios externos y transición a nuevos intereses ("proyecto de país").

---

### 2. Transcripción Verbatim

#### Archivo: PTT-20251031-WA0004.opus

**[Hablante 1]**
Y un poco eso de la KL M distance también indica cuándo es al pedo calcular la evidencia, ¿no? porque si la KL M distance por ejemplo con mucha cantidad de datos termina siendo qué sé yo, el número de parámetros dividido dos, no sé una cosa así, que que sería como como la aproximación de de Euler o no sé de Poisson, no me acuerdo, hay una una fórmula de la evidencia, eh, eso sería interesante también ver esos límites, ¿no? o sea cuándo ya tenés tantos datos que la los cálculos quizás se se simplifican, ¿no? Eh. Eso también es una una buena pregunta a hacerse.

---

#### Archivo: PTT-20251031-WA0005.opus

**[Hablante 1]**
Claro una una de las hipótesis es que si vos tenés priors más anchos, vos podrías no diferenciar este, o sea vos estarías castigando eh más firmemente eh modelos con más parámetros. Eh, ese ese sería un poco una de las eh hipótesis a a testear. O sea yo creo que casi casi que eso podría ser un un paper digamos este que analice todas estas cosas a fondo, ¿no? o sea. No sé si realmente no da para un un tercer paper, es un paper extra. Este, nada tendría que eso eso pensarlo un poco. Porque realmente da, da para bastante y y digamos eh quizás hacer un trabajo sistemático no esté mal. Pero no estaría mal la digamos dejar planteado el tema y como para después eh por ahí que otra gente lo haga esto, chau, eso me estio sistemático después simplemente, bueno, tirar la piedra.

---

#### Archivo: PTT-20251031-WA0008.opus

**[Hablante 1]**
Bueno, creo que hice como tres audios de diez minutos en Ideas hablando de la multiplicación de matrices y y eso como como a ver de de ahí salga algo alguna idea que permita avanzar algo en a en en algo, ¿no? o sea. En fin. Y bueno y llegué a algo que es aplicable aquí y ahora eh para MacroIR. Y lo siguiente, digamos, o sea ya me había olvidado pero eh uno de los resultados teóricos más importantes que saqué en estos últimos tiempos es darme cuenta que la diferencia entre la expected eh la creo que es la expected posterior likelihood, log no, la sí, la posterior log likelihood y la prior log likelihood o algo así, eh, es la eh KL distance entre el prior y el posterior.

**[Hablante 1]**
Entonces, eso qué implica, que de cualquier experimento que hagamos con un modelo vamos a tener dos variables, o sea la likelihood la evidencia y la KL divergence. Y eso lo podemos ver como dos parámetros dos este sí parámetros o variables lo afectan. Uno es el ancho y la ubicación, la distancia entre el prior y el y el digamos y no el posterior sino el el verdadero valor este, cómo eso me afecta y me eh eh y eso y en función de del número de datos, ¿no? O sea como para ver un poco el tema de cuánto pesa el prior, ¿no? Cómo hacer para que que para entender, digamos, el peso del prior, cuál es un un prior realmente no informativo, o sea eh y especialmente eso en el caso de de modelos más complejos que que son los que digamos estos de de muchas eh interacciones y todo que ya ya rondan en lo en lo muy dudoso.

**[Hablante 1]**
Eh, yo creo que ese tema eh sería excelente como digamos este para cerrar el paper que estoy haciendo, que digamos para presentar el método, porque eh realmente eh cierra una duda que yo tengo este y que no y este que de alguna manera digamos lograr algún tipo de de relación este entre digamos eh digamos cuán errado está el prior eh cuán uninformative es y y cuántos datos tenés y eh eso afecta la digamos el el verdadero valor que vos encontrás, etcétera, ¿no? o sea. Eh me parece que eso eso sería lo digamos el el tercer punto, ¿no?

**[Hablante 1]**
Porque el paper habíamos dicho que era la estimación de la que que el FIM sea eh o sea validar el el la likelihood, el otro es la la cross correlation, pero claro, la tercera sería ver un poco este agarrar un un este un modelo y y claro y ver digamos una especie de claro de validación cruzada pero ya con con priors este ridículos, ¿no? priors apartados de la realidad, o sea priors este serían voy a encontrar un término, ¿no? Priors errados, serían falsos priors o no sé cómo llamarlo. O sea, habría que encontrar un término que que indique realmente de de cuando vos estás meando fuera del tarro muy mal. Eh, sí. Claro, sí, sí, eso digamos este que la idea es hacer dos priors posiblemente uno con mucha más varianza que el otro y chau. Este, sí. Eh, hm. Ta, eso sería.

---

#### Archivo: PTT-20251031-WA0009.opus

**[Hablante 1]**
Bueno. Última parada con MacroIR. Eh. ¿En qué estado está MacroIR? Bueno, está en un estado extraño porque es una especie de mo programa monstruoso que que gran parte de lo los eh de las funciones más importantes son funciones estáticas de una clase que no tiene miembros, una cosa medio extraña. Eh no me acuerdo por qué hice eso. Eh, la cuestión es que un poco la política que tomé es la de tocar lo menos posible, eh asegurarme que funcione, o sea poner test que indiquen que funcione, tratar de establecer algún tipo de de línea de comandos que sea estable y bueno y este y y hacer exactamente eso, o sea los los test que indican que el algoritmo es confiable.

**[Hablante 1]**
Y luego una vez que tenga todo eso, bueno tratar de de este reorganizar el programa para que sea un poco eh más tratable, o sea que tenga una estructura un poquito este mejor. Eh para que digamos pueda este continuar existiendo o lo que sea. Eh. Eh, claramente a ver en qué estado estamos con MacroIR. Bueno, logramos sacar un paper, lo cual no es poco. Este igualmente yo todavía no hice los testeos de del algoritmo que me gustaría hacer, que eso es el segundo paper que voy a hacer, después voy a hacer un tercer paper con toda esta serie de estudios que yo había planteado respecto de cómo afecta el prior y el número de datos la evidencia y la la contracción del del prior, o sea la la ganancia de información, no sabría qué ponerle un nombre, que sería eh digamos lo que sería lo que uno aprende acerca de los de los parámetros, ¿no? O sea eh que que es la esa diferencia sí, la contracción, o sea es la KL distance entre ey [audio se corta abruptamente].

---

#### Archivo: PTT-20251031-WA0010.opus

**[Hablante 1]**
Bueno, a ver. Entonces, ¿cuál es la situación de MacroIR? Bueno, tengo ese programa así medio inentendible y complejo. Entonces, la idea es bueno, más o menos este tratar de así con unos clavos, digamos, yo me imagino como unos clavos que son estos test, que fijan comportamientos, cosa que digamos yo pueda después empezar a cambiarlo, eventualmente, si quiero. Eh. Pero bueno un poco yo ya veo que esto va a llegar a un punto donde no va quedar ahí. Va a quedar fosilizado en algo que que terminaré con estos este este paper y quizás el paper con Cecilia Bouzat y algún otro paper más con con Gustavo. Y después bueno, que se arreglen los la gente que lo quiera usar. Yo después de eso, la idea mía es cerrar ese boliche y meterme con eh con la fundación me sale decir, que sería no sería la fundación sería a ver bueno, eso voy a hablar ahora en otro lado, a ver. Qué es lo que voy a hacer.

---

#### Archivo: PTT-20251031-WA0017.opus

**[Hablante 1]**
Bueno, más o menos, eh quedó qué es lo que voy a hacer ahora que termine MacroIR, cierre esta etapa. Una etapa que que ya tengo que planificar la el cierre de la empresa, como me acuerdo que una vez hice un hicimos un curso de de startup o qué sé yo. Y había uno que que se puso a charlar que era un alguien que había hecho una empresa y todo nos contó su experiencia. Y una cosa que me quedó es que bueno, que el chabón había hecho una empresa y que bueno vino el 2001 y le dijeron los socios dijeron no, no nos queremos ir, hay que cerrar. Y él pataleó y la socia le dijo, mirá, tenés que cerrar, es así, o sea si se quieren ir, tenés que cerrar. No hay tu tía, y entonces se la pasó un año entero cerrando la empresa esa, qué sé yo.

**[Hablante 1]**
Y bueno, qué sé yo, o sea, yo creo que Macro R, MacroIR es un proyecto que que voy a tener que cerrar, entonces tengo que dejarlo digamos empaquetado de una manera decente, que la gente no me odie demasiado después. Eh. Qué sé yo. Entonces, bueno, esa sería un poco poco la idea. Eh, también puede ser que algo lo usen y que eso me dé alguna guita. Digamos, eso no estaría mal después de haber puesto tanto esfuerzo en eso. Si bien, digamos, un poco la idea es que eso me sirvió eh de entrenamiento como para hacer el proyecto este de proyecto de país que que es el que realmente me interesa hacer. Eh que es como me encanta así pensarlo como un modelo total en que todo pase por ahí. Eso me me fascina, es como como todo es para un guion. Pero bueno. En fin.

**[Hablante 1]**
Eh, claro es el tratado del mundo de mi papá, claro, el tratado del mundo. Sí. Es eso. Es un proyecto que me, no, es alucinante, no sé. Sí. Bueno. Entonces, bueno, ¿cómo cerramos MacroIR? Bueno, entonces lo que hay que tener es básicamente comandos, que la gente pueda usar, que la gente pueda entender lo que hacen. Que sean fáciles de entender, que fase sean fáciles de comprobar que digamos que lo que hacen sea correcto. Y después bueno, adentro de eso es una maraña de cables que agarrate si lo entendés, pero bueno, eh pero bueno es este eso después de que yo haga todo todo el cableado de todas las eh los comandos y y cómo se hace funcionan cada uno, después eso eh lo voy a ir este ajustando un poquito.

**[Hablante 1]**
Pero bueno, no quiero perderme hace tiempo porque yo sé que si me pongo a ajustar voy a tardar mucho. Este y la idea es esto irlo cerrando, porque yo creo que ya cumplió su ciclo, o sea tuvo muchos años con esto, bueno, me sirvió para algunas cosas. Para fundamentalmente estar dormido y no pensar, pero bueno. Eh es hora de de poner mi mis habilidades al servicio de la patria como diría Candela. Y sí, qué sé yo. No sé. Ver qué hacer algo un poco más interesante o qué sé yo, no sé, interactuar con más gente, no sé. Pero bueno, la cuestión es que eh nada, eh esa es la idea, ¿no? O sea la idea es a ver eh tener estos cosa, porque claro, si yo lo hago al revés digo bueno, no, acá lo que el motor, la guía, el carro que guía es terminar un nuevo paper, eso es más fácil.

**[Hablante 1]**
Eh digamos es un objetivo bien claro, terminar el paper, pum. Entonces qué tengo que hacer para terminar el paper. Y bueno, y entre otras cosas unas cosas que tengo que hacer para terminar el paper es tener un un cli más o menos decente que que la gente pueda replicarlo. Digamos, yo creo que eso está bien, eso es lo que tiene que ser y eso es lo que va a ser. Eh entonces sí, fundamentalmente y y que hacer los gráficos sea digamos un algo digamos que se puede escribir en algún script o lenguaje, ¿no? o sea. Bueno. En principio eh podría ser en en R o o un makefile o o algo así. Eso eso es una cosa que que bueno, igual todavía no llegué a ese punto pero pero no estaría mal tener digamos alguna especie de de makefile que que instale todo más o menos automáticamente.

**[Hablante 1]**
Claro si yo hago como una especie de docker entonces todo eso es más o menos automático y rápido. Sí. Igual eso creo que es un paso más allá que digamos ya ya ir hacia eso de entrada es un poco demasiado, pero bueno. Pero sí la idea de bueno, tenerlo este el programa que eh bueno, sí. Este que responda preguntas, más o menos rápido. Eh, bueno tengo que sacar estos dos papers y el tercero que sería el el de eso. Claro, pasa que el tercero bueno yo ya ahí me asusto de que no me va a dar los tiempos o lo que sea. Bueno, ahora tengo tengo más velocidad de cómputo. Eh acceso a cómputo y nada, eso tendría que ver qué tan rápido es. Este bien.

**[Hablante 1]**
Y bueno y además tengo que ver de de optimizar un poco Macro de R. Eh claro, sí, o sea, bueno ahí lo que tendría que hacer es hacerlo correr en en GPUs, pero eh es un poco mucho, no sé, o o tratar de optimizar el algoritmo, no sé. Eso qué sé yo, no sé. Me atrae bastante, o sea sería lindo hacerlo. No sé si lo voy a cerrar todo totalmente eh porque bueno esas cosas no, pero creo que lo mejor es es este si querés cerrarlo ahí. Liberar mis neuronas para para otras cosas, ¿no? Sí, para hacer eh sí modelos eh quizás de otro estilo, no tan detallistas, no tan eh digamos zarpados, sino por ahí un poco más más globales, más eh más rápidos, digamos, ¿no? o sea que no sean unas cosas tan lentas como esto que hago que es un una carreta.
