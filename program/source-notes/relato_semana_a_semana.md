# MacroIR: relato semana a semana del paper

Reconstruido leyendo los 274 audios de los grupos "MacroIR" 1 a 13, del 26 de agosto de 2025 al 2 de agosto de 2026.

**Cobertura y advertencias.**
- No hay audios de MacroIR entre el 16 de enero y el 22 de marzo de 2026. Buscado en todo el disco: no existe un grupo "MacroIR 11" en ninguna parte, ni exportado ni en el disco viejo. En esa ventana sí hay unos 190 audios grabados, pero en `Ideas 6` a `Ideas 9` y en `plan 2026`, y son de otros temas. El único que habla del paper es el del 6 de marzo (Ideas 8), y explica el hueco: "Macro IR me quedó totalmente fuera de mi mente por alguna razón". En ese mismo audio ya está formulada la pregunta del paper tal como aparece el 23 de marzo, incluida la descomposición de la matriz de distorsión en componente de correlación y componente de sample. Conclusión: no es material perdido, son dos meses en los que el paper estuvo parado.
- Un solo audio no tiene transcripción: `PTT-20250926-WA0011.opus`.
- Las transcripciones de agosto a diciembre de 2025 las hizo Gemini en tandas. Casi todas son verbatim, pero un bloque (27 ago WA0013/0018/0019, y los del 2, 4 y 8 de septiembre) sale mucho más corto y más pulido que el resto, con un estilo que no es el de Luciano hablando. Los archivos de audio existen, así que no son inventados, pero conviene volver al audio antes de citar esos textualmente. De diciembre en adelante son transcripciones con whisper.cpp y con marcas de tiempo, ruidosas pero fieles.
- Archivos de trabajo con las transcripciones concatenadas en orden de fecha: `tmp/audio_relato/P5..P8*.txt`.

---

## Etapa 0. Antes de que hubiera un paper (26 ago 2025 – 29 sep 2025)

### Semana del 25 de agosto (26 y 27 de agosto)
Arranca caminando por la 9 de Julio, sin saber bien qué quiere decir, planificando qué puede publicar en los próximos meses. Tiene dos colaboraciones en marcha, Gustavo Pierdominici y Cecilia Bouzat, y aparece por primera vez la idea de un trabajo propio: "una colaboración conmigo mismo", publicar todos los controles y pruebas que verifican que el algoritmo es correcto. Lo describe como "un paper largo y tedioso" para bioRxiv.

En el mismo bloque queda fijado el test que va a sostener todo el resto del año: la esperanza del score tiene que ser cero y la covarianza del score tiene que ser igual a la esperanza de la matriz de información de Fisher. También queda listada la arquitectura de comandos (Evidence, Likelihood, Simulation, Sampling, Idealization) y la decisión operativa que va a ordenar el trabajo: dos versiones del código, la vieja para largar corridas ya en el cluster y la nueva para el trabajo de caracterización.

Hay una admisión que vale la pena guardar, porque explica meses de demora posteriores: "yo necesitaba recrear todo MacroR como sentirme que hacía algo nuevo; hay una cosa ahí que la parte mía emocional va en contra de lo práctico".

### Semanas del 1 y 8 de septiembre
Vuelve después de unos días. La versión vieja ya corre en el cluster. Encuentra y corrige un error en la derivada de la likelihood con estados redundantes, y con eso el test de Fisher da bien en modelos simples. Aparece el primer problema estructural: al aumentar la complejidad del modelo, la FIM se vuelve muy mal condicionada.

El 8 de septiembre escribe el primer esqueleto del paper de caracterización, en cinco secciones: framework, teoría de los tests de consistencia, precisión contra velocidad, selección de modelos y priors, y un caso de uso alostérico contra secuencial. Es un paper distinto del que va a terminar escribiendo, pero el punto 2 sobrevive hasta el final.

### Semanas del 15 y 22 de septiembre
Dos semanas casi enteramente de software. Rediseña el DSL: tipos, funciones y variables, el environment guardado en JSON, el help como caracterización completa del programa. Aparece la formulación que le va a servir después: un test es una tripleta (función, postcondición, dominio de inputs), y lo que hay que encontrar no es si la función anda sino en qué región del espacio de parámetros anda, con qué precisión y a qué velocidad.

El 20 de septiembre, sábado a la mañana, con los pajaritos de fondo, comenta al pasar que ya envió la versión final del manuscrito anterior. Queda libre para MacroR.

El 24 y 25 tiene sesiones intensas con Codex y consigue integración continua en GitHub y los primeros comandos que exponen el interior del programa (`patch state`, autovalores, probabilidad de transición). El 26 se mete de lleno en Homotopy Type Theory: qué es un dominio, qué es un prior, qué es un significado, qué es un objeto. Es el punto más lejano del paper en todo el año.

### Semana del 29 de septiembre: primer punto de inflexión
Se queda sin tokens de Codex y se frena en seco. El audio del 29 es el primero donde dice en voz alta que trabaja mucho y no avanza, y donde toma la decisión que reorienta todo:

> "Mi objetivo primario es tener confianza en mis corridas. Entonces tendría que sí o sí testear la likelihood con la FIM y el score. [...] Si yo tengo test de Fisher Information Matrix, del sampling y de la evidencia, con esas tres cosas yo podría plantear un paper."

Se pone una fecha, el 15 de octubre, y manda la hotificación a segundo plano. A partir de acá hay un paper.

---

## Etapa 1. El paper existe, pero todavía es otro (2 oct – 31 oct 2025)

### Semana del 29 de septiembre (2 de octubre)
Define la secuencia concreta: correr muchas simulaciones, calcular likelihood con score y FIM en cada una, y sobre esa muestra hacer el test. Descarta explícitamente dos tentaciones: recuperar el Levenberg-Marquardt y recuperar cumulative evidence.

### Semana del 6 de octubre (7 de octubre)
Descubre que su sistema de derivadas tiene punteros colgando en el delta X, con matrices de tamaño absurdo apareciendo de la nada. Lo tapa pero no lo arregla. Además le da una diferencia entre la log-likelihood directa y la calculada, y no sabe de dónde sale.

### Semana del 13 de octubre (17 y 20 de octubre)
Andrea le dice que es importante publicar un segundo paper. Lleva dos semanas atascado con las derivadas. El diagnóstico llega el 17: no se puede derivar autovectores porque no son una función sino un espacio. La salida es no derivarlos y derivar directamente las magnitudes que sí están bien definidas, la matriz de probabilidad de transición y las conductancias condicionales al estado inicial y final. Con ChatGPT saca la fórmula de la exponencial matriz por bloques con Padé y le hace escribir un LaTeX con los tres métodos. Considera y descarta backpropagation: "no es la idea de que sea más rápido, la idea es tener algo robusto y confiable".

El plan del paper todavía es modesto: E[score]=0 y Var(score)=FIM, "una nota más o menos corta en Biophysical Journal o algo así y listo, con eso salvo el año".

### Semana del 20 de octubre (22, 23, 28 de octubre): el paper se define por primera vez
El 22 de octubre aparece la primera formulación explícita del contenido:

> "Los dos conceptos que quiero imponer en el paper. Uno: MacroIR es mejor que MacroR, MacroNR y MacroINR. Dos: MacroIR permite diferenciar esquemas cinéticos. [...] El paper son dos figuras: los tests de FIM para esos cuatro algoritmos, y la confusion matrix para cuatro esquemas."

También aparece el problema del forzamiento: qué hay que hacer para que el algoritmo no se vaya a la mierda (forzar la probabilidad a mantenerse estable), y la idea de que ese dato es en sí mismo parte del resultado. Ese forzamiento va a reaparecer ocho meses después como el bug más difícil del año.

El 23 resuelve la derivada de Q(t) haciendo una trampa razonable: saca del test la derivada de la conductancia condicional (que se divide por una P_ij chiquísima) y testea en cambio las dos magnitudes con menos error. De ahí sale, de paso, la reflexión sobre que necesita un álgebra del error de las variables derivadas.

### Semana del 27 de octubre (28 y 31 de octubre)
El 31 saca un resultado teórico que él mismo califica como de los más importantes de los últimos tiempos: la diferencia entre la esperanza de la log-likelihood posterior y la prior es la distancia KL entre prior y posterior. Lo propone como tercer punto del paper (cuánto pesa el prior, qué es un prior realmente no informativo) y en el mismo audio sospecha que da para un paper aparte.

En el mismo día hace un diagnóstico brutal del estado del código ("un programa monstruoso"), define la estrategia de los clavos (tests que fijan comportamiento para poder refactorizar después) y habla de MacroIR como una empresa que hay que cerrar bien, "dejarlo empaquetado de una manera decente, que la gente no me odie demasiado después".

---

## Etapa 2. El algoritmo pasa a ser el paper (5 nov – 19 dic 2025)

### Semana del 3 de noviembre (5, 6, 7 de noviembre): segundo punto de inflexión
El 5 piensa figuras: una infografía del algoritmo y un test gráfico por cada parte, con la exigencia de que el test muestre que es un test, es decir que distinga modelos válidos de inválidos. Cierra con la nota práctica: "no me tengo que enredar, tengo que sacar el paper lo más rápido posible; con la confusion matrix y el FIM, y un par de modelos de juguete, alcanza".

El 6 llega el giro:

> "Ayer tuve un momento de cierto brillo intelectual: me di cuenta de que puedo prescindir completamente de la Monte Carlo Markov Chain y publicar el algoritmo MacroIR basado pura y exclusivamente en la likelihood de procesos estocásticos promediados."

La validación por score y FIM no necesita MCMC. Eso saca del paper toda la maquinaria de evidencia y sampling, y de un plumazo lo vuelve escribible. La confusion matrix, que hasta acá era la mitad del paper, se cae.

El 7 empieza a escribir y se traba con algo revelador: ChatGPT le convierte el planteo en un filtro de Kalman y le vuela los meta-estados. Decide escribirlo pedagógicamente alrededor del boundary state, el estado definido por el par (estado inicial, estado final) del intervalo, "porque los biólogos que quieren entender esto necesitan las explicaciones".

### Semana del 10 de noviembre (10, 12, 14 de noviembre)
El 10 arma la narrativa en cuatro puntos, que es en esencia la que sobrevive hasta agosto: (1) qué problema resuelve, que las mediciones siempre están integradas en el tiempo y eso vale para cualquier proceso markoviano, no solo canales; (2) cuál fue la estrategia, el meta-estado inicial-final; (3) cómo se prueba, el atajo del score contra la FIM en lugar del MCMC; (4) en qué condiciones el algoritmo funciona, que es "la carne del paper".

Ese día también enumera por primera vez el roster completo de algoritmos (R, NR, IR, INR, más las variantes promediadas tipo Münch) y nace la idea de micro R / micro IR como gold standard, con el comentario de que complicaría el trabajo pero daría la verdad de referencia.

El 12 pasa algo importante para el paper siguiente y para la confianza en este: le pregunta a Gemini si poner un filtro de Bessel después de un Kalman multiplica los estados por cuatro, y la respuesta es que suma cuatro, no multiplica. Se da cuenta de que se había engañado durante años confundiendo ruido coloreado en paralelo con filtrado en serie. Y de paso confirma que su aporte no es trivial: su algoritmo no sale de agregarle una línea a un Kalman, se apoya en la propiedad de Markov, y Münch no resolvió el problema del intervalo.

El 14 planifica los datos de la figura 1 (un experimento de tres puntos, la trayectoria del número de canales abiertos dentro del intervalo) y se mete en el DSL a resolver cómo construir vectores y tuplas. Decide reemplazar ATP por agonist. Aparece también la idea de los vectores indexados, que va a reaparecer en abril.

### Semanas del 17 y 24 de noviembre
El 17 y 18 resuelve técnicamente la figura 1 con template metaprogramming (los simulation tags), que le permite decidir en compile time qué guardar de la simulación, y usar la misma maquinaria para retener prior y posterior del algoritmo. El 18 estima que podría tener manuscrito a principios de diciembre.

El 23 y 25 llega el segundo error de cuentas del año: encuentra, trabajando con DeepSeek, que su corrección de Taylor estaba mal. La fórmula vieja daba drift de la media aun cuando lo medido coincidía con lo esperado, algo que siempre le había llamado la atención. La nueva es compacta e interpretable, con un vector V que suma las contribuciones de la media y de la varianza de la conductancia. Reescribe toda la teoría de MacroIR, MacroTaylor y MacroTaylorIR.

En el mismo audio del 25 hace el cálculo político: son demasiados avances conceptuales para un solo paper, y en el CV cuenta como uno solo. Habla explícitamente de salami slicing, de un techo de cristal en eLife, y de cuatro papers (IR, Taylor, Taylor-IR, Bessel).

El 26 lo formula como "los tres papers del buen humor, Matrix 1, 2 y 3", y admite que los plantearía abiertamente como una serie.

El 27 y 28 aparecen dos ideas grandes que quedan flotando: los intervalos crecientes exponencialmente (posteriors de posteriors sobre varianzas a distintas escalas) y el mecanismo de seguridad del simplex, un factor que infla la varianza de medición cuando la probabilidad se quiere salir de [0,1], que él mismo llama ad hoc "pero con su elegancia dentro de lo ad hoc". El roster se reduce a cinco algoritmos.

### Semanas del 1 y 15 de diciembre
Poca actividad. El 3 de diciembre reorganiza clases y complica el `patch state`. El 19, ya con los datos para la figura 1, escribe el audio más honesto del período: "estoy totalmente bloqueado, no sé bien qué poner". Se desatasca hablando: define qué mostrar (probabilidades de los estados, corriente esperada contra medida, prior y posterior), y descubre dos cosas concretas. Sin la corrección de varianza el gráfico de MacroIR sale cualquier cosa, así que la corrección es obligatoria. Y para que la comparación sea justa hay que centrar MacroR y MacroNR en el medio del intervalo, no al principio.

---

## Etapa 3. De presentar el algoritmo a medir la distorsión (20 dic 2025 – 15 ene 2026)

### Semana del 22 de diciembre
El 20 decide no incluir micro IR, por una razón práctica (necesitaría la varianza condicionada al estado inicial y final, que es lo del paper de Taylor) y una argumental (MacroR ya está justificado contra micro R en el paper de 2007). En el mismo día encuentra el argumento pedagógico más limpio que va a tener:

> "Cuando partís de que todos los canales están cerrados, la medición durante ese primer intervalo no te puede cambiar el prior, porque tu prior ya es único. La única manera de actualizar tu conocimiento es con el estado al final de la medición."

Y el converso, con la corriente cayendo. De ahí que lo ideal sea condicionar a los dos extremos.

El 21 aparece la honestidad metodológica: si presenta MacroR prediciendo al principio del intervalo está construyendo un hombre de paja, así que lo pone en el medio, "en condiciones de combate razonables".

El 23 diseña la figura 2 en el espacio (probabilidad de estar abierto al inicio, probabilidad al final): prior, corriente predicha, likelihood, posterior, aproximación gaussiana.

### Semana del 22 de diciembre, parte 2 (25 y 26 de diciembre)
El 25, en un audio de 16 minutos, tira abajo la figura 2. El argumento es que para justificarla bien necesitaría el microscópico recursivo estocástico, que es abrir una caja de Pandora, y que la figura abre más preguntas de las que responde. Decide reducir la superficie de ataque y pasarse al gradiente. La pregunta del paper se reformula:

> "¿Cuál es el aporte de este paper? En qué situaciones conviene usar MacroIR y en qué situaciones MacroMR. [...] La eficacia es más importante que la eficiencia: primero hay que determinar si el algoritmo da la respuesta correcta o no."

El 26 agrega repeticiones al programa (n simulaciones con seeds independientes) y propone una trazabilidad basada en data frames en vez de scripts.

### Semana del 29 de diciembre: la crisis
El 30 de diciembre, en cuatro minutos, se le cae el plan:

> "Estamos con una crisis. La crisis consiste en que no hay mucha diferencia entre macro MR y macro IR. Es más, hasta parece que macro MR es mejor porque tiene mayor varianza en el gradiente. Entonces voy a tener que optar por la opción nuclear, que es calcular la Fisher Information Matrix."

El test del gradiente solo no distingue los algoritmos. Ese fracaso es lo que empuja el paper hacia su forma definitiva.

### Semanas del 5 y 12 de enero
El 11 de enero, discutiendo con ChatGPT, ordena las dos estimaciones de la FIM (covarianza del score y esperanza del Hessiano) y encuentra la distinción que va a ser el motor del análisis: la covarianza de la suma de los scores contra la suma de las covarianzas de los scores individuales. La diferencia entre ambas mide correlación temporal, es decir cuánta información no estás deconvolucionando.

El 14 y 15 nace el objeto central del paper. Lo describe como un factor de expansión de la varianza, producto de un estimador de Fisher por la inversa del otro, y lo acompaña de dos indicadores más simples: los residuos estandarizados (media cero, varianza uno, sin correlación entre sí) y el bias del score. En el mismo audio observa que ese factor traduce directo a lo que uno quiere estimar, los posteriores de los parámetros y la evidencia. Cierra la semana definiendo qué variables hay que registrar por medición: residuo normalizado, gradiente de la log-likelihood, gradiente de la media, gradiente de la varianza y la varianza.

**Parate: 16 de enero al 22 de marzo de 2026.** No hay audios de MacroIR porque no hubo trabajo sostenido en MacroIR. En esos dos meses graba en otros grupos (Ideas, plan 2026) sobre política, el Valle Plurinacional y una consultoría de network meta-analysis. El 6 de marzo, en Ideas 8, lo dice al pasar: "Macro IR me quedó totalmente fuera de mi mente por alguna razón, que la razón no sé cuál es; tengo que terminar, tengo que ponerme". Y sin embargo, en ese mismo audio la pregunta del paper ya está formulada como va a quedar, con la descomposición de la matriz de distorsión en componente de correlación y componente de sample incluida. Es decir que el planteo del grupo 12 no nace el 23 de marzo, viene armado de antes.

---

## Etapa 4. La maquinaria de validación (23 mar – 1 may 2026)

### Semana del 23 de marzo
Abre el grupo 12 con un roadmap seco: teoría hecha, código hecho, teoría de validación hecha, código de validación hecho; falta correr, analizar, hacer figuras. El 26 aparece una idea lateral que va a durar (fitear tramos cortos con 2 a 4 estados hiperconectados, como un microscopio de cinética) y las variables adimensionales del estudio: τ sobre Δt, ruido de gating sobre ruido instrumental, número de canales, número de estados.

El 26 y 27 pelea con la cross-covarianza inicio-final para poder graficar la distribución del boundary state, y consigue las elipses de prior y posterior. El resultado es tibio: en IR las elipses son más redondeadas y en MR más alargadas, lo que muestra que el update está más balanceado entre los dos extremos, pero no es una diferencia que convenza a nadie por sí sola. Además la figura se puede hacer también para MR, que no usa el boundary state, lo cual confunde.

El 27 hay un audio de 23 minutos que empieza con "volví a casa y no pude hacer nada, me tiré en el sillón, miré videos, no sé bien qué pasa, estoy atascado". Se destraba reformulando la pregunta: en qué ámbitos MacroIR es mejor que MacroMR, y con qué confianza se puede afirmar eso.

### Semanas del 30 de marzo y 6 de abril
Trabajo de infraestructura: los tipos indexados y los ejes en el DSL, que le permiten correr los diagnósticos para varios modelos y condiciones a la vez. El 4 de abril habilita el algoritmo Taylor porque el gradiente le da bias en algunos parámetros.

### Semana del 6 de abril (11 de abril): el primer hallazgo real
> "Ayer fue un día excepcional para MacroIR porque finalmente pude analizar los datos. Me puse a pensar que en realidad todo dependía del número de transiciones que ocurren en un intervalo."

Tres resultados de una sentada: el ruido instrumental no afecta ninguno de los diagnósticos; el bias depende del número de canales; y la distorsión de la likelihood, sorprendentemente, no depende del número de canales ni del ruido, solo del largo del intervalo de integración. Los parámetros cinéticos y los dimensionales (número de canales, conductancia) siguen el mismo patrón; los instrumentales (ruido, línea de base) no se alteran.

Aparece también la formulación que da la intuición: es como si tuvieras el doble de muestras de las que realmente tenés, pseudo-réplicas.

### Semana del 13 de abril
El 14 se pasa el fin de semana angustiado por una matriz de 600×600×1000 para el bootstrap de la cross-correlación, decide bajarse del bootstrap, y al rato se da cuenta de que la matriz se instancia por sample y que el cálculo es barato. En el medio formula por primera vez la paradoja: "cómo puede ser que mi resolución en parámetros baje al aumentar la resolución en corriente; tengo más información, cómo puedo tener menos información".

El 16, en un audio de 23 minutos, llega el gráfico que estaba buscando: la cross-correlación del residuo estandarizado y de la log-likelihood para distintos lags, por algoritmo. La correlación decae más rápido en Taylor, después IR, después MR, y NR es un desastre. Ese gráfico dice dónde está la eficacia de los algoritmos: en reducir la correlación temporal. En el mismo audio aparece la duda sobre la descomposición de la distorsión en componente de sample y componente de correlación, que le da un número menor que uno y no logra interpretar.

El 17, en varias tandas, resuelve cómo traducir la matriz de distorsión en una corrección de la evidencia usando la aproximación de Laplace, y arma la lista de lo que falta calcular: determinantes de las matrices de distorsión, K efectivos, cross-correlación con lag máximo configurable.

### Semana del 20 de abril: la paradoja confirmada
El 21, en un audio de 14 minutos:

> "Para mi sorpresa absoluta, el error de los parámetros corregido por la covarianza del score me daba constante. Vos tomabas una medición cada un tau o cien, el error es prácticamente el mismo. [...] Es como que hay un mecanismo de compensación perfecto, pero que no pude formular matemáticamente."

Encuentra la analogía que lo hace tolerable: en un proceso de Poisson, saber cómo se distribuyeron las cuentas dentro del intervalo no agrega nada, solo importa el acumulado. Y toma una decisión que le da carácter al paper: contarlo aunque no lo entienda, "honestidad intelectual ante todo", "tenemos que hacer un paper más o menos humano, mostrar que uno no entiende eso me parece que está bueno y llama la atención".

El 21 también decide quedarse en dos estados, con el argumento de que la carga mental de dos constantes cinéticas interactuando disipa el mensaje, y con el argumento oblicuo de que se guarda cartas para "la segunda temporada".

El 22 aparece la crisis de sentido: "estaba pensando si MacroIR fue un fracaso, porque la mejoría respecto de MacroR es pequeña". La respuesta que se da es la que va a sostener el paper: no todos los descubrimientos pueden ser enormes, el trabajo es un trabajo normal y bueno, y el hecho de que quede un problema abierto es en sí mismo interesante en una revista donde los lectores pueden comentar.

### Semana del 27 de abril: micro IR cierra el argumento
El 28 está implementando micro IR con Claude, peleándose con la arquitectura. Encuentra que no puede estimar la FIM con la aproximación gaussiana y que tiene que usar la derivada numérica del gradiente.

El 29, en un audio de 13 minutos, llega el cierre conceptual:

> "Fue una excelente decisión incluir a micro IR en el análisis. Micro IR efectivamente no tiene correlación temporal, las medidas son completamente independientes en el residuo estándar, y no hay distorsión en la matriz de información. Lo cual me indica que los cambios que sí veo en macro son producto de la aproximación normal del espacio de probabilidades."

Es decir: el diagnóstico distingue un algoritmo bueno de uno aproximado, porque en el exacto da la identidad. Y en el mismo audio arma la línea histórica: "esto cierra 20 años de un algoritmo", el paper de 2007 presentaba micro R, macro R y macro NR, y este presenta los mismos tres con la I.

El 30 y el 1 de mayo resuelve la construcción de la Q microscópica combinando canales de a pares y duplicando (la idea tipo pirámide de Pascal, que llama Tartaglia), porque los autovalores se rompen con 50 canales. Y pone el límite del paper: la comparación camino a camino entre el update microscópico y el macroscópico queda para un trabajo posterior.

---

## Etapa 5. Escribir, y las crisis del Hessiano (19 may – 26 jun 2026)

### Semana del 18 de mayo
Abre el grupo 13 el 19 "porque empiezo a escribir el manuscrito". El primer problema es el abstract: el que hay describe MacroIR y no dice nada nuevo respecto de Communications Biology. La reformulación es que este paper aporta la derivación, la génesis y la prueba de validez.

El 20 juega con plantear todo como una adjunción entre sampling y likelihood, con teoría de categorías, y él mismo se lo desarma: sería una herramienta perfecta para que lo rechacen por unanimidad porque nadie conocería todos los conceptos. Queda como anclaje conceptual mencionado al pasar, no como marco.

### Semanas del 25 de mayo y 1 de junio
El 25 nota un hueco argumental del paper anterior: mostró que el boundary state anda bárbaro, pero nunca mostró que el estado sin boundary no anda. El 27 pasa a modo cluster (Dirac) y a paralelizar dentro de cada combinación de condiciones.

El 31, después de un domingo familiar, se destraba de lo que lo tenía frenado hacía días: la Fisher Information Matrix se volvía singular, o peor, indefinida, con autovalores negativos, cosa perfectamente posible para un Hessiano medido fuera del máximo. Sale por dos lados: medir el Hessiano en el óptimo (donde por definición tiene que ser definido positivo) y duplicar el análisis en versión likelihood y versión posterior.

El 2 y 3 de junio la crisis vuelve por otro lado: gran dispersión de la matriz de distorsión con 10.000 canales y intervalo 0.01, condiciones donde no debería pasar nada. Sospecha de un bug de paralelismo. Decide dejar de mirar bootstraps y mirar réplicas individuales para cazar outliers.

### Semana del 8 de junio: el bug
El 9 encuentra el salto: está en la función trust coefficient, el mecanismo de seguridad del simplex que venía de noviembre. Lo reescribe desde primeros principios como un mínimo suave, un log-sum-exp diferenciable, y queda contento con la fórmula.

El 10, en un audio de 12 minutos, cuenta el bug completo, que era más sutil que eso:

> "El salto se transportaba inmediatamente a la matriz de covarianza pero no a la matriz de P_min. [...] Si estabas en una región donde justo coincidía lo esperado con lo encontrado, el D va a ser cero, entonces una pequeña variabilidad te mueve para arriba o para abajo y calculás la derivada en un régimen o en el otro. [...] El error fue que se usaba el mismo alfa para P_min y para P_cov."

Sacada la corrección de la covarianza, la variabilidad desaparece. En el mismo audio hay un párrafo sobre la situación personal (dinero, buscar trabajo) que explica algo de la urgencia de los meses siguientes.

### Semanas del 15 y 22 de junio
El 15 y el 18 agrega la optimización por máxima likelihood en grupos de réplicas, para comparar la covarianza empírica de los parámetros recuperados contra la predicha por el sandwich. Con pocos canales la empírica sale mucho más grande, y muchos parámetros quedan indeterminados, así que tiene que agrupar. Termina eligiendo grupos de 10, 100 y 1000 con 10.000 réplicas.

El 22 hay dos cosas. Primero, la constatación incómoda: no hay casi diferencia entre MacroIR y MacroR en la covarianza de los parámetros; la ventaja de IR está en no tener bias y en estimar bien el error, no en achicarlo. Segundo, el mapa de regímenes que se vuelve la columna vertebral del paper:

> "Me imagino un gráfico con dos regiones: la región multinomial con pocos canales, la región poissoniana con intervalos muy cortos, y después la zona gaussiana, con muchos canales e intervalos no tan cortos. Ahí es donde este algoritmo es ideal."

Son las dos aproximaciones que hace el algoritmo, la distribución de estados por una normal multivariada y la distribución de la corriente media por una normal, cada una con su borde.

El 23, en un audio de 32 minutos, dicta el paper entero de punta a punta: procesos markovianos como herramienta universal, la tensión entre estados instantáneos y mediciones integradas, MacroR y su limitación reconocida desde 2007, MacroIR y el boundary state, los tres estimadores de validez (residuos, score, igualdad de Bartlett), el estudio en dos estados con tres variables, y los resultados. La conclusión es dura con el resto: los no recursivos inflan la covarianza tantas veces como mediciones por intervalo tengas, MacroR se queda en factor 2 o 3, MacroIR llega como mucho a 1.3.

El 26 agrega dos cosas. Que MacroIR sería isomorfo a un Kalman aumentado e integrado (lo deja como sospecha, no como demostración). Y la explicación de por qué MacroMR falla: cuenta dos veces la varianza entre estados finales, una en la likelihood y otra en el estado donde termina. MacroMR pasa a ser, en sus palabras, un método falso, útil para mostrar que la idea simple no funciona.

---

## Etapa 6. El giro final: comparar contra lo que la gente usa (2 jul – 2 ago 2026)

### Semanas del 29 de junio y 6 de julio
El 2 y 3 de julio tiene cinco figuras y cinco suplementarias, y el problema de cómo mostrar la exploración completa del espacio de condiciones. La respuesta es heat maps con contornos, más algunos puntos elegidos a mano para ver las curvas. Decide también definir la distorsión con la Fisher gaussiana en vez de la numérica, porque es más estable y tiene menos error, y volver a correr todo.

El 6 y 8 de julio trabaja la descomposición de la distorsión en componente de sample y componente de correlación, que es la misma historia que cuentan el r² estandarizado y la autocorrelación pero dicha de manera más precisa. Delimita el "cubo crítico": ruido entre 0.05 y 1, entre 10 y 100 canales, todo el rango de intervalos.

### Semana del 6 de julio (11 de julio)
Audio de 28 minutos donde fija el alcance y saca la conclusión más linda del año. El alcance: dos estados, no estacionario, sin micro IR, sin datos experimentales, y todo eso dicho explícitamente en el paper "para maximizar el quantum cognitivo" y para dejarle la puerta abierta a un revisor que pida más.

El hallazgo:

> "Cómo la matriz de Fisher registra el hecho de que vos no extraés más información acerca del número de canales una vez que deja de subir el número de canales abiertos. [...] Dado que sabés cuántos hay ahora, no ganás más información sobre cuántos había originalmente. Eso me voló la cabeza."

Y el resultado práctico que ordena las recomendaciones: aumentar la resolución temporal no mejora la resolución de las constantes cinéticas, pero sí la del número de canales y la conductancia.

### Semana del 20 de julio
El 20 incorpora cuadrados mínimos no lineales (el Levenberg-Marquardt de su propio trabajo de 2007) como término de comparación. La justificación es que es importante ver, para los que usan estos métodos, qué tan confiables son. Él mismo registra que esto le atrasa la escritura y lo embola.

El 23 y 24 tiene el paper prácticamente escrito en la cabeza y las cinco figuras armadas: esquema del método, el problema de la covarianza mal estimada y cómo IR lo arregla, el mecanismo, dónde falla IR, y si vale la pena medir con más resolución (no).

El 24 cuenta que estuvo a punto de dividir el paper en tres y lo volvió a juntar, ahora con LSE adentro:

> "Analizo las estrategias más comunes para tratar con corrientes macroscópicas, y lo que se ve es que básicamente son todas horribles menos MacroIR. Con less squares no tenés bias en los parámetros, pero el error que te da está subestimado hasta cien veces. Si tenés datos con ruido de gating, no podés usar least squares. Es un mensaje fuerte, es un mensaje para eLife."

### Semana del 27 de julio: el giro editorial decisivo
El 28, en el resumen del período, explica por qué se juntó todo:

> "Estaba escribiendo el primer paper, pero me di cuenta de que era invendible para eLife. Porque si ves el uso de MacroR en la literatura es mínimo, no lo usa nadie. Todo el mundo usa least squares. Si vos planteás 'tengo este algoritmo que no usa nadie', no tiene mucho sentido. En cambio si comparo el algoritmo que usan todos con el nuevo y muestro que es mejor, ahí el paper tiene más gancho."

En el mismo audio aparece la figura que cierra la argumentación: un mapa de cinco regiones (el ruido tapa todo; least squares empata; least squares distorsiona pero no hay información extra de corriente; MacroIR gana y recupera conductancia; MacroIR distorsiona), con las condiciones experimentales reales dibujadas encima. Patch escindido, célula entera y oocitos caen todos en la región donde conviene MacroIR, con los oocitos en el borde. Y menciona que extrajo, con ayuda de IA, una implementación de MacroIR usable desde R y Python, contrastada contra la suya.

El 28 dicta también la introducción completa: la independencia de los residuos como supuesto fundamental de la estadística clásica, las cadenas de Markov como manera natural de modelar dependencia temporal manteniendo parámetros universales al registro, la escalera de métodos (least squares, los que agregan varianza de gating, MacroR, MacroIR), y la razón por la que nadie usa los sofisticados: no está caracterizada su validez y no tienen la garantía visual que da superponer predicción y datos.

### Semana del 27 de julio, parte 2 (29 de julio y 1 de agosto)
El 29 el mensaje se termina de ordenar en una sola palabra: autocorrelación. MacroIR resuelve la autocorrelación de procesos modelables con cadenas markovianas. Las tres ventajas quedan enumeradas: más gente puede acceder a la corriente unitaria y el número de canales; con una likelihood sana se pueden comparar modelos; y se puede distinguir variabilidad individual de variabilidad estocástica.

El 1 de agosto corrige el último error: MacroMNR no incluía G_bar_i, la varianza de la conductancia media. Restaurado el código publicado en Communications Biology, lo renombra MacroINR. Y de esa corrección sale el resultado más limpio del paper:

> "La corrección de intervalo te elimina el bias, y la recursión te elimina la inflación de la varianza. Para restaurar la media alcanza con la corrección de intervalo; para restaurar la varianza necesitás la recursión, y la recursión tiene que considerar ambos extremos del intervalo."

En el mismo audio ordena las dos fronteras que quería tocar y que efectivamente tocó: dónde deja de haber ventaja sobre least squares, y dónde MacroIR deja de ser confiable.

### Semana del 3 de agosto (1 y 2 de agosto)
El 1 hace un raconto del paper entero (que es, sin saberlo, el pedido de este documento) y repasa las decisiones que hay que justificar por escrito: por qué dejó afuera el prior, por qué terminó necesitando optimizar cuando el plan era simular y medir sin optimizar, y qué pasa cuando la optimización no converge.

El 2 agrega el octavo algoritmo (cuadrados mínimos sin promediar), verifica que la matriz de distorsión da distinto evaluada en el punto de simulación que en el punto pool para todos los recursivos que no son de intervalo, y confirma que para least squares hay que usar la Fisher numérica porque la gaussiana da mal, cosa que lo sorprendió.

El último audio, con el micrófono abierto por si se le ocurre algo más, deja anotado el pendiente: la anisotropía de la matriz de distorsión de least squares. Si la distorsión fuera isotrópica alcanzaría con corregir el N efectivo por la integral de la autocorrelación; si no lo es, distintas constantes se distorsionan distinto, y eso hay que mostrarlo.

---

## Los seis giros que hicieron el paper

1. **29 sep 2025.** Sin tokens y sin avanzar, decide que los tests de validez son el paper y manda la refundación del software a segundo plano.
2. **6 nov 2025.** Se da cuenta de que puede prescindir del MCMC. Eso saca la mitad del contenido planeado y vuelve el trabajo escribible.
3. **30 dic 2025.** El test del gradiente no distingue MR de IR. La crisis empuja hacia la Fisher Information Matrix y, en enero, hacia la matriz de distorsión.
4. **21 abr 2026.** La resolución de los parámetros no mejora con el sample rate. Decide publicar lo que no entiende.
5. **29 abr 2026.** Micro IR da identidad exacta. El diagnóstico queda validado por construcción.
6. **28 jul 2026.** MacroR no lo usa nadie: hay que comparar contra cuadrados mínimos. Todo se junta en un solo paper y aparece el mapa de regiones de uso.

## Cosas que se abrieron y quedaron para después

Cumulative evidence, el filtro de Bessel (desde el 12 de noviembre), MacroTaylor y MacroTaylorIR, micro IR como paper propio, el régimen estacionario, más de dos estados, los intervalos espaciados exponencialmente, la corrección de la evidencia por la matriz de distorsión, el microscópico estocástico con muestreo de configuraciones, y Luthier.
