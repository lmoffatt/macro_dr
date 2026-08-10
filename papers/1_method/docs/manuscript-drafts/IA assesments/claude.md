alisis critico. Como lo ves para eLife?

Leí el paper completo. Verifiqué la álgebra de las Ecs. 11–19 y A4 (la covarianza del par de frontera, el colapso a $K\times K$, la identidad MR/IR de varianza total): todo cierra. La crítica que sigue no es sobre corrección sino sobre alcance, fuerza de evidencia y presentación.

Veredicto corto sobre eLife
Vale la pena, pero tal como está apostaría a un assessment de "valuable" + "solid", no de "important" + "convincing". Y en el modelo actual de eLife eso queda publicado y permanente. El riesgo real no es el rechazo (ya no existe después de revisión), es el desk reject y la evaluación tibia. Hay precedente directo a favor: Münch et al. 2022 es eLife y está citado, así que el editor no puede argumentar que el tema no pertenece.

Dos agregados concretos, que creo que están a tu alcance con la maquinaria que ya tenés, mueven eso a "important/convincing". Los pongo primero porque son la decisión.

Los cuatro golpes que vas a recibir
1. K=2 y nada más. Decís que el esquema de dos estados es "the most favourable case for a Gaussian treatment of the occupancy and therefore a lower bound on the error for anything richer". Eso está afirmado, no medido. Para el cierre de ocupancia el argumento es plausible (más estados, menos canales por celda, peor TCL). Para el cierre de señal-intervalo no es obvio en absoluto: con conductancias agregadas y $\mathbf{Q}$ no lumpable, $A(\Delta)$ tiene una ley distinta y el reuso de la marginal final como prior del intervalo siguiente (que vos mismo admitís que no es exacto) puede degradarse de otra forma. Y el golpe es evitable: en Data availability decís que el core fue chequeado "on eight cells, covering four kinetic schemes". Un referee lee eso y pregunta por qué no corriste ninguno. Un esquema K=3 o 4 no lumpable, en una subgrilla rala (3 conteos de canal × 3 niveles de ruido × 2 intervalos), convierte una extrapolación en una medición. Es el ítem de mayor retorno de toda la lista.

2. $P_{open}=0.5$ fijo. Vos mismo decís que por eso no podés nombrar un umbral en canales, porque conteo de canales y conteo de canales abiertos difieren por un factor fijo de dos. Pero el problema es peor que una unidad ambigua: el eje $N_{ch}$ de la Figura 6 es en realidad $N_{open}$ al pico, y la Figura 6 es lo que le vendés al experimentador. Un receptor con $P_{open}=0.05$ está en otro lugar del mapa y el lector no tiene cómo trasladarse. Barrer tres valores de $P_{open}$ arregla la unidad y hace transferible el mapa.

3. El kernel Bessel. Lo tratás honestamente en Discussion, pero la honestidad no salva la aplicabilidad. La ventaja estructural de IR es exactamente el confinamiento del peso dentro de un intervalo, y el Bessel lo rompe justo ahí. Decir "no lo medimos" sobre la premisa central de la recomendación es la clase de hueco que un referee convierte en revisión mayor. Simular datos con kernel Bessel y puntuarlos con la verosimilitud de ventana uniforme, en dos o tres celdas, acota el costo. Si la distorsión se queda cerca de uno, desactivás la objeción entera con una figura suplementaria; si explota, es mejor saberlo vos.

4. El brazo clásico es débil, y la defensa ya está en tus datos. Vas a recibir: "nadie usa mínimos cuadrados homocedástico, la reparación estándar es WLS con varianza siguiendo la media". El paper ya reconoce que el factor 15 mezcla heterocedasticidad con correlación. Lo que no hacés, y deberías, es señalar que NR es esencialmente esa reparación heterocedástica y no ayuda: 13.63 contra 14.96. Eso está medido, es una respuesta de una frase, y sin ella el titular "no hay región donde las fluctuaciones informen y el intervalo clásico sea confiable" queda expuesto.

Huecos evitables
No reportás ningún costo. El paper se organiza alrededor de una "cost ladder", contesta a Del Core y Mirams que argumentan contra el filtrado por costo, y dice "the cost is real, and what it buys is now quantified" sin cuantificar el costo en ninguna parte. Una tabla de tiempo de pared por evaluación de verosimilitud, por miembro, a K=2 y K=5, cierra el pedido antes de que lo hagan.

No invocás la propiedad de martingala. Para una factorización prequencial correctamente especificada, los incrementos del score son una diferencia de martingala y los términos cruzados se anulan idénticamente. Toda tu Figura 3 (identidad por intervalo que se sostiene, identidad acumulada que falla) es esa afirmación. Una oración lo convierte de observación empírica en consecuencia teórica.

No distinguís dos fallas que tu descomposición mezcla. En NR la factorización no es prequencial: es una verosimilitud compuesta de bloques independientes, y su distorsión de información es la ineficiencia conocida de Godambe, esperada y calculable, no un defecto. En R e IR la falla residual es error de cierre. Son animales distintos y ambos aparecen como $\mathbf{R}\neq\mathbf{I}$. Un estadístico lo va a marcar. Enmarcarlo vos es más fuerte: "medimos el tamaño de la ineficiencia composite en este dominio".

Definí $\mathbf{H}$ sin ambigüedad para los miembros recursivos. La Ec. 24 lleva $\mathbb{E}[\cdot]$, pero $\sigma_t^2$ y $\nabla\mu_t$ de un filtro dependen de los datos. Métodos no dice si se evalúa por grabación y se promedia entre réplicas o de otra forma. Es una línea y evita una pregunta incómoda.

La Introducción promete más de lo que los Resultados sostienen. Decís que el conteo de canales y el ruido entran en la respuesta "through that one statistic and through nothing else", y después la Fig. 4–supp 3 muestra que la distorsión es anisotrópica y que ningún reescalado escalar sirve. Las dos frases conviven mal. Ajustá la Introducción a lo que medís: el número dice cuánto está mal el intervalo, no cómo corregirlo.

El repositorio se confiesa en Métodos. "The tables hold stale values the figures never used", "the production scripts do pass a number_of_substeps argument, and on this branch it has no effect", la semilla resuelta nunca registrada. Lo primero y lo segundo se arreglan en el repo, no se documentan. Lo de la semilla es un error no forzado en un paper cuyo tema es el rigor metodológico; es defendible por equivalencia estadística, pero preferiría no tener que defenderlo.

Pendientes de las fronteras sin incertidumbre. 1.01 y 1.03 sobre cuatro cruces cada uno, y el titular es "ruido proporcional al número de canales". Poné un intervalo sobre la pendiente. Y corré las diez celdas que fijan el vértice inferior de la Figura 6: decís explícitamente que no se corrieron, y son baratas contra todo lo demás que ya corriste.

Sacá VR y MR del cuerpo. Vos mismo decís que el resultado de VR es de suplemento y que el cuerpo no depende de él. Con seis columnas en vez de ocho, las Figuras 2 y 3 respiran y no se pierde nada central. La Figura 4 sigue siendo el problema: dos mitades × seis miembros × dos parámetros × cuatro filas, con dos códigos de línea y dos de trazo. Es ilegible en pantalla. Necesita un panel resumen de una línea por miembro al lado, o partirse.

"Seven-rung ladder" en el abstract contra ocho miembros en el cuerpo. Menor, pero es lo primero que lee el editor.

Justificá la exclusión de MRT/IRT en una oración. Existen en el código, el marco es una escalera, y no decir por qué no están invita la pregunta.

Sobre el encuadre para eLife
Tu reclamo de novedad más fuerte es la medición, no el método: IR ya está publicado y coincide con un Kalman de medición integrada de 1988 a $10^{-8}$. Lo concedés limpiamente, y eso te honra, pero también es lo que baja el techo del assessment. La carta de presentación debería vender la lección general, no el nicho: cualquier modelo de estado con observaciones promediadas en ventana, cualquier verosimilitud aproximada para un proceso que se puede simular exacto pero no escribir. Esa frase está en tu Discussion y merece estar en el título o el abstract.

El título actual es preciso y plano, y "distort" hace trabajo técnico que el lector todavía no tiene. Consideraría uno que diga el hallazgo accionable: que no existe grabación que dé a la vez la corriente unitaria y una barra de error clásica confiable.

Si no querés correr más simulaciones
Entonces eLife no es la mejor apuesta por valor esperado. Biophysical Journal es la casa de Milescu 2005, Celentano y Hawkes 2004, Moffatt 2007 y Sigworth 1981, con revisores que entienden el argumento sin traducción. PLoS Computational Biology está pegado a IonBench. Phil Trans R Soc A es donde vive la comunidad de discrepancia y calibración (Lei 2020, Del Core y Mirams 2025) y donde la crítica al brazo clásico se recibiría como conversación en vez de como ataque. Cualquiera de los tres da un resultado más predecible.

Mi recomendación: los ítems 1 y 2 (un esquema K≥3 y un barrido de $P_{open}$), más el Bessel acotado en dos celdas y la tabla de costos, y ahí sí eLife con expectativa razonable de "important". Sin ellos, mandalo a Biophysical Journal.