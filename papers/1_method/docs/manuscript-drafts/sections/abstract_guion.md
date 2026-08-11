# Guion del abstract, sintetizado desde Intro + Results + Discussion (2026-08-11)

**El trabajo del abstract:** que un Senior Editor de eLife con quince minutos escriba dos palabras,
"valuable/important" en significancia y "compelling" en evidencia, sin tener que abrir el paper.
Presupuesto 200-220 palabras. Cada beat lleva su costo, para que las decisiones sean de suma cero.

---

## La espina: nueve beats

### 1. LA APUESTA. Lo que una aproximación se equivoca sobre la información pasa al mecanismo que elige. 51-53 palabras
Fuente: `05_discussion.tex:228-229`, que es la forma comprimida de `01_introduction.tex:254-263`.
El movimiento tiene dos tiempos y los dos son necesarios: hubo una época en que un mecanismo dejaba
una marca que se señalaba en el registro sin ajustar nada (el flip state, un retardo antes de que
suba la corriente), y las afirmaciones finas de hoy no dejan ninguna y descansan en razones de
evidencia entre esquemas. La consecuencia es la frase tesis de todo el manuscrito:
*whatever an approximation gets wrong about the information in a record passes into the mechanism it selects.*
**Para el editor:** convierte el paper de "barras de error" en "cómo se dirimen los mecanismos".
Es el beat que decide entre valuable e important.

### 2. POR QUÉ NO SE VE. 20-35 palabras
Fuente: `04_results.tex:98-100` y `01_introduction.tex:344-348`.
Un ajuste por cuadrados mínimos predice cada punto desde los parámetros, así que un parámetro
equivocado desplaza toda la traza y se ve; una verosimilitud recursiva se re-ancla en cada muestra y
sigue al registro reporte bien o mal su incertidumbre. Versión barata: "a recursive likelihood
follows the record whether or not the uncertainty it reports is the one it delivers" (18).
**Para el editor:** explica veinte años de no-detección sin acusar a nadie.

### 3. POR QUÉ ES CONTESTABLE AHORA. 18-30 palabras
Fuente: `01_introduction.tex:350-356`, más `339-342` para el hueco.
Ninguna estadística sobre una grabación lo decide (una grabación es un sorteo de la distribución
muestral, y el modelo generador tampoco se conoce). Hace falta preguntarlo desde afuera, sobre un
ensemble con verdad conocida, y para esta clase se puede:
*this process simulates exactly where its likelihood cannot be evaluated exactly.*
**Para el editor:** es la frase que compra "compelling". Sin ella el estudio es "simulaciones".

### 4. LOS INSTRUMENTOS SON CLÁSICOS. 26 palabras
Fuente: `01_introduction.tex:358-361`.
Score con media cero en los parámetros verdaderos, y su covarianza entre grabaciones igual a la
información de Fisher reportada. La glosa "the log-likelihood's gradient" es obligatoria (test de
paráfrasis del brief). **Para el editor y el referee metodológico:** contesta de antemano "¿esto es
nuevo o es la identidad de White?" (regla 6) y compra al lector Mirams/Münch/Del Core.
Ubicado DESPUÉS del beat 3, "the instruments" tiene antecedente y se cierra el seam abierto en el brief.

### 5. QUÉ SE MIDIÓ. 28-45 palabras
Fuente: `01_introduction.tex:371-380`.
Ocho miembros, de cuadrados mínimos sobre la corriente media a un filtro que condiciona cada
promedio de intervalo en sus dos extremos; plano de número de canales, ruido e intervalo; diez mil
grabaciones por celda. Y el escudo: dos estados es **el caso más favorable** para la clausura
gaussiana, así que lo reportado es cota inferior. **Para el editor:** los dos números son evidencia
barata; el escudo desactiva la objeción que apareció en tres rondas externas seguidas.

### 6. RESULTADO 1: los estimadores casi no se distinguen; los reportes sí. 26-33 palabras
Fuente: `04_results.tex:414-419` y `:295-304`, `05_discussion.tex:90`.
Por intervalo cada miembro reporta bien su información; la falla aparece en cómo se acumulan los
intervalos, y la lleva la correlación temporal del score. Un orden de magnitud donde el residuo
conserva memoria (la condición NO es opcional). **Para el editor:** dice que el defecto es
estructural y no un descuido de implementación, y explica por qué es invisible.

### 7. RESULTADO 2: no es un escalar, así que no hay factor de corrección. 25-35 palabras
Fuente: `04_results.tex:592-602` (factor 5-6 con dispersión 3.2, direcciones erradas por veinte y
direcciones casi bien) y `:217-220` (el sandwich, que es matricial, sí lo repara en los miembros que
modelan el gating, y deja corto al de cuadrados mínimos).
**Para el editor:** mata la primera objeción que se le va a ocurrir a él y a cualquier referee
("multiplicá las barras por una constante"). Es además lo más nuestro del abstract.

### 8. RESULTADO 3: un miembro calibrado casi en todo el plano, y sus propias fallas localizadas. 26-40 palabras
Fuente: `04_results.tex:433-435` (93 % / 70 % / 31 % de puntos dentro de ±15 %) y
`05_discussion.tex:298-300` (la desviación es de dos lados: 0.645 a 1.701).
La frase que compra confianza es de `01_introduction.tex:425-426`:
*a method whose limits are unmapped has been endorsed rather than characterised.*
**Para el editor:** es la prueba de que el paper caracteriza y no promociona (regla 5). Un editor que
sospecha promoción de algoritmo se desarma acá.

### 9. RESULTADO 4 y CIERRE. 18 + 20-24 palabras
Ordenamiento, libre de criterio (`04_results.tex:734-737`): ninguna región del plano medido da a la
vez información sobre la corriente unitaria y un intervalo clásico confiable.
Cierre, dos candidatos:
 - **P** (`01_introduction.tex:430-433`, `04_results.tex:526-540`): hacer la medición requiere un
   ensemble, usarla no; la autocorrelación integrada del residuo estandarizado cuenta los intervalos
   efectivos y sale gratis. Concreto, accionable. **Trampa:** lee el RÉGIMEN, no certifica
   (`04_results.tex:572-574`: NR se equivoca hasta 2217× con residuo indistinguible de blanco). El
   abstract no puede insinuar certificado.
 - **AA** (`05_discussion.tex:294-296`): lo que el lector no puede medir por su cuenta, y es lo que
   se entrega acá, es si el intervalo reportado es el intervalo que se entrega. Más abstracto, más
   seguro, y es la respuesta directa al OPEN del brief sobre la frase de significancia.

### 10. LA SEGUNDA MITAD DEL ENCUADRE VIEJO: entregar la maquinaria. 25-35 palabras
Luciano, 2026-08-11, recuperando el frame de `00_abstract.tex` seccion 5: la familia recursiva tiene
veinte anios y no se usa, y hay DOS explicaciones plausibles (no objeciones publicadas). Una, nadie
sabe cuando fallan, y la respuesta es el mapa, que son los beats 7 a 9. Dos, son dificiles de
implementar, y la respuesta es el codigo en R y Python. La extension nueva es que no se entrega solo
el MLE: se entrega **el aparato de medir bias y matriz de distorsion en el modelo y el punto que el
lector quiera**, de modo que un usuario ajuste con Levenberg-Marquardt y despues estudie si ahi
distorsiona.

**LA MAQUINARIA EXISTE DESDE 2026-08-11, macroir `911828d`.** El inventario de esa manana decia que
estaban las cuatro piezas (esquema arbitrario, simulador exacto con semilla, `loglik` con gradiente,
ajuste Levenberg-Marquardt) y NO el lazo que las junta. El lazo se escribio ese mismo dia:
`core/include/macroir/distortion.hpp`, expuesto como `mi.distortion(...)` en Python y
`macroir_distortion(...)` en R, con la misma firma (esquema, theta, experimento, n_grabaciones,
semilla) y devolviendo la matriz de distorsion, sus autovalores, y el sesgo de primer orden.
**El sesgo NO cuesta ajustar**: como la covarianza del score esta centrada, la media del score sale
aparte y el sesgo es G^(-1)·s̄, una pasada por grabacion. Medido: 0,52 ms por grabacion en la celda
del fixture.

Verificado antes de que la frase se pueda escribir: once suites en siete segundos, el checkpoint
contra el binario congelado en 8 celdas y 0 off, y los tres lenguajes coincidiendo bit a bit en la
diagonal, los autovalores y el sesgo. Contra macro_dr, con las diez mil grabaciones de la bateria,
las diagonales dan 0,900 / 1,276 / 0,867 / 0,937 / 0,957 / 0,919 contra 0,907 / 1,252 / 0,864 /
0,940 / 0,916 / 0,923.

**LO QUE LA FRASE TODAVIA NO PUEDE DECIR.** macroir implementa el miembro calibrado y nada mas
(README: "The R and NR variants of the filter, and LSE, are not written"), asi que el lector audita
SU punto con ESE miembro. Comparar peldanos en su propio punto necesita R y NR como enum en runtime,
que es el item 2 de `macroir/docs/next.md`.

**DOS ESCALONES, Y NO SE PUEDEN CONFUNDIR.** El barato lee el REGIMEN desde una sola grabacion y no
certifica: `04_results.tex:572-574`, NR se equivoca hasta 2217x con residuo indistinguible de blanco.
El caro simula en el punto ajustado y mide. El abstract promete uno solo.

**LIMITE DEL ESCALON CARO**, ya escrito en `05_discussion.tex:224`: simular desde el modelo ajustado
sigue simulando desde el modelo que uno queria verificar, asi que contesta si la aproximacion
distorsiona BAJO ESE MODELO, no si el modelo es correcto.

**DOBLE TRABAJO.** Es tambien la respuesta a la objecion de dos estados, que aparecio en tres rondas
externas seguidas y hoy vive solo en `05_discussion.tex:430`: la biblioteca es mas general que el
estudio, el lector corre el diagnostico en su esquema en vez de extrapolar el mapa.

**REGLA 5.** Lo que se corto en 2026-08-07 (nota 1f(iv)) fue terminar en un producto con nombre
despues de decir que la familia esta rota, que arma el arco "aca esta el mio". Este encuadre no lo
arma, porque lo entregado es el aparato de auditoria y no el algoritmo. Dejar el nombre MacroIR
afuera del abstract y en Code Availability. Y no nombrar "rmacroir" ni "pymacroir": no existen, el
paquete se llama `macroir` en los dos lenguajes.

**COMPITE CON EL CIERRE P DEL BEAT 9**, no se suma. Los dos son "que hace el lector"; juntos son
cincuenta y cinco palabras.

---

## Aritmética

Núcleo obligatorio (beats 1, 3, 4, 5, 6, y el ordenamiento del 9): ~170 palabras. Quedan ~50 y hay
cuatro candidatos que piden 25-35 cada uno: el beat 2 (por qué no se ve), el 7 (no es escalar), el 8
(miembro calibrado y auditado) y el par de cierres que se excluyen, P (gratis, una grabación) contra
10 (la maquinaria en el punto del lector). **Entran dos.** Esa es la única decisión real del abstract.

Mi orden: 8 primero, porque sin él el paper no tiene resultado positivo y el mapa no existe, y porque
ahí va la frase que desarma la sospecha de promoción (*a method whose limits are unmapped has been
endorsed rather than characterised*). Después el cierre, y entre los dos elijo **10 sobre P**: hace
doble trabajo (significancia y defensa del alcance de dos estados), mientras que P es más barato pero
más débil, porque lee el régimen y no certifica nada. Quedan afuera el 7 y el 2, y el 7 es el primero
que vuelve si el techo sube de 220.

**CORRECCIÓN, 2026-08-11, después de armarlo de verdad.** La estimación de 235-240 que estaba acá
era baja. Ensamblado beat por beat, con gancho, escudo de cota inferior, instrumentos clásicos y
entrega, el abstract sale en **250** (`tmp/v11_full.txt`). El techo de 220 se alcanza
(`tmp/v11_tight.txt`) y cuesta exactamente dos cosas: el escudo de dos estados (18 palabras) y la
localización de la falla del miembro calibrado ("toward few channels and low noise", 6). El resto es
elisión sin pérdida. Esas treinta palabras son la decisión, y no hay una tercera opción de redacción
que las evite.

---

## El elevator speech, y lo que enseñó

Sintetizado desde los beats el 2026-08-11 y conservado porque el ejercicio cambió el presupuesto.

**Sesenta segundos.** When we said P2X2 passes through a flip state before it opens, you could point
at the delay in the trace. Anyone could see it. Twenty years later the claims we make about the same
receptor come out of model comparison, and there is nothing to point at. So when two analyses of the
same recordings put different mechanisms first, and that happened to us, the difference is in the
likelihood and not in the data. The trouble is you cannot check a likelihood on a real recording. You
would need to know the answer. This is one of the rare problems where you can simulate the process
exactly even though you cannot evaluate its likelihood exactly, so we ran the classical test: at the
true parameters the score should average to zero, and its spread should equal the information the
likelihood claims. Eight approximations, over channel number, noise and sampling interval, ten
thousand recordings per condition. The estimates barely move. The error bars are off by an order of
magnitude, and no correction factor repairs them, because the distortion has directions. One of the
eight is honest over almost the whole plane, and we mapped where it is not, its own corner included.
And the machinery ships in R and Python, so you run the same test on your own model at your own
operating point instead of trusting our map.

**Quince segundos.** A likelihood that misreports its own uncertainty does not just widen an error
bar, it changes which mechanism wins. We measured how far eight of them are off, against a process
that can be simulated exactly, and mapped where each can be trusted.

**TRES COSAS QUE ENSEÑÓ.**
1. Los beats 6 y 7 se fusionan solos al decirlos en voz alta: "the error bars are off by an order of
   magnitude, and no correction factor repairs them, because the distortion has directions", 20
   palabras, contra 56 presupuestadas por separado. El beat 7 dejó de ser el primer candidato a
   quedar afuera: entra pegado al 6.
2. El beat 2 no existe como beat hablado, aparece como subordinada dentro del 3 ("you cannot check a
   likelihood on a real recording, you would need to know the answer"). Confirma que era el más
   prescindible como unidad, y que la mitad que importa cuesta seis palabras.
3. Lo que se cae solo del discurso es la atribución clásica del test: hablando nadie pregunta si la
   identidad es de White, leyendo el referee metodológico sí. Se queda por razón defensiva, no porque
   sea parte del argumento natural, y por lo tanto es candidata a salir ANTES que cualquier resultado
   si alguna vez hay que cortar.

## Deliberadamente afuera, y de dónde vienen si alguien los pide

- El 2217× con residuo blanco (`04_results.tex:572-574`). Es la advertencia sobre el cierre P y
  necesita treinta palabras para no leerse como que el método propio también falla.
- Las cinco regiones del mapa y los pisos de alcanzabilidad, 17 y 52 canales (`04_results.tex:742-744`).
- Que el número de canales y el ruido entran sólo por su cociente (`05_discussion.tex:92`). Es
  elegante y es lo que hace transferible el diagnóstico, pero cuesta veinte palabras.
- La concesión: cuadrados mínimos recupera las tasas con sesgo cero, lo que erra es cuán bien dice
  que las recuperó (`04_results.tex:410-412`). Trece palabras, y es el mejor candidato a volver.
- El puente a la evidencia bayesiana (volumen + tamaño muestral efectivo, `05_discussion.tex:235`).
  Es lo que conecta el beat 1 con lo medido, y el paper NO lo calcula.
- Dónde la ruta clásica es fuerte: muchos sweeps intercambiables, NSFA (`05_discussion.tex:226`).
- El costo computacional y la comparación con Münch (`05_discussion.tex:138`, `:183`).
