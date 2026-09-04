# Síntesis de las dos revisiones (orden de secciones + economía de figuras) — 2026-09-03

Fuentes: workflow section-order (4 lentes: experto JGP, generalista, cazador de artefactos,
mapeador; salida completa tasks/wzwctp30l.output) + workflow figure-reorg (6 familias;
tasks/wqclyx0hi.output; digest en figure_reorg_findings.md). Veredicto global: el orden
intro→teoría→methods→results→discussion es sano y ninguna figura sobra; el daño es local y
enumerable, casi todo artefacto de las dos mudanzas de hoy.

## Arreglado ya (2026-09-03)
- BLOCKER: 06_methods.tex:547 "Figures~4 and~6" hardcodeado (numeración vieja) → \ref a clouds +
  plane--memory-s2, con [L-REVIEW] (pendiente: config de la columna ILSE de la familia temporal).

## Lote A — texto/orden, sin decisión de fondo (ejecuto en la próxima pasada)
1. Fósiles de taxonomía eLife en 6 captions: "the body figure" (Figs 3, 8, 9), "the preceding
   supplement" (Fig 6), "the body's twenty" (Fig 9), etc. → número de figura explícito.
2. Ídem en Methods: "the supplements that restrict themselves to converged fits say so" (599-601)
   → nombrar figuras concretas o eliminar el contraste; frases forward-leaning (543-546, 684-729)
   → "as the Results show".
3. Párrafo-resumen del material suplementario (requisito JGP) al CIERRE de Methods. La
   declaración de IA SE QUEDA en Methods (política JGP verificada: "documented in the Materials
   and methods section" — los agentes propusieron moverla a end matter y están equivocados);
   queda antes del párrafo-resumen.
4. Resultados, orden: mover R5 ("Sampling a hundred times faster", lee Fig 11) a continuación de
   R3 (bloque del plano), antes del trío τ9. Discharge de: figura 11 citada tras 12-14.
5. Familia limit: reordenar trazos → esquina → ley (cero R; mover dos entornos + ~3 oraciones del
   opener). Con eso 15/16/17 quedan en el orden en que la prosa los consume.
6. Memoria: mover la oración introductoria de Fig 12 al frente de su párrafo (12 antes que 13);
   en mi oración nueva del colapso, quitar el paréntesis a fig:limit-s1 y decirlo en palabras
   (evita otra inversión).
7. Recablear punteros forward restantes: Methods cita Figs 4/5/6 antes de que exista la 2 →
   nombrar corridas, no figuras ("the time-resolved run"); R1 cita 7A temprano; R2 cita 8 antes
   de 5-6; R3 cita 17. Tras esto + (4) + (5), el residuo de inversiones es chico; renumeración
   total NO (dispersaría los grupos padre-hijo).
8. Discusión: subir "The calibrated member and the control select different mechanisms" (el pago
   P2X2) a después de los dos párrafos de apertura; separar el bloque prior-art (Kalman/Münch/
   costo/NSFA) del heading "Grouping samples" con heading propio.
9. Costura Methods→Results: recortar redefiniciones duplicadas (unidades, celda, dos brazos).
10. Oración ilegible de números cabecera (p.24: "five per cent... twenty-nine and seventy") →
    conteos explícitos por miembro.
11. Backmatter: retitular "Data availability" (JGP) — ya estaba en el plan.

## Lote B — figuras (R; para firma de L., costos por figura en figure_reorg_findings.md)
1. FUSIONAR time-s1+s2 en una figura de 4 bandas ~7.0×8.8 in (la pregunta pendiente de L. queda
   respondida con evidencia: NUNCA se citan por separado, y el contraste del washout hoy queda
   partido entre los dos PDFs). figure_3_S1_S2.Rmd ya dibuja ambas en un loop.
2. En esa misma figura: columna rotulada "LSE" dibuja el digest de ILSE (DTOK figure_3.Rmd:143)
   → RELABEL a ILSE. Cuestión de verdad, no cosmética. VERIFICAR el DTOK antes.
3. QQ (Fig 3): terminación de promoción: quitar título horneado, recortar ~0.9 in de lienzo
   muerto, y reconciliar "noise 0.05" (etiqueta ×10) vs S̃=0.005 de la caption. Opcional: small
   multiples por miembro (el propio Rmd lo pre-planifica).
4. Plano padre (Fig 7): apagar el overlay de contornos del SE corregido (se_lines=FALSE) ahora
   que el campo es la Fig 11 de primera clase; una oración en la caption lo remite.
5. Vestigios: tag 'A' huérfano en figure_2.Rmd:524; "N_ch" literal en la key de figure_7.Rmd:1243.
6. Opcionales con firma: recortar/agrisar el outlier NR de Fig 14; NO tocar fila C del padre
   temporal (lean keep).
7. En revisión (no ahora): re-render familia plano a medida JGP (fuentes primero: el piso de 7 pt
   se presupuestó con el estirón eLife 6.1%; en JGP es ~1.2%).

## Los dos forks para L.
FORK 1 — ¿18 figuras o democión parcial? El generalista pide demover 8 (5,6,8,9,10,13,14,17): el
tramo 7-11 son cinco mapas de plano casi idénticos seguidos y ahogan el skim. Los seis agentes de
familia + la filosofía escrita de JGP dicen keep. RECOMENDACIÓN: mantener 18 y atacar el skim con
A.4/A.5 (orden), B.4 (diferencia 7 de 11), A.1 (captions), y la fusión B.1 (18→17). Camino medio
si L. prefiere: las marginales identificadas son plane-s3 (Fig 10, una sola cita) y el QQ (Fig 3).


FORK 2 — ¿Framework? Experto: se queda, renombrar "Theory" (JGP-nativo). Generalista: fundirlo en
Methods y mandar la construcción Eq.2/3 al Apéndice 1 (que ya la duplica). RECOMENDACIÓN: renombrar
"Theory" ya; mover a apéndice las dos salvaguardas numéricas de Methods (~1 página) y evaluar la
mudanza de Eq.2/3 con calma; la fusión completa es cirugía mayor y no la pide la venue.


## EJECUTADO 2026-09-03 (Lote A completo + fork 2 parcial)
- Fósiles de taxonomía eLife: 0 en el PDF ("body figure", "preceding supplement", "figure supplement",
  "Supplementary File" todos en cero; 15+ sitios reescritos, captions y Methods).
- Párrafo "Online supplemental material" agregado al cierre de Methods; documento suplementario
  retitulado "Supplemental material"; backmatter → "Data availability".
- R5 movida al bloque del plano (su cierre ahora introduce τ9); familia limit reordenada
  trazos→esquina→ley con opener reescrito [L-REVIEW]; Fig 12 introducida antes que 13; P2X2 subida
  en la Discusión + 2 headings nuevos [L-REVIEW]; Theory renombrada [L-REVIEW]; costura
  Methods→Results recortada; oración de números cabecera reescrita.
- ORDEN DE PRIMERA MENCIÓN: 1..18 EXACTO, cero violaciones (recableo de ~10 punteros forward a
  palabras + cita nueva de plane-s1 en su sitio del plano). El ítem de etapa B quedó cumplido HOY.
- SALTEADO CON CAUSA: mover salvaguardas numéricas a apéndice — el fuente registra decisión
  vinculante previa ("every DISCLOSURE stays here", 06_methods comment, ajuste 2026-08-30);
  ídem números del differenced-Fisher (son la divulgación; la tabla completa ya está en S1.5).
  La sugerencia aceptada quedó sin efecto por conflicto con decisión registrada.
- Build: 65 páginas, 0 "??". ~9 sitios [L-REVIEW] pendientes de lectura de L.
- PENDIENTE Lote B (R): fusión time-s1+s2, relabel ILSE, terminación QQ, se_lines off, vestigios.
