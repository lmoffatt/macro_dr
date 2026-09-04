# Plan de adaptación a JGP — CONSOLIDADO 2026-09-03

Decisión de venue: JGP primera opción (fees $0 + format-neutral + scope resuelto con precedente de
género + linaje flip/JGP + "specialized community" del desk de eLife 2026-09-03). BJ = fallback
(plan en ../bj-draft/, compuerta = waiver ANTES de enviar). Tipo: **Article** (no M&A: 16/18 M&A
2019-2026 son técnica húmeda). Fuentes: HTML guardados en esta carpeta (Submission Guidelines,
Fees, Editorial Policies) + snapshots About 2023/2026 + 2 workflows (censo 718 registros PMC;
medición del repo) + 5 ejemplares medidos por texto completo.

## 1. Verificado que gobierna (todo resuelto; cita textual en los HTML de esta carpeta)

| ítem | JGP | nosotros |
|---|---|---|
| costo | Green $0, sin page/color charges; "inability to pay... will not affect publication". Gold $2.000, waived para low/middle income (¿Argentina? lista RUP pendiente) | Green por defecto; preprint cubre el embargo de 12 meses |
| 1ª submission | format-neutral: PDF único, páginas numeradas ✓ (cls ya numera), leyenda bajo su figura | PDF actual casi sirve |
| límites | "No limits are imposed on the number of words or figures in Articles" | 66 pág OK; rango medido: RyR 2026 = 15.1k palabras + 30 figs + 3 apéndices; Benndorf 2023 (género exacto) = 8.4k + **22 figs** |
| supplemental | desalentado ("thus, it should be possible to read and understand the article without consulting supplemental materials"); sólo blots/datasets/videos(≤10)/lengthy derivations; párrafo-resumen al final de Methods | figsupps → PROMOVER a cuerpo; apéndices quedan |
| citas | Harvard autor-año alfabético, sin tope | **ya somos autor-año alfabético** (bbl verificado); formato fino de lista → revisión |
| abstract | ≤250, sin citas | 188, sin citas ✓ |
| título | <100 chars | 93 ✓ |
| running title | <50 chars, requerido | ESCRIBIR |
| Summary eTOC | ~40 palabras, 3ª persona (existe como elemento en el XML de producción) | ESCRIBIR |
| orden | Intro→Materials and methods→Results→Discussion (aplica al formatear; 1ª submission libre) | mover Methods (ver A3) |
| estadística | p exactos, N, centro/dispersión, software+versión | esencia cumplida |
| IA | documentar en Methods contenido generado con IA (correctores exentos) | DECISIÓN L. (texto) |
| código | público obligatorio ("publicly available database or as supplemental") | Zenodo snapshot desacopla del repo (ver §3.8) |
| aceptación | DOCX + figuras editables individuales + ORCID + CRediT | costo diferido (pandoc+math) |
| edición | online = registro, publicación diaria; impresa existe (ISSN 0022-1295, números a pedido) | figuras pensadas para 85 mm (revisión) |

## 2. Etapa A — hasta la primera submission (lo único bloqueante)

A1. **Tabla de promoción figura por figura** (la preparo, L. veta): los 12 figsupp con su mención
    en prosa (23 menciones: 16 Results, 6 Discussion, 1 Methods), qué argumento carga cada uno,
    voto promover/suplementario. Referencia: hasta ~19 figs de cuerpo < Benndorf 22 < RyR 30.
A2. Ejecutar promoción: \figsupp → entornos figure normales tras su figura madre; renumeración
    automática por \ref; reescribir las 23 menciones compuestas "Figure~\ref{...}--figure
    supplement N" → \ref simple; leyendas quedan bajo cada figura. Absorber
    supplementary_file_1 (1.960 palabras + 3 tablas → cuerpo o Tables S; decidir en A1).
A3. Methods tras la Introducción: mover \input, leer costuras (refs adelante→atrás), párrafo-
    resumen de supplemental al final de Methods si queda algo suplementario.
A4. Neutralización cosmética: fuera eLifeMediumGrey (cls+bst), abstract env neutro, [lineno] ON,
    pie "page X of Y" queda en negro. Sin tocar contenido.
A5. Textos nuevos: (i) running title <50; (ii) Summary ~40 palabras; (iii) carta JGP: conceptual
    advance + precedente ("similar in kind to Benndorf & Schulz 2023") + pisar "theoretical
    research... grounded on established experimental evidence" y "helps design new experiments"
    (mapa Fig 7) + preprint + COI none + datos accesibles para revisores + editor sugerido +
    reviewers (Qin, Benndorf, Mirams, Sivilotti, Plested, Kinz-Thompson; Milescu a criterio) +
    exclusiones (¿ninguna?); (iv) Data availability formato JGP (arriba de Acknowledgments) con
    DOIs; (v) acknowledgment del preprint + funding con nombres exactos de agencia (text mining).
A6. Checks de build: pdftotext | grep "??"; swallowed_sentences.py; figuras citadas en orden
    numérico de primera mención (requisito JGP; extender check.sh item 11); leyendas bajo figura;
    páginas numeradas.
A7. Form EM: ORCID, preprint DOI, tipo Article, editor sugerido, reviewers, sin related papers.
A8. Verificaciones externas — CORREGIDO 2026-09-03: los depósitos Zenodo YA EXISTEN con DOI
    (código 10.5281/zenodo.22167744, datos 22168409, macroir 22168263, P2X2-2025 17162475) y la
    sección "Data and code availability" del backmatter ya los cita; el error "sin DOI" fue mío
    (warnings de bibtex por tipo de entrada + captura vieja del form). Queda: bioRxiv (subir, con
    el PDF actual sirve), ORCID en el sistema, lista países RUP (¿Gold gratis?), y para JGP sólo
    RE-UBICAR la sección de datos arriba de Acknowledgments con título "Data availability".

### Estado de ejecución (2026-09-03)
- PASO 0 ✓ línea de base: 66 pág, 0 "??" (los 7 iniciales eran bbl sin regenerar), 46 "figure supplement".
- PASO 1 (A2) ✓: 11 promovidas (plane-s1..4 y memory-s1..2 en fullwidth a tamaño de diseño,
  deshaciendo el 0.81 que elife.cls les imponía; resto a ancho de texto); Fig3-S3 eliminada y su
  cláusula única ("noise level keeps a white score") caída CON MARCADOR % [L-REVIEW]; 22 menciones
  reescritas; oración nueva del colapso insertada CON MARCADOR; compila 65 pág, 0 "??",
  0 "figure supplement", Figs 1-18 completas. Diff: +131/−35 en 3 archivos.
- ABIERTO (etapa B, no bloquea la 1ª submission format-neutral): orden de primera mención ≠
  numeración en 10 casos (secuencia real 1,2,4,7,3,8,5,6,9,17,10,13,16,12,14,11,15,18); la prosa
  cruza referencias entre secciones, herencia eLife donde el orden no importaba. Resolver en
  revisión (renumerar/reubicar o tocar prosa) si JGP lo exige.
- PASO 2 (A3) ✓ 2026-09-03: Methods movido entre framework y Results (marcador [L-REVIEW] en el
  main); costuras limpias (framework cierra en unidades de distorsión; Methods abre "Minimal model
  and simulation"; Results abre "A recursive likelihood cannot be checked by eye"). HALLAZGO: la
  declaración de IA YA EXISTÍA al final de Methods (decisión #5 CERRADA; cumple la política JGP tal
  cual). Main renombrado: elife_paper.tex → jgp_biorxiv_moffatt.tex (pedido L.); elife_paper.pdf
  restaurado a HEAD = artefacto enviado a eLife, intocado.
- xelatex devuelve rc=1 sin "!" en el log desde la línea de base; PDF completo; vigilar.

## 3. Etapa B — si invitan a revisar
Referencias a formato JGP fino (et al.>10, abreviaturas, Reference Guidelines); figuras a 85/180 mm
y ≤1 página c/u; estadística en leyendas; SciScore (Methods al sistema); punto por punto.

## 4. Etapa C — aceptación
DOCX final (pandoc + pasada manual de math), figuras editables individuales, CRediT, license,
proofs en 48 h; tapa opcional (300 dpi, 8.75×11.25 in).

## 5. Decisiones de L. pendientes
1. Veto de la tabla A1 (bloquea A2).
2. Green vs Gold (tras lista de países RUP).
3. Editor sugerido (mirar rupress.org/jgp/pages/editors-and-staff; Faraldo-Gómez fue el AE
   computacional; Eisner EIC hasta dato en contrario).
4. ¿Consulta presubmission opcional a jgp@rockefeller.edu? (ya no bloqueante; puede citar el
   precedente).
5. ~~Texto de la declaración de IA en Methods~~ RESUELTA: ya existía al final de Methods, conforme a la política JGP.
6. Desk de eLife: NO se declara (sólo "related or competing papers in press or consideration
   elsewhere"); el preprint SÍ.

## 6. Esfuerzo estimado
A1 hoy (tabla para veto). A2-A6: 1-2 jornadas con vetos resueltos, compilando por paso y mirando
el PDF. Etapas B/C sólo con buenas noticias.
