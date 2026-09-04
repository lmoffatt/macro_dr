# Plan de adaptación a Biophysical Journal — 2026-09-03

Fuentes: HTML vivo guardado por L. 2026-09-03 (`Information for authors_ Biophysical Journal.html`,
completo, 98k chars de texto) + workflow wf_908d97df (Wayback 2025-07-08 coincide; 4 ejemplares PMC;
medición del repo). Extracto verbatim: `tmp/bj/authors_page_extract.md`. Salida completa del
workflow: `~/.claude/.../tool-results/b1qz4tsyf.txt`. OJO: `TBP_Author_Guidelines.docx` = The
Biophysicist (revista educativa BPS), NO sirve, descartar.

## 1. Requisitos que gobiernan (verificados)

| requisito | BJ | nuestro estado |
|---|---|---|
| tipo | Research article; sin límite de páginas ("no page limitations for research articles") | ok |
| abstract | ≤300 palabras, sin citas, no-especialista | 188 (def. wordcount.py) — OK |
| Significance | ≤120 palabras, bajo el abstract, obligatorio, público general | NO EXISTE → escribir |
| título | ≤150 chars con espacios | 93 — OK |
| orden | intro→materials and methods→results→discussion(→conclusion); "STAR methods acceptable" | Methods está al FINAL; precedente publicado con Methods tras Discussion: Pinto-Anwandter BJ 124:2500 (2025). DECISIÓN L. |
| apéndices | slot propio EN el manuscrito (Main text→Appendices→Data availability→...) | los 4 apéndices (8.357 palabras) SE QUEDAN |
| referencias | numeradas (1,2)/(3–5), et al. tras 10, abrev. Biosis, "Supporting citations" al final | hoy vancouver-elife.bst → swap + auditoría \citet |
| figuras | 3.25 / 6.75 in, máx 1 página c/u, TIFF/PDF, 300dpi (línea 1000), Arial, RGB, ≤20MB | eLife 7.0→6.75 in (−3,6%; ¿7pt→6.75pt bajo piso?); Helvetica≈Arial verificar; alturas Fig 3/4 vs 1 página |
| supporting | UN PDF "Document S1" <10MB; "Figure S1…"; "Supplemental", NO "Supplementary" | armar |
| estadística | test, N exacto, centro/dispersión; software+versión en Methods | ya cumplido (hash macro_dr) |
| data availability | statement obligatorio, fórmula "[Type] data have been deposited at..." | reusar el de eLife; completar <URL>/<DOI> Zenodo (¡pendiente de antes!) |
| cover letter | OBLIGATORIA, confidencial: contenido+significancia, competing interests, apariciones previas, relacionados, estadístico si hubo | rehacer venue (§5) |
| preprint | bioRxiv OK; statement en acknowledgments; submit directo B2J posible | VERIFICAR si el depósito bioRxiv se hizo |
| IA | política Elsevier: declaración al final si se usó en la escritura | decisión L. |
| peer review | ~1 mes a decisión; UNA revisión mayor; 3 meses para revisar | — |

## 2. Ruta LaTeX

- Inicial format-neutral de facto: sólo manuscrito+carta; Word/LaTeX/PDF/.txt; figuras embebidas
  preferidas. PERO no mandar el build vestido de eLife → puerto liviano a clase neutra ANTES.
- Base recomendada: template LaTeX oficial Cell Press v1.10 (article 12pt + geometry + numbered.bst,
  derivado de model6-num-names; NO elsarticle). Overleaf "Biophysical Journal Template": 7–10 años,
  clase biophys(-new), EVITAR. CSL existe (biophysical-journal.csl, num., actualizado 2025-12) por si
  hiciera falta cotejar formato.
- Reemplazos (censo medido): \figsupp ×12 (todo en 04_results; Fig2:1 Fig3:3 Fig4:4 Fig5:2 Fig6:2)
  → S-figs; fullwidth ×8 (2 framework + 6 results) → figure*/6.75in; appendixbox ×5 + \appendix →
  secciones planas (refs todas vía \ref, nada hardcodeado); abstract env (elife.cls lo redefine);
  \author/\affil/\corr → bloque BJ; lineno: ACTIVAR (elife lo trae sin activar); stix ya resuelto.
  "eLife" literal en texto vivo: 0 (sólo \graphicspath).
- Final files (post-aceptación): .tex una columna + .bbl + PDF con refs; figuras "Figure 1.tif";
  producción convierte LaTeX→Word (+3 días hábiles).

## 3. Números medidos (workflow; def. wordcount.py = texto corrido, captions aparte)

- MAIN TEXT def. eLife (intro+framework+results+discussion): 11.596 + 4.289 captions ≈ Münch 11.515.
- Cuerpo total (con methods+abstract+backmatter): 16.992. Apéndices: 8.357. Supp file: 1.960 + 3 tablas, 0 figs.
- PDF actual: 66 pág. eLife + 5 supp. Figuras: 7 cuerpo + 12 figsupp.
- Ejemplares BJ: abstracts 210–263, Significance 93–122, Methods antes de Results en 3/4 (después en
  Pinto-Anwandter 2025), hasta 16 figs de cuerpo (Sigg 2025) → sin presión de tamaño.

## 4. Adaptar vs mudar — veredicto

NADA se muda por obligación. Adaptación = FORMA:
(a) Significance nuevo (semilla: impact statement 30 palabras + cierre del plato);
(b) refs numeradas; (c) 12 figsupp → Document S1 Figs S1–S12 + mapa de renumeración + reescribir
"Figure N—figure supplement M" → "Fig. SK" en texto; fusionar supplementary_file_1 (→ Tables S1–S3);
(d) resize figuras; (e) front/back matter BJ (Data availability, Author contributions,
Declaration of interests "The authors declare no competing interests.", acknowledgment preprint).
Opcional: Conclusion corta (3/4 ejemplares la tienen); promover 1–2 figsupp a cuerpo si conviene.

## 5. Carta y ruteo

- Reescribir el párrafo de venue: en BJ el argumento es que ESTA es la comunidad que usa estas
  verosimilitudes (Moffatt 2007, Milescu 2005, Celentano 2004, Clerx 2019 = todos BJ). Pisar el
  criterio positivo del scope ("significant methodological or technological advances... open new
  areas") y esquivar el negativo ("improvements in accuracy or speed of existing methods... not
  suitable"). Declarar: preprint, sin relacionados en revisión, no competing interests.
- Handler (Associate Editor): Valeria Vásquez (channels, exp+comp) vs Jeremy Smith (computational).
  Recomendación: Vásquez. DECISIÓN L.
- Reviewers al form EM: Qin, Benndorf, Mirams, Sivilotti, Plested, Kinz-Thompson (+Milescu si L.
  acepta el riesgo de que decline).

## 6. Tarifas y waiver — GATE previo a todo envío

- Post-2024 tarifa plana: no-OA $1.200 / $950 miembro BPS; OA $3.000 / $2.400. Sin page/color charges.
  (El "encourage ≤12 pages + higher charges" es del régimen viejo; asterisco de "None*" irresuelto.)
- WAIVER (paste L. 2026-09-03): sólo circunstancias excepcionales, SOLO no-OA, pedir a
  bj@biophysics.org y tener APROBACIÓN ANTES de enviar; enviar sin waiver = compromiso de pago.
  Ejemplo dado: país Research4Life (Argentina probablemente NO califica, ingreso medio-alto; el
  ejemplo no es taxativo). Si funder exige OA: descuento sobre OA = tarifa no-OA.
- Membresía BPS ahorra $250 en no-OA. ¿L. es miembro? DECISIÓN L.

## 7. Secuencia

0. GATE: OA/no-OA + membresía + waiver (si va: email y ESPERAR respuesta). Verificar bioRxiv y Zenodo.
1. Decisiones de forma: Methods al final (sí/no), Conclusion (sí/no), promociones de figsupp.
2. Puerto mecánico: clase, refs numeradas, figsupp→S, appendices, front matter, lineno.
3. Textos nuevos: Significance, carta BJ, acknowledgment preprint, data availability con DOIs.
4. Figuras: resize + specs; ensamblar Document S1.
5. Build + checks: swallowed_sentences.py, pdftotext | grep "??", refs, orden de citas numeradas.
6. EM: form (reviewers, AE), submit (o B2J desde bioRxiv).

Pendientes heredados: asterisco "None*"; si EM pregunta rechazos previos (transfer intra-BPS dice
que no hace falta mencionar reviews previas; eLife no es BPS → desk reject de eLife no se declara
salvo que el form lo pida explícitamente — la carta sí debe declarar el PREPRINT, no el rechazo).
