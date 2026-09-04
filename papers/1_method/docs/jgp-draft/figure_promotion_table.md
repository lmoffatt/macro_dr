# A1: promoción de los 12 figure supplements — tabla para veto — 2026-09-03

Regla JGP: sin límite de figuras; supplemental desalentado ("it should be possible to read and
understand the article without consulting supplemental materials"). Referencia de género:
Benndorf & Schulz 2023 = 22 figuras de cuerpo. Voto = mío; veto = L.

| # | supp | título corto | usa la prosa | carga | voto |
|---|---|---|---|---|---|
| 1 | Fig2–S1 | Sandwich vs distribución empírica (QQ Mahalanobis) | Results ×1 ("in all six directions at once") | única validación directa de la reparación sandwich | **PROMOVER** |
| 2 | Fig3–S1 | Información por paso vs varianza del score (k_on, k_off) | Results (washout) + Methods (anchor singular) | evidencia de la subsección de Discusión "channel number stops informing" | **PROMOVER** (fusionable con S2 en una sola) |
| 3 | Fig3–S2 | Ídem (i, N_ch) | ídem | ídem | **PROMOVER** (o fusionar) |
| 4 | Fig3–S3 | Residuo y score en los 3 parámetros no coloreados | ninguna mención directa | completitud | **ELIMINAR** (veto L. 2026-09-03: "puede obviarse"; verificar cero refs antes de borrar) |
| 5 | Fig4–S1 | Segundo momento factorizado: por-muestra × correlación | Results ×2 + Discusión ×2 (Milescu "now measured"; Münch) | LA figura de mecanismo; sostiene el 20× de la carta | **PROMOVER** |
| 6 | Fig4–S2 | Las cuatro rungs recursivas (MR/VR) | Results + Discusión ("neither partial correction reaches") | evidencia de una oración del ABSTRACT | **PROMOVER** |
| 7 | Fig4–S3 | Magnitud vs anisotropía (m, a) | Results (los dos brazos LSE se separan) | evidencia de "no single rescaling repairs it" (ABSTRACT) | **PROMOVER** |
| 8 | Fig4–S4 | Error estándar corregido por distorsión | Results + Discusión ("the repair is the sandwich") | qué entrega cada miembro; cierre práctico | **PROMOVER** |
| 9 | Fig5–S1 | Ruido y N_ch actúan sólo por su cociente | Results ×3 (colapso en r) | el colapso GRADÚA con el miembro: gobierna a los no recursivos, R intermedio (23%), IR plano | **PROMOVER + ORACIÓN** (veto L. 2026-09-03; borrador de la oración en el chat: cociente S̃/N para el fallo de los no-recursivos vs producto N·S̃ para el límite de IR; verificar paneles dibujados antes de afirmar miembros) |
| 10 | Fig5–S2 | Cada miembro contra la memoria de su propio residuo | Results ×2 ("a white residual does not certify calibration") | evidencia de una oración del ABSTRACT | **PROMOVER** |
| 11 | Fig6–S1 | La ley del límite de IR, con su dispersión | Results + Discusión + caption de Fig 7 (extrapolación) | el límite que ahora está en la CARTA | **PROMOVER** |
| 12 | Fig6–S2 | El registro a ambos lados del límite (trazos) | Results ("what these cells look like") | la foto física de la frontera macro/micro | **PROMOVER** |

Saldo con vetos de L. (2026-09-03): 7 + 11 promovidas − 1 eliminada = **18 figuras de cuerpo**
(17 si Fig3-S1/S2 se fusionan; L. no se expidió → default NO fusionar), **cero figuras
suplementarias**; supplementary_file_1 queda como documento suplementario (re-anclaje, 3 tablas;
default, L. no se expidió) → párrafo-resumen al final de Methods obligatorio.

ALTERNATIVA "todo adentro": promover las 12 → 19 figuras, cero suplementario salvo el
supplementary_file_1 (o absorberlo también y quedar en cero). Máxima ortodoxia JGP, costo:
dos figuras de completitud diluyen el cuerpo.

Nota de ejecución (A2): la renumeración es automática vía \ref al convertir \figsupp en figure;
lo manual son las 23 menciones compuestas "Figure~\ref{...}--figure supplement N" y el orden de
primera mención = orden numérico (requisito JGP; chequear tras mover Methods).


## Fusión time-s1+s2: RECHAZADA por la prueba de página (2026-09-03, gate de L.)
Aritmética: cada mitad 5.6×5.6 in (2 bandas c/u, ggsave figure_3_S1_S2.Rmd:367); 4 bandas fusionadas
≈ 10.3 in de arte (ahorro de un header y una leyenda incluido) + caption fusionada ~210 palabras
≈ 2.4 in = ~12.7 in, contra ~9.4 in de página útil. Comprimir filas 40% viola fuentes-primero.
Quedan DOS figuras. El costo aceptado: el contraste del washout sigue partido entre ambas; lo
carga la prosa (cita conjunta) como hasta ahora.

## Relabel LSE→ILSE en Figs 5-6: EJECUTADO 2026-09-03
Verificado contra DTOK (figure_3.Rmd:143): digest "LSE" = av=1 = ILSE; el LSE verdadero (LSE_av0)
no se dibuja en estas figuras. Editado ALAB/FULL en figure_3_S1_S2.Rmd (token de archivo intacto),
comentario de procedencia agregado; la caption del tex YA decía ILSE (04_results.tex:576), o sea
que el conflicto era arte-vs-caption y muere con el re-knit.
