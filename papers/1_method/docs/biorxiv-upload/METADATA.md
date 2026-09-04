# DEPOSITADO 2026-09-04 20:03 — MS ID BIORXIV/2026/749531 — CC-BY 4.0 — en screening; DOI pendiente.
# La conversion de bioRxiv fue passthrough: su PDF es byte-identico al subido (md5 e667f101...).

# Paquete bioRxiv — listo para subir — 2026-09-04

## Archivos (en esta carpeta)
1. jgp_biorxiv_moffatt.pdf — el manuscrito completo, 66 pag, figuras embebidas (build 48c8d06c,
   rc=0, tipografia reparada). bioRxiv acepta PDF unico directamente.
2. supplemental_material.pdf — 5 pag (retitulado "Supplemental material: ..."), subir como
   Supplementary Material.

## Campos para pegar en el formulario

**Title**
Likelihood approximations distort the ion channel kinetic information in macroscopic currents

**Author**
Luciano Moffatt — lmoffatt@qi.fcen.uba.ar — [tu ORCID] — corresponding.
Affiliation: Instituto de Quimica Fisica de los Materiales, Medio Ambiente y Energia (INQUIMAE),
CONICET; Facultad de Ciencias Exactas y Naturales, Universidad de Buenos Aires, Ciudad de Buenos
Aires, C1428EHA, Argentina

**Abstract** (188 palabras; los simbolos Δ̃/S̃ son Unicode, pegan bien)
Kinetic inference from macroscopic ion channel currents rests on likelihoods that approximate an intractable process, and what the approximation costs the uncertainty they report had not been measured. The diagnostics are classical: for a correctly specified likelihood the score, the log-likelihood's gradient, has mean zero at the generating parameters and its covariance equals the Fisher information it reports. We measured both for eight likelihoods in a two-state scheme, over channel number, instrumental noise and acquisition interval, with ten thousand exact simulations per condition. Those leaving the correlated gating fluctuation unmodelled, among them the control that reversed our evidence for the asymmetric activation of P2X2, understate their error bars two- to threefold wherever the residual retains memory, in directions no rescaling repairs, and a white residual does not certify calibration. Two partial corrections reduce the bias and worsen the error bar. A filter conditioning each interval average on both endpoints stays calibrated except where single openings become resolvable, and there its departure follows one power law. Least squares becomes calibrated about two decades of noise above where the fluctuations stop giving the unitary current, so no condition gives both.

**Subject Area**: Biophysics

**Licencia** — DECISION TUYA en el formulario; las opciones y su letra chica:
- CC-BY: maximo reuso; si algun dia fueras a Gold OA en JGP (CC-BY) queda todo coherente.
- CC-BY-NC-ND: la mas restrictiva de las CC; comun en preprints que van a revista de suscripcion.
- "All rights reserved / no reuse": tambien valida para JGP (solo pide que el preprint exista).
Ninguna afecta la submission a JGP (Green): eleccion libre.

**Declaraciones que el formulario pide:**
- Competing interests: The author declares no competing interests.
- Funding: [PENDIENTE L.: agencia y numero de subsidio exactos, o "self-funded"; el mismo texto
  ira a los Acknowledgments de JGP, que hacen text-mining de nombres de agencia]
- Data/code availability (por si hay campo): All data and code are archived at Zenodo:
  code 10.5281/zenodo.22167744, data 10.5281/zenodo.22168409, macroir library
  10.5281/zenodo.22168263. (Ya esta tambien dentro del PDF, seccion Data availability.)

## Al terminar el deposito (destraba la submission a JGP)
1. Anotar DOI y fecha que asigne bioRxiv.
2. Yo inserto en el manuscrito (Acknowledgments): "A preliminary version of this work, [DOI],
   was deposited in bioRxiv on [date]." y recompilo.
3. Relleno el [PENDIENTE] de la carta (cover_letter_jgp.md) y compilo su PDF.
4. En el form de JGP: declarar el preprint con su DOI.
5. CONFIRMADO en el Submission Help (2026-09-04): 'The Journal of General Physiology' ESTA en la
   lista B2J de bioRxiv ('Submit bioRxiv Preprint to a Journal'). Tras el screening, la submission
   a JGP puede dispararse desde el Author Area con los mismos archivos, sin re-subir nada.

## Nota
El PDF principal es el MISMO que ira a JGP (format-neutral): un solo artefacto, dos destinos.
