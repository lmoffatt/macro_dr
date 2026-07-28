# Recording configurations: currents, channel counts and noise

Collected 2026-07-26 to place real preparations on Figure 6's design plane. Each number below is
either **sourced** (with the link) or **derived** (with the arithmetic). Nothing here is a guess left
unlabelled.

UPDATED 2026-07-27: the noise numbers are now stated on a single basis, **rms in a 1 kHz bandwidth**,
for TEVC oocyte / HEK293 whole cell / excised outside-out patch. See the next section. Bib entries
for every source added to `papers/1_method/docs/manuscript-drafts/biblio_full.bib`. Two corrections
are flagged there; one of them (the macropatch S convention) is NOT applied, pending Luciano.

## The conversion, verified from the engine

`legacy/qmodel.h:3740`: `e = Current_Noise * fs / number_of_samples`, so the variance of one recorded
point is `Current_Noise / Δ`. **`Current_Noise` is therefore a noise power spectral density**, and
with the dispatcher's `label = 1000 · Current_Noise` and the dimensionless group `ν = S·k_off/i²`:

    label = 10 · k_off · σ_raw² / (fs · i²)        general
    label = σ_Δ² / i²                              at Δ = 0.1 τ = 1 ms, k_off = 100 /s

where σ_Δ is the RMS of a point averaged over Δ. For white noise of RMS σ_B measured in a bandwidth
B, σ_Δ² = σ_B² / (2BΔ). So **√label is the RMS noise of a 1 ms point in units of the unitary
current** — the y axis has a direct physical reading.

Channel counts follow from the macroscopic current: N_ch = I_max / (i · P_open), and with the paper's
model (i = 1 pA, P_open = 0.5), N_ch = 2 · I_max / pA.

## Noise per preparation, all on ONE basis: rms in a 1 kHz bandwidth

Rebuilt 2026-07-27. Everything below is rms current noise in a bandwidth from DC to B = 1 kHz.
One basis, because the published figures come at four different bandwidths and are not comparable
until carried to a common one. Bib keys refer to `papers/1_method/docs/manuscript-drafts/biblio_full.bib`.

| preparation | rms @ 1 kHz | range | how obtained | source |
|---|---|---|---|---|
| oocyte, conventional two-electrode (TEVC) | **3.3 nA** | 1.5–6 nA | measured, native 1 kHz, no conversion | `nagel2025p2x4` Fig 8b |
| oocyte, cut-open vaseline gap (the LOW-noise oocyte method) | **0.1–0.5 nA** | | carried from 1.2 nA @ 5 kHz | `stefani1998cutopen` |
| HEK293 whole cell | **~1 pA** | 0.3–1.5 pA | computed from Axon Guide eq. 16 at 1 kHz | `axonguide1993` + traces below |
| excised outside-out patch, routine | **~0.05 pA** | 0.05–0.2 pA | component budget at 1 kHz (below) | `axonguide1993` |
| excised outside-out patch, exceptional | **0.01–0.04 pA** | | carried from 0.083 pA @ 5 kHz | `rae1992exceptionally`, `levis1993quartz` |
| outside-out MACROpatch (big tip, many channels) | **0.07–0.7 pA** | | carried from 2.28 pA @ 10 kHz | Moffatt & Hume 2007 own record |
| capacitive-feedback headstage, open circuit | **0.01–0.03 pA** | | carried from 0.06 pA @ 5 kHz | `axonguide1993` (Axopatch 200A) |

### The carry factors, and why the ranges are wide

Carrying an rms figure from bandwidth B₀ down to 1 kHz needs the spectral shape, because rms is an
integral over the band. Three regimes, and they differ by up to 10× over a decade of bandwidth:

| spectrum | physical origin | σ scales as | ×factor 10 kHz→1 kHz | 5 kHz→1 kHz | 2 kHz→1 kHz |
|---|---|---|---|---|---|
| white, S const | seal resistance, feedback resistor | B^0.5 | 0.316 | 0.447 | 0.707 |
| S ∝ f | pipette + holder dielectric | B^1.0 | 0.100 | 0.200 | 0.500 |
| S ∝ f² | voltage noise across a capacitance (eₙ·2πf·C) | B^1.5 | 0.032 | 0.089 | 0.354 |

Real records are a mixture, so a single carried number is a range, not a value. **This is why the
1 kHz basis is worth the trouble: the only entry that needs no conversion at all is the TEVC one.**

### Whole cell, computed rather than carried

Axon Guide ch. 12, "Noise in Whole-Cell Voltage Clamping", eq. 16:

    S_wc(f) = 4π²f²·e_s²·C_m² / (1 + 4π²f²·τ_sr²),   e_s² = 4kT·R_s,   τ_sr = R_sr·C_m

with R_s the series resistance, R_sr the residual (uncompensated) part, C_m the cell capacitance,
k Boltzmann's constant, T absolute temperature. Integrated 0 to 1 kHz, at HEK293-typical values:

| R_s | C_m | compensation | rms @ 1 kHz |
|---|---|---|---|
| 5 MΩ | 10 pF | 0% | 0.32 pA |
| 5 MΩ | 20 pF | 0% | 0.59 pA |
| 10 MΩ | 20 pF | 0% | 0.68 pA |
| 10 MΩ | 20 pF | 70% | 0.89 pA |
| 10 MΩ | 33 pF | 0% | 0.86 pA (the Guide's own worked example) |

The trace-measured whole-cell median of §"Whole-cell noise, MEASURED" below is 1.9 pA at 2 kHz. The
same equation's own 1 kHz : 2 kHz ratio for these parameters is ~0.58, so that median back-converts
to ~1.1 pA at 1 kHz. Computed and measured agree. Note compensation makes whole-cell noise WORSE,
not better: it restores bandwidth that the R_s·C_m pole was filtering away.

### Excised patch, component budget at 1 kHz

Rather than carry a 5 or 10 kHz figure down through an unknown spectral mixture, build it at 1 kHz:

| term | value @ 1 kHz | basis |
|---|---|---|
| seal thermal noise, 4kT/R_seal, R_seal = 10 GΩ | 0.040 pA | white |
| headstage (capacitive feedback) | ~0.025 pA | 0.06 pA @ 5 kHz carried |
| pipette dielectric, Sylgard-coated borosilicate | ~0.015 pA | 0.15 pA @ 10 kHz, S ∝ f |
| **quadrature sum** | **0.049 pA** | |

A 5 GΩ seal alone raises the first term to 0.057 pA. Routine patches land 2–4× above this budget
(imperfect seals, larger tips, 1/f, mains), hence the 0.05–0.2 pA range in the main table.

### The oocyte TEVC number, verified three ways (2026-07-27)

The TEVC figure is the one most worth doubting, so it was attacked from three sides and survived.

**1. It needs no bandwidth conversion.** `nagel2025p2x4` Methods: currents were *"filtered at 1 kHz,
and digitized at 1 kHz"*. The 3.3 nA read off Fig 8b is therefore natively an rms in a 1 kHz band.
Every other entry in the table above is a carried number; this one is not.

**2. First-principles floor, 2.9 nA.** The voltage electrode's thermal noise is imposed across the
oocyte capacitance by the clamp:

    S_I(f) = 4kT·R_e·(2πf·C_m)²      σ = sqrt( 4kT·R_e·(2πC_m)²·B³/3 )

with R_e the voltage-electrode resistance. This is Axon Guide eq. 16 in the limit of full series-
resistance compensation, so it is not an ad hoc model. At R_e = 1 MΩ, C_m = 200 nF, B = 1 kHz it
gives **2.92 nA**, against 3.3 nA measured. The measurement sits just above its own thermal floor,
which is where a real record should sit.

  - C_m is sourced, not assumed: `baumgartner1999tevc` measured specific membrane capacitance
    45.6 ± 5 mF/m², which on a 1.2 mm oocyte is 206 nF.
  - Sensitivity: R_e 0.5–5 MΩ and C_m 150–250 nF span 1.55–8.16 nA. Hence the 1.5–6 nA range.
  - A finite clamp-loop corner lowers it (corner 1 kHz gives 2.34 nA, 500 Hz gives 1.69 nA), so
    2.9 nA is an upper bound on the thermal-only part.

**3. Ordering against cut-open holds.** `stefani1998cutopen` report 1.2 nA rms at 5 kHz (time
constant 24 µs) for the cut-open clamp. Carried to 1 kHz that is 0.1–0.5 nA, so conventional TEVC
sits 6–30× above the low-noise oocyte method. Direction correct. The factor is 6–30× and not the
"33×" claimed in the older text below, because that comparison was made in S rather than in rms and
S is bandwidth-invariant only for white noise, which this is not.

**Bandwidth sensitivity is extreme here and must always be quoted with the number.** While the f²
term dominates, σ ∝ B^1.5, so the SAME rig gives ~0.09 nA at 100 Hz, which is where TEVC is
normally actually filtered. "3 nA" without "at 1 kHz" is meaningless-to-misleading.

### Two corrections to the older text below

1. **`baumgartner1999tevc` contains NO noise analysis.** It is "voltage errors and compensation for
   local current flow". Full text checked 2026-07-27. It is a good citation for oocyte membrane
   capacitance and for nothing else here. It was previously listed under oocyte links in a way that
   implied it backed a noise figure.

2. **The macropatch S = 1.04e-4 pA²/Hz mixes conventions with the whole-cell S values.** NOT CHANGED,
   flagged only, pending Luciano. It is σ²/fs = 2.28²/50000, i.e. the ENGINE's convention, in which
   `Current_Noise = σ_raw²/fs` (`qmodel.h:3740`: variance per point = `Current_Noise/Δ`, presuming
   independent raw samples). That is internally correct and is what `figure_6.Rmd` consumes. But the
   whole-cell entries in the measured table below are σ²/B_filter with B_filter = 2 or 10 kHz. The
   two are not the same quantity: the raw record is a 10 kHz Bessel sampled at 50 kHz, so adjacent
   samples are correlated and the physical in-band PSD is 2.28²/10 kHz = 5.2e-4 pA²/Hz, 5× larger.
   Consequence if confirmed: the "whole cell sits 1.4× to 86× above it, median 17×" ordering check
   compares two different conventions and the true span is ~0.3–17×, median ~3.5×, with the quietest
   whole-cell records NOT above the macropatch. Also, a first-difference estimator underestimates σ
   when adjacent samples are correlated (σ_fd = σ·sqrt(1−ρ)), so the true 10 kHz σ exceeds 2.28 pA.
   The 1 kHz table above sidesteps the whole issue by never using S.

## Noise, as originally published (provenance for the table above)

| configuration | RMS noise | bandwidth | source |
|---|---|---|---|
| capacitive-feedback headstage, open circuit | 0.06 pA | 5 kHz | Axon Guide ch. 12 (Axopatch 200A) |
| resistive-feedback headstage | 0.25–0.30 pA | 10 kHz | Axon Guide ch. 12 |
| pipette dielectric noise, ordinary glass | 0.15–0.44 pA | 10 kHz | Axon Guide ch. 12 |
| excised patch, exceptional (quartz + integrating) | 0.083 pA | 5 kHz | `rae1992exceptionally` |
| excised patch, real cell patches, quartz | ≤0.10 pA, best 0.06 pA | 5 kHz | `levis1993quartz` |
| whole cell, integrated amplifier | 8 pA | 10 kHz | patch-clamp-on-a-chip, PMC2978236 |
| oocyte, cut-open / vaseline gap | 1.2 nA | 5 kHz | `stefani1998cutopen` |
| oocyte, cut-open / vaseline gap | ~1 nA | 3 kHz | Pantazis & Olcese, PMC11549981 |
| oocyte, conventional two-electrode (TEVC) | 3.3 nA | 1 kHz | `nagel2025p2x4` Fig 8b, measured |

Whole-cell noise is dominated by the series resistance in series with the membrane capacitance
(Axon Guide ch. 12, "Noise in Whole-Cell Voltage Clamping"; their worked example is Cm = 33 pF,
Rs = 10 MΩ). Standard two-electrode oocyte clamp is **worse** than the cut-open figure above, which
is quoted as the low-noise technique. That gap is now measured rather than assumed, so the oocyte
box is drawn from ~3 nA upward, not from 1 nA.

Note the two cut-open figures are mutually inconsistent under any single spectral exponent (1.2 nA
at 5 kHz and 1 nA at 3 kHz imply σ ∝ B^0.36, flatter than white). Different rigs, so do not fit an
exponent to them; they bound the cut-open value rather than determining it.

## Currents, sourced

| configuration | macroscopic current | source |
|---|---|---|
| whole cell, modest expression (HEK) | ~300 pA | PMC / J Neurophysiol, transfected HEK 293 |
| whole cell, robust expression (HEK, Na) | > 2 nA, several nA | same |
| oocyte, two-electrode | up to several tens of µA without losing clamp control | ScienceDirect, TEVC review |

Macropatch: no single citable range found in this pass; drawn from the patch-clamp current range
(unitary currents < 10 pA resolved, so a macropatch of tens to hundreds of channels gives tens to
hundreds of pA). **This is the weakest of the three boxes and should be replaced with a sourced
figure.**

> **SUPERSEDED 2026-07-27 by `EXPRESSION_LEVELS.md` beside this file.** The thin table above is replaced
> by a sourced sweep of six configurations, with the macropatch box now carried by counted numbers.

## Whole-cell noise, MEASURED off published traces (2026-07-27)

Because no whole-cell rms noise measured with a cell attached exists in the literature (§2 of the
verified sweep: every published number is an amplifier floor, a model cell, or theory), it was read
off published traces with the recipe in `TRACE_WIDTH_METHOD.md`. Eight papers attempted, 24 panel
segments measured, 11 returned a number, 13 hit the ink floor or had no usable segment.

| Citation | Figure / panel | σ (pA) | B (Hz) | S = σ²/B (pA²/Hz) |
|---|---|---|---|---|
| Li 2019 eLife 8:e47060, hP2X3, HEK293 | Fig 3B Mg+Ca, pre-ATP | 0.55 | 2000 | 1.53e-4 |
| Sokolov 2013 Front Pharmacol 4:78, Na_V1.5, HEK293 | Fig 6A inset | 1.38 | 10000 | 1.90e-4 |
| Li 2019 | Fig 3B MgCl₂, pre-ATP | 0.78 | 2000 | 3.06e-4 |
| Karasawa 2016 eLife 5:e22153, pdP2X7, HEK293 | Fig 1C pre-ATP | 1.65 | 2000 | 1.36e-3 |
| Li 2019 | Fig 3B EDTA, pre-ATP | 1.84 | 2000 | 1.70e-3 |
| Karasawa 2016 | Fig 1C deactivation tail | 1.92 | 2000 | 1.85e-3 |
| Karasawa 2016 | Fig 1F inter-pulse, 2nd cell | 2.03 | 2000 | 2.07e-3 |
| Li 2019 | Fig 3B CaCl₂, pre-ATP | 2.96 | 2000 | 4.38e-3 |
| Harnau 2025 Sci Rep 15:32686, GlyR α3L, HEK293T | Fig 1C pre-glycine | 4.17 | 2800 | 6.21e-3 |
| Dallas 2021 Sci Rep 11:8194, Kv2.1, HEK293 | Fig 1A zero-current sweeps | 4.21 | 2000 | 8.87e-3 |

**Range S = 1.5e-4 to 8.9e-3 pA²/Hz, median 1.8e-3**, over five independent papers, different labs,
four channel systems, three amplifier families. Ten segments across about seven cells, so not ten
independent draws. The factor-30 spread is cell-to-cell, as expected when C_m and R_s vary
several-fold between HEK cells; within one figure Li spans 0.55 to 2.96 pA across four cells.

**Ordering check passes.** The outside-out macropatch of Moffatt & Hume 2007 measures S = 1.04e-4
pA²/Hz. Whole cell sits 1.4× to 86× above it, median 17×.

**These are lower bounds on typical noise**: published traces are the best records of a set.

One measurement is excluded from the table because it is not baseline: Karasawa Fig 1C at the peak of
the 1 mM ATP response gives σ = 5.49 pA. Its excess over the same trace's own baseline,
5.49² − 1.65² = 27.4 pA² on 289 pA of mean current, is **P2X7 gating noise**, not instrumental. Listed
so it is not mistaken for a floor — and it is the same phenomenon found in Moffatt & Hume's own
baseline, where 99 % of the variance of a 13.7 ms point is gating rather than amplifier.

## Oocyte TEVC noise, MEASURED off published traces (2026-07-27)

Four papers reached the measurement stage. One returned a number that survives every check; the rest
either floor out or land in a physically impossible place.

| Citation | Figure / panel | σ | B (Hz) | S = σ²/B (pA²/Hz) | verdict |
|---|---|---|---|---|---|
| Nagel 2025 Nat Commun 16:10367, hP2X4-E307T | Fig 8b left, pre-ATP baseline of the black trace | 3.30 nA | 1000 | 1.09e4 | **use this** |
| Ito 2020 Sci Rep 10:13999, zGlyR α1+βb | Fig 3c, six inter-application baselines | 0.49 nA | 10 | 2.44e4 | upper bound only |
| Jackson 2025 IJMS 26:9506, hα3β2 nAChR | Fig 2b, post-wash and pre-ACh baselines | 10.9 pA | 2000 | 5.9e-2 | **rejected, see below** |
| Jackson 2025 | Fig 2a, post-wash plateau | 26.3 pA | 2000 | 3.4e-1 | rejected, same reason |

**The anchor is Nagel: σ = 3.3 nA rms at 1 kHz, S = 1.1e4 pA²/Hz.** It sits 33× above the cut-open
figure, which is the direction the ordering demands, and two independent physical routes agree with
it. Scaling the whole-cell median (1.9 pA at 2 kHz) by the capacitance ratio of an oocyte to a HEK
cell (200 nF / 33 pF = 6000) and the electrode-noise ratio (284 nV/√Hz over 400) predicts 8.2 nA at
2 kHz; Nagel's own value carried from 1 kHz to 2 kHz on the e_n·2πf·C_m spectrum (σ ∝ B^1.5) gives
9.3 nA. The measurement was also validated on synthetic renders, and the panel's occlusion was
handled with a top-edge estimator that uses neither the ink width nor the occluded lower edge.

> SUPERSEDED IN PART, 2026-07-27, by §"The oocyte TEVC number, verified three ways" above. The
> measured 3.3 nA stands and is now known to need NO bandwidth conversion (Nagel's Methods: filtered
> at 1 kHz AND digitized at 1 kHz). Two claims in the paragraph above are weakened:
> (a) the "33×" over cut-open was computed in S, which is bandwidth-invariant only for white noise
>     and this spectrum is f²; done in rms at a common 1 kHz the factor is 6-30×.
> (b) the capacitance-ratio route used e_n = 284 nV/√Hz, which implies a ~5 MΩ electrode; oocyte TEVC
>     electrodes are 0.5-2 MΩ. The cleaner check is the direct thermal floor,
>     σ = sqrt(4kT·R_e·(2πC_m)²·B³/3) = 2.92 nA at R_e = 1 MΩ, C_m = 200 nF, B = 1 kHz, against
>     3.3 nA measured. C_m from `baumgartner1999tevc` (45.6 ± 5 mF/m² on a 1.2 mm oocyte = 206 nF).
> The conclusion is unchanged and now rests on a better argument.

**The two Jackson numbers are rejected.** They fall 1000× to 5700× BELOW the cut-open value, which is
impossible: cut-open is the low-noise oocyte method and conventional TEVC must sit above it. They
also put an oocyte within 33× of a HEK cell that has 6000× less membrane capacitance. Either the
displayed trace is not a 2 kHz record (decimation or smoothing in figure preparation, in which case
both the bandwidth and the assumed 610 samples per pixel column are wrong), or the bandwidth
attribution is wrong. Nothing in the PDF settles which. Do not use these numbers.

**Ito is an upper bound, not a value.** The by-the-book binary step returns extent − ink = 0; the
0.49 nA comes from subpixel variants at 6.8 nA per pixel. Its S also cannot be compared with the
others: S = σ²/B is bandwidth-invariant only for white noise, and the dominant oocyte TEVC term rises
as f², so σ² grows like B³ and S like B². Nagel's spectrum carried down to 10 Hz predicts 3.3 pA,
150× below what Ito can resolve, which is consistent with an unresolved trace.

**These are lower bounds**: published traces are the best records of a set, and Nagel's is explicitly
a representative sweep.

The oocyte box on the figure should now be drawn from **~3 nA rms at 1 kHz (S ≈ 1e4 pA²/Hz)**, not
from the ~1 nA cut-open number, which is the low-noise method and underestimates conventional TEVC by
about 30×.

## Links

Bib keys are in `papers/1_method/docs/manuscript-drafts/biblio_full.bib`, section
"Recording-configuration noise", all metadata verified against CrossRef 2026-07-27.

- `axonguide1993` Axon Guide ch. 12, Noise in Electrophysiological Measurements. PDF beside this file.
  https://www.ifsc.usp.br/~reynaldo/eletronica_instrumentacao/axon-guide/GuideCh12.pdf
- `rae1992exceptionally` **Rae & Levis** (NOT Levis & Rae, author order corrected 2026-07-27 from
  CrossRef), A method for exceptionally low noise single channel recordings, Pflügers Arch 420:618-620.
  https://doi.org/10.1007/BF00374642
- `levis1993quartz` Levis & Rae, The use of quartz patch pipettes for low noise single channel
  recording, Biophys J 65:1666-1677. https://doi.org/10.1016/S0006-3495(93)81224-4
- `levis1998lownoise` Levis & Rae, Low-noise patch-clamp techniques, Methods Enzymol 293:218-266.
  https://doi.org/10.1016/S0076-6879(98)93017-8
- `stefani1998cutopen` Stefani & Bezanilla, Cut-open oocyte voltage-clamp technique, Methods Enzymol
  293:300-318. Source of 1.2 nA rms at 5 kHz. https://doi.org/10.1016/S0076-6879(98)93020-8
- `baumgartner1999tevc` Baumgartner, Islas & Sigworth, Two-microelectrode voltage clamp of Xenopus
  oocytes: **voltage errors and compensation for local current flow**, Biophys J 77:1980-1991.
  **CITE FOR CAPACITANCE ONLY (45.6 ± 5 mF/m²), it contains no noise analysis.**
  https://pmc.ncbi.nlm.nih.gov/articles/PMC1300479/
- `nagel2025p2x4` Nagel et al., Nat Commun 16:10367. PDF beside this file. Fig 8b + Methods
  ("filtered at 1 kHz, and digitized at 1 kHz"). https://doi.org/10.1038/s41467-025-66244-3
- Patch-clamp amplifiers on a chip — https://pmc.ncbi.nlm.nih.gov/articles/PMC2978236/
- Pantazis & Olcese, cut-open oocyte technique (~1 nA at 3 kHz) —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11549981/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC4145744/
- Whole-cell patch-clamp recording and parameters — https://pmc.ncbi.nlm.nih.gov/articles/PMC10133435/

### PDFs that could NOT be retrieved (2026-07-27)

Attempted and blocked, so only the verified metadata and the quoted numbers are held here, no PDF:
`levis1993quartz`, `baumgartner1999tevc` (PMC serves these free to a browser but blocks scripted
download via Cloudflare; not in the PMC Open Access subset, so no bulk route either);
`stefani1998cutopen` and `levis1998lownoise` (Methods in Enzymology, paywalled);
Pantazis & Olcese chapter (host timed out). Publisher endpoints at cell.com returned 403.
Retrieve by hand from a browser if the PDFs are wanted beside the others.
