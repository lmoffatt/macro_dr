# Expression levels per recording configuration: macroscopic current and channel number

Collected 2026-07-27, seven parallel literature lanes plus one gap-fill pass, every claim re-verified
against the PDF by a second reader. Same rule as `SOURCES.md`: each number below is either **sourced**
(with the locator) or **derived** (with the arithmetic shown). Nothing here is a guess left unlabelled.

This file **supersedes the thin table** in `SOURCES.md` §"Currents, sourced", which had three rows and
two of them unattributed, and it answers that section's own admission that the macropatch box "is the
weakest of the three boxes and should be replaced with a sourced figure."

**Symbols.** `N_ch` channel number in the recorded membrane; `i` unitary (single-channel) current;
`γ` single-channel conductance, so `i = γ · V_drive`; `P_open` open probability; `I` macroscopic
current; `C_m` membrane capacitance; `A` membrane area; `Q_max` saturating gating charge; `e`
elementary charge, 1.602177e-19 C; `EC50` agonist concentration giving half-maximal response;
`n` number of cells or patches; SEM standard error of the mean, SD standard deviation. Two abbreviations
earn their place because they are the field's own keys into this literature: **TEVC** two-electrode
voltage clamp, and **NSFA** non-stationary fluctuation analysis (variance against mean current fitted
to `σ² = i·I − I²/N_ch + σ_B²`, which returns `i` and `N_ch` together with a baseline variance `σ_B²`).

**The forward and reverse maps**, used throughout:

    I = N_ch · i · P_open              N_ch = I / (i · P_open)

so only the product `i · P_open` matters for converting a current into a count. The figure's model
value is `1 pA × 0.5 = 0.500 pA per channel`. The one measured value for P2X2 in the same preparation
is `0.80 pA × 0.695 = 0.556 pA per channel` (Sattler 2020), an 11 % agreement, which is worth stating
because the whole reverse map runs through that product.

**Two counts, never merged.** A count from fluctuation analysis, gating charge, antibody binding or
direct single-channel resolution is evidence. `I/(i·P_open)` is arithmetic on top of an assumed `i` and
`P_open`. They are kept in separate columns in §2 and separated again in §3, and the two do not measure
the same thing even when both are right: see §4.5 flag 4.

**Exclusion, absolute.** The project's own P2X2 work supplies no number anywhere in this file. See §5b.

---

## 1. Headline table

One row per configuration. Currents are macroscopic per record, `N_ch` per recorded membrane.

| configuration | current, low end | current, high end | N_ch, low end | N_ch, high end |
|---|---|---|---|---|
| excised patch, outside-out | **5 pA** | **2 nA** | **10** | **300** (largest single patch ~860) |
| excised patch, inside-out | single-channel level, tens of pA | no macroscopic ligand-gated record found | **1** | **25** |
| MACROpatch (excised, 0.5–3 MΩ tip) | **6.5 pA** | **12 nA** (analysed band 2.5–6 nA) | **~90** | **2e4** |
| whole cell (HEK293, CHO, tsA201) | **100 pA** | **5–7 nA** clamp-controlled, 31 nA achieved | **~1e3** | **2e4** (6e5 for a sub-pS channel) |
| Xenopus oocyte, TEVC | **3 nA** | **65 µA** (distorted above ~20 µA) | **8e5** | **1.2e8** functional |
| Xenopus oocyte, cut-open vaseline gap | ionic 0.27–0.35 µA | ionic amplitudes essentially unreported | **1.5e9** | **4e9** (surface count only) |

### Basis of each end

| edge | value | which source, which method |
|---|---|---|
| outside-out, `I` low | 5 pA | Khanra 2026 Methods, stated acceptance floor: patches below 5 pA excluded |
| outside-out, `I` high | 2 nA | Ivica 2022 *J Physiol* measured mean 1160 ± 630 pA (n=6, ~1.8 nA at +1 SD); Schoppa & Sigworth reject records above 2.5 nA for series resistance |
| outside-out, `N_ch` low | 10 | Han 2024 NSFA, source's own "~10–50"; Mørkve & Hartveit lowest patch 10.3, NSFA. 1–3 channels reachable by design (Ding & Sachs 2002) |
| outside-out, `N_ch` high | 300 | Rook 2020 NSFA per-patch mean ~390 (individual patches ~150–860); Ivica derived 159. 300 is a mean-of-patches edge, ~860 a single-patch edge |
| inside-out, `N_ch` low/high | 1 / 25 | Geng 2023, 1 channel by construction (8–12 MΩ pipettes); Mauerer 1998, source's own "2 to ~25 K+ channels" per seal over 559 seals |
| macropatch, `I` low | 6.5 pA | Csánády 2005, CFTR in an oocyte macropatch, read against the panel's own 2 pA bar |
| macropatch, `I` high | 12 nA | Geng 2023 Fig 3A read against its own 10 nA bar, +240 mV. Analysed band from three labs' own rejection rules: >2.5 nA (Schoppa), ≤5 nA (Wang & Brenner), <15 mV error at 1.0–2.5 MΩ = 6 nA (Ma). Horrigan reports >20 nA as "too large to measure" |
| macropatch, `N_ch` low | ~90 | Alvarez 2002, NSFA on an oocyte macropatch, N = 94 ± 2.8 |
| macropatch, `N_ch` high | 2e4 | Liu 2016, patch-clamp fluorometry plus fluctuation analysis, "up to 20,000 per macropatch"; Alvarez Shaker 6000 ± 291 |
| whole cell, `I` low | 100 pA | Shi 2019 selects cell lines at −100 to −250 pA on purpose; Montnach 2021 lowest of 52 cells 391 pA; Li 2017 acceptance floor 500 pA over 42,272 wells |
| whole cell, `I` high | 5–7 nA | Wang & Brenner ≤5 nA, Eltokhi <6 nA, Montnach "larger than 7 nA should be prevented or discarded". Wand 2025 records 31 nA at −140 mV and states the voltage error "was not neglectable" |
| whole cell, `N_ch` low | ~1e3 | set by `I_min`, not by expression: 500 pA / 0.5 pA = 1000. The lowest sourced count is 1180 (Ivica derived) |
| whole cell, `N_ch` high | 2e4 | Stelmashenko 796 ± 149 pA/pF × 10–25 pF, derived 8.7e3–2.3e4; Ivica 8.6 nA derived 1.7e4; Khadra + Sattler 5036. Del Core & Mirams count 8.8e4–6.3e5 hERG per CHO cell, which needs γ = 0.92 pS and so sits outside the 1–10 pA family |
| oocyte TEVC, `I` low | 3 nA | Jackson 2025 Methods, stated range 3–200 nA across oocytes for human α3β2 nAChR, a poor expresser |
| oocyte TEVC, `I` high | 65 µA | Fujiwara & Kubo 2004 Fig 3C, top of 102 oocytes. Baumgartner: records above ~20 µA "are likely to be distorted"; 100 µA attainable |
| oocyte TEVC, `N_ch` low | 8e5 | Fujiwara & Kubo bottom of their own 1:1 to 1:1000 cRNA dilution series, derived at the ramp's −75 mV |
| oocyte TEVC, `N_ch` high | 1.2e8 | Fujiwara top, derived. Firsov counts 1.6–1.8e8 **surface** channels by antibody binding, which is 3.4e6–1.4e7 functional at Firsov's own P_open |
| cut-open, `I` | 0.27–0.35 µA | Savalli 2021 Fig 1A ionic trace at +20 mV, still rising at the end of the sweep |
| cut-open, `N_ch` | 1.5e9–4e9 | Rodriguez 1998, gating charge: Q_max = 3.2 nC over 4.83 e per channel gives 4.1e9; over the conventional 13 e for Shaker gives 1.5e9 |

### Drop-in form for the PREP table (`figure_6.Rmd:412`)

Only the four columns this file owns. `σ` and `τ` are `SOURCES.md`'s business, not this file's.

    #  name                          n_lo    n_hi    Imin(pA)  Imax(pA)
    "excised outside-out patch",       10,    300,        5,      2e3
    "whole cell",                    1e3,    2e4,      5e2,      6e3
    "oocyte, two-electrode",         8e5,    1e8,      3e3,      2e7

Three differences from the file as it stands, all in the direction of a sourced edge replacing a
placeholder, and none applied here (this file does not edit `figure_6.Rmd`):

1. **patch `Imin` 10 pA → 5 pA.** The 10 is flagged unsigned in the code. Khanra 2026's outside-out
   acceptance floor is 5 pA, stated in Methods. Consequence: `nmin_star` falls from 20 to 10, which
   equals `n_lo`, so the top-left diagonal vanishes, `i` returns to 1 pA at the left edge, and the
   patch's worst `r` returns from 0.63 to 2.5. This runs against the argument the figure makes.
2. **whole cell `Imin` 100 pA → 500 pA, `Imax` 5 nA → 6 nA.** Li 2017's >500 pA acceptance is measured
   over 42,272 wells; Eltokhi's 6 nA is the analysis ceiling at <5 mV residual error after 90 %
   compensation. Consequence: `nmin_star` rises to 1000, so `i` is forced to 3.3 pA at `n_lo` and the
   top-left corner drops about 11-fold. This helps the argument.
3. **oocyte `n_lo` 1e5 → 8e5.** 1e5 is unsourced and unreachable at `i` = 1 pA (1e5 × 0.5 pA = 50 nA,
   below the row's own `I_min`). 8e5 is Fujiwara's own dilution floor. Consequence: the top-left corner
   cut disappears and the oocyte's worst `r` goes from 22 to about 88.

Two stale strings in `figure_6.Rmd`, flagged only. Line 359 points at this file, which did not exist
until now. Line 358 says "patch 10-300, whole cell 300-5000, oocyte 1e5-5e7", while the table at line
414 already reads `20, 1e4` for whole cell; the comment was not updated with the table.

---

## 2. Detail tables

Only CONFIRMED claims. Everything that failed verification is in §2h. "Sat?" is whether the agonist
was saturating, which decides whether a derived `N_ch` is a floor or an estimate. Voltages are as the
paper states them.

### 2a. Excised patch, outside-out

| channel | expression system | current | N_ch, and the method | V | Sat? | n | citation |
|---|---|---|---|---|---|---|---|
| Drosophila NMJ iGluR + Neto (6 subunit/splice combinations) | HEK293T, transient, 3 d at 30 °C | fitted peak 52.4 ± 41.1 to 222.0 ± 307.8 pA (SD > mean in 5 of 6 groups) | **~10–50, NSFA**, per-patch values plotted in Fig 2E (~10 to ~45); γ 149–164 pS by NSFA, 169–175 pS directly; P_open 0.1–0.6 | −60 mV | 10 mM glutamate, **not declared saturating** | 7–8 per construct; 3–5 for NSFA | Han 2024 *J Physiol* 602:7043 |
| zebrafish α1 GlyR; rat α1β GlyR | HEK293, transient, 5 % receptor plasmid + 75 % empty pcDNA3 to hold expression down | 1160 ± 630 pA (α1, symmetric 131 mM Cl⁻); 348 ± 200 pA (α1β) | none in source. γ = 73.8 ± 6.4 pS measured in outside-out patches, same solutions | −100 nominal, −102 mV junction-corrected | 3 mM glycine, 2 ms jump, ~16× EC50 | 6 (α1), 7 (α1β); 4 for γ | Ivica 2022 *J Physiol* 600:333 |
| chicken ASIC1 (cASIC1), WT and L414A | HEK293, piezo-driven fast perfusion | ~200 pA (WT), ~500 pA (L414A), read off the Fig 2B axes; average of ~70 sweeps | **~390 WT (individual ~150–860), ~305 L414A, NSFA**, plotted in Fig 2C on a 0–1000 axis; peak P_open 0.86 ± 0.02, γ 10 ± 1 pS | −60 mV | pH 5, saturating; P_open measured, 0.86 | 5 per construct | Rook 2020 *eLife* 9:e51111 |
| GluA2(Q) AMPA ± TARP γ8, γ2 | HEK293T/17, transient | tens of pA (Fig 2A inset axes to 50 and 30 pA); 30–200 sweeps averaged | none printed. γ 17.3 ± 4.4 to 32.1 ± 7.6 pS; P_open,peak 0.33 ± 0.22 to 0.77 ± 0.15, both NSFA | −60 mV | 10 mM glutamate, 200 ms, saturating | 8–13 for γ and P_open | Coombs 2022 *Mol Pharmacol* 101:343 |
| adult mouse muscle nAChR (αβδε) | tsA 201, transient calcium phosphate | ~131–133 pA, **average of 25 sweeps** | **~70 (65–90), maximum-likelihood fit of the whole trace with N_ch free**; `i` = 3.4 pA assumed, not fitted | −40 mV | 1 mM carbachol, piezo, <50 µs exchange | 2 patches | Milescu 2005 *Biophys J* 88:2494 |
| GluK2, GluK2/GluK5 kainate receptors | HEK293, transient | **acceptance floor: outside-out patches below 5 pA excluded** (lifted whole cells below 20 pA); up to 10 sweeps averaged | none | −70 mV | dose-response 0.3–300 µM glutamate | not resolvable for the patch subset | Khanra 2026 *Nat Commun* 17:3789 |
| human α3β4 nAChR, 1:9 and 9:1 subunit ratios | HEK293, transient, ratios skewed on purpose | single-channel: 2.6 ± 0.1 pA (26 pS) and 3.9 ± 0.04 pA (39 pS) | 1 channel resolved, by design | −100 mV | **5 µM ACh, far below EC50 (138–309 µM)**; chosen low to isolate openings | 7 and 5 | Krashia 2010 *PLoS ONE* 5:e13611 |
| rat P2X2 (rP2X2-3T construct) | HEK293, transient; 7–20 MΩ pipettes | single-channel: −3.51 ± 0.2 pA | at least 3 (O2, O3 levels resolved); no count | −120 mV | **1 µM ATP, subsaturating on purpose** | 20 | Gasparri 2019 *J Gen Physiol* 151:e201912347 |
| rat P2X2 | Xenopus oocyte, cRNA cut from 50 to 25 ng and temperature lowered to reduce expression | up to 11.9 pA (three levels open) | **1–3, direct counting of simultaneous open levels** plus a binomial check; 4.0–4.1 pA per level | −120 mV | **10 µM ATP ≈ EC50 (9.8 ± 0.8 µM)** | 7 multi-channel + 2 single-channel patches (15 in the lifetime analysis) | Ding & Sachs 2002 *BMC Neurosci* 3:17 |
| GABA_A, somatic extrasynaptic, A17 amacrine cell | native, rat retinal slice | 154 ± 119 pA (range 24–344), average of 5–30 repetitions | **134 (range 24–291), NSFA**, "available channels"; γ 24.7 ± 1.3 pS; max P_open 0.65 ± 0.05 | not stated for patches; recovered as −60 mV from the paper's own 1.38 pA / 23.0 pS | 3 mM GABA, 2–3 ms; P_open measured, 0.65 | 9 | Beltrán-Matas 2022 *Eur J Neurosci* |
| glycine receptor, rod bipolar axon terminal | native, rat retinal slice | **93.4 ± 9.6 pA (range 44–171)** | **28.1 (range 10.3–70.3), NSFA**; γ 64.2 ± 2.5 pS; max P_open 0.88 ± 0.01 | −60 mV | 3 mM glycine, ~1 ms; P_open measured, 0.88 | 22 for N_ch, 13 for the peak | Mørkve & Hartveit 2009 *J Physiol* 587:3813 |
| Na_V, basket cell axon and soma | native, rat hippocampal slice | not printed. Conductance density 31.9 ± 3.7 (soma), 310.7 ± 32.7 (proximal axon), 574.3 ± 120.8 pS/µm² (distal axon); traces are averages of 5–100 | **density 2.6 / 25.0 / 46.1 channels/µm², NSFA**; γ 12.5 ± 1.0 pS, max P_open 0.66 ± 0.05. Derived N_ch = density × area, 5.5–756 over 4.9–32 MΩ | test pulse to 0 mV from −50 mV; reversal 70 mV | voltage-gated; the density is already P_open-corrected | 5 for NSFA, 7 axonal + 9 somatic for the area calibration | Hu & Jonas 2014 *Nat Neurosci* 17:686 |
| GABA_A, AII amacrine cell (**nucleated** patch, a large somatic patch) | native, rat retinal slice | 47.3 ± 19.7 pA (range 19.9–81.1), average of ~36 trials | **68.2 ± 29.0 (range 26.5–110), NSFA**; γ 23.2 ± 2.8 pS; max P_open 0.56 ± 0.06 | −60 mV, E_Cl ≈ 0 | 3 mM GABA, ~2 ms; the authors attribute the low P_open to slower exchange in a larger patch | 7 for N_ch, 11 for the peak | Beltrán-Matas 2023 *Front Ophthalmol* 3:1134765 |
| Piezo1 (**cell-attached**, not excised; kept for the density) | native endogenous, Neuro2A | I_max = N·i ≈ 3.4 pA average | **3.5 ± 3.1 channels/patch (n = 281 patches, 33 with zero), direct counting** as I_max/i with i = 0.98 ± 0.04 pA from 35 one-channel patches; **density 1.75 channels/µm²** on a 2 µm² imaged dome, ~45/µm² when Piezo1 is overexpressed, ~100/µm² in the most extreme patch | +60 mV (outward, chosen because inactivation precludes I_max at negative potentials) | pressure-clamp, ramp peak at saturating stimulus | 281 patches | Lewis & Grandl 2021 *eLife* 10:e70988 |
| voltage-gated Na⁺ and K⁺ (**nucleated** somatic macropatch) | native adult neuron, brain slice | a few hundred pA to ~0.65 nA, read against the 300 pA bar | **none stated anywhere in the chapter**; it works in pA/pF and pA/µm² throughout | V_hold = −90 mV, 10 mV/30 ms steps | voltage-gated | not stated | Tamagnini 2020 *Methods Mol Biol* 2188:229 |

### 2b. Excised patch, inside-out

Thin, and the reason is structural rather than accidental. Inside-out exposes the cytoplasmic face, so
a fast extracellular agonist jump is not available to it. The macroscopic inside-out literature is
voltage-gated and calcium-gated, and it uses large tips, so it appears in §2c instead.

| channel | expression system | current | N_ch, and the method | V | Sat? | n | citation |
|---|---|---|---|---|---|---|---|
| inwardly rectifying ATP-sensitive K⁺, basolateral membrane | native, dissociated axolotl proximal tubule | single-channel level; whole-cell reference for the same cells ~240 pA at 0 mV, 6 nS | **2 to ~25 per seal, direct counting of resolved levels** over 559 seals (98.6 % had activity; single-channel patches only 0.7 %). Authors separately derive ~4 per patch from density × area (2.23/µm² × 1.77 µm²) | outward conductance extrapolated to 0 mV | not agonist-gated | 559 seals | Mauerer 1998 *J Gen Physiol* 111:139 |
| human Slo1 / BK (KCa1.1) | Xenopus oocyte, 0.1–18 ng cRNA; 8–12 MΩ pipettes | single-channel | **1 per patch by construction**; γ 312 ± 4 pS (WT), 245 ± 6 (hybrid), 190 ± 14 (mutant), all at +100 mV in symmetric 160 mM KCl | +100 mV | voltage- and Ca-gated | 7 / 25 / 8 channels | Geng 2023 *J Gen Physiol* 155:e202213302 |

### 2c. MACROpatch, called out separately

Tips 0.5–3 MΩ, 3–20 µm diameter, against 4–32 MΩ and 0.5–2 µm for a conventional patch. That is a
20–200× step in membrane area (§4.1), which maps one-to-one onto `N_ch`, so the two must never be
pooled. Excision direction is stated per row because several of these are cell-attached.

| channel | expression system | current | N_ch, and the method | V | Sat? | n | citation |
|---|---|---|---|---|---|---|---|
| human Slo1 / BK, WT and 1:1 G375R+WT | Xenopus oocyte, 0.1–18 ng cRNA, inside-out excised | **12.25 nA** (WT) and 11.1–11.4 nA (1:1), read against the panel's own 10 nA bar; "mean response" over an unstated number of sweeps; minus P/4 subtraction | source's own: **"many hundreds to thousands"**, no method given. Derived: 164 at P_open = 1, 178 at the paper's own G/G_max = 0.920, using its own 312 pS | −80 mV holding, steps to +240 mV; symmetric 160 mM KCl, free Ca²⁺ <0.01 µM | voltage-gated, G/G_max = 0.92 at +240 mV | 6 (WT), 12 (1:1) patches | Geng 2023 *J Gen Physiol* 155:e202213302 |
| BK (KCa) in oocyte | Xenopus oocyte, macropatch | plateau ~1.6 nA, **mean of 256 traces** | **94 ± 2.8, NSFA**; fitted i = 21 ± 4.6 pA (175 pS at 120 mV driving force); printed P_open,max 0.8 | +120 mV from 0 mV holding | 100 nM internal Ca²⁺, **not saturating for BK**; activation is voltage-driven | 1 patch (a worked tutorial example) | Alvarez 2002 *Adv Physiol Educ* 26:327 |
| Shaker H4 Δ(6–46) in oocyte | Xenopus oocyte, macropatch | plateau ~7.9–8.1 nA, **mean of 256 traces** | **6000 ± 291, NSFA**; fitted i = 1.4 ± 0.03 pA (12 pS); printed P_open,max 0.8 | +120 mV from −100 mV | voltage-gated | 1 patch | Alvarez 2002 *Adv Physiol Educ* 26:327 |
| Shaker H4ir (**cell-attached**) | Xenopus oocyte, 1.5–3 MΩ | not printed numerically; Fig 2A carries a 60 pA bar, measured ~90 pA at +26 mV; gating currents 20 sweeps averaged | **2250 in one patch, NSFA** from 50–60 traces at +70 mV, cross-checked against idealised single-channel sweeps in the same patch. P_max = 0.79 measured | +70 mV for NSFA; −80 mV holding | voltage-gated; P_max measured | 1 patch quantified | Islas & Sigworth 1999 *J Gen Physiol* 114:723 |
| Kv2.1 (**cell-attached**) | Xenopus oocyte | Fig 2A bar 1 nA | **967, 1250, and 1500 (legend) or ~1650 (text) for the Fig 6 patch: three patches, NSFA**; P_max = 0.69 | +70 mV | voltage-gated | 3 patches | Islas & Sigworth 1999 *J Gen Physiol* 114:723 |
| Shaker H4ir (**cell-attached**, 5 µm pipette) | Xenopus oocyte, 50 nl of 2–3 mg/ml cRNA | not printed; **derived I = N·i·P_open = 4.10 nA**, corroborated by the Fig 6 inset x-axis reaching ~4200 pA | **11,800 in one patch, NSFA** from 257 traces; i = 0.46 pA, γ 13.8 pS, P_open 0.7551, all fitted | −20 mV test pulse | voltage-gated; P_open,max 0.71 ± 0.02 measured | 1 patch for the fit; 12/8/7 for γ at three temperatures | Rodriguez 1998 *J Gen Physiol* 112:223 |
| mSlo1 / BK, gating charge | Xenopus oocyte AND HEK 293 macropatches, pooled by the paper; ~50 ng cRNA | gating currents only, no macroscopic ionic amplitude anywhere in the paper; signal-averaged over ≥8 pulses | source's own **Q_max = 1–30 fC per patch** and its own counting formula `N = Q_Tfast/4z_J`, but it never prints a resulting count. Derived: **1.4e3–7.9e4** using 2.36–4.4 e per channel. Patch C_m 0.25–1 pF, so A = 25–100 µm² | 1 s ramp −160 to +200 mV | 0 Ca²⁺; Q-V saturated | 15 patches for the C_g-V fits | Horrigan & Aldrich 1999 *J Gen Physiol* 114:305 |
| mSlo1 / BK, ionic | Xenopus oocyte, 0.5–5 ng cRNA, inside-out | **>20 nA "too large to measure"** at maximally activating voltages; one example tail <5 nA; averages of 4–8 pulses | source's own **"hundreds"**, twice, and separately "hundreds of channels are present and functional"; method `nP_omax = G_Kmax/g_K` plus a Poisson fit to the amplitude histogram of the same patch | −80 mV holding; tails −150 to −360 mV | 0 Ca²⁺ (0.8 nM); voltage-activated | not stated for the "hundreds" | Horrigan, Cui & Aldrich 1999 *J Gen Physiol* 114:277 |
| Shaker 29-4 | Xenopus oocyte, ~3 ng cRNA, inside-out | **records above 2.5 nA rejected** to avoid series-resistance error; 10–100 sweeps averaged | none stated. Tips **3–10 µm (0.5–3.0 MΩ)** against 1–2 µm (4–10 MΩ) for single channels in the same paper | −93 mV holding | voltage-gated | rejection rule applied across the study | Schoppa & Sigworth 1998 *J Gen Physiol* 111:271 |
| mSlo1 / BK | Xenopus oocyte, inside-out | peak 5.1–7.0 nA over five holding-potential blocks of one patch; P/5 subtraction, no sweep averaging stated | none; NP_open only, by Gaussian fits to the total amplitude histogram. **Gating-current macropatches 10–20 µm diameter** | 1.6 ms steps to +160 mV | **300 µM Ca²⁺, saturating for BK**, so close to N·i·P_omax | 1 patch | Zhou & Lingle 2014 *J Gen Physiol* 144:415 |
| mSlo1 / BK | Xenopus oocyte, 0.5–2 MΩ electrodes | not printed (records taken at very negative voltages, so deliberately tiny) | source's own **"hundreds of channels"**, with typically one open at a time; R210C patches **selected to contain 2–10 channels**. Method: NP_open from summed open time over the resolved levels | −50 to −200 mV | 0 and 10 mM Mg²⁺, 0 Ca²⁺ | 6–14 patches per point | Chen, Geng & Magleby 2011 *J Gen Physiol* 138:593 |
| mSlo1 / BK | Xenopus oocyte, excised macropatch | not printed | source's own **"excised macropatches containing hundreds of channels"**, by limiting slope down to P_open ≈ 1e-6 combined with the macroscopic G-V in the same patch | limiting slope at strongly negative V | various [Ca²⁺] | not stated | Sun & Horrigan 2022 *Sci Adv* 8:eabq5772 |
| mSlo1 / BK, S1–S4 mutants | Xenopus oocyte, 0.5–50 ng cRNA | not printed; Fig 2A read ~9–10 nA at +240 mV. **Series-resistance error <15 mV at 1.0–2.5 MΩ, i.e. 6 nA** | source's own "hundreds of channels are present and functional"; macroscopic G_K plus single-channel currents in the same patch | tails at −80 mV after 30 ms pulses | 50 µM Ca²⁺ as the saturating reference | not stated | Ma, Lou & Horrigan 2006 *J Gen Physiol* 127:309 |
| BK α and α+β4 (**HEK293, not oocyte**) | HEK293, excised inside-out | **"currents 5 nA or less were used for steady-state G-V"** | none | various | various [Ca²⁺]ᵢ | not stated | Wang, Rothberg & Brenner 2006 *J Gen Physiol* 127:449 |
| mHCN2-EGFP | Xenopus oocyte, 40–50 ng cRNA | not printed. **Derived I = N·i·P_open = 1.96 nA** (activation) and 474 pA (tail); traces are single sweeps out of 100 | **8715 (activation) and 4707 (tail) in the SAME patch**, fluctuation analysis with variance from neighbouring-trace differences; i = 0.298 pA at −130 mV and 0.114 pA at −40 mV; P_open 75.5 % and 88.4 %, all fitted. Separately **"up to 20,000 per macropatch"** | −130 mV (activation), −40 mV (tail) | saturating cAMP (3 µM); P_open measured | 1 patch for the pair; 3–4 per construct | Liu 2016 *J Gen Physiol* 148:65 |
| HCN1, HCN2, HCN4 | Xenopus oocyte macropatch | not stated | source's own field-norm statement: macropatches "containing typically hundreds of channels" | not applicable | not applicable | review-style statement | Benndorf 2025 *PNAS* 122:e2422533122 |
| CFTR Cl⁻, WT and severed constructs | Xenopus oocyte, ~2 MΩ macropatch pipettes | **~6.5 pA** in saturating ATP, read against the Fig 3 "2 pA" bar; continuous records at 1 kHz | source's own **"tens or hundreds"** of channels; no method | pipette +80 mV, i.e. V_m = −80 mV | 2 mM MgATP, saturating; 300 nM PKA to phosphorylate | not stated for the amplitude | Csánády 2005 *J Gen Physiol* 125:43 |
| spHCN, gating currents | Xenopus oocyte, **giant** inside-out patch, 180–250 kΩ pipettes | not stated; P/4 subtraction | none. **All Q-V curves normalised, so no Q_max in coulombs and no count is reconstructible** | 5–10 kHz filter | 100 µM cAMP, saturating | not stated | Ryu & Yellen 2012 *J Gen Physiol* 140:469 |
| not a channel study; giant-patch geometry | rat RBL secretory line | not applicable | **patch size 2–4 pF at 10–15 µm pipette diameter**, measured as capacitance. The only measured giant-patch area in the set, and the source of the dome factor in §4.1 | not applicable | not applicable | not applicable | Wang & Hilgemann 2008 *J Gen Physiol* 132:51 |
| generic, instrumentation manual | not applicable | worked example, 100 pA macropatch current against 100 MΩ cell input resistance gives 10 mV of cell-voltage change | definitional: **"patches of membrane containing tens to hundreds of ion channels"**; macropatch tips 5–10 µm against 0.5–2 µm for single channels; giant patch 12–40 µm | not applicable | not applicable | not applicable | Axon Guide ch. 5 |

### 2d. Whole cell

| channel | expression system | current | N_ch, and the method | V | Sat? | n | citation |
|---|---|---|---|---|---|---|---|
| rat P2X2a; P2X2b | HEK293, transient Lipofectamine 2000, C_m ~10 pF | 3.0–3.5 nA (P2X2a) and 4.5–6.6 nA (P2X2b), four and five successive applications to ONE cell each; continuous single-sweep records | none in source | −60 mV | 100 µM ATP, called saturating for wild-type P2X2 by Sattler 2020 on the same platform | 1 cell per receptor for the amplitudes | Khadra 2012 *J Gen Physiol* 139:333 |
| rat P2X2, wild type | Flp-In-T-REx 293 stable inducible, induced 2–6 h to keep expression LOW | absolute I_max not reported (all responses normalised; the string "nA" appears nowhere) | none. **`i` = 0.80 ± 0.03 pA measured cell-attached, and P_open,max = 0.695 ± 0.009 by stationary fluctuation analysis at the peak of a saturating response.** This is the one measured `i · P_open` for the project's channel from an independent lab | −50 mV both configurations | **100 µM ATP stated saturating** | 6 for both; 13 pooled with H319K | Sattler 2020 *Sci Rep* 10:21751 |
| rat P2X2, wild-type monomer; concatenated trimer | HEK293, transient | **current DENSITY 796 ± 149 pA/pF (monomer), 341 ± 39 pA/pF (trimer)**, not an absolute current | none. Own γ = 22 ± 1.8 pS at −120 mV. Derived at C_m 10–25 pF: I = 8.0–15.9 nA, N_ch = 8.7e3–2.3e4 | −60 mV | ATP concentration for the monomer row **not stated**; the paper's own EC50 is 16 ± 2.0 µM | 4 (monomer), 23 (trimer) | Stelmashenko 2012 *Mol Pharmacol* 82:760 |
| human P2X2; P2X2/3 heteromer | Flp-In-T-REx 293 stable inducible, 16–24 h after induction | **31.0 ± 8.6 nA** and 22.7 ± 4.0 nA; G353R 18.5 ± 2.8 and 14.1 ± 2.1 nA. 142 mM KCl on both faces, which is part of why they are this large | none | −140 mV (0 mV holding, 150 ms steps) | 100 µM ATP | not stated for wild type; 8 and 9 for G353R | Wand 2025 *Cells* 14:510 |
| human α3β4 nAChR | HEK293, transient calcium phosphate, 14–48 h | **I_max = 6.67 ± 1.63 nA**, a per-cell Hill asymptote then averaged | none in source; γ 26 and 39 pS from its own outside-out patches. Derived at 39 pS and −30 mV: 5.7e3 at P_open = 1, 1.14e4 at 0.5 | −30 mV | yes by construction (fitted asymptote); EC50 91.1 ± 10.7 µM | 8 | Krashia 2010 *PLoS ONE* 5:e13611 |
| zebrafish α1 GlyR | HEK293A, 2 % receptor plasmid + 78 % empty pcDNA3 | **I_max 4.3 ± 1.3 nA** (glycine); 5.8 ± 1.8 (AMS), 3.5 ± 0.3 (β-alanine), 1.1 ± 0.5 (taurine, a partial agonist). Hill asymptote per cell, mean ± SD | none. **max P_open = 0.96 ± 0.06 measured cell-attached in the same study**, so this row is one number short of a count | −40 nominal, −50 mV corrected | 100 mM glycine explicitly saturating | 8, 9, 6, 6; 8 patches / 92 clusters for P_open | Ivica 2022 *eLife* |
| zebrafish α1 GlyR (same lab, same cells as the outside-out row in §2a) | HEK293, 5 % receptor plasmid | **I_max = 8.6 ± 1.2 nA** | none. Derived at the same paper's 73.8 pS and −50 mV: **1180–2400** | −40 nominal, ~−50 mV corrected | 10 or 30 mM glycine, explicitly saturating | 8 | Ivica 2022 *J Physiol* 600:333 |
| human GlyR α3L185L and six point mutants | HEK293T, FuGENE HD, 1 µg DNA, ~1 d | 1.595 ± 0.225 nA at 1000 µM glycine; 1.621 ± 0.300 nA at 0.1 mM in a second set; mutants 0.563–1.595 nA | none | −50 mV | plateau over 100–1000 µM (80.7 / 90.9 / 93.0 % of each cell's own maximum) | 6, 12, and 5–7 per mutant | Harnau 2025 *Sci Rep* 15:32686 |
| α1β3γ2L GABA_A; binary α1β3 | HEK293T, transient, 48 h, lifted cells, C_m 8–12 pF | **7038 ± 302 pA** (ternary), 1204 ± 235 pA (binary), 931 ± 58 pA (γ2L P302L). Cells pre-selected at ~5 nA | none. Own γ 25.34 ± 1.77 pS main level (cell-attached, +80 mV). A derived floor is possible but the ohmic extrapolation from +100 mV to −20 mV is refuted by the paper's own voltage series, so it is reported in §2h instead | −20 mV | 1 mM GABA, ~130× the EC50 of 7.50 ± 0.77 µM | 5–6 | Hernandez 2017 *eNeuro* 4:e0251-16 |
| α1β2 and α3β2 GABA_A | HEK293, transient | 1.6 ± 1.0 nA (α1β2), 1.3 ± 0.2 nA (α3β2); one sweep per cell | none | −60 mV | 30 mM GABA, 28 s | 17, 15 | Olander 2020 *PLoS ONE* 15:e0234080 |
| human P2X3, slow-desensitising mutant | HEK293, transient FuGENE6, 18–30 h | ~0.93–0.95 nA | none | **not stated anywhere in the paper**; "mV" occurs nowhere in the text | 10 µM free ATP, explicitly saturating | 1 cell | Li 2019 *eLife* 8:e47060 |
| panda P2X7 | HEK293, transient, 1 µg | ~2.75–2.85 nA at 1 mM ATP. Separately, **the 27th of 27 successive 100 µM applications to ONE patch reaches ~3.3 nA against ~0.8 nA for the first**, a 4.2-fold run-up over minutes | none | −60 mV | 1 mM ATP ≈ 8× the paper's EC50 of 122 µM | 1 cell each | Karasawa 2016 *eLife* 5:e22153 |
| GluA3i-G AMPA | HEK293 stable lines, expression deliberately kept low | **selection window −100 to −250 pA** at saturating agonist, chosen so that patches from the same line hold one channel ~40 % of the time | the calibration itself: ~200 pA whole cell corresponds to ~1 channel per patch | −60 mV | 10 mM glutamate + cyclothiazide, saturating | stable lines | Shi 2019 *J Gen Physiol* 151:156 |
| human Na_V1.7 | CHO-S15 stable line, automated planar patch (SyncroPatch 768PE) | **mean 1.66 ± 0.01 nA over 42,272 recordings; 82 % above 0.5 nA; acceptance floor >500 pA.** C_m 25 ± 0.1 pF, 95 % between 10 and 40 pF | none | step to −10 mV from −120 mV | voltage-gated | **42,272 recordings over 16 chips** | Li 2017 *PLoS ONE* 12:e0180154 |
| Na_V1.5 | COS-7, transient (manual); HEK293 stable (automated) | **391 pA to 17.8 nA across 52 transiently transfected cells**, a 46-fold spread in one lab with one construct. Ceilings: "current amplitudes larger than 7 nA should be prevented or discarded" (manual, residual R_s 2.3 ± 0.2 MΩ); 500 pA–3.5 nA window (automated); 10 nA for outward K⁺ at R_s ≤5 MΩ | none | −20 mV, holding −100 mV | voltage-gated | 52 cells | Montnach 2021 *Sci Rep* 11:3282 |
| human Na_V1.2 + β1 | HEK293, BacMam transduction with thymidine G1/S arrest | **3000 pA/pF peak density**, a mean of ~19–20 cells; ~30 nA at the paper's own 10 pF for untreated cells. **Analysis restricted to <6 nA** "to assure the fidelity of the voltage clamp", and virus was down-titrated to stay there. Only cells with <5 mV residual error after 90 % compensation kept | none. Whole-cell ON gating charge measured but used only as a relative surface-expression index (5-fold), never converted to a count | −100 to +50 mV, holding −120 mV | voltage-gated | ~19–20 cells | Eltokhi 2023 *Cell Rep Methods* 3:100559 |
| rat Na_V1.4 | HEK293 stable, IonFlux Mercury **20-cell ensemble plates** | **9.54 ± 0.51 nA per ensemble, which is a SUM over up to 20 cells**, i.e. ~0.48 nA per cell; exclusion floor 2 nA per ensemble ≈ 100 pA per cell. The paper states occupancy per hole cannot be verified | none | 17-pulse protocol | voltage-gated | 60 ensembles, up to 1200 cells | Lukacs 2021 *Front Pharmacol* 12:738460 |
| human Na_V1.5 (hH1α) | HEK293t, transient calcium phosphate | **74 ± 8 pA/pF**, and C_m 43 ± 1 pF, so ~3.2 nA. Superimposed single sweeps, capacitance and leak subtracted on line | none | test to −30 mV from −150 mV | voltage-gated | 27 for the density, **121 for C_m** | Xiao 1998 *PNAS* 95:2680 |
| human Na_V1.5 + rat β1 | HEK293, transient Polyfect, 24–36 h | ~3.55 nA representative; −P/4 leak subtraction | none | holding −130 mV | voltage-gated | 1 cell | Sokolov 2013 *Front Pharmacol* 4:78 |
| rat Kv2.1 | HEK293 stable | **415–418 pA/pF** (n = 17) and 3.85–3.87 nA in the representative cell, implying C_m = 9.3 pF | none | +80 mV from −70 mV | voltage-gated, G/G_max = 1.0 | 17; 1 for the trace | Dallas 2021 *Sci Rep* 11:8194 |
| hERG1a (Kv11.1) | CHO, stable, room temperature | 0.4–3 nA across the nine cells, and the ordering tracks the fitted count | **8.8e4 to 6.3e5 channels per cell (median 4.0e5), maximum-likelihood fit of a five-state stochastic gating model** with `η` free and `g_s` box-constrained from independent single-channel work to 0.919 pS. The paper warns that ignoring the non-conducting open state makes the classic `η = g/g_s` route recover only ~27 % of the true count | Beattie sinusoidal protocol, [K]ₒ = 4 mM | voltage-gated | 9 cells | Del Core & Mirams 2025 *Phil Trans R Soc A* 383:20240224 |
| Ca_V1.4 | HEK-293 | tail currents to ~1800 pA; **mean current of 250–500 traces** for the fluctuation analysis. Peak density 25.99 ± 4.18 pA/pF | **density 0.97 ± 0.153 channels/µm², NSFA**, with i = 4.56 ± 1.178 pA and P_open 0.80 ± 0.052. A count per cell is NOT derivable: see §2h | tail at −49 mV after +41 mV | voltage-gated, with 3 µM BayK8644 | 3 for the fluctuation analysis, 29 for the density | Heigl 2023 *Channels* 17:2192360 |
| HyNaC2/3/5, peptide-gated | HEK 293T and COS-7, native versus codon-optimised cDNA | **113 ± 9 pA/pF (HEK, optimised) against 40 ± 6 (native); COS-7 17 ± 9 against 2 ± 3.** At the paper's own ~20 pF that is ~2.3 nA | none; no γ stated, so no count | −70 mV | 1 µM Hydra-RFamide II, "a concentration that elicits a maximal response in oocytes" | 15, 12, 10, 12 | Bachmann 2021 *Channels* 15:198 |
| **endogenous** voltage-gated K⁺, the background a heterologous current must beat | HEK293, **untransfected** | **up to 16.8 ± 2.0 pA/pF at +80 mV**, i.e. ~230–290 pA in a 14–17 pF cell with no transfection at all. Rises with passage number (11.8 ± 1.7 at passage 20 to 20.8 ± 1.4 at passage 70) | none | +80 mV from −120 mV | not applicable | C_m 14.0 ± 0.3 (high plating density) and 17 ± 0.5 pF | Ponce 2018 *Physiol Rep* 6:e13663 |
| hERG / KCNH2 | Flp-In HEK293 stable, SyncroPatch 384PE | **no current minimum in the quality-control set at all**; the four criteria are seal >300 MΩ, **C_m between 5 and 50 pF**, R_s <20 MΩ, leak within ±40 pA | none | −120 mV leak step | voltage-gated | 2304 wells screened, 52.1 % pass | Ng 2021 *Biol Methods Protoc* 6:bpab003 |
| **native** TRPV1 and outward K⁺, mouse DRG nociceptors | acutely dissociated mouse lumbar DRG | capsaicin density mean ~−22 pA/pF with single cells to −151; representative trace ~2.65 nA. Outward K⁺ 108 and 177 pA/pF, traces ~3.9 and 6.2 nA | none | −70 mV | 100 nM capsaicin is **sub-maximal** | 25 and 70; 10 and 17 | Defaye 2024 *J Clin Invest* 134:e176474 |
| **native** K⁺ in rat hippocampal neurons | primary culture, 7–14 days in vitro | 182 pA/pF at +60 mV, representative trace ~2.31 nA, implying C_m ≈ 12.7 pF | none | +60 mV from −70 mV, with a 30 ms −50 mV prepulse | voltage-gated | 8 | Dallas 2021 *Sci Rep* 11:8194 |

### 2e. Xenopus oocyte, two-electrode voltage clamp

| channel | expression system | current | N_ch, and the method | V | Sat? | n | citation |
|---|---|---|---|---|---|---|---|
| rat P2X2 | Xenopus stage V oocyte, 50 ng cRNA, 2–5 d | **3150 ± 260 nA** (mean ± SEM across oocytes of a single 10 s application peak). Rat P2X1 in the same study 1600 ± 160 nA | none. Derived at Sattler's measured `i·P_open`: 5.7e6 | −60 mV | **100 µM ATP stated to evoke maximal currents** for both P2X1 and P2X2 | **47** | Werner 1996 *PNAS* 93:15485 |
| rat P2X1 | same | **142 to 3467 nA, a stated RANGE over 26 oocytes at one fixed 50 ng dose**, a 24-fold spread. The best single datum on how wide a fixed injection leaves the level | none; no γ in source, and the project's `i` is a P2X2 value | −60 mV | as above | 26 | Werner 1996 *PNAS* 93:15485 |
| rat P2X2, expression **titrated** | Xenopus stage V oocyte, ~50 nl of cRNA at serial dilutions **1:1 to 1:1000**, 1–2 d at 17 °C | **~0.8 to ~65 µA across 102 oocytes**; three representative ramps at 2.7, 21.8, 42.2 µA. Reproduced in 11 batches. Clamp-quality rule: "data with an error of over 2 mV from the command potentials were discarded" | none. Derived: 8e5 to 5.2e7 at P_open = 1, taking the ramp's own −75 mV | Fig 3A ramps −75 to +75 mV; the Fig 3C axis voltage is **not stated** | **100 µM ATP stated saturating** | 102 (Fig 3C), 68, 41 | Fujiwara & Kubo 2004 *J Physiol* 558:31 |
| human α3β2 nAChR | Xenopus oocyte (Ecocyte), 50.6 ng total mRNA, 6–9 d at 17–19 °C | **stated range 3 to 200 nA** across oocytes and ACh concentrations. The only stated across-oocyte range found in the oocyte lane, and the lowest expression in it | none anywhere in the paper | −60 mV | range spans sub- to fully saturating; ACh tested 100 nM to 100 mM, EC50 12.2 ± 1.7 and 264 ± 1.6 µM | not stated for the range; 14 and 12 for the dose-response | Jackson 2025 *Int J Mol Sci* 26:9506 |
| zebrafish α1 GlyR; α1/βb heteromer | Xenopus oocyte, **5 fmol** cRNA, 24–48 h at 18 °C | 2.25–2.27 µA at 1 mM glycine (saturated, the 500→1000 µM step adds 4 %); 2.15–2.20 µA at 200 µM; heteromer 1.83 µA at 200 µM. Single continuous records | none anywhere; no γ | −70 mV | **1 mM saturating** (EC50 112 ± 10 µM, Hill 2.7); 200 µM is ~0.8 of maximum | 1 oocyte per trace; the curves are n = 10 but normalised per oocyte | Ito 2020 *Sci Rep* 10:13999 |
| human P2X4-E307T | Xenopus oocyte, 25 ng cRNA, ≥36 h at 16 °C | **~0.46 µA** (one oocyte) and ~1.10 µA (a second), a 2.4-fold spread at the same concentration | none | −60 mV | **NO: 10 µM ATP is the paper's own EC50**, so a maximal response is roughly twice this | 1 per trace; 3–8 per dose-response point | Nagel 2025 *Nat Commun* 16:10367 |
| human α1β2γ2L GABA_A | Xenopus oocyte; cRNA amount and incubation **not stated**, deferred to a previous paper | **4.24–4.27 µA at the paper's own maximal-activation reference** (1 mM GABA + 50 µM propofol). The same cell gives only 0.61 µA to saturating 1 mM P4S, because P4S has a maximal P_open of 0.18 | none. The paper names non-stationary noise analysis as the alternative route to P_open,max and does not use it | −60 mV | **yes, this is the maximum-activation reference**; note P_open = 1 is an inherited assumption, not measured here | 1 cell for the traces, 5 for the curves | Germann 2025 *J Gen Physiol* 157:e202413644 |
| human Slo1 / BK, WT and 1:1 G375R+WT | Xenopus oocyte, 0.5–150 ng cRNA, 2–5 d at 18 °C; ND96 + 1 mM DIDS | **29.8–30.9 µA** (WT, +100 mV) and 34.4–34.6 µA (1:1, +60 mV), read against the panel's own 10 µA bar | none. A derived count is **rejected** here (§2h): the paper's own 312 pS is a symmetric-KCl number and the whole-oocyte record is in 2 mM external K⁺ | −80 mV holding; steps to +100 and +60 mV | voltage-gated, and **far from maximal**: the paper's own relative conductance is ~0.18 at the most positive step | 4 | Geng 2023 *J Gen Physiol* 155:e202213302 |
| human α3β4 nAChR (same lab and cDNAs as the whole-cell row in §2d) | Xenopus oocyte, 1:1, 1:9 and 9:1 cRNA | **2050 ± 301 nA** (1:1), 1270 ± 316 (1:9), 686 ± 213 (9:1); per-oocyte Hill asymptote then averaged | none in source. Derived at its own 26 pS and −70 mV: 1.13e6 at P_open = 1 | −70 mV | yes by construction | 5, 4, 4 | Krashia 2010 *PLoS ONE* 5:e13611 |
| α1β2 GABA_A, two expression vectors | Xenopus oocyte, 27.4 nl at 50 pg/nl per subunit (1.37 ng each), 18 °C, 2–14 d | **10.3 ± 1.1 µA (pGH19, day 5) and 14.7 ± 1.0 µA (pUNIV); by day 13, 7.5 ± 1.1 and 13.0 ± 0.9 µA.** Two knobs at a fixed dose: days after injection (peak at day 5, then a 27 % or 12 % decline) and the vector (1.4× at day 5, 1.7× at day 13) | none. Radioligand B_max was measured for the mammalian arm but reported only as a fold-change | −80 mV | 10 mM GABA, maximal | printed inside Fig 4c only | Venkatachalan 2007 *Pflügers Arch* 454:155 |
| rat Na_V1.4 α + Na_Vβ1 | Xenopus oocyte, 10 ng + 2.5 ng in 24–46 nl, 2–4 d | **1, 17, 25 and 43 µA: four single-oocyte examples, explicitly "not leak-subtracted"**; the population is split at 20 µA. The same paper states its own large TEVC values "are estimates, with the largest values possessing the greatest error" | none | peak inward under TEVC; ~−50 mV for the action-potential protocol | voltage-gated | 4 examples | Corbin-Leftwich 2018 *J Gen Physiol* 150:1583 |
| rat ENaC αβγ, FLAG-tagged | Xenopus oocyte, binding assayed 20–24 h after injection, then the same oocytes clamped | amiloride-sensitive I_Na; slopes 4.82 µA/fmol (low Na⁺) and 1.13 µA/fmol (high Na⁺), so ~1.4–4.3 µA and ~0.3–1.0 µA over the stated 0.3–0.9 fmol binding range. A steady current, no sweeps | **B_max = 4.9e8 antibody sites per oocyte by radioligand binding, and the paper's own division by three subunits gives 1.6–1.8e8 SURFACE channels per oocyte.** Independent of the current, which is the point. Its own summary: 125–550 sites, i.e. 42–167 channels, per pA. γ 5.5 pS in Na⁺; i = 0.6 pA at −100 mV | −100 mV | not agonist-gated. P_open is very LOW here (their 0.010–0.040 under their own stoichiometry), which is the whole finding | 12 oocytes per isotherm; 4 experiments in Fig 6A | Firsov 1996 *PNAS* 93:15370 |
| Shaker V478W and W434F, non-conducting | Xenopus oocyte, TEVC, K⁺-free external | **ionic <50 nA** for the non-conducting mutant, and a stated counterfactual that wild-type Shaker at the same density would be "in the low mA range" | **Q_max = 22 nC in one oocyte, which the paper converts to ~1e9 channels** at 13.6 e per channel, plus ">10^8 channels per cell" as the routine level. **Its own arithmetic gives 1.0e10, so it disagrees with itself 10-fold** (§4.5 flag 3) | −100 mV holding, steps to 0 mV; P/−4 subtraction, Q the average of the ON and OFF integrals | voltage-gated, Q saturated | 1 example oocyte | Kitaguchi 2004 *J Gen Physiol* 124:319 |
| not channel-specific: the clamp ceiling, and the oocyte's own area | not applicable | **currents of "100 µA or more" are attainable; "recordings of currents larger than ~20 µA are likely to be distorted".** The cause is potential differences inside the 1 mm cell, which series-resistance compensation does NOT fix. Bath series resistance ~100 Ω | not applicable. **C_m = 45.6 ± 5 mF/m² measured, "6.5 times larger than expected for a lipid bilayer"**, radius ~0.55 mm | not applicable | not applicable | theory plus measurements on suspensions and single cells | Baumgartner 1999 *Biophys J* 77:1980 |
| not channel-specific | not applicable | a virtual-ground bath electrode "allows currents as large as **several tens of microamps** to be recorded without losing control of the cell's membrane potential" | not applicable; the review states the relation as I = (E_m − E_rev)·N·P_open·γ | not applicable | not applicable | review | Papke & Smith-Maxwell 2009 *Comb Chem High Throughput Screen* 12:38 |
| not channel-specific; **oocyte area, second lab, second method** | native, uninjected oocytes | not applicable | **C_m = 166.4 ± 24.7 nF measured against 31.4 nF predicted for a smooth 1 mm sphere at 1 µF/cm², a 5.3-fold amplification**; "oocytes may have 5–10 times the area needed to enclose their geometric volume" | not applicable | not applicable | 11 oocytes, 2 frogs | Zhang & Hamill 2000 *J Physiol* 523:101 |
| Kv7.2/7.3 (KCNQ2/3), the worked example of a TEVC protocol chapter | Xenobus oocyte, 46 nl = 3 ng mRNA, 1–2 d at 18 °C | **NO current amplitude anywhere in the chapter.** No numeral followed by pA, nA, µA or mA in eleven pages, verified by exhaustive search; its one current figure has no current scale bar and no current axis. It discusses "large currents" three times and quantifies the threshold zero times | none; no fluctuation analysis, no charge, no counting | resting potential ~−50 mV; Fig 3 steps −80 to +40 mV | voltage-gated | protocol chapter, no data set | Guan, Chen & Zhang 2013 *Methods Mol Biol* 998:79 |

### 2f. Xenopus oocyte, cut-open vaseline gap, called out separately

| channel | expression system | current | N_ch, and the method | V | Sat? | n | citation |
|---|---|---|---|---|---|---|---|
| Shaker H4ir (ShH4ir), and W434F non-conducting | Xenopus oocyte, cut-open after 0.3 % saponin permeabilisation; series resistance <0.4 kΩ | **Q_off saturates at ~3.2 nC** at 22.3 and 11.7 °C, ~2.8 nC at 4.6 °C, read off the paper's own "Q_off (nC)" axis. Unsubtracted records, 100–300 ms test pulses | **n = 4 × 10^9 channels in the clamped membrane, gating charge**, from a two-step Boltzmann in which n is an explicit free parameter and the charge per channel is z1 + z2 = 4.83 e. Arithmetic check: 3.2 nC / (4.83 × 1.602e-19) = 4.13e9. **Their 4.83 e is well below the 12–13 e usually accepted for Shaker; at 13 e the same charge gives 1.5e9** | Q-V from about −125 to +42 mV, holding −90 mV | voltage-gated, Q-V saturated, so this is not a floor | 1 oocyte for the Q-V; the fitted z values rest on two temperatures of one experiment, not two oocytes | Rodriguez 1998 *J Gen Physiol* 112:223 |
| Shaker W434F, gating currents | Xenopus oocyte, 50 ng cRNA, 3–6 d at 16 °C; Dagan CA-1B | **peak ON gating current ~4.30 µA**, read against the panel's own 2 µA bar. **Composition stated and it is an average of four traces**, so the plotted amplitude is one sweep's | **none. Every Q-V is normalised; "nC" and "fC" appear nowhere, so no count is reconstructible** | steps to about +40 mV; K⁺-free external | voltage-gated, Q-V saturates near 0 to +20 mV | 1 for the trace; 3 for this construct's Q-V | Priest, Lee & Bezanilla 2021 *eLife* 10:e58148 |
| human CaV1.1 with β1a, α2δ1, STAC3 | Xenopus oocyte, 50 nl at 0.1–0.5 µg/µl, 4–5 d | **peak ON gating ~2.1 µA; ionic ~0.27–0.35 µA and still deepening at the end of a ~59 ms step**, so not a steady state (the paper's own text: pore activation "takes several tens of milliseconds to reach steady-state open probability") | **none. All Q-V, G-V and F-V normalised; no Q_max, no count, no γ** | step from −90 to +20 mV | voltage-gated; +20 mV is near the peak of the I-V | 1 oocyte for the panel; 3–7 for the curves | Savalli 2021 *J Gen Physiol* 153:e202112915 |
| not channel-specific; amplifier characterisation on a model cell | bench, no oocyte in the loop | not applicable | not applicable. **Current noise 1.0–3.0 nA rms at 3 kHz and 3.6–5.4 nA at 10 kHz across five amplifier configurations**, so the low-noise oocyte method is itself a 3-fold range at fixed bandwidth. The paper quotes prior work at 1 nA @ 3 kHz (a Dagan CA-1B manual figure) and **1.2 nA @ 5 kHz (Stefani & Bezanilla 1998)**, which corroborates `SOURCES.md`'s figure from a retrievable PDF | not applicable | not applicable | 630 ms records, five configurations | Koerner 2024 *Biophys Rep* 4:100185 |

### 2g. Not an expression level: the excised-patch noise floor, now independently sourced

Included here because it was the one gap the critic named as fatal, and the gap-fill pass closed it.
`SOURCES.md`'s excised-patch noise row (0.05–0.2 pA at 1 kHz) is a computed component budget, and the
0.1–0.5 pA in `figure_6.Rmd` brackets one measurement from the excluded record. Below are seven
measurements or acceptance ceilings from six independent laboratories, all with a patch attached, all
with the bandwidth stated, none from this project. Carried to a 1 kHz band under two bounding
assumptions, `σ ∝ B` (rising spectrum, the figure's p = 1 patch case) and `σ ∝ B^0.5` (white).

| channel and preparation | σ as published | band | → 1 kHz, σ ∝ B | → 1 kHz, white | n | citation |
|---|---|---|---|---|---|---|
| native AMPA, mouse cerebellar granule cell, outside-out | **132 ± 3 fA (121–144)**, a measured mean | 2 kHz | 0.066 pA | 0.093 pA | 8 patches, seals 70–770 GΩ | Smith, Wang & Howe 2000 *J Neurosci* 20:2073 |
| recombinant AMPA, HEK293, outside-out | ≤300 fA (acceptance ceiling) | ~3.7 kHz | 0.081 pA | 0.156 pA | 10–40 MΩ pipettes | Prieto & Wollmuth 2010 *J Neurosci* 30:4449 |
| recombinant AMPA, HEK293, outside-out | <500 fA (acceptance ceiling) | ~4.5 kHz | 0.111 pA | 0.236 pA | same | Prieto & Wollmuth 2010 |
| native NMDA, rat dentate granule cell, outside-out | <300 fA (acceptance ceiling, checked before sealing) | 5 kHz | 0.060 pA | 0.134 pA | 20–30 MΩ pipettes | Rycroft & Gibb 2002 *J Neurosci* 22:8860 |
| **rat P2X2**, HEK293 stable, outside-out | **0.24 pA**, shut-level SD of the all-points histogram (the legend prints 0.20; 0.24 is the self-consistent value, since it alone reproduces the paper's own excess open noise of 0.92 pA) | 5 kHz | 0.048 pA | 0.107 pA | one worked record; i = 3.2 pA, 32 pS | Ding & Sachs 1999 *J Gen Physiol* 113:695 |
| recombinant NMDA, oocyte-excised outside-out | **0.27–0.62 pA**, shut-level SD across six subunit patterns | 2 kHz | 0.135–0.31 pA | 0.19–0.44 pA | 44 patches (Table III) plus 13–15 (Table II) | Premkumar & Auerbach 1997 *J Gen Physiol* 110:485 |
| GluA2(Q)/GSG1L AMPA, HEK293, outside-out | **0.21–0.35 pA**, natively at a 1 kHz filter, no conversion needed | **1 kHz** | 0.21–0.35 pA | 0.21–0.35 pA | 5 per condition; 6–8 MΩ pipettes | McGee 2025 *J Neurosci* 45:e1930242025 |

The band this lands in is **0.05 to 0.35 pA at 1 kHz** under either assumption, which contains the
0.1–0.5 pA the figure uses and contains the 0.049 pA component budget at its bottom. So the patch box's
vertical position survives on independent evidence, and the correct citation for it is now Smith 2000
plus McGee 2025 (the only record already at 1 kHz) rather than the excluded one. Two caveats. Published
traces are the best records of a set, so these are lower bounds on typical noise. And `σ ∝ B` is a
convention, not a bound: a patch dominated by pipette capacitance times amplifier voltage noise has a
spectrum going as f², hence `σ ∝ B^1.5`, which would put the 1 kHz values below 0.066 pA. The bracket
above is therefore one-sided at the bottom.

The NSFA route the critic proposed for this, extracting the fitted baseline variance `σ_B²` from papers
that report a count for the same patch, does not work. Four papers fit and subtract a measured baseline
variance at a stated bandwidth (Rook, Coombs, Han, Beltrán-Matas) and **none of them prints its value**.
Only Coombs plots variance against mean on a labelled pA² axis, so an intercept could in principle be
digitised. Route recorded as closed.

### 2h. Rejected, and corrected

Knowing which number NOT to use has already proved worth its space in this folder (see the Jackson
rejection in `SOURCES.md` §"Oocyte TEVC noise"). REJECTED means do not use. CORRECTED means the claim
survives with a changed value or a changed label.

**Rejected**

1. **Jackson 2025 Fig 2a and 2b trace amplitudes (12.5 nA and 5.7 nA).** Both misread. Re-measurement
   gives ~15 nA and ~4.4 nA; the 5.7 nA reading implied a trough 250 px below the panel, inside the
   cartoon under the trace. Worse, the corrected pair makes the 1:5 example three times LARGER than the
   5:1 example, opposite to the paper's own statement that 5:1 injections gave larger currents. Single
   illustrative cells. Jackson's **stated** 3–200 nA range is CONFIRMED and is used in §1; the trace
   readings are not. (Jackson's noise numbers were already rejected in `SOURCES.md`, for a different
   and independent reason.)
2. **Geng 2023 whole-oocyte derived floor N ≥ 9.6e5.** Wrong configuration. `i` = 31.2 pA is 312 pS ×
   100 mV in symmetric 160 mM KCl, where E_K = 0. The Fig 2 record is a whole oocyte in ND96 with 2 mM
   external K⁺, where the paper's own Fig 2B legend puts E_K near −80 mV, so the driving force at a
   +100 mV step is ~180 mV and the unitary current is ~45–56 pA. The quotient is then ~5.4–6.7e5, so
   the stated floor does not bound N from below. The 30.9 µA current itself stands.
3. **Hernandez 2017 derived N = 11,970.** The 29.38 pS input is the paper's ZINC-experiment value, not
   its plain wild-type conductance; the primary figure is 25.34 ± 1.77 pS. And the ohmic extrapolation
   from +100 mV to −20 mV is refuted by the paper's own voltage series for the same receptor (24.84 pS
   at +80, 29.38 at +100, 35.83 at +120 mV, a 44 % rise over 40 mV). With the primary conductance the
   floor is 13,880, and even that inherits the bad extrapolation. Report as a range or not at all.
4. **Heigl 2023 count per cell (970–1940).** The paper's own three quantities do not close with its own
   tail currents under any C_m: inverting its relation gives 493 channels for a 1800 pA tail, which at
   0.97/µm² implies a 5 pF cell, and a 5 pF cell at 20.64 pA/pF cannot produce 1800 pA. Keep the
   density, discard the count.
5. **Kitaguchi 1.0e10 channels per oocyte.** Its own printed value is 1e9 and its own text says >10^8.
   The density ceiling arbitrates: 1e9 is 48.5/µm², at the measured ceiling, while 1e10 is 485/µm²,
   5–10× above anything ever counted in a patch. Use 1e8–1e9, and only as a surface count.
6. **Horrigan macropatch top of 7.2e4 channels as a functional count.** It needs 380 channels/µm² in a
   188 µm² macropatch, 4–8× above the empirical density ceiling, and it is a gating-charge count in a
   heavily injected oocyte, so it is subject to the surface-versus-functional discount of §4.5 flag 4.
7. **Lukacs 2021's 9.54 nA as one cell's current.** It is an ensemble sum over up to 20 cells, ~0.48 nA
   per cell, and the paper says occupancy per hole cannot be verified. Listed so it cannot be misread.
8. **Nekouzadeh & Rudy 2007's "macroscopic current" as an expression level.** It is the SUM of 100 real
   single-channel sweeps, so its channel count is the number of records summed, not a property of a
   preparation. Its amplitude is bounded at 100 × 19.5 pA = 1.95 nA by construction.
9. **Wang 2012's "N_c = 250 pS" as a channel count.** Stated in siemens. It is a lumped maximal
   conductance. Its figures separately set N_c = 1 as a normalisation.
10. **Firsov's "125–550 sites per pA" reported as wrong by a factor 1000.** That flag is itself the
    error, and it would have moved every ENaC density by three decades. Redone: 1 fmol = 6.022e8
    molecules; the 4.82 µA/fmol slope gives 6.022e8 / 4.82e6 pA = 125 sites per pA exactly, and 1.13
    µA/fmol gives 547. The paper is internally consistent. Its own P_open range, however, does use the
    site count where its own sentence assumes three subunits per channel: under its stated stoichiometry
    the numbers are 0.010–0.040, not 0.004–0.014, a factor 3. The conclusion survives either way.
11. **Sigworth's CMP 610b course notes as an independent source for the TEVC ceiling.** Unrefereed, and
    the same author as Baumgartner 1999, so consistent with it rather than independent of it. Kept in
    §4 only because it is the one place the R_s·I arithmetic is written out.
12. **`Zhang_Hamill_2000_...pdf` in this folder.** A 57-byte JSON error stub, not a PDF. Its quotes were
    read from the PMC HTML at `https://pmc.ncbi.nlm.nih.gov/articles/PMC2269787/`, which does serve to
    a plain request with a browser user-agent.

**Corrected, and the corrections add data**

13. **Han 2024 "no peak amplitude reported".** False. Table 1's fourth data column is headed "A (pA)"
    and holds the fitted deactivation amplitudes, 52.4 to 222.0 pA across six HEK293 constructs. Cite
    the measured amplitudes rather than deriving them. (The related remark that SDs exceed means in
    three of six groups should read five of six.)
14. **Rook 2020 "N plotted but no numeric value".** N per patch IS readable: Fig 2C's middle panel plots
    "Number of Channels" on a labelled 0–1000 axis, wild-type mean ~390 with individual patches at
    ~150, ~185, ~370 and ~860. Internally consistent with the Fig 2B ~200 pA peak at i = 0.6 pA and
    P_open 0.86. Outside-out peak currents are likewise readable, ~200 pA and ~500 pA.
15. **Mørkve & Hartveit "peak amplitude not reported".** It is: 93.4 ± 9.6 pA (range 44–171), 13 patches,
    which lands almost exactly on the 95 pA that had been derived from their N, γ and P_open. Also their
    equation (6) as transcribed is identically zero; the correct form is `σ² = i·I − I²/N + σ_B²`.
16. **Krashia 2010 outside-out agonist "5 mM ACh, saturating".** It is 5 **µM**, roughly 1/28 to 1/60 of
    the paper's own EC50 (138–309 µM), chosen low on purpose to isolate single openings. Single-channel
    amplitude is concentration-independent, so the 2.6 and 3.9 pA survive; the saturating label does not.
17. **Mauerer 1998 "approximately four channels per patch" as the observed count.** That is the authors'
    back-calculation from whole-cell conductance. Their direct observation, in the same paper, is
    "Seals on the BLM typically contained from 2 to ~25 K⁺ channels" over 559 seals, with one-channel
    patches in only 0.7 %.
18. **Premkumar & Auerbach n = 113 and n = 17 patches.** Those are per-pattern occurrence counts, and a
    patch is counted once per pattern it showed. The paper's Methods give the real figures: 44 patches
    (NR1 mixtures, Table III) and 13, or 15 by its own inconsistent Results, for Table II.
19. **Defaye 2024 absolute currents.** Fig 5F "a few hundred pA" is ~2.65 nA (the trace is ~13 scale
    bars deep, not ~2), and Fig 5J "1–2 nA" is ~3.9 and ~6.2 nA. The corrected values are the
    self-consistent ones: divided by the panels' own densities they give 35–36 pF, ordinary for a DRG
    neuron, where the low readings would imply ~10 pF.
20. **Savalli 2021 ionic "steady-state 270 nA".** The red trace never plateaus in the displayed window.
    Report ~0.27–0.35 µA and still rising.
21. **Del Core & Mirams 2025 title.** The paper is "Parameter inference for stochastic reaction models
    of ion channel gating from whole-cell voltage-clamp data". The title first reported for it does not
    exist.
22. **Rook 2020 citation.** *eLife* 9:e51111 is Rook, Williamson, Lueck, Musgaard & MacLean, "β11-12
    linker isomerization governs acid-sensing ion channel desensitization and recovery". The channel is
    chicken ASIC1, and the paper says explicitly that "Tachyphylaxis of mammalian ASIC1a in patches
    precludes using NSFA", so attributing the analysis to ASIC1a asserts the one thing it rules out.
    Cell line HEK293/HEK293T, not CHO-K1.
23. **`Pantazis_Olcese_2024_...pdf` authorship.** The paper is Koerner, Delgadillo Bonequi, Shogren,
    Stroschein, Haag & Boland. Pantazis and Olcese are a reference inside it. The filename misleads.
    Also, its 1.0 nA at 3 kHz is its own measurement of its own amplifier, and only coincidentally
    equal to the Dagan CA-1B manual figure it quotes.
24. **Alvarez 2002 method description.** Its Eq. 4 is `σ² = i⟨I⟩ − ⟨I⟩²/N` with no background-variance
    term; it returns i and N only, and P_open,max is computed afterwards from Eq. 6. Also, the Shaker
    panel is NOT internally consistent: 6000 × 1.4 pA = 8.4 nA against an 8.1 nA plateau implies
    P_open = 0.96, while the panel prints 0.8. The KCa panel is consistent (0.79 against a printed 0.8).
25. **Fujiwara & Kubo voltage.** The Fig 3A ramps are recorded from −75 to +75 mV, so the troughs are
    the current at −75 mV, and the Fig 3C x-axis voltage is never stated. Any count derived from them
    carries a 1.25-fold voltage uncertainty. Fig 4 separately uses the inward current at −80 mV, so
    "the expression index throughout is the inward current at −60 mV" is not the paper's convention.
26. **Islas & Sigworth Kv2.1 "four independent counts".** Three. The Fig 6 legend's 1500 and the text's
    ~1650 are the same patch, reported inconsistently.
27. **Csánády 2005 "~2 pA".** That is the scale bar. The record sits ~6.5 pA below the zero-current
    level in saturating ATP, so the gap to a BK macropatch is 3 to 3.5 decades, not four.
28. **Shi 2019 conductance levels.** Five values against five labels: C = 0.053 pA, O1 = 1.5, O2 = 3.0,
    O3 = 4.6, O4 = 5.9 pA. The version that dropped C makes O1 read as 3.0 pA, a factor-2 error in any
    `i` taken from it.
29. **Karasawa 2016 panel B against panel C.** ~9-fold, not 20-fold. The 20 came from comparing scale-bar
    labels rather than trace depths.
30. **Werner 1996's 3150 ± 260 nA described as "the figure `SOURCES.md` currently carries unattributed".**
    `SOURCES.md` carries no such figure. It is new to this file.

---

## 3. Source's own N kept apart from derived N

### 3a. Counted. Fluctuation analysis, gating charge, antibody binding, direct resolution of levels.

| N_ch | configuration | channel | method | citation |
|---|---|---|---|---|
| 1 | outside-out, oocyte | rat P2X2 | direct counting of open levels, single-channel patches | Ding & Sachs 2002 |
| 1–3 | outside-out, oocyte | rat P2X2 | direct counting of simultaneous levels + binomial check | Ding & Sachs 2002 |
| 1 | outside-out, HEK293 | human α3β4 nAChR | single-channel amplitude histograms | Krashia 2010 |
| ~1 in 40 % of patches | cell-attached, HEK293 | GluA3i-G AMPA | resolution of individual openings, calibrated against a ~200 pA whole cell | Shi 2019 |
| 2–25 per seal | inside-out + cell-attached, native | axolotl K_ATP | direct counting of resolved levels, 559 seals | Mauerer 1998 |
| 3.5 ± 3.1 (n = 281) | cell-attached, Neuro2A | Piezo1 | I_max/i with i from 35 one-channel patches | Lewis & Grandl 2021 |
| ~10–50 | outside-out, HEK293T | Drosophila iGluR + Neto | NSFA | Han 2024 |
| ~15–30 | outside-out, HEK293T/17 | GluA2(Q)/TARP | NSFA, from the plotted variance-mean parabola | Coombs 2022 |
| 28.1 (10.3–70.3) | outside-out, native | rat GlyR, rod bipolar terminal | NSFA | Mørkve & Hartveit 2009 |
| ~70 (65–90) | outside-out, tsA 201 | mouse muscle nAChR | maximum-likelihood fit of the trace, N_ch free | Milescu 2005 |
| 68.2 (26.5–110) | nucleated patch, native | rat GABA_A, AII amacrine | NSFA | Beltrán-Matas 2023 |
| 94 ± 2.8 | macropatch, oocyte | BK (KCa) | NSFA | Alvarez 2002 |
| 134 (24–291) | outside-out, native | rat GABA_A, A17 amacrine | NSFA | Beltrán-Matas 2022 |
| ~390 (150–860) | outside-out, HEK293 | chicken ASIC1 | NSFA, per-patch values plotted | Rook 2020 |
| "tens to hundreds" | macropatch, definitional | generic | technique definition | Axon Guide ch. 5 |
| "tens or hundreds" | macropatch, oocyte | CFTR | not stated | Csánády 2005 |
| "hundreds" | macropatch, oocyte | mSlo1 BK | G_Kmax/g_K plus a Poisson fit in the same patch | Horrigan, Cui & Aldrich 1999; Chen 2011; Sun & Horrigan 2022; Ma 2006 |
| "many hundreds to thousands" | macropatch and whole cell | human BK | **no method given** | Geng 2023 |
| "typically hundreds" | macropatch, oocyte | HCN1/2/4 | field-norm statement | Benndorf 2025 |
| 967, 1250, ~1500–1650 | cell-attached macropatch, oocyte | Kv2.1 | NSFA at +70 mV | Islas & Sigworth 1999 |
| 2250 | cell-attached macropatch, oocyte | Shaker | NSFA, cross-checked against single channels in the same patch | Islas & Sigworth 1999 |
| 6000 ± 291 | macropatch, oocyte | Shaker | NSFA | Alvarez 2002 |
| 4707 and 8715, same patch | macropatch, oocyte | mHCN2 | fluctuation analysis with neighbouring-trace variance | Liu 2016 |
| 11,800 | cell-attached macropatch, oocyte | Shaker | NSFA, 257 traces | Rodriguez 1998 |
| up to 20,000 | macropatch, oocyte | mHCN2 / CNG | patch-clamp fluorometry calibrated against an electrical count | Liu 2016 |
| 8.8e4 to 6.3e5 (median 4.0e5) | whole cell, CHO | hERG1a | maximum-likelihood fit of a stochastic gating model, η free, g_s constrained | Del Core & Mirams 2025 |
| 1.6–1.8e8 **surface** | oocyte TEVC | rat ENaC | radioligand antibody binding on the same oocytes then clamped | Firsov 1996 |
| ~1e9 **surface** (>1e8 routine) | oocyte TEVC | Shaker V478W | gating charge, Q_max/13.6 e | Kitaguchi 2004 |
| 4e9 **surface** (1.5e9 at 13 e) | oocyte cut-open | Shaker | gating charge, n a free parameter of a two-step Boltzmann | Rodriguez 1998 |
| **densities**: 2.6 / 25.0 / 46.1 per µm² | outside-out, native | rat axonal Na_V | NSFA, already P_open-corrected | Hu & Jonas 2014 |
| **density**: 1.75 per µm² native, ~45 overexpressed, ~100 extreme | cell-attached | Piezo1 | I_max/i over an imaged dome area | Lewis & Grandl 2021 |
| **density**: 0.97 ± 0.153 per µm² | whole cell, HEK-293 | Ca_V1.4 | NSFA normalised to C_m at 1 µF/cm² | Heigl 2023 |
| **density**: 2.23 per µm² | native tubule | axolotl K_ATP | whole-cell G / γ / P_open over a 10,000 µm² surface | Mauerer 1998 |
| **charge**: Q_max = 1–30 fC per patch | macropatch, oocyte + HEK | mSlo1 BK | admittance analysis; the paper prints the charge and its own formula, never a count | Horrigan & Aldrich 1999 |
| **charge**: Q_max = 22 nC per oocyte | oocyte TEVC | Shaker V478W | integration of gating current | Kitaguchi 2004 |
| **charge**: Q_off = 3.2 nC | oocyte cut-open | Shaker | integration, unsubtracted records | Rodriguez 1998 |

### 3b. Derived. `N_ch = I / (i · P_open)`, arithmetic on top of a stated current and an assumed or borrowed `i` and `P_open`.

Every row states which `i` and which `P_open`, and whether the result is a floor.

| N_ch | I | i, and where it comes from | P_open | floor? | configuration and citation |
|---|---|---|---|---|---|
| 154–159 | 1160 pA | 7.53 pA, from the same paper's 73.8 pS at its own −102 mV in matched symmetric Cl⁻ | 0.97, the same paper's max cluster P_open (cell-attached, so a voltage mismatch) | yes, twice over: a 2 ms pulse need not reach the equilibrium maximum, and 73.8 pS is a main-conductance figure | outside-out, Ivica 2022 *J Physiol* |
| 46–48 | 348 pA | same | same | yes | outside-out, α1β GlyR, Ivica 2022 |
| 164–178 | 12.25 nA | 74.9 pA, the same paper's 312 pS taken as ohmic to +240 mV in matched symmetric KCl | 0.920, from the paper's own Boltzmann | soft: if the unitary I-V is sublinear at +240 mV the true count is larger (i = 50 pA gives 265), which is the direction that reconciles this with the authors' own "hundreds" | macropatch, Geng 2023 |
| 1.4e3–7.9e4 | Q = 1–30 fC | 2.36–4.4 e per channel (the paper's own 2.36, and 2.6–4.4 from the literature it quotes) | not applicable, charge | no; Q-V saturated. But see §2h item 6 | macropatch, Horrigan & Aldrich 1999 |
| 5.5–756 | not used | density × area | already included in the density | no | outside-out, Hu & Jonas 2014, area from its own capacitance regression |
| 5036 | 2800 pA | **0.80 pA measured** for the same channel in the same cell background | **0.695 measured** | no, both factors measured | whole cell, Khadra 2012 current + Sattler 2020 `i·P_open` |
| 5.7e3 (1.14e4 at P = 0.5) | 6670 pA | 1.17 pA, from the paper's own 39 pS at its own −30 mV | 1 | yes | whole cell, Krashia 2010 |
| 1180–2400 | 8600 pA | 3.69 pA, the paper's own 73.8 pS at ~−50 mV | 0.97 | yes | whole cell, Ivica 2022 *J Physiol* |
| 8.7e3–2.3e4 | 796 ± 149 pA/pF × 10–25 pF | 1.32 pA from the paper's own 22 pS taken ohmically to −60 mV (or 1 pA on the figure's model) | 0.695 | yes, and the ATP concentration is unstated | whole cell, Stelmashenko 2012 |
| 1.13e6 | 2050 nA | 1.82 pA, the paper's own 26 pS at −70 mV | 1 | yes | oocyte TEVC, Krashia 2010 |
| 5.7e6 | 3150 nA | 0.556 pA/channel, from Sattler's measured pair | included | no | oocyte TEVC, Werner 1996 + Sattler 2020 |
| 8e5–5.2e7 | 0.8–65 µA | 1 pA at the ramp's −75 mV | 1 | yes | oocyte TEVC, Fujiwara & Kubo 2004 |
| 3.4e6–1.4e7 **functional** | ~1–4 µA | 0.6 pA, the paper's own | 0.010–0.040, the paper's own under its own stoichiometry | this is the functional count implied by a **surface** count of 1.6–1.8e8, a 25–100× discount | oocyte TEVC, Firsov 1996 |

The cross-configuration ratio worth having, because it is one paper, one lab, one receptor, one agonist
and one fitting procedure: **Krashia 2010** gives oocyte TEVC 2050 ± 301 nA against HEK293 whole cell
6.67 ± 1.63 nA. Raw 307×; corrected to equal driving force (−70 against −30 mV, reversal 0 mV), **132×**.
The same comparison assembled across labs (Venkatachalan's oocyte 10.3 µA against Olander's HEK 1.6 nA)
gives 4.8e3 after the same correction, 36× larger. That gap is the argument for preferring single-source
ratios, and it belongs on the figure rather than hidden. A second single-source pair, weaker because the
solutions and voltages differ: **Geng 2023** whole oocyte ~30.9 µA against an excised macropatch of the
same oocytes ~12.25 nA, a factor 2.5e3.

---

## 4. Consistency check

Everything in this section is either an input with its source or a computation. Script:
`bridge.py`, `exp.py` in the session scratchpad; every figure below reproduces.

### 4.1 The unit bridge, and one convention that has to be fixed first

Specific membrane capacitance 1 µF/cm², stated as the assumption by three of the sources (Hu & Jonas
supplement, Mauerer, Heigl):

    1 µF/cm² = 0.01 pF/µm²   ⇒   A(µm²) = 100 · C(pF)

"Patch area" appears in two inequivalent senses. Mauerer states his: a 1.5 µm inner tip gives "the
**minimal** membrane area of a patch is 1.77 µm²", which is exactly the projected mouth π(0.75)², and
he adds that a larger real area makes his conductance smaller. Lewis & Grandl use the same projected
convention (an imaged 0.8 µm dome radius gives 2.01 µm², and 89.2/2.01 = 44.4, their stated ~45/µm²).
Capacitance-measured areas are larger by a dome factor, and Wang & Hilgemann measure it:

| pipette diameter | A_disc = πd²/4 | area from measured capacitance | factor |
|---|---|---|---|
| 10 µm | 78.5 µm² | 2 pF = 200 µm² | **2.55** |
| 15 µm | 176.7 µm² | 4 pF = 400 µm² | **2.26** |

Matched ends give 2.3–2.6, slightly deeper than a hemisphere (2.0). Taken as unmatched ranges instead
they span 1.1–5.1, so the factor is good to about a factor 2. Used below as **A_cap = 2.4 · A_disc**.

| tip diameter | A_disc | A_cap | C_patch |
|---|---|---|---|
| 0.5 µm | 0.196 µm² | 0.47 µm² | 0.0047 pF |
| 1.5 µm | 1.77 | 4.24 | 0.042 |
| 2 µm | 3.14 | 7.54 | 0.075 |
| 5 µm | 19.6 | 47.1 | 0.47 |
| 10 µm | 78.5 | 188 | 1.88 |
| 15 µm | 176.7 | 424 | 4.24 |

So conventional patch 0.5–7.5 µm², macropatch 47–190 µm², giant patch 190–420 µm². Area goes as d², so
the conventional-to-macropatch step is exactly 100× for a 0.5 µm against a 5 µm tip, which is the
20–200× span of realistic pairings. The oocyte, from Baumgartner's measured 45.6 ± 5 mF/m² on a 1.2 mm
sphere: geometric area 4.52e6 µm², folding factor 4.56 ± 0.50, true area 2.06e7 µm², whole-cell
capacitance **206 nF**, which reproduces `SOURCES.md`'s own figure and sits beside the 200 nF that its
thermal-floor derivation already uses. Zhang & Hamill get 5.3 by a different method in a different lab,
so the folding factor is 4.6–5.3 and the honest excess-area range is 5–10.

### 4.2 Does density × area reproduce the counted N?

Yes, and with no free parameter left over. Empirical densities from patches counted directly: Hu & Jonas
2.6, 25.0, 46.1 per µm² (3.9, 37.9, 69.8 if divided by their own measured P_open of 0.66, which they
do not do); Lewis & Grandl 1.75 native, ~45 overexpressed, ~100 in the most extreme patch; Mauerer 2.23.

| configuration | A_cap | N at 1/µm² | at 45/µm² | at 100/µm² |
|---|---|---|---|---|
| conventional patch | 0.47–7.5 µm² | 0.5–7.5 | **21–339** | **47–754** |
| macropatch | 47–188 µm² | 47–188 | 2120–8460 | 4710–18,800 |
| giant patch | 188–424 µm² | 188–424 | 8460–19,100 | 18,800–42,400 |

The sourced outside-out ligand-gated counts (Han 10–50, Coombs ~15–30, Milescu ~70, Ivica 159,
Rook ~390, and 1–760 including the density route) are reproduced by a 0.5–2 µm tip at 45–100 per µm².
Inverting for the figure's box, `N_ch` = 10–300 in a mid-size 4 µm² patch needs 2.5–75 per µm², inside
the measured band. The macropatch column likewise contains Alvarez's 94 and 6000, Islas's 2250, and
Liu's 8715 and 20,000. **Packing is not the binding constraint**: a 10 nm channel footprint is 7.9e-5
µm², so 100 per µm² occupies 0.8 % of the membrane and even 485 per µm² only 3.8 %. The ceiling of
~100 per µm² is empirical, not steric.

### 4.3 Is the ordering across configurations monotone?

In current and in count, yes.

| configuration | C_m | current | N_ch |
|---|---|---|---|
| excised patch | ~0.04 pF | 5 pA – 2 nA | 1 – 760 |
| macropatch | 0.5–4 pF | 6.5 pA – 12 nA | 90 – 2e4 |
| whole cell | 10–40 pF (25 pF mean over 42,272 wells) | 100 pA – 31 nA | 1e3 – 6e5 |
| oocyte TEVC | 206 nF | 3 nA – 65 µA | 8e5 – 1.2e8 |

In surface density, **no, and this is the place to say it plainly.** Channels per pF needs no folding
factor and no diameter, so it is the only comparison that survives the geometry, and the same channel
at the same expression level must give the same value in all three configurations. It does not:

| box, as `figure_6.Rmd` stands | N_ch | C_m used | channels/pF | channels/µm² |
|---|---|---|---|---|
| excised patch (A = 4 µm²) | 10–300 | 0.04 pF | **250 – 7500** | 2.5 – 75 |
| whole cell | 20 – 1e4 | 25 pF | **0.8 – 400** | 0.008 – 4.0 |
| oocyte | 1e5 – 5e7 | 2.06e5 pF | **0.49 – 243** | 0.005 – 2.4 |

The patch box's top edge is 19× above the whole-cell box's and 31× above the oocyte's, and its bottom
edge is 300× above the whole-cell bottom. The three overlap only on a sliver, 250–400 channels/pF.
Excised patches come off these same cells, so they cannot all be right. See flag 1.

### 4.4 What hangs together

Four routes agree without adjustment, and they are the load-bearing ones.

1. **The patch count axis.** 0.5–2 µm tips at the measured 45–100 per µm² give 21–754 channels; the
   sourced outside-out set spans 1–760; and the figure's `n_hi` = 300 equals the density ceiling times
   a typical 4 µm² patch (200–400). This supports the count edge by a route independent of the
   outside-out sorting.
2. **The dome factor**, 2.26–2.55 from two matched endpoints of one paper, consistent with Mauerer
   explicitly calling the projected disc a minimum.
3. **The oocyte capacitance**, 206 nF derived here against 166 nF measured by Zhang & Hamill and 200 nF
   assumed in `SOURCES.md`'s own thermal-floor calculation, three routes inside 20 %.
4. **`i · P_open`**, 0.500 pA/channel on the figure's model against 0.556 pA/channel measured for P2X2
   itself by Sattler, an 11 % agreement. Since the whole reverse map runs through that product, this is
   the input most worth having validated, and it is validated.

### 4.5 Flagged disagreements: routes to the same count differing by more than 3×

**Flag 1, the largest, and it is the figure's own numbers.** The patch and whole-cell boxes are disjoint
in surface density by 19–300× (§4.3). Read one way, whole cell `n_hi` = 1e4 implies 16 channels in a
4 µm² patch against a sourced outside-out range of 10–760. Read the other, patch `n_hi` = 300 implies
1.9e5 channels in a 25 pF cell (4e4 to 2.6e6 over the plausible area and capacitance ranges). The likely
wrong input is the whole-cell count edge, for three reasons: 400 channels/pF is 4 per µm², the density
of a **native** somatic membrane (Hu & Jonas' soma is 2.6 per µm²) rather than of a transfected
overexpressing cell; the independently sourced whole-cell counts already run past it (Stelmashenko
2.3e4, Del Core & Mirams 8.8e4–6.3e5, the latter overlapping the patch band in channels/pF); and the
density ceiling allows 1.25e5–2.5e5 channels in a 25 pF cell. The reason current-derived whole-cell
counts look small is the clamp, not the biology: 1.25e5 channels at 0.5 pA each is 62 nA, ten times the
ceiling. So a whole-cell count edge set from current is a clamp limit wearing an expression limit's
clothes, and the figure already imposes the clamp limit through its `Imax` diagonal. The comment block
at `figure_6.Rmd:406` has already moved this row onto exactly that footing (`20 = 100 pA/(10 pA × 0.5)`
and `1e4 = 5 nA/(1 pA × 0.5)`), which is the right fix; what remains is that the same row is then not
an expression range and should not be read as one, and that the stale comment at line 358 still calls it
one.

**Flag 2.** Two measured routes to patch area at the same pipette resistance differ by 9.2×. Hu & Jonas'
capacitance regression `A = 0.08271·g_P − 0.47526 µm²` at 4.35 MΩ gives 18.5 µm²; Lewis & Grandl imaged
their own 3–6.5 MΩ pipettes (mean 4.35 ± 0.8 MΩ) and state ~2 µm². Correcting the projected disc to a
capacitance-equivalent dome leaves 3.8×, still flagged. At a stated density this maps one-to-one onto
`N_ch`. Neither measurement is obviously wrong; the wrong input is **pipette resistance used as a proxy
for tip diameter across labs**, since resistance depends on taper as well as tip and brain-slice
pipettes have long tapers. A 2.5 µm tip at 2.4× gives 11.8 µm², close to Hu & Jonas' 16.4 µm² at
4.9 MΩ, which accounts for the sign and roughly the size of the gap. **Anchor the count axis on tip
diameter or on measured patch capacitance, never on pipette resistance**, and do not apply Hu & Jonas'
regression to another lab's resistances.

**Flag 3.** Kitaguchi disagrees with itself 10-fold and the recomputation lands on the high end:
22 nC / (13.6 × 1.602177e-19 C) = **1.01e10**, against its printed 1e9 and its own text's >1e8. The
charge that would give 1e9 is 2.18 nC. The density ceiling arbitrates (§2h item 5): use 1e8–1e9.

**Flag 4, the one that reverses a recommendation.** Surface count and electrically effective count
differ by 25–100× in one sourced paper, and that paper says so. Firsov reports 125–550 binding sites per
pA, three subunits per channel, so 42–167 channels per pA, while the current-derived count at
i = 0.6 pA is 1.67 channels per pA. If the gating subset has P_open ≈ 0.5, only 2–8 % of surface
channels are electrically effective. The figure's `N_ch` enters as `I = N_ch·i·P_open`, so it is the
**functional** count, and every binding-derived or gating-charge-derived count is an upper bound on it,
sometimes by two decades. Firsov's 1.6–1.8e8 surface channels are 3.4e6–1.4e7 functional ones, which puts
the figure's oocyte `n_hi` = 5e7 **4–15× above** Firsov rather than 3.4× below it. So treating Firsov as
sourced support for raising the oocyte top edge does not survive the arithmetic, and the direction of the
discrepancy reverses. The same caution does not apply to Fujiwara's 8e5, which is a functional count and
therefore directly comparable, though it is a deliberate 1:1000 dilution sitting at 3.9 channels/pF,
12× below the lowest whole-cell density in the set.

**Flag 5, genuinely open.** Whether the oocyte's microvillar folding survives the suction into a
macropatch is a 4.56× question nobody has answered. A 10 µm tip on an oocyte is 78.5 µm² projected;
flattened villi give 188 µm² (1.9 pF), surviving villi give 860 µm² (8.6 pF). At 8.3 channels/µm² that
is 1560 against 7100 channels, and the sourced macropatch counts straddle both, so they cannot
discriminate. The route that avoids it is channels per pF, and it needs a **measured** macropatch
capacitance; the only measured giant-patch capacitance in the folder is Wang & Hilgemann's, on RBL cells
rather than oocytes.

**Flag 6, soft.** Horrigan's macropatch top of 7.2e4 needs 383 channels/µm² in a 188 µm² macropatch,
4–8× above the empirical ceiling, and in the same direction as flag 4's surface discount. Rejected as a
functional count (§2h item 6).

---

## 5. What is still missing

A count for a ligand-gated channel in an excised **macropatch** exists nowhere: every hard macropatch
count in this sweep is a voltage-gated or calcium-gated channel counted by fluctuation analysis or
gating charge, because those methods force a count while the ligand-gated macropatch literature has no
need of one. No paper measures patch **area** except Hu & Jonas, and its regression is calibrated over
31–204 nS, so applying it at a macropatch's 500–2000 nS is a 2.5–10× extrapolation. Nobody reports
channels per square micrometre for an oocyte, so every oocyte density here is built from a count and a
separately sourced area. Between-frog and between-batch spread is asserted qualitatively three times and
quantified never, while within-batch spread is measured once (Werner, 142–3467 nA over 26 oocytes at a
fixed dose). And no cut-open paper states a clamped membrane area, an aperture diameter or a clamped
capacitance, so the cut-open row has an absolute count and an absolute charge and no density at all. One
gap the critic named as fatal has closed: the excised-patch noise floor is now carried by six
independent laboratories with the patch attached and the bandwidth stated (§2g).

### 5b. THE PROJECT'S OWN P2X2 STUDY IS DELIBERATELY ABSENT FROM EVERY TABLE HERE

No number in §1 through §4 comes from Moffatt & Hume 2007 *J Gen Physiol*, Moffatt 2007 *Biophysical
Journal*, or Moffatt & Pierdominici-Sottile 2025 *Communications Biology*. Two reasons, both from the
author, 2026-07-27, and both bearing on the numbers rather than on provenance. The records **sum three
traces per response**, so the 1777 pA peak now carried in `SOURCES.md` and in `figure_6.Rmd` is not one
patch's current: one patch would be about 590 pA, with σ smaller by √3 so that the noise power falls
3-fold as well. And the preparation shows **unliganded gating**, with 99 % of the variance of a 13.7 ms
point coming from the channel rather than the instrument, so it is not a typical excised patch either.

**Consequence, stated and not acted on.** Two things still rest on that number. The excised-patch box of
the PREP table takes its noise edges from it (`sig_lo`/`sig_hi` = 0.1–0.5 pA bracket
S = 1.04e-4 pA²/Hz, which is Moffatt & Hume's own record, `SOURCES.md:41` and `:210`), and the EXP block
that draws the same record as "one REAL experiment" is built from it (already switched off at
`SHOW_EXP <- FALSE`, and the code's own withdrawal note gives the corrected figures). If the summing is
confirmed, **the box moves left by about a factor of three** and down by about the same in noise power.
The ordering check at `SOURCES.md:210-211` ("whole cell sits 1.4× to 86× above it, median 17×") is
computed against the same excluded S and inherits the same correction. `figure_6.Rmd` is **not edited
here**. Note also that its comment at line 360 says the project's own study "supplies NO number here",
which is true of the count axis and false of the noise axis eight lines above. Flagged only.

---

## 6. Links

Bib keys, where they exist, are in `papers/1_method/docs/manuscript-drafts/biblio_full.bib`. PDFs are
beside this file in `docs/bibliography/recording_configurations/` unless noted.

**Excised patch, outside-out and inside-out**

- Han, Vicidomini, Ramos, Mayer & Serpe, *J Physiol* 602:7043-7064 (2024) — https://doi.org/10.1113/JP287331
- Ivica, Lape & Sivilotti, *J Physiol* 600:333-347 (2022) — https://doi.org/10.1113/JP282171
- Ivica, Kaur, Harvey, Sivilotti et al., *eLife* (2022), AMS and the glycine receptor
- Rook, Williamson, Lueck, Musgaard & MacLean, *eLife* 9:e51111 (2020) — https://doi.org/10.7554/eLife.51111
- Coombs, Sexton, Cull-Candy & Farrant, *Mol Pharmacol* 101:343-356 (2022) — https://doi.org/10.1124/molpharm.121.000473
- Milescu, Akk & Sachs, *Biophys J* 88:2494-2515 (2005) — https://doi.org/10.1529/biophysj.104.053256 (PDF in `docs/bibliography/`)
- Khanra, Strauss, Moreno Wasielewski, Lenze, Meyerson, Reiner & Levitz, *Nat Commun* 17:3789 (2026) — https://doi.org/10.1038/s41467-026-72226-w
- Krashia, Moroni, Broadbent, Hofmann, Kracun, Beato, Groot-Kormelink & Sivilotti, *PLoS ONE* 5:e13611 (2010) — https://doi.org/10.1371/journal.pone.0013611
- Gasparri, Wengel, Grutter & Pless, *J Gen Physiol* 151:e201912347 (2019) — https://doi.org/10.1085/jgp.201912347
- Ding & Sachs, *BMC Neurosci* 3:17 (2002) — https://doi.org/10.1186/1471-2202-3-17
- Beltrán-Matas, Castilho, Veruki & Hartveit, *Eur J Neurosci* (2022) — https://doi.org/10.1111/ejn.15634
- Beltrán-Matas, Hartveit & Veruki, *Front Ophthalmol* 3:1134765 (2023) — https://doi.org/10.3389/fopht.2023.1134765
- Mørkve & Hartveit, *J Physiol* 587:3813-3830 (2009) — https://pmc.ncbi.nlm.nih.gov/articles/PMC2746612/ **NO PDF**; the PMC HTML does serve to a plain request with a browser user-agent, and the quotes were read there
- Hu & Jonas, *Nat Neurosci* 17:686-693 (2014) — https://doi.org/10.1038/nn.3678, with the patch-area regression in Supplementary Fig. 4
- Mauerer, Boulpaep & Segal, *J Gen Physiol* 111:139-160 (1998)
- Lewis & Grandl, *eLife* 10:e70988 (2021)
- Tamagnini, in *Patch Clamp Electrophysiology*, Methods Mol Biol 2188:229-242 (2020) — https://doi.org/10.1007/978-1-0716-0818-0_11 (accepted version from the Reading repository)

**Macropatch**

- Geng, Li, Butler, Wang, Salkoff & Magleby, *J Gen Physiol* 155:e202213302 (2023) — https://doi.org/10.1085/jgp.202213302
- Alvarez, González & Latorre, *Adv Physiol Educ* 26:327-341 (2002)
- Islas & Sigworth, *J Gen Physiol* 114:723-742 (1999)
- Rodriguez, Sigg & Bezanilla, *J Gen Physiol* 112:223-242 (1998)
- Horrigan & Aldrich, *J Gen Physiol* 114:305-336 (1999)
- Horrigan, Cui & Aldrich, *J Gen Physiol* 114:277-304 (1999)
- Schoppa & Sigworth, *J Gen Physiol* 111:271-294 (1998)
- Zhou & Lingle, *J Gen Physiol* 144:415-440 (2014)
- Chen, Geng & Magleby, *J Gen Physiol* 138:593-607 (2011)
- Sun & Horrigan, *Sci Adv* 8:eabq5772 (2022)
- Ma, Lou & Horrigan, *J Gen Physiol* 127:309-328 (2006)
- Wang, Rothberg & Brenner, *J Gen Physiol* 127:449-465 (2006)
- Liu, Xie, Grant, Su, Gao, Liu & Zhou, *J Gen Physiol* 148:65-78 (2016)
- Benndorf, Enke, Tewari, Kusch, Liu, Sun, Schmauder & Sattler, *PNAS* 122:e2422533122 (2025)
- Csanády, Chan, Nairn & Gadsby, *J Gen Physiol* 125:43-55 (2005)
- Ryu & Yellen, *J Gen Physiol* 140:469-479 (2012)
- Wang & Hilgemann, *J Gen Physiol* 132:51-65 (2008)
- Axon Guide ch. 5, Advanced Methods in Electrophysiology — PDF beside this file

**Whole cell**

- Khadra, Tomić, Yan, Zemková, Sherman & Stojilkovic, *J Gen Physiol* 139:333-348 (2012) — https://doi.org/10.1085/jgp.201110716
- Sattler, Eick, Hummert, Schulz, Schmauder, Schweinitz, Unzeitig, Schwede & Benndorf, *Sci Rep* 10:21751 (2020) — https://doi.org/10.1038/s41598-020-78672-w
- Stelmashenko, Lalo, Yang, Bragg, North & Compan, *Mol Pharmacol* 82:760-766 (2012) — https://doi.org/10.1124/mol.112.080903. **NO PDF** (not open access; PMC blocks scripted download; the publisher returns 403). Verified full text saved as `Stelmashenko_2012_..._current_density.fulltext.txt`
- Wand, Bruenings, Tewari, Reuter, Mrowka, Benndorf, Zimmer & Sattler, *Cells* 14:510 (2025) — https://doi.org/10.3390/cells14070510
- Harnau, Zeller, Fricke & Meier, *Sci Rep* 15:32686 (2025) — https://doi.org/10.1038/s41598-025-20807-y
- Hernandez, Kong, Hu, Zhang, Shen, Jackson, Liu, Jiang & Macdonald, *eNeuro* 4:e0251-16.2017 — https://doi.org/10.1523/ENEURO.0251-16.2017
- Olander, Janzen, Villmann & Jensen, *PLoS ONE* 15:e0234080 (2020) — https://doi.org/10.1371/journal.pone.0234080
- Li, Wang, Banerjee, Marinelli, Silberberg, Faraldo-Gómez, Hattori & Swartz, *eLife* 8:e47060 (2019)
- Karasawa & Kawate, *eLife* 5:e22153 (2016)
- Shi, Yuan, Sipple, Srinivasan, Ptak, Oswald & Nowak, *J Gen Physiol* 151:156-173 (2019) — https://doi.org/10.1085/jgp.201812209
- Li, Lu, Chiang, Chernov-Rogan, Grogan & Chen, *PLoS ONE* 12:e0180154 (2017) — https://doi.org/10.1371/journal.pone.0180154
- Montnach, Lorenzini, Lesage, Simon, Nicolas, Moreau, Marionneau, Baró, De Waard & Loussouarn, *Sci Rep* 11:3282 (2021) — https://doi.org/10.1038/s41598-021-82077-8
- Eltokhi et al., *Cell Rep Methods* 3:100559 (2023) — https://doi.org/10.1016/j.crmeth.2023.100559
- Lukacs, Mille, Snyder, Neubauer, Ilyes, Alsaloum, Dib-Hajj, Waxman & Mike, *Front Pharmacol* 12:738460 (2021) — https://doi.org/10.3389/fphar.2021.738460
- Xiao, Kang, Morgan & Leaf, *PNAS* 95:2680-2685 (1998)
- Sokolov, Peters, Rajamani & Ruben, *Front Pharmacol* 4:78 (2013)
- Dallas, Al-Owais, Hettiarachchi, Vandiver, Jarosz-Griffiths, Scragg, Boyle, Steele & Peers, *Sci Rep* 11:8194 (2021) — https://doi.org/10.1038/s41598-021-87198-8
- Del Core & Mirams, "Parameter inference for stochastic reaction models of ion channel gating from whole-cell voltage-clamp data", *Phil Trans R Soc A* 383:20240224 (2025) — https://doi.org/10.1098/rsta.2024.0224 (PDF in `docs/bibliography/`)
- Heigl et al., *Channels* 17:2192360 (2023) — https://doi.org/10.1080/19336950.2023.2192360
- Bachmann, Dürr, Gräschus, Assmann, Kirchner, Ehrhardt, Gräf & Gründer, *Channels* 15:198-210 (2021) — https://doi.org/10.1080/19336950.2021.1877594
- Ponce, Castillo, Hinojosa, Martinez-Rendon & Cereijido, *Physiol Rep* 6:e13663 (2018) — https://doi.org/10.14814/phy2.13663
- Ng, Farr, Young, Windley, Perry, Vandenberg & Hill, *Biol Methods Protoc* 6:bpab003 (2021) — https://doi.org/10.1093/biomethods/bpab003
- Defaye, Bradaia, Abdullah et al. (Altier), *J Clin Invest* 134:e176474 (2024) — https://doi.org/10.1172/JCI176474

**Xenopus oocyte, TEVC and cut-open**

- Werner, Seward, Buell & North, *PNAS* 93:15485-15490 (1996)
- Fujiwara & Kubo, *J Physiol* 558:31-43 (2004)
- Jackson, Hall & Sudweeks, *Int J Mol Sci* 26:9506 (2025)
- Ito, Kawazoe, Sato, Uesugi & Hirata, *Sci Rep* 10:13999 (2020)
- Nagel et al., *Nat Commun* 16:10367 (2025) — https://doi.org/10.1038/s41467-025-66244-3
- Germann, Shin, Steinbach & Akk, *J Gen Physiol* 157:e202413644 (2025)
- Venkatachalan, Bushman, Mercado, Sancar, Christopherson & Boileau, *Pflügers Arch* 454:155-163 (2007) — https://doi.org/10.1007/s00424-006-0183-1 (PMC author manuscript)
- Corbin-Leftwich, Small, Robinson, Villalba-Galea & Boland, *J Gen Physiol* 150:1583-1593 (2018)
- Firsov, Schild, Gautschi, Mérillat, Schneeberger & Rossier, *PNAS* 93:15370-15375 (1996)
- Kitaguchi, Sukhareva & Swartz, *J Gen Physiol* 124:319-332 (2004)
- Baumgartner, Islas & Sigworth, *Biophys J* 77:1980-1991 (1999). **No longer unretrievable**: the PDF is now beside this file, via `https://europepmc.org/articles/PMC1300479?pdf=render`. `SOURCES.md:290` says "CITE FOR CAPACITANCE ONLY", which is right about the absence of noise analysis and understates the paper: it is also the best source for the TEVC current ceiling and for why series-resistance compensation does not remove the error
- Papke & Smith-Maxwell, *Comb Chem High Throughput Screen* 12:38-50 (2009) — PMC author manuscript
- Zhang & Hamill, *J Physiol* 523:101-115 (2000) — https://pmc.ncbi.nlm.nih.gov/articles/PMC2269787/ **NO PDF** (the file of that name here is a 57-byte error stub)
- Guan, Chen & Zhang, "Two-Electrode Voltage Clamp", Methods Mol Biol 998:79-89 (2013) — https://doi.org/10.1007/978-1-62703-351-0_6 (PDF at `docs/bibliography/Two-ElectrodeVoltageClamp.pdf`)
- Sigworth, CMP 610b Lecture 3, "Two-Microelectrode Voltage Clamp" (Yale course notes). **Not peer reviewed**, and by a co-author of Baumgartner 1999, so not independent of it
- Rodriguez, Sigg & Bezanilla, *J Gen Physiol* 112:223-242 (1998), cut-open
- Priest, Lee & Bezanilla, *eLife* 10:e58148 (2021)
- Savalli, Angelini, Steccanella et al. (Olcese), *J Gen Physiol* 153:e202112915 (2021)
- Koerner, Delgadillo Bonequi, Shogren, Stroschein, Haag & Boland, *Biophys Rep* 4:100185 (2024) — PMC11549981. **The file here is named `Pantazis_Olcese_2024_...`, which is wrong**; Pantazis and Olcese are a reference inside it. It corroborates `stefani1998cutopen`'s 1.2 nA at 5 kHz verbatim, which unblocks a figure `SOURCES.md` holds from a paywalled chapter

**Excised-patch noise (§2g)**

- Smith, Wang & Howe, *J Neurosci* 20:2073-2085 (2000) — PMC6772487
- Prieto & Wollmuth, *J Neurosci* 30:4449-4459 (2010) — PMC2857311
- Rycroft & Gibb, *J Neurosci* 22:8860-8868 (2002) — PMC6757676
- Ding & Sachs, *J Gen Physiol* 113:695-719 (1999) — PMC2222910
- Premkumar & Auerbach, *J Gen Physiol* 110:485-502 (1997) — PMC2229386
- McGee, Bats, Farrant & Cull-Candy, *J Neurosci* 45:e1930242025 (2025) — https://doi.org/10.1523/JNEUROSCI.1930-24.2025

### PDFs newly saved into this folder by this sweep

Each verified present with `ls` before being listed. Not listed are the files already in the folder from
the earlier noise and recording-configuration passes.

| file | lane |
|---|---|
| `Khanra_2026_NatCommun_kainate_receptor_outsideout_patch_amplitude_criterion.pdf` | outside-out |
| `Coombs_2022_MolPharmacol_AMPA_TARP_gamma8_NSFA_conductance_Popeak.pdf` | outside-out |
| `Rook_2020_eLife_ASIC_outsideout_patch_NSFA_open_probability.pdf` | outside-out |
| `Ivica_2022_JPhysiol_GlyR_outsideout_patch_peak_currents_acidic_pH.pdf` | outside-out |
| `Shi_2019_JGP_AMPA_wholecell_amplitude_vs_channels_per_patch.pdf` | outside-out |
| `Benndorf_2025_PNAS_HCN_single_channel_conductance_macropatch.pdf` | macropatch |
| `Chen_Magleby_2011_JGP_BK_Mg_oocyte_macropatch_channel_counts.pdf` | macropatch |
| `Csanady_2005_JGP_CFTR_oocyte_macropatch_current_amplitude.pdf` | macropatch |
| `Sun_Horrigan_2022_SciAdv_BK_gating_lever_oocyte_macropatch.pdf` | macropatch |
| `Wang_Brenner_2006_JGP_BK_beta4_macropatch_series_resistance_limit.pdf` | macropatch |
| `Corbin-Leftwich_2018_JGP_oocyte_action_potentials_Nav_current_range_TEVC.pdf` | oocyte |
| `Kitaguchi_2004_JGP_Shaker_V478W_gating_charge_per_oocyte_TEVC.pdf` | oocyte |
| `Baumgartner_1999_BiophysJ_two_microelectrode_voltage_clamp_oocyte_errors.pdf` | oocyte |
| `Li_2017_PLoSONE_SyncroPatch768PE_CHO_Nav_current_amplitude_QC.pdf` | whole cell |
| `Stelmashenko_2012_MolPharmacol_..._current_density.fulltext.txt` (text, no PDF available) | whole cell |
| `AxonGuide_Ch5_Advanced_Methods_in_Electrophysiology_macropatch_channel_counts.pdf` | cross-cutting |
| `Olander_2020_PLoSONE_GABAA_alpha1beta2_alpha3beta2_HEK293_wholecell_amplitudes.pdf` | cross-cutting |
| `Tamagnini_2020_MethodsMolBiol_nucleated_somatic_macropatch_native_neurons.pdf` | cross-cutting |
| `Ding_Sachs_1999_JGP_P2X2_single_channel_outsideout_patch_noise.pdf` | noise gap-fill |
| `McGee_2025_JNeurosci_GSG1L_AMPAR_outsideout_patch_rms_noise.pdf` | noise gap-fill |
| `Premkumar_Auerbach_1997_JGP_NMDA_stoichiometry_outsideout_patch_oocyte.pdf` | noise gap-fill |
| `Prieto_Wollmuth_2010_JNeurosci_AMPA_gating_modes_outsideout_patch_rms_noise.pdf` | noise gap-fill |
| `Pantazis_Olcese_2024_BiophysRep_digital_amplifier_cutopen_oocyte_noise.pdf` (misnamed, see above) | noise gap-fill |
| `Rycroft_Gibb_2002_JNeurosci_NMDA_calmodulin_outsideout_patch_rms_noise.pdf` | noise gap-fill |
| `Smith_Wang_Howe_2000_JNeurosci_native_AMPA_outsideout_patch_rms_noise.pdf` | noise gap-fill |

Retrieval routes, so nobody repeats the failures. What works when PMC, Wiley, rupress, cell.com and
ASPET all return 403 or HTML: `curl -sL --http1.1 "https://europepmc.org/articles/PMC<id>?pdf=render"`.
The `--http1.1` flag matters, since HTTP/2 gets the stream reset. Also working: publisher direct at
`nature.com/articles/<DOI>.pdf`, `cdn.elifesciences.org/articles/<N>/elife-<N>-vN.pdf`, and
institutional repositories. Dead: the `oa_pdf` paths under `ftp.ncbi.nlm.nih.gov/pub/pmc/`, the AWS
`pmc-oa-opendata` bucket for these accessions, and the Europe PMC `fullTextPDF` REST endpoint (zero
bytes). This unblocks several of the PDFs that `SOURCES.md:299-306` lists as unretrievable.
