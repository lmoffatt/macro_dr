# Cross-paper decision log

> Updated: 2026-07-28. Split from `macroir-elife-2025/02_decision_log.md`, which logged a single
> paper. **Only decisions that bind more than one paper live here.** A decision that binds one paper
> lives in that paper's folder; §5 lists what was left behind and where it went.
>
> Ledger of **settled** decisions, so we can rewind via git history. Open decisions live where their
> owner does; the program-level ones are in `program.md` §9.

## 1. The program

- **TWO papers, split on the macro/micro boundary** (2026-07-23, superseding the three-paper split of
  2026-07-20). The former papers 1 and 2 are **merged into one broad macro paper** aimed at eLife;
  the micro paper (the multinomial boundary) stays separate. The map and the N_ch partition are in
  `program.md`; do not restate them here.
- **Why the merge, since it reverses a decision three days old.** The alternative was macro-method to
  eLife and a near-identical usage-map paper after it. Desk rejections at one journal are **not
  independent draws**: same triage editor, same topic, and a real salami perception, so two shots at
  eLife is one poisoned shot. Two near-identical papers is also redundant publication, which is a
  problem prior to strategy. The merged paper's breadth lives in the **combination** (least squares
  convicted, the ladder, and the map), not in either half alone.
- **What unblocked it technically:** eLife has no display-item limit, so the R/MR/VR ladder detail
  goes to a supplement; and **faceting by N_ch dissolves the two-axis problem**, because at fixed N_ch
  the single-channel noise scale and the fraction-of-total are the same variable relabelled. Methods
  must state one convention and derive the other per facet.
- **The comparison anchor is least squares** (2026-07-28). MacroR has essentially no uptake, so a
  paper comparing one unused algorithm against another is unsellable. The paper compares the method
  the field actually uses against the new one. LSE is a **full body citizen** across the figures, not
  a cited background baseline.
- **The validation machinery is written once**, in the macro paper, and cited by the micro paper.
- **The micro paper stays a stub in `_program/` until it starts drafting.** Structure built ahead of
  content rots; `01_workboard.md` retired with all 25 of its checkboxes unticked, including the ones
  whose work had been done.
- **Citation runs one way:** a paper may cite `_program/`; `_program/` never cites a paper; papers do
  not cite each other in the planning layer (`00_index.md` rule 2).
- **Venue: eLife is the single shot**, with Biophysical Journal as a clean fallback. Münch et al. 2022
  was published in eLife under the older, harder model, which is evidence the topic clears triage and
  simultaneously raises the novelty bar: **the abstract must position against it.** What is offered is
  the test of whether such filters tell the truth about their own uncertainty, plus the map of where
  they fail, not another filter.

## 2. Model, methods, naming

- Minimal **two-state** model (`scheme_CO`), single K_on/K_off, **non-stationary** protocol
  (single concentration jump). Both papers.
- **Five methods on two levels.** Off the lattice: classical nonlinear least squares on the mean
  current, data key `nonlinearsqr`, display `LSE`, engine flag `family_approximation = 2`. On the
  lattice: `NR`, `R`, `MR`, `IR`, plus `VR`. Which of them sit in the body is in `program.md` §1.
- **LSE is not a rung of the family** (2026-07-20). In the dispatcher it carries the same two knob
  settings as the dropped `NMR` (`recursive=false, averaging=1`) and differs only by the third flag.
  The "one object with two knobs" framing is retired; the structure is a root question with the ladder
  hanging from it. **This does not make LSE peripheral** (§1): it is the root of the ladder and the
  paper's anchor.
- **`NMR`: the drop is REOPENED (2026-07-29), and the reason it was dropped is false.**
  - The 2026-07-20 reason was "no literature attribution, no mechanistic role". **The first clause is
    wrong and has been wrong the whole time.** `NMR` **is** `MacroINR`, and `MacroINR` is the
    **published control** of Moffatt & Pierdominici-Sottile 2025 (Comm Biol), the method whose failure
    carries that paper's central methodological claim: *"Model ranking was sensitive to likelihood
    approximation: the control method (MacroINR) systematically underestimated evidence for schemes
    with conformational intermediates"*, and *"leading to systematically different evidence values
    (Supplementary Table S1)"*. That is attribution of the strongest kind: not a citation, a published
    demonstration.
  - **The published-name bridge `NMR = MacroINR` is CORRECT** (Luciano, 2026-07-29), closing the open
    item in `nomenclature.md`. `MacroINR` parses as **I**nterval (the averaged conductance, `av = 1`)
    + **N**on-**R**ecursive, which is exactly `NMR`. The Comm Biol Introduction's phrase "ignores time
    averaging" is loose prose for "does not do IR's boundary-conditioned interval treatment"; it is not
    a statement about the `av` flag, and it must not be read as one.
  - **What survives is the measured reason**: recomputed on the freeze, NMR is **numerically
    indistinguishable from NR** (envelope 1.23 to 3.6 × 10⁴ against NR's 1.21 to 3.5 × 10⁴; 0 of 336
    points within ±15% for both). Evidence in `../1_method/decisions/D-4_ranking_verdict.md` §5.
  - **That redundancy is itself a result, not only a reason to cut.** If the interval-mean conductance
    buys nothing without recursion, then recursion is the step that matters at the bottom of the ladder
    and interval conditioning is the step that matters at the top. That is the ladder's shape, measured.
  - **Open, and it is a scope decision, not a naming one:** body column (no: it duplicates NR), cut
    entirely (costs the bridge to the paper's only real-data demonstration), or **named with its
    published attribution and measured once in a supplement beside NR, where the near-identity is the
    finding**. The third is the standing recommendation.
  - **This entry supersedes the "six methods" listing that stood two bullets above it until 2026-07-28**
    and contradicted it from within the same section; the count is five or six depending on how the
    item above resolves.
- **`VR` keeps its letter** (Luciano, 2026-07-28), closing the open item in `program.md` §9. It earned
  it: VR was predicted to come out over-confident and more so than MR, and it did. The `V`/Taylor
  collision warning in `nomenclature.md` still stands and Methods must carry the one sentence.
- Naming standardized on **NMR** (scripts have used MNR) for the historical record only.
  Published-name bridge: IR = MacroIR.
- `nonlinearsqr` must appear **verbatim** end to end (`.macroir` label → CSV `algorithm` cell → R
  `ALGOS` entry). A mismatch silently drops rows, the same failure class as the MNR/NMR bug.

## 3. Method and anchor

- **Gaussian Fisher** is the distortion anchor; the numerical finite-difference Fisher only gauges how
  good the Gaussian one is.
- Distortion and bias evaluated at the **optimum / θ_pool**; θ_sim used to expose the bias.
- Diagnostics: residual mean/variance/whiteness; score bias; Var[score] against the Gaussian Fisher →
  the distortion matrix as a symmetric sandwich, decomposed into correlation and sample/geometric
  parts; plus the direct empirical-vs-sandwich covariance test. Definitions, sign conventions and
  thresholds: `machinery.md`.
- **Likelihood-only.** MLE / Gauss-Newton local maximum kept only to obtain the empirical parameter
  covariance. The posterior information-distortion framework and full model-comparison results are a
  later program component; the likelihood-side evidence correction stays as **motivation**, derivation
  deferred.

## 4. Data and provenance

- **D-0 (2026-07-15):** freeze at `1c2ae6f`; **multi-commit provenance accepted**, each CSV
  self-stamping its engine hash. `433ed13` kept as the numerical-Fisher equivalence demo. E-1…E-5
  decoupled to `main` as code hygiene.
- **No paper number may be quoted from `433ed13`** (2026-07-28). It is the demo, not the basis. The
  ranking verdict was computed there and has been **recomputed on `1c2ae6f` + `87889e6` against the
  Gaussian Fisher**; the verdict survives and the magnitudes move by five to ten percent
  (`../1_method/decisions/D-4_ranking_verdict.md`, scripts committed under `decisions/recompute/`).
  That the two anchors agree is exactly what `433ed13` was run to show: state it once in Methods and
  quote the freeze everywhere else.
- **The non-IR noise columns are on `87889e6`, not on the freeze.** `1c2ae6f` carries R, MR and NR at
  noise 0.1 only (plus R at 100); noise 1 and 10 for those three landed with the D-0 fill on
  `87889e6`. Any grid-wide statement must read both directories or it silently reports a
  single-noise slice as the whole plane.
- **`seed = 0` means random.** It is the sentinel for `std::random_device` and the resolved value was
  never logged, so every simulated ensemble is statistically equivalent but **not bit-reproducible**,
  and cannot be fixed retroactively. Methods must say so plainly in both papers.
- **Do not pool cells across n_sims.** Every scalar summary of the distortion matrix carries a Jensen
  bias in n_sims. The grid is ragged across the program: band-A cells at 10⁴, the 2026-07-20 fill at
  1000, `433ed13` also holding 200, and the micro cells at 100/1000/10⁴. Hold n_sims fixed within any
  panel or use the debiased quadratic. **This is the likeliest way for a wrong result to reach print**
  and it now sits on more than one paper's headline figure.
- Audio sources and their transcripts are both tracked (author preference).

## 5. Reopened by the split, or left behind

**Reopened.**
- **"One paper = one repo."** Settled when there was one paper: this work carves out to a dedicated
  repo at code freeze, with `macro_dr` referenced by pinned tag. With two papers sharing one engine,
  one machinery and one data tree, the question is now whether that is one program repo or two, and
  it is **not settled**. Owner: `carve_plan.md`.
- **Venue.** Was one question; is now three (`program.md` §9).

**Left with paper 1**, in `1_method/decisions.md`: its thesis and scope sentence; the band-A results
table and the two contested cells of `1_method/decisions/D-4_ranking_verdict.md`; its figure arc and
figure count; **Q-2** (MR and **VR** in main text or supplement, itself reopened by the reframe since
"strawman" is a ranking word with no meaning on a map); **Q-3** (where the Fisher-to-zero result
goes). Those two were `D-1` and `D-3` until 2026-07-21, when paper 1's open questions were relabelled
`Q-n` to stop colliding with the `D-n` decision briefs; and it read "MR and NMR" until the same date,
which was stale — NMR is dropped from the program (§2) and VR is the method in question.

**Working plan, not a decision — the Comm Biol erratum (gvar_i).** Program-level, so it stays here.
Intent: disclose, done properly. **Decouple** the erratum (re-run at the same fidelity to isolate the
gvar_i fix's effect on the Bayes factors) from any Bessel-filter high-fidelity redo, since bundling
them confounds attribution. **Triage first**, cheaply, via the distortion correction or by reweighting
the deposited MCMC samples, to learn whether the fix moves the ranking before committing to a re-run.
Low external urgency, but load-bearing: the program uses Comm Biol as its demonstration, and a
programme whose whole machinery exists to catch this class of error is the worst place for it to go
unmentioned.

## 6. Superseded (kept for rewind)

- **One paper** → three (2026-07-20) → **two** (2026-07-23, §1). The 1|2 cut was the artificial one
  (mid against high N_ch, single-channel noise scale against fraction-of-total) and is deleted. The
  surviving split is the natural one, the Gaussian macro closure against the multinomial micro one.
- **`paper-2.md`** → archived 2026-07-28 under `archive/`, with a migration ledger. It is the ancestor
  of both the current body roster and the region map, so it is archived rather than deleted.
- **A-strict** ("the non-recursive members are named once in Theory and measured in no figure") →
  **DEAD**. NR is measured, and it is in the body.
- **"Paper 1's own novelty is the VR mechanism"** → demoted. The headline is the machinery plus the
  validity map; VR is a supplement that confirms a prediction.
- **"The temporal correlation of macroscopic currents is unused"** → **false, do not write it**
  (2026-07-28). It has been used since 1973: the closing rate from the Lorentzian power spectrum
  (Katz & Miledi 1970/1972; Anderson & Stevens 1973), the exact nonstationary two-time covariance
  (Conti et al. 1980; Sigworth 1981), and a covariance likelihood three years before MacroR
  (Celentano & Hawkes 2004). The defensible gap is that **nobody characterised when these methods are
  valid**. Sources and the replacement paragraph:
  `docs/bibliography/temporal_correlation_and_AR_errors_2026-07-28.md` and §A.10 of
  `docs/bibliography/MacroIR_prior_art_map.md`.
- **"10 to 16" as the non-recursive overconfidence factor** → **10 to 15**, once NMR is dropped and
  the anchor moves to the freeze (`D-4` §3). "14 to 21" was never the ellipse-area factor at all.
- **"Five algorithms" as the closed roster** → six methods on two levels (2026-07-20). LSE was
  previously present only as cited background describing what the field does; it is now a measured arm,
  in paper 2.
- **The ranking as the deliverable** → a usage map per paper (2026-07-20). Retired phrasings, so they
  are not re-copied out of the older documents: *"IR sole survivor"*, *"only MacroIR stays calibrated
  across the practical regime"*, *"MR strawman"*. The first two state a one-band result globally; the
  third is a ranking word with no meaning on a map, where a method has a domain, possibly empty.
- **The noise axis as a nuisance third dimension** → the organizing axis, because it carries the two
  crossovers that decide which method is needed (`axes.md`).
- Validity-map axes N_ch / Δτ / Noise → briefly **N_ch × K_off** → corrected 2026-07-15 back to
  **N_ch × noise** with interval internal. K_off was never swept (fixed at 100 on disk) and the
  K_off framing never had data behind it.
- LID↔Evidence as a standalone finding (Δlog Z = ½ log det C) → the later component's **motivation**,
  likelihood-side, derivation deferred.
- "Keep MicroIR out to reduce attack surface" → it is **the micro paper**.
- `elife-macroir-merged.tex` as manuscript source of truth → `elife_paper.tex`, which belongs to
  the macro paper.
