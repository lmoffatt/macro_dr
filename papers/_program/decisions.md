# Cross-paper decision log

> Updated: 2026-07-20. Split from `macroir-elife-2025/02_decision_log.md`, which logged a single
> paper. **Only decisions that bind more than one paper live here.** A decision that binds one paper
> lives in that paper's folder; §5 lists what was left behind and where it went.
>
> Ledger of **settled** decisions, so we can rewind via git history. Open decisions live where their
> owner does; the program-level ones are in `program.md` §9.

## 1. The program

- **Three papers on one axis** (2026-07-20). Two Gaussian closures plus the question prior to both;
  one paper each. The map, the N_ch partition and the publication order are in `program.md`; do not
  restate them here.
- **Publication order 1 → 2, 3** with 2 and 3 back to back. Valid structurally, not merely convenient,
  **because paper 2's floor sits at 10² channels**; that floor is therefore a decision, not a grid
  parameter (`program.md` §3).
- **The validation machinery is written once**, in paper 1, and cited by papers 2 and 3. This is what
  keeps the three from becoming three copies of the same Methods section.
- **Papers 2 and 3 stay as stubs in `_program/` until they start drafting.** Structure built ahead of
  content rots; `01_workboard.md` retired with all 25 of its checkboxes unticked, including the ones
  whose work had been done.
- **Citation runs one way:** a paper may cite `_program/`; `_program/` never cites a paper; papers do
  not cite each other in the planning layer (`00_index.md` rule 2).

## 2. Model, methods, naming

- Minimal **two-state** model (`scheme_CO`), single K_on/K_off, **non-stationary** protocol
  (single concentration jump). All three papers.
- **Six methods on two levels.** Off the lattice: classical nonlinear least squares on the mean
  current, data key `nonlinearsqr`, display `LSE`, engine flag `family_approximation = 2`. On the
  lattice: `NR`, `NMR`, `R`, `MR`, `IR`, plus `VR` pending its name. Which paper carries which is in
  `program.md` §1.
- **LSE is not a rung of the family** (2026-07-20). In the dispatcher it carries the same two knob
  settings as NMR (`recursive=false, averaging=1`) and differs only by the third flag. The "one object
  with two knobs" framing is retired; the structure is a root question with the ladder hanging from it.
- **`NMR` is dropped from the program**, not relocated: no literature attribution, no mechanistic
  role. It existed to complete a grid that is no longer the argument.
- Naming standardized on **NMR** (scripts have used MNR). Published-name bridge: IR = MacroIR,
  NMR = MacroINR.
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
- **`seed = 0` means random.** It is the sentinel for `std::random_device` and the resolved value was
  never logged, so every simulated ensemble is statistically equivalent but **not bit-reproducible**,
  and cannot be fixed retroactively. Methods must say so plainly in all three papers.
- **Do not pool cells across n_sims.** Every scalar summary of the distortion matrix carries a Jensen
  bias in n_sims. The grid is ragged across the program: band-A cells at 10⁴, the 2026-07-20 fill at
  1000, `433ed13` also holding 200, and the micro cells at 100/1000/10⁴. Hold n_sims fixed within any
  panel or use the debiased quadratic. **This is the likeliest way for a wrong result to reach print**
  and it now sits on more than one paper's headline figure.
- Audio sources and their transcripts are both tracked (author preference).

## 5. Reopened by the split, or left behind

**Reopened.**
- **"One paper = one repo."** Settled when there was one paper: this work carves out to a dedicated
  repo at code freeze, with `macro_dr` referenced by pinned tag. With three papers sharing one engine,
  one machinery and one data tree, the question is now whether that is one program repo or three, and
  it is **not settled**. Owner: `carve_plan.md`.
- **Venue.** Was one question; is now three (`program.md` §9).

**Left with paper 1**, in `1_method/decisions.md`: its thesis and scope sentence; the band-A results
table and D-4's two contested cells; its figure arc and figure count; D-1 (MR and NMR in main text or
supplement, itself reopened by the reframe since "strawman" is a ranking word with no meaning on a
map); D-3 (where the Fisher-to-zero result goes).

**Working plan, not a decision — the Comm Biol erratum (gvar_i).** Program-level, so it stays here.
Intent: disclose, done properly. **Decouple** the erratum (re-run at the same fidelity to isolate the
gvar_i fix's effect on the Bayes factors) from any Bessel-filter high-fidelity redo, since bundling
them confounds attribution. **Triage first**, cheaply, via the distortion correction or by reweighting
the deposited MCMC samples, to learn whether the fix moves the ranking before committing to a re-run.
Low external urgency, but load-bearing: the program uses Comm Biol as its demonstration, and a
programme whose whole machinery exists to catch this class of error is the worst place for it to go
unmentioned.

## 6. Superseded (kept for rewind)

- **One paper** → three (2026-07-20).
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
- "Keep MicroIR out to reduce attack surface" → it is **paper 3**.
- `elife-macroir-merged.tex` as manuscript source of truth → `elife_paper.tex`, which belongs to
  paper 1.
