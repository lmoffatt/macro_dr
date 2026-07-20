# Discussion — what it must do, and a paragraph plan

> Working doc, same genre as `abstract.md`. Opened 2026-07-14. Comment and correct inline.
> Depends on: `introduction.md` (the two open questions it must pay off), `results.md` (the numbers it may claim), `00_plan.md` §4 (the ranking).

## The job

The Results say what happened. The Discussion has to say **what it means for someone who is about to fit a macroscopic recording next week**, and it has to close the two loops the Introduction opened. Everything else is optional.

Three things it must deliver, in this order of importance:

1. **A decision rule.** Not "MacroIR is best" but the conditions: how many channels, what acquisition interval relative to the relaxation, how much instrumental noise, and what the cheap approximations cost you in that corner. This is the paper's utility and it is what gets cited.
2. **The two answers.** Milescu 2005 asked whether non-filtering estimates are intrinsically biased. Moffatt 2007 said a time-averaged formulation was needed. Both get answered here, explicitly, with the papers named.
3. **The honest perimeter.** Two states, non-stationary, simulation. Plus the one perimeter the reviewer will find if we do not state it: the Gaussian Fisher is itself the model's own object, so the anchor is internal (see below, it is a real issue and it has a real answer).

## The trap

The Discussion is where a paper drifts back into being about its method. Every sentence that begins "MacroIR is..." is a candidate for deletion. The subject of this Discussion is *the distortion*, and MacroIR is where the distortion happens to be small. Keep the grammar that way and the paper stays what it claims to be.

## Paragraph plan

### D1 — The verdict, as a decision rule

Lead with the finding, in the register of the abstract's last sentence: the choice of likelihood is now a matter of evidence.

Then the ranking, with the *mechanism* attached to each entry, because a ranking without mechanism is a benchmark table and a ranking with mechanism is an explanation. The current table (`00_plan.md` §4, to be re-checked against the final Gaussian rerun):

| Algorithm | In the literature | Bias | What it does to the reported information | Verdict |
|---|---|---|---|---|
| `IR` (MacroIR) | Moffatt & Pierdominici-Sottile 2025 | none | ≈1 in the Gaussian regime, up to ~1.3 in the corners | calibrated across the practical regime |
| `R` | Moffatt 2007, Münch 2022 | none | residual factor ~2 | usable; correlation distortion at short Δ |
| `MR` | (intermediate) | none | overestimates variance | drops the boundary cross-covariance N·γᵀΣγ that IR keeps |
| `NR` | ≈ Milescu 2005 | biased | inflation ∝ N (hundreds) | often does not even reach the MLE |
| `NMR` | (intermediate) | none | underestimates | unbiased; possible speed niche (non-recursive, parallelizable) |

**Sign discrepancy, unresolved.** `00_plan.md` says MR *overestimates* variance (under-confident); the `Figure_S3` caption reports MR as *over-confident* (empirical/Fisher 1.5 to 2.1). One of the two has the sign backwards, and the bidirectionality claim in the abstract currently rests on this. **Resolve against the Gaussian-Fisher data before this table is written.** (Same flag as `abstract.md`, "Verify before submission".)

The mechanism sentence that makes the table an explanation rather than a list: every member makes the same two Gaussian approximations (occupancy multinomial → Gaussian; single-interval signal → Gaussian), and the differences among them are entirely in **how much of the acquisition interval the hidden state is conditioned on**: no endpoints, one endpoint, both endpoints (the endpoint ladder, `../_program/nomenclature.md`). The distortion falls monotonically along that ladder. That is why the verdict feels structural rather than empirical.

### D2 — Milescu's question, answered

Milescu 2005: *"It remains an open question whether all estimates should be intrinsically biased if obtained by any method that is not a Bayesian filtering algorithm."*

The answer we can support is more interesting than a yes:

**Recursion is not what makes the estimator unbiased. It is what makes the reported uncertainty correct.** `NMR` is non-recursive and unbiased; `NR` is non-recursive and biased. So filtering is not the dividing line for bias. But every non-recursive member misreports the information it claims to carry (NR inflating it, NMR deflating it), because with no conditioning the likelihood mistakes the sequential structure of the data for independent evidence. The two questions Milescu ran together, "is the estimate biased" and "is the estimate trustworthy", separate cleanly under the distortion matrix, and only the second one is really about filtering.

And his own diagnosis lands: Milescu named his dominant error source as *"the local time correlation of the current"*. That is the **correlation distortion** of our decomposition, now measured rather than suspected. Say so, and cite him. This is the paper's most gratifying single connection and it costs two sentences.

### D3 — The averaging question, answered, without re-announcing MacroIR

Moffatt 2007 said a time-averaged formulation was *"fundamental to the widespread use"* of the recursive algorithms. That formulation exists (Comm Biol 2025). **What this paper adds is not the formulation but the measurement of what its absence costs**, and the phrasing has to keep that distinction sharp or a reviewer will say the paper re-presents a published method.

The honest form of the claim: treating each sample as an instantaneous observation costs a distortion that grows with the number of process samples per acquisition interval, continuously, and it is largest exactly where modern acquisition sits (fast gating, finite bandwidth). It does not "break down at" a threshold.

**Concede the device, here, once.** The prior-art map (`docs/bibliography/MacroIR_prior_art_map.md`) establishes that MacroIR is, formally, an integrated-measurement augmented Kalman filter (verified to ~1e-8), a device known in control and in target tracking. The novel part is the realization for an exact continuous-time Markov chain of channel occupancies, its first use in this field, and now its characterization. Conceding this in the Discussion is strictly better than having a reviewer find it: it is a small concession, and volunteering it buys credibility for the claims we do not concede.

**[Q]** Does the concession go in the Discussion or the Theory section? Recommendation: Discussion. In Theory it interrupts a derivation; in the Discussion it is a scholarly note.

### D4 — Answering the cost argument, which is live

Del Core and Mirams 2025 argue against filtering on cost: the moment integration and the correction step run between every pair of consecutive measurements. That is true, and the paper should say it is true. The question they could not answer, and this paper can, is **what the saving buys and what it costs**, in units of misreported uncertainty.

The shape of the answer (fill from the final maps): in the Gaussian corner, with many channels and intervals not much shorter than the relaxation, the cheap approximation's distortion is small enough that its confidence intervals are usable, and the filter's cost is not obviously worth paying. Toward the multinomial and telegraphic corners the distortion grows to where a non-recursive fit reports intervals that are wrong by more than an order of magnitude in the information, which no amount of compute saving justifies. **State both halves.** A paper that concludes "always use the expensive one" is read as advocacy; a paper that says exactly when the cheap one is fine is read as evidence, and it is also more likely to be true.

`NMR` is the interesting entry here: unbiased, non-recursive, parallelizable. If the speed niche is real it belongs in this paragraph, and it is the most useful thing the paper can offer to the group that is arguing against filters.

### D5 — Where the information lives, and the identifiability limit

The Fisher-to-zero result. Once the open-channel population stops rising and relaxes, the per-step information about the number of channels, the association rate and the unitary current falls to zero, while the dissociation rate stays informative through the decay.

Two reasons this is a Discussion item and not a footnote:

1. It is a statement about **the macroscopic observable**, not about any algorithm. It survives every approximation in the family, so it is the one result in the paper that no future better likelihood can overturn.
2. **Del Core and Mirams 2025 declare exactly this problem open.** Landing a concrete answer on a 2025 open problem is worth a paragraph.

Practical consequence, which is what the electrophysiologist takes home: the informative part of a macroscopic recording is the *transient*, and prolonging the plateau adds data without adding information about N, k_on or the unitary current. That is a design statement about experiments, and it is the closest this paper gets to advice about the bench.

Author's own reaction in the 11 July audio, kept as a note on placement: *"me voló la cabeza"*, but *"no creo que vaya el abstract"*. Agreed on both. Main text, Discussion, prominent. Not the abstract.

### D6 — What the distortion matrix is for, beyond error bars

One paragraph, forward-pointing, no derivation.

Near the maximum, the same matrix that corrects the error bars propagates into the Bayesian evidence, through a volume term (½ log det C) and an effective-sample rescaling (α⋆ = p / tr C). So the distortion is not only a statement about confidence intervals: it is the thing that decides whether model comparison under an approximate likelihood ranks mechanisms correctly. That is the next study (`00_plan.md` §3; the program's bridge 3), and this paper supplies its input, which is why the matrix and not a scalar summary is the object we report.

**Restraint required.** Do not state a Bayes-factor correction as a result. It is motivation, and the derivation is deferred. One paragraph, no equations beyond the two scalars, and an explicit "derived elsewhere".

### D7 — The perimeter, stated before the reviewer states it

Four limits, each with the reason it is a *choice* and where it goes next.

- **Two states.** The point of the minimal model is to isolate the algorithm's own error, with no confounding from scheme misspecification. Richer schemes introduce interacting timescales, and the natural expectation, which we do not test, is that the distortion grows where the interval straddles more than one relaxation.
- **Non-stationary only.** The protocol is a single concentration jump. The stationary regime is a different information structure (no transient, so by D5 the informative content differs qualitatively) and is a separate study.
- **Simulation only, no experimental data.** Deliberate: the ground truth has to be known for any of these diagnostics to mean anything. Biophysical relevance of the likelihood family is already established elsewhere (Comm Biol 2025). The map is a statement about algorithms, and algorithms are the kind of thing you can only audit against a known truth.
- **The few-channel corner is where the whole family degrades**, because the multinomial-to-Gaussian step is the first approximation to fail. That corner needs a different solver (the microscopic-recursive line of work), and the map is what says how few is too few.

### D8 — The generalization, kept modest

The machinery (simulate exactly, test the score, sandwich the Fisher, read off the distortion) is not specific to ion channels. It applies to any approximate likelihood for a process that can be simulated exactly and whose likelihood cannot be computed exactly, which is most of stochastic biophysics and much of systems biology.

Say it in **two sentences and stop**. It is true, it is worth saying, and it is exactly the kind of claim that, over-extended, reads as promotional and invites a reviewer to ask for a second application domain that the paper does not have.

## The anchor problem (a reviewer will raise it; have the answer ready)

The distortion matrix is anchored on the model's **own** Gaussian Fisher, so a sceptic can say: you are comparing the approximation to itself. The answer, which must be *in the paper*, not just in our heads:

- The Gaussian Fisher is the information the likelihood **claims** to carry; the score covariance is the information it **actually** delivers under data from the true process. The two are different objects and the mismatch between them is exactly what misspecification means. The anchor being internal is not a circularity; it is the definition.
- Independently, the numerical finite-difference Fisher was computed to check how faithful the Gaussian one is (the F-vs-G bridge, `../_program/machinery.md`), and it is *not* used as the anchor because it is itself numerically fragile in this regime (finite-difference bias makes it indefinite; see the figure-2 diagnostics).

Give this half a paragraph in the Discussion or a boxed note in Methods. **[Q] which?** Recommendation: Methods for the mechanics, one sentence in the Discussion that names the objection and points at Methods.

## Open questions

- **[Q] The Comm Biol erratum.** The published P2X2 validation was run with the buggy per-state variance (`gvar_i`), fixed since. This paper cites Comm Biol as the demonstration that licenses the component study, so the two are linked in print. Options: (a) say nothing here and handle the erratum separately; (b) one neutral sentence noting that the demonstration's numbers are being re-checked under the corrected variance. This is a strategic call, not a writing call, and it should be decided before the Discussion is drafted, not after. (See `../_program/decisions.md` §5, Comm Biol erratum.)
- **[Q] Does the research-program framing (three bridges) get a paragraph?** Recommendation: yes, compressed into two or three sentences inside D6, not a section of its own. The program is why the paper exists but the paper must stand without it.
- **[Q] Ordering of D4 and D5.** D5 (the Fisher-to-zero result) is the more striking, D4 (the cost argument) is the more useful. Current order puts utility first, which suits a Discussion that is trying to be cited by practitioners. Reversing it suits a Discussion that is trying to be remembered.

## Verify before submission

- Every number in the D1 table, against the final Gaussian rerun. The MR sign in particular.
- The claim that NMR is unbiased. It carries the "recursion is not the dividing line for bias" argument in D2, which is the Discussion's best paragraph; if it turns out NMR is biased too, D2 collapses to "yes, Milescu was right", which is still publishable but much less interesting.
- That Del Core and Mirams 2025 do in fact declare the channel-number identifiability problem open (quoted in `abstract.md`); re-read the passage in context before building D5's second reason on it.
</content>
