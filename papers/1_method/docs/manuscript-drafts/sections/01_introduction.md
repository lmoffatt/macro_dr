# Introduction — the brief for `01_introduction.tex`

> **This file holds no prose and no numbers.** That is the rule that keeps it from going stale. The
> eight section plans it replaces rotted because they carried drafts, paragraph plans and figures
> arcs, all of which the manuscript then superseded. What lives here is policy: the job, the
> constraints that bind this section, what is open or blocked, and the verify-before-submission
> list. Policy is superseded only by a decision, and decisions are logged.
>
> The prose is `01_introduction.tex`, beside this file. The provenance of every number is a `% src:` comment next
> to the claim it supports. The rules that govern **every** section are `README.md` in this
> directory; do not restate them here. Roster, figure set and everything cross-section:
> `../../../decisions.md`.

**The job.** Take a reader who fits the mean by least squares and show them, in their own vocabulary,
what the fluctuations carry, why using them takes a likelihood, and why no one can currently tell
whether such a likelihood is working.

**The five moves** (2026-07-28 voice notes): independence as the foundational assumption of least
squares, and what violating it costs; Markov chains as the way to model dependence while still yielding
record-universal constants tied to biophysical structure; the ladder of methods; why least squares still
reigns (you can see, point by point, whether it predicts, which is visually convincing, while a
recursive filter follows the data closely and holds no surprises, so a working method and a bug look
alike); and why autoregressive alternatives do not solve it.

**Constraints.**
- **Do not re-announce MacroIR.** Cite it as prior work in the same breath as Milescu 2005 and
  Münch 2022. The moment the Introduction explains how it works, the referee says "you published this".
- **The NSFA naming trap.** Non-stationary fluctuation analysis *is* widely used and *does* use the
  variance. Name it and distinguish it or an electrophysiologist will think the gap is already filled.
  Stepanyuk 2014's own words end the objection: *"the unitary current is virtually the only parameter
  that can be reliably obtained from this type of analysis"*, and *"kinetic rates have never been
  estimated for any synaptic receptors in their intrinsic environment"*.
- **On ARIMA**, four points, all derivable from machinery already in the paper: an ARMA error model is
  stationary by construction while the gating covariance tracks the mean current and restarts at every
  jump; it contains no channel count and no unitary current; its timescales are free where the Markov
  model ties them to the rates that generate the mean; and whiteness is not falsification, since enough
  ARMA terms whiten anything.

**Stays out.** The mechanics of MacroIR (that is Theory). The Fisher-to-zero result. The
research-program framing, which reads as a grant proposal and belongs in the Discussion. The Comm Biol
biology beyond one citation. Any claim about experimental data.

**Open.** How much statistics vocabulary in the "why nobody could test this" move: recommendation is
all three terms, each glossed in physical terms in the same sentence, at a cost of about forty words.

**Verify.** Del Core & Mirams 2025 and Owen & Mirams 2025 quotations against the version of record.
**[SETTLED SINCE]** The old top verify item was the "no published likelihood integrates the acquisition
window" claim; it is dead as written and survives only scoped to macroscopic-N by a scaling argument
(`decisions/D-3_novelty_claim.md`).

**[OPEN 2026-08-04] Covariance-fitting cost: square or cube?** `01_introduction.tex:15` says the cost of
`celentano2004use` "grows as the cube of the record length". `moffatt2007estimation` says in its own
Introduction (p.75, and the abstract) that covariance fitting "scales as the square of the number of
samples and therefore is limited in the number of points it can use". Both are defensible depending on
what is counted (building the T x T covariance vs factorising it), but the paper cannot carry both, and
one of the two sources is this author's. Decide, state which operation the exponent refers to, and
propagate. Low priority; blocks nothing.

**[REWRITTEN 2026-08-04]** The section was rebuilt on a single thread after a controlled comparison of
five drafts on five threads (same fact packet, same length, three critique lenses). 2095 words to 1808.
Predecessor verbatim at `../archives/01_introduction_SUPERSEDED_20260804.tex`. What changed, and the
reasons, are in the header of the `.tex` itself so they travel with the prose. Two new verify items:
the "each sample equals the last plus noise" claim is an assertion and needs either its derivation
(the first difference of a mean-square continuous process vanishes with the interval, at a rate set by
the interval axis and the noise ratio) or a naive-forecast citation; and the error-bar factor quoted
at tau_int = 3 needs deciding, since the distortion is a ratio of informations and sqrt(3) = 1.73.

**[CLOSED 2026-08-04 by the rewrite] Two attributions the lineage paragraph got wrong.** (i) The Bayesian
framing is `moffatt2007estimation`, not `munch2022bayesian`: p.76 makes the state probability vector the
observer's uncertainty rather than a population frequency, and says explicitly that this is what allows
it to be used in the deduction of the recursive likelihood. 2007 is Bayesian in the state and maximum
likelihood in the parameters; Munch carries the recursion into Bayesian inference over the parameters
and models state-dependent open-channel noise, which is also the mandatory positioning sentence
(`approach.md` section 10 item 7). (ii) `moffatt2007estimation` already names Micro R, Macro NR and
Macro R and describes Macro NR as the approximation that makes the state distribution continuous and
neglects the dependence between successive measurements, so the roster nomenclature is inherited, not
coined here. Check whether that paper's reference (19) for Macro NR is `milescu2005maximum`; if it is,
NR is the field's published member and not a construction of this study, which is worth saying.

**[OPEN 2026-08-05, Luciano] Hodgkin and Huxley do not give rate constants, they give empirical
parameters.** Against `01_introduction.tex:56-57`: *"the current follows the deterministic relaxation
of the Hodgkin and Huxley description `\citep{hodgkin1952quantitative}`, and rate constants are
obtained by fitting that relaxation, in current practice by least squares."* If the sentence is to say
**rate constants**, it needs a different example.
The conflation is between two senses of the words. H&H's own text does call alpha and beta rate
constants, but they are empirical functions of voltage inside a phenomenological m^3 h formalism, fitted
curve by curve; they are not the rate constants of a kinetic scheme, which is what this paper means
everywhere else and states at `01_introduction.tex:81` (states connected by rate constants that hold
over the whole record). The opening therefore borrows H&H's authority for a claim about a different
object.
This is the **second** H&H attribution problem in the same paragraph. Header note (b) of the `.tex`
already records that the opening no longer attributes least squares to H&H, who fitted by hand. That
fix moved the least-squares attribution off them; it left the rate-constant attribution on them.
Two ways out, decide: (i) keep H&H as the illustration of the many-channel deterministic regime only,
and move the "rate constants by fitting the relaxation" clause onto a macroscopic Markov-scheme fit
(`milescu2005maximum` is the obvious candidate and is already cited elsewhere in this section); or
(ii) drop H&H from the sentence entirely and open the macroscopic regime on the scheme-fitting practice
the paper actually benchmarks. Option (i) preserves the physiological anchor and costs one clause.
Side effect worth knowing: `hodgkin1952quantitative` is one of the four keys currently undefined at
compile time and is in **neither** `biblio.bib` nor `biblio_full.bib`. Under option (ii) that bib gap
closes by itself; under option (i) the entry still has to be created.

**[OPEN 2026-08-05, Luciano] The autoregressive family is dismissed before it is introduced.** Raised
off page 2 of the built PDF. `01_introduction.tex:86` closes the Markov paragraph with *"A correlation
function fitted freely to the residuals delivers neither of those."* That is a verdict on the AR/ARMA
route, and the AR/ARMA route is not put on the table until `:138-149`, fifty lines and one full
subsection later, where it is introduced as one of the *"two alternatives to a mechanistic likelihood
worth placing here, since both are in use"*. The reader meets the dismissal of a method they have not
been offered.
It is also a **duplicate argument, weaker copy first**. Line 86 compresses into one clause exactly what
`:146-148` then makes properly: *"Its timescales are free parameters, where the Markov model ties them
to the eigenvalues of the rate matrix that already fixed the mean, so free timescales absorb the
information the mechanism would have contributed"*, backed there by `lei2020considering` reporting a
widened posterior with no better prediction. The second statement is the one that carries evidence; the
first spends the point before the evidence exists.
Interacts with header note (f) of the `.tex`, which deliberately moved the ARMA and fluctuation-analysis
concessions to *after* the criterion because in the predecessor they cut the line at the point of
maximum tension. Line 86 looks like the residue of that move: the paragraph relocated, its verdict did
not. Read that way this is finishing (f), not reversing it.
Recommendation: **delete line 86**, do not relocate it. The paragraph's point is positive and stands
alone (the constants mean something outside the record, they are properties of the protein, comparable
across conditions and laboratories, and they connect the measurement to a mechanism); it does not need
a foil, and dropping the foil also removes a contrastive close of the kind the register rule in
`README.md` (8) disfavours. If a forward pointer is wanted at all, it should be neutral and name the
comparison as coming, not deliver the outcome.

---
