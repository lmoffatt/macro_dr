# From Molecular Mechanisms to Data and Back

> Framing note for the MacroIR paper and its successors, and the shared "why" the individual papers plug into. Intended readers: biophysicists who might use MacroIR; a distilled version goes into the paper, and this may later stand alone as a framework note. Organizing principle: matrioska. A simple global form on the outside; each piece opens into its own, more detailed study on the inside. The outer form must stay graspable in one read. Comment and correct inline.

## The engine

Biology advances by ordering competing theories according to the evidence the data give them: theories with more evidence are studied further. A biological theory of a molecular mechanism becomes a rankable hypothesis by mapping it, through three bridges, onto a Bayesian evidence:

```
biological theory
      │  bridge 1: conformational-change notation
      ▼
Markov model
      │  bridge 2: MacroIR
      ▼
likelihood (model, experiments)
      │  bridge 3: tempered evidence estimator
      ▼
Bayesian evidence  ──►  ranking of theories
```

That is the whole program on the outside. Each bridge is an inner doll.

The chain runs forward to rank competing hypotheses. Read backward, from evidence to theory, it returns not a single mechanism but the class the data cannot separate (below). Both directions matter, hence the name.

## The three bridges

1. **Conformational-change notation.** Maps a molecular mechanism to the Markov model that represents it. The bridge accepts *any* Markov model, so the familiar kinetic schemes already qualify; the conformational notation is the subunit-resolved generalization that has Castillo-Katz, MWC and KNF as corners (detailed in the bridge-1 doll below). It is a power-user authoring tool, not a required entry point.

2. **MacroIR.** Given a Markov model and a set of experiments, returns the likelihood of the time-averaged macroscopic data. **This paper verifies this bridge.**

3. **Tempered evidence estimator** (name to be decided). Affine-invariant parallel-tempering MCMC plus thermodynamic integration over the tempered ladder, mapping (model, data, likelihood approximation) to the Bayesian evidence.

## What the program resolves, and what it does not

The theory-to-Markov-model correspondence is many-to-many, so evidence is computed per model and a theory's evidence is the marginal over the Markov models it admits, weighted by a within-theory prior. Theories are therefore rankable, with one exception the data cannot fix: theories whose models are observationally equivalent tie, and how tightly they tie is itself something bridge 3 measures. Two independent losses of information sit on the chain:

- **Bridge 2 (aggregation).** The macroscopic observable does not see the whole rate matrix: distinct matrices can produce identical current statistics (aggregated-Markov equivalence).
- **Bridge 1 (mechanism).** A rate matrix does not encode what its states *mean* structurally, so it is the image of an equivalence class of theories.

Inverting the chain (data to theory) therefore returns a *class* of theories, not one: "the data show two open states" constrains the model but is realized by many mechanisms. A single theory can likewise map to several Markov models when it is structurally underspecified, and that one-to-many is productive: it is what generates the candidate set the program ranks (the demonstration's nine schemes are the Markov realizations of a handful of hypotheses). The program's regime softens the aggregation loss: macroscopic, non-stationary observation reshapes identifiability relative to the single-channel bounds, breaking degeneracies unbreakable at equilibrium (Milescu et al. 2005) and adding state-count and topology information from the current covariance (Celentano and Hawkes 2004); the sourced treatment is in `docs/bibliography/identifiability/`. The mechanistic loss remains, so collapsing the theory class needs a second, independent bridge whose degeneracies differ. In the demonstration, single-ATP molecular dynamics broke the left-right coupling ambiguity that electrophysiology alone could not.

## Two phases

- **Phase 1, done.** The P2X2 study (Moffatt and Pierdominici-Sottile, Comm Biol 2025) ran all three bridges end to end on real outside-out patch data and produced biology: asymmetric, sequential activation, Scheme IX favored by Bayes factor > 5000, the classical flip state reinterpreted as an obligatory intermediate. The demonstration is what licenses studying each bridge in isolation: once the assembled chain is shown to deliver on a real mechanism, pinning down a single bridge in a minimal model is clearly worth doing, and its failure in isolation matters because the whole generates interpretable results on real data.

- **Phase 2, now.** Characterize the robustness of each bridge, one at a time, in the cleanest setting that still exercises it.

## Where this paper sits

Bridge 2 (MacroIR), minimalist context: two states, non-stationary. Goal: its validity limits and an approximation to its validity over channel count, interval length and instrumental noise. It backs, rigorously, a claim the demonstration already relied on: that the recursive interval likelihood gives different, and better-calibrated, evidence versus non-recursive and single-point approximations once integration windows exceed channel timescales (Comm Biol Fig 1d vs 1e). It also exposes a sharp identifiability limit of the macroscopic observable itself: the information about the original number of channels vanishes once the open population stops rising and relaxes. The published biology rests on bridge 2's fidelity, so nailing it down is load-bearing, not a toy.

## The other dolls (later steps)

- **Bridge 2, richer regimes:** interaction of kinetic constants across timescales (more states), the stationary regime, and MicroIR for the few-channel multinomial regime.
- **Bridge 3:** two questions. First, whether the evidence is a calibrated model-discriminator. Bayesian model selection is calibrated by construction under an exact likelihood and exact integration, so the real test is whether the approximations break that calibration: simulate data from alternative models and check whether the evidence's confidence matches the empirical frequency of confusing them (that confusion frequency is the quantitative width of the theory equivalence class from the caveat above). The evidence correction (volume ½ log det C, effective-sample α⋆ = p / tr C) is the theory of how the likelihood distortion breaks the calibration, stated as motivation in the bridge-2 paper and derived in this bridge-3 study; the simulation is its empirical check. Second and in parallel, which method computes the evidence most efficiently. The move is self-similar: bridge 2 validates the likelihood against simulation, bridge 3 validates the evidence against discriminability, one level up.
- **Bridge 1:** robustness of the conformational-change notation, and the theory-to-model correspondence itself (the caveat above). The notation generalizes the familiar schemes along two knobs they leave fixed, subunit resolution and a continuous binding-to-gating coupling (Castillo-Katz is one lumped subunit with obligatory coupling, gating only when bound; MWC is concerted, KNF sequential); in the demonstration it produced schemes VI-IX (Scheme IX: 24 states). That continuous coupling is data-forced, not aesthetic: Castillo-Katz forbids unliganded opening, yet the demonstration records spontaneous partial activation without ligand, which only a continuous coupling reproduces.
- **Experimental and sampling design** wraps all three (interval spacing, channel pooling, noise budget).

## Two names to reconcile

- Bridge 2 family: this paper uses NR, NMR, R, MR, IR; the published Comm Biol uses MacroIR and MacroINR. IR = MacroIR, NMR = MacroINR. Bridge the names so a reader coming from the published paper is not lost.
- Bridge 3 has no settled name yet.






