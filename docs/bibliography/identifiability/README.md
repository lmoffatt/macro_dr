# Two-Level Non-Identifiability in the MacroIR Program

This folder collects the sourced literature behind a theoretical caveat of the MacroIR research program. The program runs a pipeline: biological theory → Markov (rate-matrix) model → likelihood of the observed current (MacroIR) → Bayesian evidence → ranking of competing theories. The caveat is that the correspondence at each arrow is many-to-many, so inverting data back to a theory does not return one theory. It returns an equivalence class.

## The decomposition

Two distinct sources of non-identifiability sit in series, and the inversion (data → theory) composes both.

**(A) Mechanistic non-identifiability: theory ↔ rate matrix.** A rate matrix Q is a labelled graph of states and transition rates. It does not carry the molecular meaning of a state (which conformation, which subunit configuration, which binding occupancy). Different molecular stories can therefore land on the same Q, and a single Q can be narrated by many mechanistic accounts. An underspecified theory (one that fixes connectivity but not all rates, or fixes topology up to symmetry) maps to a whole family of matrices rather than one.

**(B) Aggregation non-identifiability: rate matrix → observable.** The current does not report the full state. It reports only which conductance class a channel occupies (open vs shut, or a small set of conductance levels). The observable is an aggregated Markov process: states are lumped into classes, and only the class trajectory is seen. Distinct rate matrices can produce identical observable statistics. This is the classical aggregated-Markov equivalence result.

Because the inversion travels data → observable → Q → theory, it inherits the fibres of both maps. The honest output of fitting is an equivalence class of theories consistent with the data, whose width is set by (A) composed with (B).

## What IS identifiable

The aggregation literature is not only negative. It counts and canonicalises what the data do fix. Notation: let n_O be the number of open (conducting) states, n_C the number of shut states, and R the rank of the open-to-shut transition submatrix of Q.

- **Fredkin, Montal & Rice (1985)** and the companion Fredkin & Rice (1986): single-channel open/shut dwell-time data identify at most **2·n_O·n_C** functions of the rates, not the individual rate constants. This is the canonical identifiability bound.
- **Kienker (1989)**: the rate matrices giving identical single-channel statistics are exactly those related by similarity transforms that preserve the conductance classes. The equivalence class is a group orbit, and its identifiable dimension is again the Fredkin count.
- **Bruno, Yang & Pearson (2005)**: a canonical "manifest interconductance rank" (MIR) form; single-channel dwell-time data identify exactly **2R(n_O + n_C − R)** rate constants, sharpening the Fredkin bound and giving an explicit parameterisation of the identifiable subspace.
- **Larget (1998)**: a unique minimal canonical representative for each equivalence class of observationally identical aggregated Markov processes. This is the constructive answer to "what does the data actually pin down": the class, named by its canonical form.
- **Colquhoun & Hawkes (1981, 1982)**: the forward theory. Dwell-time densities and burst/cluster statistics are given by the spectra of the Q-blocks partitioned by conductance class. This is the generative map whose non-injectivity the identifiability results characterise.

The positive content is: data fix the equivalence class (Larget's canonical form) and a counted number of parameters (Fredkin, Bruno-Yang-Pearson), and nothing finer.

## "One theory → many Markov models" in three senses

The (A) side deserves separation into three distinct meanings, because only one of them is scientifically productive.

1. **Parametric (trivial).** A theory that leaves rate values free maps to a manifold of matrices indexed by those values. This is ordinary parameter estimation and not the concern here.
2. **Structural underspecification (the productive one).** A theory constrains topology only up to choices it does not resolve: which state is the entry point, whether two pathways are ordered or parallel, whether a desensitised state hangs off the open or the shut side. Each resolution is a distinct candidate scheme. This is the generator of the candidate set that the Bayesian ranking is meant to adjudicate. The non-identifiability here is the reason the ranking step exists.
3. **Discretisation granularity.** A continuous or high-dimensional conformational landscape can be coarse-grained into Markov schemes at different resolutions. Coarser and finer schemes can be observationally indistinguishable, so the "number of states" is itself only identifiable up to the aggregation bound (see Celentano & Hawkes below on inferring state number from macroscopic covariance).

## The macroscopic-current angle

The results above (Fredkin, Kienker, Bruno-Yang-Pearson, Larget, Flomenbom-Silbey) are stated for **single-channel dwell-time trajectories**. The MacroIR program observes something different: many channels in parallel, with the current integrated over each sampling interval. That changes what is identifiable, and the shift is the reason macroscopic observation is worth modelling explicitly.

- **Flomenbom & Silbey (2006)** make the single-channel ceiling sharp: a two-state (open/shut) trajectory fixes only an equivalence class of schemes, never the rate matrix. This is the limiting case the macroscopic regime tries to escape.
- **Milescu, Akk & Sachs (2005)**: maximum-likelihood estimation from the **mean and variance of macroscopic currents** under time-varying (driving) stimuli. Non-stationary protocols break degeneracies that are unbreakable at equilibrium. Concretely, the equilibrium C-C-O versus C-O-C ambiguity (a Kienker-type equivalence) is resolved by a relaxation protocol. Macroscopic observation plus nonstationary driving enlarges the identifiable set beyond the single-channel bound.
- **Celentano & Hawkes (2004)**: fitting the **covariance matrix** of macroscopic currents identifies the number of states and the topological position of desensitised states relative to the gating core, though it still cannot distinguish linear from branched arrangements of those states. So the macroscopic second moment adds identifiability (state count, desensitisation topology) but does not eliminate the class.

The Fisher-information / information analysis in the accompanying MacroIR paper is precisely the quantitative version of this: it measures how much the macroscopic, interval-averaged current constrains the parameters, which is exactly the width of the (B) aggregation-level non-identifiability in the program's actual observation regime. Where the Fisher information is near-singular in a direction, that direction is an unresolved fibre of the aggregation map; where it is well-conditioned, the macroscopic protocol has broken a degeneracy that single-channel data could not.

## Resolution strategy

Non-identifiability inside a modality is broken by importing a second modality whose own non-identifiability lies in different directions. When two maps with different fibres are intersected, the surviving class can collapse to a point. The concrete instance is the authors' own P2X2 work (Moffatt & Pierdominici-Sottile, Communications Biology 2025): a left-right coupling ambiguity that was macroscopically unidentifiable, showing up as a bimodal posterior, was resolved by importing molecular-dynamics evidence. The MD bridge constrains the mechanistic (A) side directly (it speaks about conformations, not just conductance classes), which the current alone cannot reach. This is the operational answer to the caveat: report the class honestly, then collapse it with an independent bridge.

## Honesty note on verification

The **citations** in the index below are verified (bibliographic details, DOIs, venues). The one-line **claims** attached to each work are my reading of that work's result as it bears on (A) and (B); they are faithful summaries of the sourced results but are compressed, and the mapping of each result onto the MacroIR program's specific macroscopic, interval-averaged, many-channel regime is this program's synthesis rather than a claim made verbatim by the original authors. In particular: the single-channel identifiability bounds (Fredkin, Kienker, Bruno-Yang-Pearson, Larget, Flomenbom-Silbey) are established results for dwell-time data; their transfer to the macroscopic regime is only partial and is exactly what Milescu et al. and Celentano & Hawkes address empirically. No claim here should be read as asserting that macroscopic observation removes non-identifiability in general. It reshapes it.

## Index

| Work | One-line | Verified | File |
|---|---|---|---|
| Colquhoun & Hawkes (1981, 1982) | Forward theory of the aggregated open/shut observable; dwell-time densities from Q-block spectra. | Yes | [ColquhounHawkes_1981_1982_aggregated_markov_observable.md](./ColquhounHawkes_1981_1982_aggregated_markov_observable.md) |
| Fredkin, Montal & Rice (1985) | Single-channel open/shut data identify at most 2·n_O·n_C parameters, not the individual rates. | Yes | [Fredkin_Montal_Rice_1985_Aggregated_Markov_Identifiability.md](./Fredkin_Montal_Rice_1985_Aggregated_Markov_Identifiability.md) |
| Kienker (1989) | Conductance-class-preserving similarity transforms give identical single-channel statistics; identifiable dimension is 2·n_O·n_C. | Yes | [Kienker_1989_Equivalence_Aggregated_Markov.md](./Kienker_1989_Equivalence_Aggregated_Markov.md) |
| Bruno, Yang & Pearson (2005) | Canonical MIR form; dwell-time data identify exactly 2R(n_O+n_C−R) rate constants (Fredkin bound sharpened). | Yes | [Bruno_Yang_Pearson_2005.md](./Bruno_Yang_Pearson_2005.md) |
| Larget (1998) | Unique minimal canonical form per equivalence class of observationally identical aggregated Markov processes. | Yes | [Larget_1998_Canonical_Aggregated_Markov.md](./Larget_1998_Canonical_Aggregated_Markov.md) |
| Flomenbom & Silbey (2006) | A two-state open/shut trajectory fixes only an equivalence class of schemes, never the rate matrix. | Yes | [flomenbom_silbey_2006_two_state_trajectories.md](./flomenbom_silbey_2006_two_state_trajectories.md) |
| Milescu, Akk & Sachs (2005) | Macroscopic MLE from mean+variance under driving stimuli; nonstationary protocols break the equilibrium C-C-O/C-O-C degeneracy. | Yes | [Milescu_Akk_Sachs_2005_MLE_Macroscopic_Currents.md](./Milescu_Akk_Sachs_2005_MLE_Macroscopic_Currents.md) |
| Celentano & Hawkes (2004) | Macroscopic covariance fitting identifies state number and desensitised-state topology, but not linear-vs-branched arrangement. | Yes | [CelentanoHawkes_2004_Covariance_Matrix_Kinetic_Fitting.md](./CelentanoHawkes_2004_Covariance_Matrix_Kinetic_Fitting.md) |