# Larget (1998): A canonical representation for aggregated Markov processes

**Citation.** Bret Larget (1998). "A canonical representation for aggregated Markov processes." *Journal of Applied Probability* 35(2): 313–324. DOI: 10.1239/jap/1032192850 (Cambridge Core identifier S0021900200014972).

**Verification status.** Verified. Citation and the complete abstract were read from the Cambridge Core article page (two independent Cambridge URLs, byte-identical abstract). Core claims quoted from that abstract. Relations to Kienker (1989) and Fredkin-Montal-Rice / Fredkin-Rice are corroborated verbatim from Bruno, Yang & Pearson (2005, PNAS, PMC1088360) and the 2025 review arXiv:2510.26728. The full Larget paper body is paywalled and was not read, so theorem numbers and the exact algebraic construction are not sourced from the primary text.

## Core result

An *aggregated Markov process* is a deterministic function of a Markov process: you watch a coarse label (for ion channels, "open" vs "closed", a conductance class) rather than the underlying microscopic state. Many distinct generators (rate matrices) can produce the same law for that observed label sequence. Larget's paper does two things.

1. **Equivalence criterion.** It gives *necessary and sufficient* conditions for two continuous-time aggregated Markov processes to be equivalent, meaning they induce the same law for the observed (aggregated) process. This closes a gap left by Kienker (1989), who exhibited transformations that preserve the observable law but did not give a complete characterization of the equivalence relation.

2. **Canonical representation.** For both discrete and continuous time, any aggregated Markov process satisfying mild regularity conditions can be converted directly to a canonical representation that is (i) *unique* for each equivalence class and (ii) a *minimal parameterization of all that can be identified* about the underlying Markov process. Hidden Markov models on finite state spaces are framed as aggregated Markov processes by expanding the state space, so they inherit canonical representations too.

The practical content: the observable data pin down a point (the canonical form), and everything mapping to that same canonical form (the fiber, or equivalence class) is observationally indistinguishable. The canonical form is the coordinate system in which "what the data can know" is exactly the free parameters and nothing else is.

Characterization of the specific canonical form (from Bruno-Yang-Pearson 2005, who use it): Larget's canonical form has the same number of parameters as their "BKU" canonical form and, like it, minimizes the number of within-conductance-class (intraconductance) reactions; Larget's version uses real rate constants (it permits negative rates), so it never requires the complex-valued rates that competing canonical forms need when microscopic reversibility / detailed balance is violated.

## What is / is not identifiable

Larget's result is the general, model-class statement of a fact established earlier in the ion-channel literature:

- **Identifiable.** Only the equivalence-class label, encoded by the canonical form's free parameters. In the single-channel dwell-time framing, Fredkin-Montal-Rice / Fredkin-Rice showed the observable law of a model with nₒ open states and n_c closed states (nₒ, n_c = counts of aggregated states in each conductance class) contains at most 2·nₒ·n_c independent identifiable parameters. Larget's canonical parameterization is the general-purpose realization of "keep exactly the identifiable coordinates."
- **Not identifiable.** Anything varying within an equivalence class. Distinct rate matrices, and even distinct connectivity topologies, can give identical stationary observable statistics (Kienker's C-C-O vs C-O-C example). Which topologies below the parameter bound are identifiable is, per Bruno-Yang-Pearson, an open problem; models above the bound are provably non-identifiable.

So the canonical form does not rescue the individual microscopic rate constants or the wiring diagram: it certifies that those are unknowable from the observable and names what remains.

## Relevance to the theory-model-observable picture (A mechanistic vs B aggregation)

This paper is squarely and only about source **(B), aggregation** (rate matrix → observable). It operates at the level of a generator Q together with the aggregation map (the partition of states into observed classes), and asks which generators are indistinguishable through that map. Its equivalence class is exactly the object the program's caveat (B) needs: the set of rate matrices producing identical observable statistics. The canonical representation is the formal answer to "data return an equivalence class of models, not one" for the aggregation half of the non-identifiability, and it makes that class computable (map any fitted Q to its canonical form; two fits agree iff their canonical forms agree).

It says nothing about source **(A), the mechanistic map** (molecular theory → rate matrix). A canonical Q still carries no information about which physical state each index denotes, how many theories realize this Q, or whether a theory is underspecified. (A) compounds on top of Larget, it is not addressed by it.

Observation model, and a caveat on macroscopic data. Larget, Kienker, and Fredkin-Montal-Rice all analyze the *single-channel* observable (the law of the open/closed label process, equivalently the joint dwell-time distributions). The program instead observes many-channel, interval-averaged current. Larget does not treat this case. A logical consequence worth stating (my inference, not Larget's claim): because the many-channel current is a sum of independent and identically distributed single-channel copies, its stationary/equilibrium statistics are functionals of the single-channel observable law. So two generators in the same Larget equivalence class produce identical stationary macroscopic current statistics as well. Stationary macroscopic observation therefore does not, by itself, break an aggregation-equivalence class; it can only add power through non-stationary / relaxation structure or extra modalities. This is consistent with the authors' P2X2 case (Moffatt & Pierdominici-Sottile, Comm Biol 2025), where a coupling ambiguity stayed macroscopically unidentifiable (bimodal posterior) and had to be resolved by importing molecular-dynamics information. Whether that specific ambiguity is a Larget-type (B) aggregation degeneracy or an (A)-type mechanistic labeling symmetry is not settled by Larget and should not be asserted either way.

## Verbatim (sourced quotes)

Larget's abstract (Cambridge Core, verbatim):

> "A deterministic function of a Markov process is called an aggregated Markov process. We give necessary and sufficient conditions for the equivalence of continuous-time aggregated Markov processes. For both discrete- and continuous-time, we show that any aggregated Markov process which satisfies mild regularity conditions can be directly converted to a canonical representation which is unique for each class of equivalent models, and furthermore, is a minimal parameterization of all that can be identified about the underlying Markov process. Hidden Markov models on finite state spaces may be framed as aggregated Markov processes by expanding the state space and thus also have canonical representations."

On Kienker's equivalence (Bruno, Yang & Pearson 2005, PNAS, PMC1088360, verbatim):

> "Kienker showed that the two topologies C-C-O and C-O-C are 'equivalent': for every set of rates for one there is a set of rates for the other such that both have identical steady-state behavior."

On Larget's canonical form (same source, verbatim):

> "A generally applicable canonical form given by Larget is similar to BKU form. Both have the same number of parameters, and both minimize the number of intraconductance reactions."

> "Larget's form always has negative rates, whereas the complex rates in BKU form can be removed by a slight change in the form that adds no more parameters."

> "either Larget's canonical form should be fit (real rates suffice) or MIR or BKU form can be fit with complex eigenvalues."

On identifiability and its bound (same source, verbatim):

> "A topology is identifiable provided any changes to its rate constants result in changes to the steady-state dwell-time distributions."

> "Topologies with more than this many rates are not identifiable; topologies with fewer rates may or may not be identifiable. The problem of determining which topologies are identifiable remains unsolved."

Fredkin-Montal-Rice bound and Kienker transformations (review, arXiv:2510.26728, verbatim):

> "models with nₒ open and nc closed states that exceed the maximum number 2nₒnc of parameters are non-identifiable."

> "it is possible to reparametrise the model Q so that all components within Qcc and Qₒₒ except for the diagonal vanish—without changing the dynamics D!"

## Caveats

- Larget's full paper body is paywalled and was not read. Theorem/lemma labels and the exact algebraic construction of the canonical form are not sourced from the primary text; only the abstract is primary-verified.
- The Fredkin-Montal-Rice (1985) reference (Berkeley Neyman-Kiefer conference, vol. I, pp. 269–290, Wadsworth) and its exact pages come from web synthesis, not a page view. Note the companion Fredkin & Rice (1986), *J. Appl. Probab.* 23:208–214, is the peer-reviewed statement of the same identifiability line; Larget cites the Fredkin-Rice work and Kienker (1989) as the direct antecedents.
- The "BKU" and "MIR" (manifest interconductance rank) canonical forms are Bruno-Yang-Pearson's terminology; I did not verify the expansion of "BKU," so it is referenced only as a named comparison form.
- The macroscopic-observation paragraph is a logical inference about i.i.d. channel sums, clearly separated above from what Larget proves. Larget's scope is the single-channel aggregated process law only.

## Sources (URLs)

- Cambridge Core (Larget 1998, abstract + citation): https://www.cambridge.org/core/product/identifier/S0021900200014972/type/journal_article and https://www.cambridge.org/core/journals/journal-of-applied-probability/article/abs/canonical-representation-for-aggregated-markov-processes/D21B644D94227A9591ADE4D32708960E
- DOI: https://doi.org/10.1239/jap/1032192850
- Bruno, Yang & Pearson (2005), PNAS 102(18):6326–6331, full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC1088360/
- Review (2025), "Modelling ion channels with a view towards identifiability," arXiv:2510.26728: https://arxiv.org/html/2510.26728
- Kienker (1989), Proc. R. Soc. Lond. B 236(1284):269–309: https://doi.org/10.1098/rspb.1989.0024