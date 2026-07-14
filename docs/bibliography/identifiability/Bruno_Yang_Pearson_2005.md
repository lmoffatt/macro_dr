# Bruno, Yang & Pearson (2005): MIR canonical form for aggregated Markov ion-channel models

**Citation.** Bruno WJ, Yang J, Pearson JE (2005). "Using independent open-to-closed transitions to simplify aggregated Markov models of ion channel gating kinetics." *Proceedings of the National Academy of Sciences USA* 102(18):6326-6331. DOI: 10.1073/pnas.0409110102. PMID 15843461, PMCID PMC1088360.

**Verification status.** Verified. Citation and abstract confirmed against PubMed (PMID 15843461); the parameter-count formula, canonical-form definitions, and the sufficiency/dwell-time claims were each confirmed by verbatim sentence extraction from the PubMed Central full text (PMC1088360). The PNAS publisher page returned HTTP 403, so no quote is drawn from it. Confirmed that "macroscopic" and "whole-cell" do not appear in the paper; the analysis is single-channel steady-state dwell-time data throughout.

## Core result

The paper gives a canonical form for aggregated (lumped open/closed) Markov models of a binary-conductance channel and uses it to characterize which models are indistinguishable from single-channel steady-state data.

Setup: `No` = number of open states, `Nc` = number of closed states, `Qoc` = the block of the rate matrix carrying open-to-closed transitions, and the *interconductance rank* `R` = the matrix rank of `Qoc` (the number of independent open-closed transition routes). Fredkin and Rice's earlier bound is that `R` is at most the lesser of the number of open states connected to closed states and the number of closed states connected to open states, hence `R <= min(No, Nc)`.

The authors construct the *manifest interconductance rank* (MIR) form, defined by (i) all interconductance links being independent and (ii) every intraconductance link touching at least one state that has an interconductance link. This is obtained mathematically by "diagonalizing" the `Qoc` / `Qco` transition blocks (a Kienker-type similarity transformation). Two headline properties:

- The MIR form has the minimal number of rate constants for any given rank, strictly fewer than the older Bauer-Kienker uncoupled (BKU) canonical form for reduced-rank models. For four open and four closed states there are exactly four MIR schemes, one per rank R = 1, 2, 3, 4.
- Fitting a model to canonical form preserves detailed balance: if a model satisfies detailed balance, its BKU and MIR equivalents do too.

They also prove a strong collapse of model space at low rank: all rank-1 topologies with a given `No`, `Nc` are equivalent (make identical steady-state predictions), and they give a hierarchical search algorithm that fits progressively higher ranks to find the simplest adequate model.

## What is / is not identifiable

**The precise identifiable-parameter count (from the paper, single-channel steady-state data).** The MIR form has `R(No + Nc - R)` reactions, i.e.

> the Fredkin bound of `2R(No + Nc - R)` independent rate constants that can be recovered from steady-state data.

(The factor 2 is forward plus reverse rate per reaction.) This is the dimension of the identifiable parameter set for a rank-`R` binary channel observed as single-channel dwell times in steady state. For comparison the BKU uncoupled form uses `2·No·Nc` rate constants; the two counts coincide only at full rank `R = min(No, Nc)` and the MIR count is dramatically smaller at low rank (e.g. No=Nc=4: rank 1 gives 14 rate constants versus BKU's 32).

**What is not identifiable.** Everything that a Kienker-type transformation can move while holding the observable dwell-time statistics fixed. Concretely: the specific topology within a rank class is not recoverable (all rank-1 topologies are mutually equivalent), and any two models sharing the same MIR canonical representative fit the data equally well. The sufficient statistics are the open and closed dwell-time densities `foc` and `fco` (Fredkin, Colquhoun-Hawkes lineage); anything not encoded in those densities is invisible in steady state. The full rate matrix generally carries more free parameters than `2R(No + Nc - R)`, and that excess is exactly the non-identifiable manifold of kinetically equivalent models. The paper also notes that many models equivalent to a detailed-balance model themselves violate detailed balance, so imposing detailed balance is a genuine (identifiability-relevant) restriction on which canonical fit is admissible.

## Relevance to the theory-model-observable picture (A mechanistic vs B aggregation)

This paper is squarely and almost entirely about source **(B), aggregation non-identifiability**: the map from rate matrix to observable, where the observable aggregates states into open/closed conductance classes. It is one of the cleanest statements available of the classic aggregated-Markov equivalence problem. Its contributions to the program's caveat:

- It provides a *canonical representative* (MIR) of each equivalence class, so "the data return an equivalence class, not one model" becomes constructive: the class is the set of rate matrices mapping to one MIR scheme, and the identifiable content is the `2R(No + Nc - R)` MIR rate constants.
- It quantifies the size of the collapse. The >2 million topologies for four-open/four-closed states reduce, in steady state, to four rank classes, with rank-1 topologies fully interchangeable. This is a concrete measure of how many-to-one the rate-matrix-to-observable map is.
- It does not address source **(A), mechanistic non-identifiability** (theory to rate matrix). State labels carry no molecular meaning here; the paper works purely at the rate-matrix-to-observable layer. But (B)'s topology interchangeability compounds directly with (A): once a topology is unidentifiable, any molecular story that would have been distinguished by that topology is also unresolved.

**How macroscopic (vs single-channel) observation changes the picture.** This is the key caveat for importing the count. The paper's identifiability result is derived for single-channel *dwell-time* data in steady state, with `foc`, `fco` as the sufficient statistics. The words "macroscopic" and "whole-cell" do not appear. The program observes many-channel, interval-averaged currents (ensemble mean and, via Moffatt 2007, fluctuation/autocovariance), which is a different reduction of the same rate matrix: the macroscopic relaxation and its variance depend on the spectrum of `Q` projected onto conductance, not on the individual dwell-time densities. So the `2R(No + Nc - R)` figure is the single-channel dwell-time identifiable dimension and should not be quoted as the macroscopic identifiable dimension without separate derivation. The equivalence classes can differ between the two observation modes: quantities fixed under a Kienker transformation for dwell times need not be exactly the macroscopically observable ones, and macroscopic observation of relaxations under multiple stimulus protocols can in principle break some single-channel degeneracies while being blind to others. The P2X2 left-right coupling ambiguity (Moffatt & Pierdominici-Sottile, Comm Biol 2025), which was macroscopically unidentifiable and resolved only by molecular dynamics, is an instance of a degeneracy that survives macroscopic observation, whether or not it coincides with a Bruno-Yang-Pearson dwell-time equivalence class.

## Verbatim (sourced quotes)

From the PubMed Central full text (PMC1088360):

- "The total number of reactions in MIR form is _R_(_No_ + _Nc_ - _R_), which gives the Fredkin bound of 2_R_(_No_ +_Nc_ - _R_) independent rate constants that can be recovered from steady-state data."
- "In BKU form, every open state is connected to every closed state and vice versa, but there are no open-to-open or closed-to-closed connections. This form has 2_NoNc_ rate constants."
- "Fredkin and Rice noted that the lesser of the number of Os connected to Cs and the number of Cs connected to Os is an upper bound on a topology's rank."
- "Fredkin _et al._ showed that _foc_ and _fco_ are sufficient to fully characterize the steady-state data, under certain technical conditions that are almost always satisfied."
- "We will compare the canonical form introduced below to the canonical form described by Bauer _et al._ and Kienker, which we refer to as Bauer-Kienker uncoupled (BKU) form."

From the abstract (PubMed, PMID 15843461):

- "We have found a canonical form that can express all reaction schemes for binary channels. This form has the minimal number of rate constants for any rank (number of independent open-closed transitions), unlike other canonical forms such as the well established 'uncoupled' scheme. Because all of the interconductance transitions in the new form are independent, we refer to it as the manifest interconductance rank (MIR) form."
- "In the case of four open and four closed states, there are four MIR form schemes, corresponding to ranks 1-4."
- "By using the MIR form we prove that all rank 1 topologies with a given number of open and closed states make identical predictions in steady state, thus narrowing the search space for simple models. Moreover, we prove that fitting to canonical form preserves detailed balance."

## Caveats

- The identifiable count `2R(No + Nc - R)` is for single-channel steady-state dwell-time data. It is not a statement about macroscopic/ensemble-current identifiability, which the paper does not treat. Do not transplant the number to the MacroIR observation model without a separate derivation.
- "Binary channel" means two conductance classes (open/closed). Multi-conductance (subconductance) channels are outside the stated scope.
- Equivalence here is exact steady-state distributional equivalence; it is silent on finite-data statistical identifiability (a degeneracy that is exact in the model may still be only weakly informed by finite noisy data, and conversely).
- The rank-1 "all topologies equivalent" theorem and the parameter count assume the technical regularity conditions Fredkin et al. require ("almost always satisfied").
- Section labels and any theorem numbering were not independently re-verified beyond the verbatim sentences above (PNAS page was 403; PMC HTML lacks stable theorem numbers). The formula and definitions themselves are confirmed verbatim.

## Sources (URLs)

- PubMed abstract and metadata: https://pubmed.ncbi.nlm.nih.gov/15843461/
- PubMed Central full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC1088360/
- Publisher landing (metadata only; body returned HTTP 403): https://www.pnas.org/doi/10.1073/pnas.0409110102