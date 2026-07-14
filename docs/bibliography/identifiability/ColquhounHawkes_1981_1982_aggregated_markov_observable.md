# Colquhoun & Hawkes (1981, 1982): the aggregated-Markov theory of the single-channel record

**Citation.** Colquhoun, D. & Hawkes, A. G. (1981). On the stochastic properties of single ion channels. *Proc. R. Soc. Lond. B* **211**(1183): 205-235. DOI 10.1098/rspb.1981.0003. — Colquhoun, D. & Hawkes, A. G. (1982). On the stochastic properties of bursts of single ion channel openings and of clusters of bursts. *Phil. Trans. R. Soc. Lond. B* **300**(1098): 1-59. DOI 10.1098/rstb.1982.0156 (received 19 Feb 1982; published 24 Dec 1982; communicated by Sir Bernard Katz, F.R.S.). Authors: D. Colquhoun, Dept of Pharmacology, University College London; A. G. Hawkes, Dept of Statistics, University College of Swansea.

**Verification status.** 1982 paper verified against its primary full text (onemol.org.uk scanned OCR PDF, pages read directly: title page + contents, abstract, body pp. 4-6, 13-15, 51-53). 1981 paper: citation and core result verified through the 1982 paper's verbatim citation of it plus three independent secondary sources; its own full text and verbatim abstract could not be accessed (publisher and ResearchGate returned 403, PubMed captcha-walled). **Task citation correction:** the prompt swapped the two papers' year/journal pairs. "Single ion channels" is the 1981 *Proc. R. Soc. B* paper; "bursts... and clusters of bursts" is the 1982 *Phil. Trans. R. Soc. B* paper.

## Core result

These two papers are the foundational forward theory of the *observable* aggregated Markov process for a single ion channel. They do not fit data; they solve the map: given a rate matrix, what are the exact statistics of the record you can actually see.

Setup (both papers, stated most generally in 1982 §1). The channel is a continuous-time Markov process on *k* kinetically distinguishable states with a *k*×*k* transition-rate matrix **Q** (the CTMC generator: off-diagonal q_ij ≥ 0 are transition rates, diagonal q_ii = −Σ_{j≠i} q_ij; transition probabilities are constant in time). Write P_ij(t) = P(state j at time t | state i at time 0); the matrix **P**(t) = exp(**Q**t). The observer does not see the state; the observer sees only the conductance class. States are partitioned into aggregates:

- **A** = open (conducting) states, k_A of them;
- **F** = shut (non-conducting) states, k_F of them (in the 1982 burst analysis F is further split into short-lived **B**, "gaps within a burst," and long-lived **C**, "gaps between bursts," with E = A ∪ B and F = B ∪ C).

**Q** is partitioned conformably into blocks. In open/shut form:

- **Q**_AA — transitions among open states,
- **Q**_AF — open → shut,
- **Q**_FA — shut → open,
- **Q**_FF — transitions among shut states.

(The task's "QAA/QAB" is this partition with the shut class written A/F: Q_AA within-open, Q_AF open-to-shut.)

The 1981 paper's contribution: a *general* method (any mechanism, any number of open/shut states, cyclic reactions allowed) for the distribution of the length of a sojourn in any specified subset of states. Its key structural finding is that a sojourn-time distribution depends not only on the state entered but also on the state the process moves to when it leaves the subset. The 1982 paper's own words (p. 5): Colquhoun & Hawkes (1981) "gave a general method, applicable to any specified mechanism, for deriving the distribution of the length of time spent in any specified subset of states."

The dwell-time densities are matrix-exponential in the sub-blocks. The survivor function of an open sojourn is exp(**Q**_AA t); the canonical Colquhoun-Hawkes open-dwell density is g_A(t) = **φ**_A exp(**Q**_AA t)(−**Q**_AA)**u**_A, with **φ**_A the entry-probability row vector and **u**_A a column of ones (the closed-class density is the same with A↔F). Scalarised, this is a **mixture of exponentials whose rate constants are the negatives of the eigenvalues of Q_AA** — the "spectral expansion" of the dwell-time density. The 1982 paper exhibits the single-state special case directly (open density f_1(t) = (−q_11) exp(q_11 t), Laplace transform (−q_11)/(s−q_11)) and states that evaluating the general results requires only "subroutines for finding eigenvalues and eigenvectors, and for matrix inversion." The 1982 paper generalises the 1981 machinery from single dwell times to every observable of a burst and a cluster of bursts (number of openings per burst, total open time per burst, gaps within/between bursts, burst length, cluster statistics), each computed by multiplying the block sub-matrices along the allowed route (the G_AB, G_BA jump matrices and the geometric resolvent Σ_r (G_AB G_BA)^r = (I − G_AB G_BA)^{-1}).

## What is / is not identifiable

This is exactly the payload the caller cares about, and both papers are unusually explicit about it (1982 §6, "Discussion").

Identifiable / lower-bounded from the aggregated record:

- **A count of components lower-bounds a count of states.** The number of exponential components in the distributions of number-of-openings-per-burst, of open times, and of total open time per burst all equal the number of open states k_A; so the number of open states "must be at least as large as the number of components observed." Symmetrically, if bursts are resolvable there must be ≥ 3 states; if clusters of bursts are resolvable, ≥ 4 states.
- **Which blocks a given dwell distribution constrains.** For the open-sojourn distribution, the rate constants (the λ_m) depend on Q_AA (within-subset rates) and Q_AF (exit rates); the component weights additionally depend on the entry rates Q_FA and on the equilibrium occupancies of the shut states.

Not identifiable / invisible:

- **The observer cannot resolve individual states within a class** — only whether the channel is open or shut (1982 p. 4, problem 1). Distinct open states of identical conductance are merged.
- **Q_FF is invisible to open sojourns.** Verbatim (1982 p. 52): "The only transition rates not involved in the distribution are those between states that are not part of the specified subset (e.g. those in Q_FF)." So the open-dwell distribution carries no information about transitions among shut states, and vice versa. This is a concrete, first-principles instance of aggregation non-identifiability.
- **State count is only a lower bound.** "Some components might remain experimentally unresolved so that k_A could always be greater than the number of observed components" (1982 p. 51). Missed brief events (formalised in Hawkes & Colquhoun 1983 and later work) make this worse.
- The two papers establish the *forward* map Q → observable statistics and enumerate what that map manifestly discards; they do **not** prove the full characterisation of which distinct Q matrices yield identical single-channel statistics (the aggregated-Markov equivalence-class theorem — Fredkin, Rice, Kienker, Bruno et al.). That later work builds directly on this setup.

## Relevance to the theory-model-observable picture (A mechanistic vs B aggregation)

These papers are squarely and only about source (B), aggregation (rate matrix → observable). They take **Q** as given and characterise the many-to-one map from Q to observable dwell-time / current statistics; they are silent on source (A), the theory ↔ rate-matrix correspondence (nothing in Q encodes which molecular story a state represents). They are the mathematical substrate on which (B) rests, and they hand you two concrete, provable losses to cite: (i) states within a conductance class are unresolvable, and (ii) the within-shut generator Q_FF leaves no fingerprint on open dwell times (and dually). The equivalence "distinct Q → identical observable statistics" is the natural closure of these two facts.

How macroscopic (many-channel, interval-averaged) observation changes the picture — the program's regime, not these papers' regime. Both papers analyse the *single-channel* record, whose dwell-time distributions are the richest aggregated observable available. Macroscopic current is a *further* aggregation: an ensemble average over many channels. The single-channel dwell-time structure is washed out, and what survives is the relaxation (ensemble-mean current = a sum of exponentials whose rates are the eigenvalues of Q — the same spectral object) plus the noise spectrum / autocovariance. The 1982 paper flags this link explicitly (p. 4): the existence of multiple shut/open states "such that conventional macroscopic measurements of the total current flow, through a large number of ion channels, would result in relaxations that were not simple exponentials, but that could be described by a sum of several exponential terms" (citing their 1977 macroscopic-relaxation paper). Consequence for the program: macroscopic observation is at best as informative, and generically strictly less informative, about the underlying Q than single-channel recording, so the aggregation non-identifiability (B) is at least as severe. This is consistent with the authors' own P2X2 result (Moffatt & Pierdominici-Sottile, Comm Biol 2025), where a left-right coupling ambiguity was macroscopically unidentifiable (bimodal posterior) and had to be broken by importing an independent modality (molecular dynamics): the macroscopic likelihood returned an equivalence class, exactly what the aggregation of Q → observable predicts.

## Verbatim (sourced quotes, from the 1982 paper's rendered OCR pages)

- Abstract (pp. 3-4): "Characteristics of observed bursts of single channel openings were derived recently for two particular ion channel mechanisms. In this paper these methods are generalized so that the observable characteristics of bursts can be calculated directly for any mechanism that has transition probabilities that are independent of time as long as the process is at equilibrium or is maintained in a steady state by an energy supply. General expressions are given for the distributions of the open time, the number of openings per burst, the total open time per burst, the gaps within and between bursts, and so on. With the aid of these general results a single computer program can be written that will provide numerical values for such distributions for any postulated mechanism, given only the transition rates between the various states. ... The analogous theory is also given for the case where bursts of channel openings are grouped into clusters; many of the results bear a close analogy with those found for simple bursts."
- The aggregation problem (p. 4): "Although most plausible reaction mechanisms postulate several non-conducting (shut) states, and sometimes also more than one open state (all open states possibly having identical conductance), the actual observations generally show only whether the channel is conducting (open) or not (shut)."
- Information asymmetry (p. 4): "in so far as there are more shut states than open states, most of the information about mechanisms should come from measurements of shut time durations rather than of open time durations."
- Markov assumption (p. 5, §1a): "the behaviour of single ion channels is analysed in terms of a Markov process in continuous time. Thus it is assumed that the probability of transition from one state to another is a constant, independent of time."
- State partition (p. 6): "We define as k the number of kinetically distinguishable states in which the system can exist. ... (1) Subset A comprises the open states (k_A in number, say)."; "we define F = B ∪ C, so that F contains all the states in B and C, i.e. all the shut states".
- Q-block partition (p. 15): "The analysis of bursts is based on the definition of the subsets of states, A, B and C ... and on the consequent partition of the matrix of transition rates in (1.6)."
- Spectral computation (p. 51): "the only complicated parts of such a program (subroutines for finding eigenvalues and eigenvectors, and for matrix inversion) are all readily available in standard libraries."
- Components → states (p. 51): "the number of components in the distributions of the number of openings per burst, of open times and of total open time per burst are all predicted to be equal to the number of open states (k_A). Therefore the number of open states must be at least as large as the number of components observed"; "Of course some components might remain experimentally unresolved so that k_A could always be greater than the number of observed components."
- Q_FF invisible (p. 52): "The only transition rates not involved in the distribution are those between states that are not part of the specified subset (e.g. those in Q_FF)."
- 1982's characterisation of the 1981 paper (p. 5): Colquhoun & Hawkes (1981) "gave a general method, applicable to any specified mechanism, for deriving the distribution of the length of time spent in any specified subset of states."

## Caveats

- Task-prompt citation error (see Verification status): the 1981/1982 and Proc./Phil.-Trans. attributions were swapped; corrected above.
- The 1981 paper was not read in full (publisher + ResearchGate 403, PubMed captcha). Its citation and core result are corroborated by the 1982 paper (which cites it verbatim) and three independent secondary sources, but I hold no directly-sourced verbatim quotation of the 1981 abstract; the verbatim block above is 1982-only.
- The compact matrix-exponential density g_A(t) = φ_A exp(Q_AA t)(−Q_AA)u_A is the standard Colquhoun-Hawkes form (1981 / restated in their 1995 "Q-matrix cookbook"); I did not read that exact formula on a page. What I read directly and can vouch for is the single-state special case and the statement that the general evaluation is an eigenvalue/eigenvector problem — which together establish the spectral-expansion claim.
- These papers assume equilibrium or a steady state maintained by an energy supply, time-homogeneous rates, and (for the clean dwell-time interpretation) perfect time resolution. Finite bandwidth (missed brief events) is treated in the companion Hawkes & Colquhoun (1983) and later asymptotic-distribution papers, not here.
- Scope: single-channel record. The program's macroscopic, interval-averaged observable is a further aggregation; these papers supply the spectral/relaxation link (via their 1977 paper) but not the many-channel likelihood itself.

## Sources (URLs)

- Full text (1982, primary, read directly): https://www.onemol.org.uk/Colquhoun%20&%20Hawkes-1982-ocr.pdf
- 1982 publisher record + DOI: https://royalsocietypublishing.org/doi/10.1098/rstb.1982.0156 (metadata only; body 403)
- 1981 publisher record + DOI: https://royalsocietypublishing.org/doi/10.1098/rspb.1981.0003 (403; citation corroborated via Google Scholar lookup https://scholar.google.com/scholar_lookup?doi=10.1098/rspb.1981.0003 and PubMed https://pubmed.ncbi.nlm.nih.gov/6111797/)
- 1982 PubMed record: https://www.ncbi.nlm.nih.gov/pubmed/6131450
- OneMol theory-papers index (Colquhoun & Hawkes PDFs): https://onemol.org.uk/?page_id=175
- Berkeley "Stochastic Models for Ion Channels: Introduction and Bibliography": https://digitalassets.lib.berkeley.edu/sdtr/ucb/text/327.txt