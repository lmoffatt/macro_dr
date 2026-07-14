# Flomenbom & Silbey 2006: Utilizing the information content in two-state trajectories

**Citation.** Flomenbom, O. & Silbey, R. J. (2006). Utilizing the information content in two-state trajectories. *Proceedings of the National Academy of Sciences USA* 103(29):10907-10910. DOI: 10.1073/pnas.0604546103. arXiv preprint q-bio/0703013.

**Verification status.** Verified. The full arXiv preprint (identical in title, abstract, and results to the PNAS article) was read page by page; the PNAS/PMC pages confirmed the bibliographic details and abstract. Related works were confirmed by web search and cross-checked against the Siekmann (2025) review, which was also read directly. See the caveats for the two citations whose exact pagination was not independently web-verified this session.

Abbreviations used here: KS = kinetic scheme (the underlying multi-substate Markov model, equivalently the rate matrix / infinitesimal generator Q); WT-PDF = waiting-time probability density function (the distribution of a single dwell duration); RD = reduced dimensions (the canonical form this paper introduces); MIR = manifest interconductance rank (a canonical form from Bruno, Yang & Pearson 2005); BKU = Bauer-Kenyon uncoupled (a canonical form from Bauer, Bowman & Kenyon 1987); n_O, n_C = number of open and closed substates.

## Core result

The observable is a **two-state trajectory**: a single molecule (one enzyme, one ion channel) produces an alternating sequence of "on/open" and "off/closed" dwell periods. The molecule actually lives in a multi-substate Markov scheme (many open substates, many closed substates, with some connectivity), but the recording only reports which of the two aggregate classes the molecule is in at each instant. The transition times between the two classes are the data.

Flomenbom and Silbey ask what can be recovered about the KS from this trajectory in the **ideal limit** (their words: "an ideal (noiseless, infinitely long) two-state trajectory"), and they give a constructive, fitting-free answer plus an upper bound on the recoverable information.

Their machinery: the complete information in the ideal trajectory is carried by the four **joint distributions of two successive dwell times**, phi_{x,y}(t1,t2) for x,y in {on,off} (an on dwell of length t1 followed by an off dwell of length t2, and the other three combinations). Higher-order successive dwell-time distributions add nothing beyond these pairwise ones. Each phi_{x,y}(t1,t2), viewed as a kernel/matrix, has a **rank** R_{x,y} equal to its number of nonzero singular values, and this rank is obtainable directly from the data without fitting any functional form. The four ranks fix the topology of a canonical **RD form**: an on-off network with connections only between substates of different classes, where each connection carries a generally non-exponential WT-PDF. The RD form has the minimal number of substates needed to reproduce the data, and all KSs that map to the same RD form (same ranks, same connection WT-PDFs) are indistinguishable from the trajectory. The RD form is thus the canonical representative of an equivalence class, and it can be read off the data rather than searched for.

This paper is the ion-channel/enzyme information-theory statement of a classical fact: an **aggregated Markov model is identifiable only up to an equivalence class**. The companion and antecedent results make the bound quantitative (see Related work).

## What is / is not identifiable

**Identifiable from an ideal two-state trajectory:**
- The four ranks R_{on,on}, R_{on,off}, R_{off,on}, R_{off,off} (fitting-free), which give the RD canonical topology, including how many substates each class has in the RD form and the on-off connectivity.
- The number of exponential components in each marginal dwell-time distribution phi_on(t), phi_off(t), which bounds the number of open and closed substates in the simplest underlying KS.
- The **eigenvalues of the sub-generators** appear as the decay rates of the dwell-time distributions (the exponents of phi_on, phi_off are the negatives of eigenvalues of the open-block and closed-block submatrices of Q; stated explicitly in the Siekmann 2025 review and in Fredkin-Rice).
- The connection WT-PDFs of the RD form (generally non-exponential), which can be direction-dependent and can encode broken detailed balance at the on-off level.
- Whether two schemes have **different** ranks: any such pair is resolvable from the trajectory.

**Not identifiable:**
- The actual physical KS (the labeled rate matrix). Only its equivalence class is recovered. Any two KSs with the same RD form "cannot be distinguished," and this is not a limitation of the estimator: it holds for a noiseless, infinitely long trajectory.
- Substates that are neither entered nor exited on class transitions (not "initial" or "final") do not affect the RD topology at all: they are invisible to the trajectory.
- Which substate "means what" molecularly. The trajectory carries topology and rates of the canonical form, not the assignment of biophysical identity to states.
- The quantitative ceiling (from Fredkin-Montal-Rice 1985 / Fredkin-Rice 1986, restated in Siekmann 2025): a scheme with n_O open and n_C closed substates has at most **2·n_O·n_C** identifiable rate constants from the dwell-time data; any model with more parameters than this is non-identifiable and can be continuously reparametrised without changing the observable dynamics. For two open and two closed states that ceiling is 8. Even below the ceiling, distinct graph topologies can generate identical dynamics ("non-identifiability of model structure"), so what remains is a finite equivalence class of mechanisms, not a unique one.

A worked instance of the "two open states, which theories fit" question: if the data indicate two open states (of equal conductance, hence aggregated into "on"), the trajectory can pin the **count** (two exponentials in phi_on, or rank structure R_{on,off} = 2) but not the **mechanism**. Distinct arrangements of the two open states and their connections to the closed manifold that share the same RD/MIR form all fit equally well. The recording returns the equivalence class; choosing a member requires information from outside the trajectory.

## Relevance to the theory-model-observable picture (A mechanistic vs B aggregation)

This work is primarily a sharp statement about source **(B), aggregation** (rate matrix -> observable). The map KS -> two-state trajectory is many-to-one, and the fibre over any observable is exactly the RD equivalence class. Kienker 1989, Bruno-Yang-Pearson 2005, and Fredkin-Montal-Rice 1985 are the ion-channel backbone of the same (B) result (for example the C-C-O and C-O-C three-state topologies are equivalent); Flomenbom-Silbey generalise it to arbitrary schemes with irreversible connections and symmetry, and turn it into a data-driven, fitting-free procedure with an explicit information upper bound.

It also speaks to source **(A), mechanistic** non-identifiability (theory <-> rate matrix), by making the loss compound. Even granting a fixed observable-to-model correspondence, the model is recovered only up to reparametrisation; and reparametrisation can cross between mechanisms. Siekmann (2025) states this bluntly for aggregated Markov models: the dynamics of one model "can be generated by a different model that suggests a completely different mechanism." So the two sources stack. Theory -> KS is many-to-one (a rate matrix does not encode which state means what, and an underspecified theory maps to many matrices), and KS -> observable is many-to-one (aggregation). Data therefore returns an equivalence class of theories, which is exactly the caveat this reference is meant to ground.

**How macroscopic (this program's) observation changes the picture.** Flomenbom-Silbey concern the **single-molecule** trajectory, where individual dwells and, crucially, the correlations between successive dwells are observed. The companion result (Flomenbom-Klafter-Szabo 2005) shows that those **successive-dwell correlations are the only route to substructure beyond the marginal dwell-time distributions**: when successive waiting times are uncorrelated, phi_on(t) and phi_off(t) alone contain everything, and whole families of schemes collapse to indistinguishable.

This program observes the **macroscopic** current: many channels, summed and time-averaged over each sampling interval (the MacroIR observation model). That is a different, and in its stochastic content generally coarser, projection of the same aggregated Markov process:
- The many-channel sum is, by the central-limit effect over channels, close to Gaussian, so the macroscopic likelihood is effectively determined by the first two moments: the mean occupancy trajectory and the current autocovariance. Individual dwells, and therefore the higher-order successive-dwell correlations that carry the single-molecule topology information, are not observed.
- What macroscopic data does deliver well: the **relaxation eigenvalues of the full generator Q** (as the exponents of the mean response to a designed step or ramp, and as the exponents of the noise autocovariance), the conductance-weighted amplitudes of those modes, the equilibrium occupancies, and, from the variance amplitude, the single-channel current and channel number.
- Compensating strengths that matter in practice: designed **non-stationary protocols** (voltage, agonist, or concentration jumps) load many eigen-directions of Q and, by varying conditions, can break equivalences that a single stationary record cannot (Flomenbom-Silbey note explicitly that "varying some parameters, e.g. the substrate concentration" can distinguish schemes or resolve within an equivalence class); high signal-to-noise from ensemble averaging; and the MacroIR integrated-measurement likelihood, which models the within-interval time-average exactly and so extracts the available second-order content faithfully.

Net assessment: macroscopic observation does not escape either non-identifiability source. The aggregation equivalence classes (B) persist, and by discarding higher-order dwell-time correlations the macroscopic observable is, for resolving internal topology at fixed protocol, typically **more** degenerate than single-molecule dwell-time analysis, not less. Its advantage is protocol design and averaging, not richer per-record information. This is precisely why the authors' own P2X2 study (Moffatt & Pierdominici-Sottile 2025, *Communications Biology*; repo PDF `docs/bibliography/Moffatt_PierdominiciSottile_2025_Functional_Asymmetry_P2X2_CommBiol.pdf`) hit a left/right coupling ambiguity that was macroscopically unidentifiable (a bimodal posterior) and had to be resolved by importing an independent modality (molecular dynamics). Both Flomenbom-Silbey and Siekmann recommend exactly this move: resolve within the equivalence class using additional information "e.g. the crystal structure of the biopolymer" or other direct probes of the conformational dynamics.

## Verbatim (sourced quotes)

From Flomenbom & Silbey 2006, arXiv preprint q-bio/0703013 (read directly; text matches the PNAS article):
- Abstract: "The determination of the underlying KS is difficult and sometimes even impossible due to the loss of information in the mapping of the mutli-dimensional KS onto two dimensions." (the spelling "mutli-dimensional" is the arXiv preprint's; PNAS prints "multidimensional")
- Abstract: "Based on our approach, the upper bound on the information content in two-state trajectories is determined."
- Introduction: "Higher order successive WT-PDFs do not contain additional information on top of phi_{x,y}(t1,t2)."
- Introduction: "Moreover, there are KSs with the same phi_x(t)s and phi_{x,y}(t1,t2)s."
- Results (rank): "phi_{x,y}(t1,t2) is a matrix, whose rank R_{x,y} (=1,2,...), which is the number of non-zero eigenvalues (or singular values for a non square matrix) of it's decomposition, can be obtained without the need of finding the actual functional form."
- RD form: "RD form has the minimal number of substates needed to reproduce the data."
- Examples: "KSs with R_{x,y}=1 (x,y=on,off) and the same phi_on(t) and phi_off(t) are indistinguishable (assuming no additional information on the mechanism is known)."
- Examples: "Additional information can be inferred, under some physical assumptions, by analyzing different kind of measurements, e.g. the crystal structure of the biopolymer, or by analyzing two-state trajectories while varying some parameters, e.g. the substrate concentration."
- Concluding remarks: "The main effort in this paper is to utilize the information content in an ideal (noiseless, infinitely long) two-state trajectory for an efficient elucidation of a unique mechanism that can generate it."

From Siekmann 2025 review, "Modelling ion channels with a view towards identifiability," arXiv:2510.26728 (read directly; secondary source stating the general aggregated-Markov result and its lineage):
- Abstract: "models with n_O open and n_C closed states that exceed the maximum number 2 n_O n_C of parameters are non-identifiable."
- Section 3: "even if an infinite amount of data which is not perturbed by noise were available, the parameters of the model nevertheless could not be uniquely determined by fitting the model to these data."
- Section 3: "This makes it challenging to interpret aggregated Markov models as representations of a particular biophysical mechanism, simply because via reparametrisation, the dynamics of a given model can be generated by a different model that suggests a completely different mechanism."
- Section 3.1 (attributing Fredkin-Montal-Rice 1985, Fredkin-Rice 1986): "the dynamics of aggregated Markov models can be completely represented by the bivariate distributions f_OC and f_CO of the length of an open time t_O followed by a closed time t_C and vice versa."
- Section 3.1: "the parameters lambda_O and lambda_C are the negatives of eigenvalues of submatrices of the infinitesimal generator Q of the aggregated Markov model."

## Related work (the two-state-trajectory / aggregated-Markov identifiability lineage)

- **Fredkin, D. R., Montal, M. & Rice, J. A. (1985)** and **Fredkin, D. R. & Rice, J. A. (1986)**, *On aggregated Markov processes*, J. Appl. Probab. 23:208-214. Origin of the result that the dynamics are fully captured by the bivariate successive-dwell distributions, and of the 2·n_O·n_C identifiable-parameter ceiling. (Exact pagination of the 1985 conference-volume chapter not independently web-verified this session; the result and attribution are confirmed via Siekmann 2025, read directly.)
- **Bauer, R. J., Bowman, B. F. & Kenyon, J. L. (1987)**, Theory of the kinetic analysis of patch-clamp data, *Biophys. J.* 52(6):961-978, DOI 10.1016/S0006-3495(87)83289-7. Theory of dwell-time analysis for single-channel data; source of the "Bauer-Kenyon uncoupled (BKU)" canonical form that Flomenbom-Silbey generalise. Also the entry point to the missed-events / time-interval-omission problem that degrades real (finite-bandwidth) data below the ideal bound.
- **Kienker, P. (1989)**, Equivalence of aggregated Markov models of ion-channel gating, *Proc. R. Soc. Lond. B* 236(1284):269-309. The ion-channel equivalence theorem: distinct topologies (for example C-C-O and C-O-C) are equivalent, related by "Kienker transformations." (Exact page range from a secondary review, not the publisher page, this session.)
- **Bruno, W. J., Yang, J. & Pearson, J. E. (2005)**, Using independent open-to-closed transitions to simplify aggregated Markov models of ion channel gating kinetics, *PNAS* 102(18):6326-6331, DOI 10.1073/pnas.0409110102. Introduces the MIR (manifest interconductance rank) canonical form with the minimal number of rate constants for a given rank; the other canonical form Flomenbom-Silbey build on.
- **Flomenbom, O., Klafter, J. & Szabo, A. (2005)**, What can one learn from two-state single-molecule trajectories?, *Biophys. J.* 88(6):3780-3783. The companion result: when successive waiting times are uncorrelated, the marginal on and off WT-PDFs contain all the information, and the schemes that produce uncorrelated trajectories cannot be told apart "by any sophisticated analyses"; correlations between successive dwells are the sole extra information channel.
- **Siekmann, I. (2025)**, Modelling ion channels with a view towards identifiability, arXiv:2510.26728. Recent review that enumerates aggregated Markov models by Polya counting, restates the 2·n_O·n_C bound with two derivations, distinguishes parameter non-identifiability from model-structure non-identifiability, and argues (as this program does) that recovering biophysical mechanism needs data beyond steady-state single-channel records.

## Caveats

- Flomenbom-Silbey's bound is for the **idealised** single-molecule trajectory (noiseless, infinitely long). Real recordings lose more to finite length, noise, filtering, and missed brief events (time-interval omission), so identifiability in practice is strictly worse than this ceiling.
- The information bound is stated in terms of **ranks and connection WT-PDFs of the RD form**, which is generally non-Markovian. It is a statement about the equivalence class, not a recipe for a unique physical rate matrix.
- The macroscopic-vs-single-channel comparison in "Relevance" is my synthesis, grounded in (i) the sufficient-statistic results of Fredkin/Flomenbom (single-channel information lives in the bivariate dwell distributions), (ii) the second-order/Gaussian character of many-channel current and this program's mean-plus-covariance observation model, and (iii) the papers' own recommendation to import additional modalities. It is not a verbatim theorem from Flomenbom-Silbey, who do not treat many-channel or interval-averaged observation.
- Kienker (1989) exact pages and Fredkin-Montal-Rice (1985) exact pagination were taken from secondary sources this session, not the publisher pages; the substantive results were verified.
- The P2X2 example is context from the program's own work (Moffatt & Pierdominici-Sottile 2025); the specific bimodal-posterior / molecular-dynamics-resolution detail is stated from that context, not re-derived here.

## Sources (URLs)

- https://www.pnas.org/doi/10.1073/pnas.0604546103 (PNAS article page)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC1544147/ (PMC full text)
- https://arxiv.org/abs/q-bio/0703013 and https://arxiv.org/pdf/q-bio/0703013 (arXiv preprint, read directly)
- https://arxiv.org/pdf/2510.26728 (Siekmann 2025 review, read directly)
- https://www.pnas.org/doi/10.1073/pnas.0409110102 and https://pmc.ncbi.nlm.nih.gov/articles/PMC1088360/ (Bruno-Yang-Pearson 2005, MIR)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC1330095/ (Bauer-Bowman-Kenyon 1987)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC1305612/ and https://arxiv.org/pdf/q-bio/0502006 (Flomenbom-Klafter-Szabo 2005 companion)
- https://royalsocietypublishing.org/doi/10.1098/rspb.1988.0063 (Royal Society B, adjacent single-channel two-state / missed-events literature)