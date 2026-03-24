## 1. **Collaborations & Planned Papers**

* **With Gustavo Pierdominici**:
  Extend their earlier paper by building kinetic models with explicit coupling of subunit rotations, incorporating molecular dynamics and geometric considerations. New “schemes” (16, 17) formalize ternary and rotational couplings.
* **With Cecilia Boussat**:
  Possible application of allosteric models to CIS-LU receptors, oocytes, and mutants. Aim: model states like closed, rotated, open, desensitized, with or without asymmetry.
* **Solo Work**:

  * *MacroDR characterization paper*: formal validation of all core algorithms (simulation, likelihood, evidence, sampling, score, Fisher Information Matrix).
  * *Cumulative evidence paper*: showing how multiple datasets combine to refine inference.
  * *Allosteric modeling DSL paper*: present a language for generating complex kinetic models.
  * Longer-term: channel idealization and likelihood theory for single-channel experiments.

---

## 2. **Algorithmic Core of MacroDR**

Luciano repeatedly defines the **six to seven core modules**:

1. **Simulation**
2. **Likelihood**
3. **Evidence**
4. **Sampling** (posterior via MCMC)
5. **Score (derivative of likelihood)**
6. **Fisher Information Matrix (FIM)**
7. (sometimes) Idealization

Each module must have:

* Clear **inputs/outputs** (model, priors, data, algorithm).
* **Tests**: correctness, precision, speed.
* **Conditions of validity**: which parameter ranges work and which fail.

---

## 3. **Testing & Validation Strategies**

He elaborates a systematic framework:

* **Likelihood validation**:

  * Expectation of the score = 0.
  * Covariance of the score = Fisher Information Matrix.
* **Evidence validation**:

  * Use *confusion matrices* from simulated vs. fitted models.
  * See if evidence correctly identifies the generating model.
* **Sampling validation**:

  * Duality between likelihood and sampling, treated as categorical adjunction.
  * Compare distributions reconstructed from samples with the originals.
  * Use tests based on multinomial approximations, kernel methods, or tessellations.
* **Two-stage testing**:

  * **Primary validation**: correctness within a prior domain.
  * **Follow-up regression testing**: check stability across code changes.
* **Continuous Integration (CI)**:

  * Integrated with GitHub to store test outputs, detect regressions, and accumulate information about valid/invalid parameter regions.

---

## 4. **Software Architecture & DSL Design**

* **Domain Specific Language (DSL)**:

  * Commands for simulation, likelihood, evidence, sampling, model construction, experiment definition, priors.
  * Commands must declare **preconditions** and **postconditions**.
  * Types: models, experiments, distributions, priors.
* **Environments**:

  * Save/restore program state (JSON or dataframes).
  * Environments represent “contexts” of knowledge (variables + values).
* **Implementation strategy**:

  * Maintain an *old version* for cluster runs (fast, reliable).
  * Develop a *new refactored version* with modular classes, HoTT-style semantics, and cleaner organization.
* **Facilities**: testing, saving environments, cluster execution, algorithm optimization, experiment modeling.

---

## 5. **Conceptual / Theoretical Explorations**

* **Category Theory**:

  * Likelihood ↔ Sampling as adjoint functors.
  * Data–parameter duality framed as a triangular relation: model maps parameters → data, inference inverts data → parameters.
* **Homotopy Type Theory (HoTT)**:

  * Types defined by invariants/postconditions.
  * Objects as environments satisfying those invariants.
  * Goal: “HoTTify” MacroIR, embedding semantics in DSL.
* **Philosophy of Priors & Meaning**:

  * Distinguishes “distribution” (mathematical object) vs. “prior” (distribution with semantic context).
  * Meaning = conditional information / environment context.
* **Topology of inference**:

  * Multimodality in parameter spaces requires topological thinking (submanifolds, symmetries).
  * Idea of reducing complex models into interacting sums of simpler state-pairs.

---

## 6. **Practical Issues & Problems Identified**

* **Model errors**:

  * Posterior of inactivation parameter behaved oddly (cutoffs, absence of inactivation).
  * Equilibrium constants not well recovered in some runs.
* **Simulation tradeoffs**:

  * “Safe” multinomial simulation vs. faster QΔt/Taylor approximations.
  * Concern about validity when approximating conductance distributions (normal vs. Poisson).
* **Performance**:

  * Need to balance speed vs. correctness.
  * Interest in implementing Levenberg–Marquardt, Quasi-Newton (Cuevi), or Fisher-based approximations for acceleration.

---

## 7. **Operational Plans (Chronological Flow)**

* **Late August (Aug 27)**: Define collaborations, three paper tracks, CLI commands, first plan to run schemes 16/17 on clusters.
* **Early Sept (Sep 2–4)**: Theorizing categorical adjunctions; adding new models (16, 17, 8); discarding harmonic oscillator approach; identifying equilibrium recovery problems.
* **Mid Sept (Sep 8)**: Refactoring MacroDR with DSL, modularization, facilities, and testing.
* **Sep 18**: MacroDR runs evidence again; define papers with Cecilia, Gustavo, and solo MacroDR presentation; deeper dive into likelihood/evidence validation.
* **Sep 20**: Refining DSL design (types, functions, variables, environments, priors); defining postconditions; simulation strategies compared.
* **Sep 24**: Integrating “load model” but concerned about speed; major step: GitHub CI testing framework designed; two-stage testing (validation + regression).
* **Sep 26**: Achieved CI integration; added JSON environment saving; new commands (Q0, eigenvalues, transition probabilities). Preparing to implement HoTTification formally; planning documentation and factoring strategy.

---

## 8. **Overall Trajectory**

Luciano’s project evolves along two intertwined lines:

1. **Scientific/Collaborative**: building and validating new kinetic models (subunit rotation, allosteric couplings) and preparing papers with Gustavo and Cecilia.
2. **Methodological/Technical**: refactoring MacroDR into a DSL-driven, test-heavy, HoTT-informed framework with continuous integration, capable of rigorous validation of simulation, likelihood, sampling, and evidence.

The **recurring tensions** are:

* Publishing quickly vs. indulging in deep software/theory redesign.
* Speed of computation vs. correctness of models.
* Practical collaboration goals vs. philosophical/categorical explorations.


