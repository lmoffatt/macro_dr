# Project Roadmap: MacroIR / MacroDR

## 1. Immediate Execution

* **Run legacy MacroIR code**:

  * Launch cluster jobs with scheme 16/17 (rotation-coupling) using NIRSK filters.
  * Also run existing schemes (6–9) for comparison.
  * Accept that some anomalies (e.g., posterior of inactivation cutoff) may only show up after long runs; still worthwhile to proceed.
* **Primary scientific question**: Do rotation-coupled schemes significantly change predictions? Answering this has priority, even if some posterior pathologies persist.

---

## 2. Verification & Characterization

* **Tests to establish correctness**:

  * Simulation, Likelihood, Sampling, Evidence.
  * Fisher Information Matrix vs. score covariance tests.
  * Evidence tests via confusion-matrix approach (distinguish true vs alternative models).
* **Performance metrics**:

  * Precision vs. speed plots for minimal models.
  * Define conditions where each algorithm (simulation, likelihood, sampling, evidence) works correctly and efficiently.
* **Secondary exploration**: Characterize MCMC performance (diffusion speed, isotropic “expansion” tests).

---

## 3. Development of MacroDR (refactor / next-gen)

Planned improvements:

1. **Command-line integration with built-in tests**: each function (simulation, likelihood, etc.) declares its test routine.
2. **Environment persistence**: save/restore state as JSON, allowing interruption/resumption of scripts (especially evidence and sampling loops).
3. **Interruptible long functions**: resume sampling/evidence jobs mid-run.
4. **Abstract class hierarchy**: models and algorithms abstracted to reduce template explosion, cleaner type system.
5. **Cluster/cloud job generation**: MacroDR or companion R scripts can auto-generate SLURM/cloud job scripts.
6. **Transparent development process**: consistent documentation of issues, evolution steps, so MacroDR can be openly published.

---

## 4. Research Projects & Papers

* **Paper A (with Gustavo)**: rotational kinetic models (continuation of published work).
* **Paper B (characterization of MacroDR)**: theoretical framework + empirical evaluation of Simulation, Likelihood, Sampling, Evidence, Score, FIM.
* **Paper C (cumulative evidence)**: successor to characterization.
* **Paper D (oocytes + mutants)**:

  * Apply conformationalModel + Evidence to macroscopic oocyte currents (slow, asymmetric regimes).
  * Integrate mutant data as parameter perturbations; evaluate evidence WT vs mutant.
* **Future with Cecilia**: receptor CIS-LU, oocyte recordings, mutants, possibly single-channel extensions.

---

## 5. Theoretical Challenges

* Likelihood for idealized single-channel traces (open problem, Colquhoun/Sigworth vs. HMM approaches).
* Macroscopic evidence modeling (oocytes).
* Mutant integration into Bayesian evidence framework.
* Potential of membrane and desensitization as later projects.

---

## 6. Organization & Workflow

* **Branches / directories** per paper:

  * `paper-gustavo-rotational`
  * `paper-characterization`
  * `paper-cumulative`
  * `paper-oocytes-mutants`
* **Two code lines**:

  * *Legacy*: fast results for cluster runs.
  * *Refactor (MacroDR v1.0)*: clean, modular, with full tests and JSON environment handling.

---

## 7. Timeline (tentative)

| Phase       | Horizon               | Actions                                                                                                        |
| ----------- | --------------------- | -------------------------------------------------------------------------------------------------------------- |
| **Phase 1** | Immediate (0–2 weeks) | Launch scheme 16/17 cluster runs, run tests, accept reproducibility issues may appear late.                    |
| **Phase 2** | Month 1–2             | Verification suite (Simulation, Likelihood, Sampling, Evidence) on minimal models; JSON persistence prototype. |
| **Phase 3** | Month 2–4             | Characterization paper draft; Gustavo rotational model results; refactor abstract class system.                |
| **Phase 4** | Month 3–6             | Oocytes + mutants modeling; generate scripts for cloud/cluster runs.                                           |
| **Phase 5** | Month 4–8             | Cumulative evidence work; Cecilia collaborations.                                                              |
| **Phase 6** | Month 6–9             | Release MacroDR v1.0 (clean code + documentation); submit characterization & Gustavo papers.                   |

---

✅ This document now integrates your newest reflections:

* Proceed with cluster runs despite potential model pathologies.
* Focus MacroDR improvements on **test integration, environment persistence, abstraction, transparency**.
* Keep dual-track development (legacy vs refactor).
* Explicitly prioritize oocyte/macroscopic evidence and mutants as next-gen directions.

---


