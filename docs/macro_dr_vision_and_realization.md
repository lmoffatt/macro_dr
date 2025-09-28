# MacroDR: A White Paper on Proof-Oriented Inference Engines for Scientific Modeling

## Abstract
MacroDR is a next-generation scientific inference engine for biophysical modeling, built on the principle that **software is a mathematical instrument**. It unifies rigorous proof-oriented design, domain-specific language (DSL) semantics, modular architecture, and provenance-aware reproducibility into a system capable of both immediate research applications and long-term theoretical contributions. This white paper outlines the **vision**, the **conceptual foundations**, the **technical realization**, and the **research roadmap** for MacroDR.

---

## 1. Introduction
The central motivation behind MacroDR is to bridge the gap between **scientific aspiration** and **software engineering discipline**. Traditional scientific software often accumulates technical debt, ad hoc algorithms, and opaque data handling. MacroDR instead aspires to be a *living demonstration* of **Proof-Oriented Design (POD)**, where each program element embodies a guarantee, and the system’s coherence is its meaning (see *algebraic_CLI.md*).

The research context is ion channel kinetics and Bayesian model evidence. However, MacroDR is designed generically: any stochastic model inference task can benefit from its DSL, safety kernel, and modular structure.

---

## 2. Vision
MacroDR’s vision can be distilled into four guiding principles:

1. **Correctness by Construction**  
   Every construct encodes its invariants. A Q-matrix is not just a matrix; it is guaranteed nonnegative off-diagonals and row sums to zero (see *safety_explosion_limit.md*).

2. **Transparency and Reproducibility**  
   Provenance is tracked as a first-class property. Every simulation, parameter set, or evidence run can be traced to its origin, seed, or measurement (see *provenance_types_hott.md*).

3. **Extendability without Fragility**  
   New commands, models, or algorithms can be added without touching central registries. Invariance is enforced locally, contagiously, without brittle boilerplate (see *dispatch_decision.md*, *unique_ptr_decision.md*).

4. **Dual-Track Development**  
   Immediate cluster jobs continue to run with the legacy system, while the refactored MacroDR evolves with clean modularity and CI-driven guarantees (see *Project Roadmap MacroIR.md*).

---

## 3. Conceptual Foundations

### 3.1 Proof-Oriented and Contagiously Safe Design
- **Types as Proofs**: Constructing an object is equivalent to proving its invariants. A `ScientificSafe<Q>` is a proof that Q satisfies stochastic properties (*Proof_Oriented_Design_Manifesto.md*).
- **Contagion Principle**: Once established, invariants propagate downstream unless explicitly escaped (*contagiously_safe_design.md*).
- **Composite Domains**: To avoid exponential growth of wrapper types, multiple invariants are bundled into composite safety domains (*composite_safety_domains (1).md*).

### 3.2 HoTTification
Inspired by Homotopy Type Theory (HoTT), MacroDR treats **types as propositions**. For example:
```text
Current(Instrument=i, Time=t, Channel=c) : Type
4pA : Current(i=Rigol123, t=2025-09-20T19:45, c=3)
```
Here the *proposition* is not just “4 picoamperes” but “4 picoamperes measured on channel 3, instrument Rigol123, at time t0, with calibration Y.” Provenance becomes part of the type (*provenance_types_hott.md*).

### 3.3 Algebraic CLI
Commands are not arbitrary; they form algebraic structures. Running `simulate` followed by `infer` and verifying with `sav-check` corresponds to an approximate identity:
```text
infer(simulate(θ)) ≈ θ
```
This elevates the CLI into a mathematical contract (*algebraic_CLI.md*).

### 3.4 Lingua Franca and Boundary Contracts
- **Canonical Data Model**: JSON-like value tree at boundaries, enriched by schemas (*general_scientific_software_design_theory.md*).
- **Contracts**: pre/postconditions checked at CMD surfaces; violations return `Maybe_error` instead of propagating silently (*ADR-001-modules-and-dependencies.md*).
- **Vector View**: optional numeric projection ensures hot loops remain efficient (*modules.md*).

### 3.5 Provenance as First-Class
- **Writer Monad Semantics**: functions remain pure, provenance accumulates externally (*provenance_types_hott.md*).
- **Compression Strategies**: retain provenance at the level of human decisions (model choice, parameters, seeds) rather than raw multiplications (*provenance_types_hott.md*).
- **Experiment Layers**: Protocol + Recording = Experiment; transformations yield TExperiment; provenance-preserving views (TExperimentView) allow multiple post-processed perspectives (*experiment_design_layers.md*).

---

## 4. Technical Realization

### 4.1 Architecture & Modules
- **DSL**: Parser, AST, registry; no domain logic (*dsl_structure_summary.md*).
- **CLI**: Thin entrypoint; orchestrates script ingestion (*cli.md*).
- **CMD**: Command surface with arguments, pre/postconditions, orchestration (*ADR-001-modules-and-dependencies.md*).
- **Core**: Implementations behind CMD; orchestrates Models, Inference, IO, Math (*modules.md*).
- **Domain Entities**: Canonical types (Experiment, Recording, Parameters) (*experiment_design_layers.md*).
- **IO**: JSON/CSV serialization, schema validation (*modules.md*).
- **Models/Inference/Probability/Math**: Deterministic forward simulation, Bayesian likelihood, MCMC, linear algebra, distributions (*modules.md*).
- **Utils**: Error handling (`Maybe_error`), metaprogramming, memoization (*modules.md*).

### 4.2 Runtime & Type System
- **Grammar**: assignments, function calls, literals (*grammar.bnf.md*).
- **Type Deduction**: `compile_expression` maps untyped to typed (double, string, identifiers, function calls) (*type_rules.md*).
- **Runtime Rules**: Execution via `typed_expression::run`, assignment binding in `Environment` (*run_rules.md*).
- **Dispatch**: dynamic_cast for AST downcasts, postponed visitor pattern (*dispatch_decision.md*).
- **Ownership**: unique_ptr for IR nodes, NodePtr borrow deferred (*unique_ptr_decision.md*).

### 4.3 CLI & DSL Functions
- **CLI Verbs**: `run`, `eval`, `compile`, `describe`, `help` (*command_syntax.spec.md*).
- **DSL Functions**: simulation, likelihood, thermo evidence, experiment loaders, prior construction, utility helpers (*functions.md*).

### 4.4 Testing & CI
- **Unit Tests**: param→Q→param round-trips, expm vs spectral, likelihood finite (*run_rules.md*).
- **Property Tests**: score = 0 in expectation; Cov(score) = FIM (*Project Roadmap MacroIR.md*).
- **Regression**: DSL regression tests, CLI smoke tests (*testing.md*).
- **CI Integration**: build times, template instantiations tracked (*testing.md*).

### 4.5 Documentation & ADRs
- **spec/**: CLI syntax, DSL grammar, type rules, runtime rules (*spec/ directory*).
- **ADR-001**: modules & dependencies (*ADR-001-modules-and-dependencies.md*).
- **ADR-002**: lingua franca (*modules.md*).
- **Dispatch decision**: dynamic_cast (*dispatch_decision.md*).
- **Ownership decision**: unique_ptr (*unique_ptr_decision.md*).
- **Roadmap**: M0–M5 phases (*implementation_roadmap.md*).

---

## 5. Research Context and Papers
- **Rotational kinetic models** (with Gustavo): extend published work using rotation-coupling (*Project Roadmap MacroIR.md*).
- **Characterization of MacroDR**: present validation (Simulation, Likelihood, Sampling, Evidence, FIM, Score) (*Project Roadmap MacroIR.md*).
- **Cumulative evidence**: theoretical and empirical framework for combining datasets (*Project Roadmap MacroIR.md*).
- **Oocytes + mutants**: macroscopic evidence modeling; mutant perturbations (*Project Roadmap MacroIR.md*).
- **Future**: receptor CIS-LU collaborations, single-channel idealization (*Project Roadmap MacroIR.md*).

---

## 6. Roadmap
- **M1 Registry Consolidation**: commands extracted from CLI into CMD (*implementation_roadmap.md*).
- **M2 Core Surfaces**: `simulate` and `load_experiment` stabilized (*implementation_roadmap.md*).
- **M3 Lingua Franca Pilot**: Recording + Experiment JSON schemas (*implementation_roadmap.md*).
- **M4 CMD Migration**: likelihood, thermo evidence moved from legacy (*implementation_roadmap.md*).
- **M5 Performance**: vectorization, CMake target separation (*implementation_roadmap.md*).

---

## 7. Conclusion
MacroDR is more than a codebase: it is a **proof-oriented research program**. It embodies the principle that scientific computation must be:
- **Safe** (invariants guaranteed),
- **Transparent** (provenance explicit),
- **Reproducible** (environments persistable),
- **Evolvable** (modular architecture),
- **Mathematically principled** (HoTTification, algebraic CLI).

The union of audios and documents reveals a coherent trajectory: from philosophical reflection to engineering scaffolding. The system is already serving immediate scientific goals while charting a course toward a general theory of proof-oriented design in computational science.

