# Provenance, Types, and HoTT — Conversation Transcript (Translated)

## Executive Summary
This document explores the integration of **provenance tracking** with **Homotopy Type Theory (HoTT)** and its implications for the design of a DSL (MacroIR/MacroDR). The key insights are:

- **Provenance as part of types:** Instead of treating provenance as external metadata, it can serve as indices in dependent types, making propositions about real-world measurements explicit.
- **Functions remain pure:** Provenance is carried alongside computations (like a Writer monad), preserving mathematical semantics while allowing auditability.
- **Granularity control:** Provenance should be tracked at the level of human decisions (model choice, parameters, seeds) rather than raw operations, preventing unmanageable histories.
- **Industrial alignment:** This approach parallels practices in ML pipelines (lineage graphs), databases (provenance semirings), and security (taint tracking), but extends them into type theory.
- **Compression strategies:** Depth, width, horizon, and influence thresholds determine what provenance to retain, balancing reproducibility with tractability.

The novelty lies in making provenance **first-class** in the type system, bridging curatorial practice with theoretical rigor. This could enable a DSL where each value is not just a number but a *number-in-context*, tied to time, instrument, simulation, or derivation.

---

## Context
This document is a translated and curated record of a reflective discussion on provenance, types, and their relationship with Homotopy Type Theory (HoTT) and DSL design for MacroIR/MacroDR.

---

## Part 1: Types and Functions
**Original idea:** If we define a type like `Likelihood`, then naturally we need a function `calculate_likelihood`. If we define `TransitionProbability`, we need a `calculate_transition_probability`. The challenge is how to define the semantics of these types.

In HoTT or type theory, a type is defined by its invariant — the values that belong to it. That’s the computational object, which is all we use for calculation. But beyond that, there’s the link between the computational object and the real-world object: a measured channel, a real system, a simulation, or a random number generator. Every object has a tie to a specific spacetime context. This abstraction — the relationship between program and real world — is what we’d like to capture.

---

## Part 2: Provenance as Witness
**Observation:** Provenance is like the way we track packages or money — the origin matters. Money can be “clean” or “black”; packages have tracking histories. Similarly, every datum should carry a *witness of origin*. This is a general abstraction:

```text
Data = { value : Computational,
         provenance : ProvenanceTag }
```

Where `ProvenanceTag` might be:
- `Measured(channel, time, location)`
- `Simulated(model, seed)`
- `Derived(function, parents)`

Provenance composes: if two inputs have provenance, the output inherits a combined provenance. Contamination is possible, just like mixing black and white money.

---

## Part 3: Industrial Analogies
This aligns with practices already found in industry:

1. **Data provenance in databases:** why/where provenance (Green et al., provenance semirings).
2. **Taint tracking in security:** flags propagate with computations.
3. **Lineage in ML pipelines:** Spark, TFX, MLflow track transformations and datasets.

The novelty for our DSL would be to integrate provenance into the *type system*, not as an external add-on.

---

## Part 4: The Limits of Provenance
If every operation contributed to provenance, the history would be intractable (“Borges’ map the size of the territory”). Example: a transition matrix depends on model, parameters, time step, concentration. Full provenance would include all of that. The key question: **where do we stop recording?**

Solution: policies of retention and compression:
- **Depth k:** how many ancestors to keep.
- **Width w:** how many inputs per level.
- **Horizon t:** time range.
- **Threshold ε:** keep only ancestors above influence threshold.

This yields manageable provenance graphs.

---

## Part 5: Human Decision Rate Constraint
Sabine Hossenfelder notes: humans only generate ~10 bits/sec of decisions. Thus, human choices (model selection, parameters, running a loop) are sparse. The provenance graph is thin at the human level. Internally, billions of likelihood multiplications don’t need provenance. Instead:
- Capture provenance at *scope* level (evidence loop, MCMC run).
- Store reproducibility metadata (code version, seed, dataset hash).
- Attach provenance only at the end of major blocks.

This avoids explosion: provenance size is bounded by human interventions, not raw operations.

---

## Part 6: Function Semantics vs Provenance Semantics
Mathematical functions are pure: `likelihood(x|θ)` gives the same result regardless of whether `x` was simulated or measured. Provenance passes *around* the function, not through it. This is exactly the pattern of a **Writer monad**: the value is computed as usual, while the log (provenance) is accumulated on the side.

Thus, provenance is not central to the algorithm but crucial for curation, reproducibility, and audit.

---

## Part 7: Types as Propositions, Extended
In HoTT, `Type = Proposition`. A type is the proposition that some value inhabits it. But in practice, the proposition is not just “4 picoamperes” — it is “4 picoamperes on channel 3, instrument X, at time t0, with calibration Y”.

This unites provenance with HoTT: the provenance attributes become **indices of the type**. For example:

```text
Current(Instrument=i, Time=t, Channel=c) : Type
4pA : Current(i=Rigol123, t=2025-09-20T19:45, c=3)
```

Thus, the proposition is about the external world, not just an abstract number.

In HoTT, paths express equivalences. With provenance indices, paths could express equivalence relations between different provenances — e.g. whether two currents measured with different instruments can be treated as equal.

---

## Conclusion
- **Provenance = indices of types.**
- **Functions stay pure.** Provenance flows externally, like a Writer monad.
- **Granularity is bounded**: keep human-level decisions, not machine-level operations.
- **Novelty:** Most industrial systems treat provenance as metadata; in MacroIR/MacroDR, provenance could be first-class in the type system, grounded in HoTT.

---

## Next Step
Potential formalization: a two-layer provenance system (macro vs micro) with semiring composition and policies for compression, expressed as dependent types in the DSL.

