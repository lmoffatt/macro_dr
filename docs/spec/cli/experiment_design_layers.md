# MacroDR Experimental Data Design Reference

This document captures the layered design around experiments, recordings, and their transformed counterparts within the MacroDR codebase. It outlines the semantic distinctions, data types, and evolution path.

---

## 🔹 Core Concepts

### 1. `Protocol`

- Describes the intended stimulus conditions (e.g., ATP concentrations, time steps).
- Immutable plan or configuration of a single trial.

### 2. `Recording`

- Captures the actual outcome of applying a `Protocol`, either from simulation or measurement.
- May be raw or preprocessed, but always observational.

### 3. `Experiment`

- Bundles a single `Protocol` with its corresponding `Recording`.
- Canonical representation of a single trial.

```cpp
struct Experiment {
  Protocol  protocol;
  Recording recording;
};
```

---

## 🔹 Transformed Experiments

### 4. `TProtocol`

- A modified or idealized version of a `Protocol`, typically with altered resolution or interpolated stimuli.

### 5. `TRecording`

- A transformed or processed view of a `Recording`, e.g., binned, averaged, sliced.

### 6. `TExperiment`

- Combines a `TProtocol` and a `TRecording`.
- Represents the working object for current CLI and simulation logic.
- Often standalone, with no source `Experiment`.

```cpp
struct TExperiment {
  TProtocol  protocol;
  TRecording recording;
};
```

---

## 🔹 Future: View over Full Experiment

### 7. `TExperimentView`

- A transformed view *of* a full `Experiment`.
- Retains a reference to the original for traceability.
- Enables multiple post-processed representations of a single experiment.

```cpp
struct TExperimentView {
  const Experiment* source;
  TExperiment        view;
};
```

---

## ✅ Naming Summary

| Concept                      | Type Name         |
| ---------------------------- | ----------------- |
| Raw protocol                 | `Protocol`        |
| Raw result                   | `Recording`       |
| Canonical trial              | `Experiment`      |
| Idealized protocol           | `TProtocol`       |
| Processed result             | `TRecording`      |
| Standalone transformed trial | `TExperiment`     |
| View over full experiment    | `TExperimentView` |

---

## 🧩 Notes and Decisions

- `TExperiment` is the current working type (used in `frac_simulation`, `simulate`, etc.)
- `Experiment` will be introduced when upstream pipelines are integrated.
- `TExperimentView` allows analysis to trace transformations without losing provenance.
- CLI and DSL operate primarily on `TExperiment` now; will evolve to support full tracing.

---

## 📌 Future Considerations

- Should `TExperimentView` inherit from or wrap `TExperiment`?
- Introduce CLI metadata to distinguish between `TExperiment` and `TExperimentView`?
- Migrate loaders (`load_experiment`, etc.) to support both standalone and traced modes?

---

## 📌 Transitional Decision on Naming

Currently, the system operates using an `Experiment` type that effectively encodes only protocol information, and a `Recording` that is independently managed. These serve the same roles as `X` and `y` in a standard fitting pipeline.

Although this does not reflect the final design (which will include a true `Experiment = { Protocol, Recording }` and distinct `TProtocol`, `TRecording` types), there is no immediate benefit in renaming these structures.

**Decision:**

- Do **not** rename `Experiment` to `TProtocol`, nor `Recording` to `TRecording` yet.
- These names will be updated when the full `Experiment` type is introduced.
- This avoids breaking working code and keeps current logic simple and practical.
- Comments or documentation should note the transitional role of current names.

This approach prioritizes stability and clarity while leaving a clear path for migration when the system's data model matures.

