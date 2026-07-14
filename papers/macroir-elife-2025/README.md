# MacroIR (eLife 2025) — Paper 2 Planning Pack

This directory contains **editable planning documents** for the **second paper** we’re currently targeting:

- **Paper:** MacroIR (macroscopic interval-recursive inference for time-averaged patch-clamp data)
- **Draft/manuscript sources:** `papers/macroir-elife-2025/docs/manuscript-drafts/`
- **Repro/figures pipeline:** `projects/eLife_2025/`
- **Audio notes (transcribed):** `program/source-notes/audios/`

## How to use this pack

1. **Start at `00_master_list.md`.** It says which document owns which topic, what each one's status is, and where the documents currently contradict each other. One topic, one owner: if a fact is owned elsewhere, cite it, do not restate it.
2. Read `00_master_plan_v2.md` for the thesis and the open decisions (D-1…D-7).
3. Record any settled choice in `02_decision_log.md`, and only there.

Retired, do not read as current: `00_master_plan.md`, `01_workboard.md`, `06_repro_pipeline.md`.

## Section-by-section writing docs (one per manuscript section)

Each one says what the section must do, the hard constraints on it, a structure or draft, what stays out, and what must be verified before submission. They are the layer between the plan and the `.tex`.

- Title: `title_options.md`
- Abstract: `abstract_draft.md`
- Introduction: `introduction_plan.md`
- Theory (the likelihood family): `theory_plan.md`
- Diagnostics (the validation machinery): `diagnostics_plan.md`
- Results: `results_plan.md`
- Discussion: `discussion_plan.md`
- Materials and Methods: `methods_plan.md`
- Naming of the five algorithms: `nomenclature.md`
- Pending correction (code + supplement): `correction_idm_reconstruction.md`

## Quick links

- **Who owns what:** `00_master_list.md`
- Plan and open decisions: `00_master_plan_v2.md` (v1 is a retired stub)
- Settled decisions: `02_decision_log.md`
- Metric registry: `03_metrics_diagnostics.md`
- Sweep design: `05_experiment_grid.md`
- Figure arc: `results_plan.md`; run manifest: `docs/figure_provenance.md`
- Engine backlog: `07_code_tasks.md`
- Sources (docs + audios): `08_sources_audio_notes.md`
- Repo carve / freeze: `09_carve_plan.md`

## Repo/branch intent

- We keep everything in the **same `macro_dr` repository**.
- We keep working on **`main`** unless there is a concrete reason to branch.
- For now, we focus on **documents and planning**, not on modifying C++ code.
