# Paper 2 (MacroIR / eLife 2025) — Sources & Audio Notes Index

Index of the sources this plan derives from. Treat transcripts as idea-capture, not citable text (automated transcription, occasional word errors).

## Manuscript source (head)
- `papers/macroir-elife-2025/docs/manuscript-drafts/elife_paper.tex` — the live draft to build on.
- Superseded drafts (`elife-macroir*.tex`, incl. `-merged`, `-revised*`) are earlier variants kept as history.

## Framing / plan
- `From molecular mechanisms to data back and forth PROGRAM.md` — research-program framing (this paper = bridge 2).
- `00_master_plan_v2.md`, `02_decision_log.md`, `09_carve_plan.md`.

## Prior published work
- Moffatt & Pierdominici-Sottile, Comm Biol 2025 (P2X2) — the demonstration; `docs/bibliography/Moffatt_PierdominiciSottile_2025_*`.
- Moffatt 2007, Biophys J — the origin of MacroR; `docs/bibliography/Moffatt_2007_*`.

## Prototype pipeline
- `projects/eLife_2025/ops/local/*.macroir`, `ops/slurm/`, `ops/clusters/`; `projects/eLife_2025/figures/paper/*.Rmd`.

## Audio transcripts
All under `program/source-notes/audios/`.

### MacroIR 13 (2026-05 → 2026-07) — current source of record
Folder `Chat de WhatsApp con MacroIR 13/`. Carries the current consensus. High-signal:
- **2026-05-19** — inaugurates the paper; the abstract problem (must say more than Comm Biol).
- **2026-05-31 / 06-02 / 06-03** — Information Distortion Matrix blow-ups; Fisher singular/indefinite; measure the Hessian at the optimum; posterior-vs-likelihood split (later cut to likelihood-only).
- **2026-06-09 / 06-10** — the trust-coefficient discontinuity bug found and fixed; decouple α_μ (mean) from the covariance down-date; IR canonical, IRT/Taylor cut.
- **2026-06-22 / 06-23** — the paper narrative dictated: two Gaussian approximations, three regimes (multinomial / telegraphic / Gaussian), the diagnostics, the ranking; MacroMR strawman; MacroIR ≈ time-augmented Kalman as hypothesis.
- **2026-07-02 → 07-08** — figures: heatmaps of bias/distortion, correlation/sample decomposition, the shift to the Gaussian-Fisher anchor for the definitive figures.
- **2026-07-11** (18.21 + 20.34/35/36) — the scope lock: 2-state, non-stationary, macroscopic-only; MacroIR sole survivor; the Fisher-to-0 result (info about the original channel number vanishes on relaxation); evidence correction as discussion; micro / >2-states / stationary / experimental-data as explicit open doors.

### MacroIR 10 (2025-12 → 2026-01) — earlier, still valid pointers
- **2025-12-20** — 4-figure narrative; avoid MicroIR to reduce attack surface.
- **2025-12-25** — score/FIM efficacy first; partial vs total likelihood gradients.
- **2026-01-11** — Cov(Σ score) vs Σ Cov(score) diagnostic.
- **2026-01-14** — variance inflation factor; residual tests; inference/evidence distortion.
