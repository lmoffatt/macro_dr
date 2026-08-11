# Sources & Audio Notes Index

> Updated: 2026-08-11. Shared across the three papers. Which audio/chat is the source of record, and
> which wins when two disagree.

Index of the sources the program derives from. Treat transcripts as idea-capture, not citable text (automated transcription, occasional word errors).

## Manuscript source (head)
- `papers/1_method/docs/manuscript-drafts/elife_paper.tex` — the live draft to build on.
- Superseded drafts (`elife-macroir*.tex`, incl. `-merged`, `-revised*`) are earlier variants kept as history.

## Framing / plan
- `program.md` — the three-paper map; `research_program.md` — the wider program framing.
- `1_method/00_plan.md`, `decisions.md` (+ `1_method/decisions.md`), `carve_plan.md`.

## Prior published work
- Moffatt & Pierdominici-Sottile, Comm Biol 2025 (P2X2) — the demonstration; `docs/bibliography/Moffatt_PierdominiciSottile_2025_*`.
- Moffatt 2007, Biophys J — the origin of MacroR; `docs/bibliography/Moffatt_2007_*`.

## Prototype pipeline
- `projects/eLife_2025/ops/local/*.macroir`, `ops/slurm/`, `ops/clusters/`; `projects/eLife_2025/figures/paper/*.Rmd`.

## Audio transcripts
All under `program/source-notes/audios/`.

### MacroIR 13 (2026-05 → 2026-08) — current source of record
Folder `Chat de WhatsApp con MacroIR 13/`. Carries the current consensus. Each `.mp3` has a `.md`
beside it with timestamps; the originals are under `raw/` as `.ogg`. High-signal:
- **2026-08-11** (14.05 → 14.15, four audios) — the venue holds and the route is bioRxiv first
  (`decisions.md` §1); automatic differentiation is an advance over the P2X2 build and must be
  stated, now in `1_method/docs/manuscript-drafts/sections/06_methods.tex` under Reproducibility;
  the deliverable must carry the validation and not only the algorithm (`carve_plan.md`, last
  section); and a pass over the figure roster. **Read that last one on the numbering from before
  2026-08-10**, confirmed by the author on 08-11: its Figure 5 and Figure 6 sit one off against
  `1_method/decisions.md` "The figure set", nothing is being moved, and the list is a list of what
  should exist. One item in it was already stale when spoken, the magnitude and anisotropy
  decomposition having been built on 08-02 as a Figure 4 supplement.
- **2026-08-10** (14.18 → 14.32, six audios) — the Introduction dictated aloud, and the strongest of
  the batch. The purpose stated as model comparison: a valid likelihood is a necessary condition for
  computing evidence, so the question is the conditions under which an undistorted likelihood is
  available. The rationale of the three ranges, which was written nowhere and is now `axes.md` §2ab.
  Usage recommendations by regime, which sit against the demarcation of `1_method/approach.md` §6.
  And a reorientation of the whole paper to MacroIR plus its validation, **which was withdrawn by the
  author on 2026-08-11 as a temporary weakness and never adopted**; the merge stands, and the
  tombstone is at `program.md` §9. Read the batch against the clock: the Introduction was rewritten
  and committed at 23:25 the same night, so several of these were answered within hours. This is the
  clearest case in the whole corpus of why a transcript is idea-capture and not a decision: the audio
  states the reorientation and carries no retraction, and the retraction is what happened.
- **2026-08-09** (08.32, 08.36) — derive it twice, from MacroR and from the Kalman side, and show the
  two routes meet, on the precedent of the 2007 paper; why MR fails, stated as a mechanism; and a
  display item for the numerical Fisher being indefinite. All three in `1_method/approach.md` §10,
  items 9 to 11.
- **2026-08-10 at 22.29.03 carries nothing.** Fifty-one seconds, aborted mid-sentence, the tail of the
  transcript is noise. Recorded here so that nobody spends time on it twice.
- The nine audios of 08-09 and 08-10 were pulled on 08-11, after the 08-11 batch had already been
  transcribed, so on disk they look newer than audios that are in fact older.
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
