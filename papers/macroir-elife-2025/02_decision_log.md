# Paper 2 (MacroIR / eLife 2025) — Decision Log

Record decisions here so we can “rewind” using git history.

## Confirmed (current defaults)

- **Repo:** single repository (`macro_dr`), no secondary repo.
- **Branch:** `main`.
- **Paper identity:** “Paper 2” = MacroIR (active draft sources in `papers/macroir-elife-2025/docs/manuscript-drafts/`).
- **Validity map axes:** Δ/τ_min vs N_ch.
- **Algorithms to include in comparisons:** `NR`, `R`, `MNR`, `MR`, `MNRV`, `MRV`, `IR`, `IRV`.
- **Keep MicroIR out of the main paper** to reduce scope/attack surface.
- **Audio archiving:** keep both audio sources and transcripts (user preference: “both”).

## Open (needs explicit decision)

### D-001 — Biological case study placement

- **Question:** Should the real case study be in main text or Supplementary?
- **Default:** Supplementary.
- **Why it matters:** main-text case study increases narrative strength but increases review surface area and risk of scope creep.

### D-002 — Audio sources in git (normal vs LFS)

- **Question:** If audio is tracked, do we want Git LFS?
- **Default:** track normally until it becomes painful, then migrate to LFS.
- **Notes:** current audio footprint lives under `program/source-notes/audios/`.

### D-003 — Thresholds for “valid” on the regime map

- **Question:** Which numeric cutoffs define “valid/invalid” for each metric?
- **Default:** start with rank-based comparisons (MacroIR best → green) + add absolute thresholds later.
