# D-0: freeze commit and re-run scope — DECIDED 2026-07-15

Decided with Luciano over the 2026-07-15 session. This header is the decision; B-0's original
three-option brief (built on a premise that turned out false) is kept below as the evidence trail,
marked superseded. The one-line record is in `../../_program/decisions.md` §4 (this was
`02_decision_log.md`, which the 2026-07-20 split retired into the cross-paper log).

## The decision

- **Freeze commit: `1c2ae6f`.** It is the home of the gaussian-Fisher macro production data.
- **Multi-commit provenance accepted.** Every CSV stamps its own engine hash in row 1 (verified:
  `1c2ae6f` files read `1c2ae6f`, `87889e6` files read `87889e6`). Provenance is self-documenting
  per file, so there is no need to homogenize everything to one commit, and no rerun for that reason.
- **`433ed13` stays as-is:** the numerical-Fisher family, kept only to demonstrate that the gaussian
  Fisher ≈ the finite-difference Fisher. It is a separate artifact with a separate purpose; a
  numerical cell costs ~400-450 CPU-h against 181 for a gaussian one, so re-running it would be
  pure waste.
- **`87889e6` is out of scope:** `micro_R` / `micro_IR` reference data for future (micro) papers,
  not this two-state macro paper.
- **IR gaussian is complete:** 10 N_ch × 6 noise = 60 cells on disk. No work.
- **E-1…E-5 are decoupled to `main`** as code hygiene; they do not touch this paper's data. The
  existing ensembles were seeded from `random_device` with `seed=0` and the resolved seed was never
  logged (`legacy/mcmc.h:37`), so they are not bit-reproducible and cannot be made so retroactively.
  The Monte-Carlo standard is statistical equivalence, not bit identity; the Data Availability
  statement says so. Landing E-1 on `main` is what stops the mistake from recurring (future papers
  born reproducible), which is the reason to do it there, not here.

## The fill (in progress on Dirac, 2026-07-15)

The only data the paper needs beyond disk: the four non-IR algorithms brought onto the gaussian
anchor at the reference noise grid.

- **Grid:** algorithms {NR, NMR, R, MR} × N_ch {10, 100, 1000, 10000} (canonical) × noise {0.1, 1, 10}.
- **Interval is internal to every cell:** each `battery_*_G.csv` carries the 7-value
  `interval_in_tau` sweep {1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01} inside it (dispatch sets
  `axis_interval`, `dispatch_figure_3_G.sh:169`; verified present in a non-IR file). So N_ch ×
  interval × noise, the paper's three map axes, all come from this one fill at no extra cost.
- **Size:** 36 cells. NR/R/MR already have their noise-0.1 column on disk, so they need noise {1, 10}
  (2 × 4 = 8 each = 24); NMR was absent at *every* noise, so it needs {0.1, 1, 10} (3 × 4 = 12).
  24 + 12 = 36.
- **Cost:** ~6,500 CPU-h at 181 CPU-h/gaussian cell (measured, see below). ~2-3 days of Dirac queue.
- **Status:** running now. The batch covers noise 1 and 10 for all four, and NMR at all three noises
  (0.1 included), confirmed with Luciano. **The fill is complete when this batch finishes; nothing
  more to launch.**

## Per-cell cost, measured (kept from B-0, still valid)

From file modification times, 32 CPUs per SLURM job (`dispatch_figure_3.sh:200`):
- **Gaussian cell** (no finite-difference Fisher): 5.67 h end-to-end median across 39 cells ⇒
  **181 CPU-h/cell**.
- **Numerical cell** (with FD Fisher): ~12-14 h ⇒ ~400-450 CPU-h/cell (the FD Fisher is a 2.3× tax).
  This is why `433ed13` is not re-run.
- **Run 3 (time-resolved dumps, Figs 3/4/S4/S5):** the whole five-algorithm ensemble regenerates in
  ~4 minutes (`runs/run-20260630-202545`). Effectively free, and E-3 does not break its notebooks
  (they collapse the doubled rows by mean; the collapse is idempotent).

## What B-0's original three-option brief got wrong (superseded)

The (a)/(b)/(c) options at 900 / 14,000 / 39,000 CPU-h all existed to satisfy one premise: that the
deposited code must be a single commit that produced all deposited data. **Once multi-commit
provenance is accepted (each CSV self-stamps its hash), that premise dissolves and the real job is
the 36-cell fill above.** Specifically:

1. **Option (a) "freeze at the current hash" was already impossible.** HEAD `829bca8` had already
   diverged from both `433ed13` and `1c2ae6f` (commit `019dbb6`, a micro-only change), so the
   code-equals-data invariant was broken before the decision was even posed. Moot under multi-commit.
2. **The "1,830 CPU-h for 4 cells" estimate was wrong twice:** the missing cells are gaussian
   (181 CPU-h, not the numerical 400-450), and NMR needed noise 0.1 too, so it was more cells at a
   lower unit price.
3. **The claimed NR-10k gap did not exist** (NR at 10,000 channels is present and complete on
   `1c2ae6f`); the only genuinely-absent algorithm was NMR.
4. **`05_experiment_grid.md` §2 calls the primary map "N_ch × K_off", but K_off was never swept.**
   On disk `k_off` is fixed at 100 and the second axis is instrumental noise (`Current_Noise`). The
   paper's real map axes are **N_ch × interval × noise**, matching the abstract. `02_decision_log.md`
   carried the same "N_ch × K_off" error and is corrected in the same commit. Whoever writes Results
   §5 must promise the noise/interval map, not a K_off map, unless a brand-new K_off sweep is
   commissioned (no option here funds one).
