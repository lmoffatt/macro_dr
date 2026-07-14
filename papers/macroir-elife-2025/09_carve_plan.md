# Carve plan: macro_dr → macroir-validity

> Frontier map for splitting this paper into its own repo. Decide the boundary now (no files move); execute the 3 steps at **code freeze** (when the Gaussian-Fisher figures are final and the Dirac runs are done). Open decisions marked **[Q]**. Comment inline.

## Why curation, not a move

`projects/eLife_2025/` is **92 GB**: `figures/data/` 73 GB, `runs/` 19 GB. `papers/macroir-elife-2025/` is 6 MB. So ~99% of the material must NOT enter git. The paper repo is a small, curated bundle; the bulk is regenerable from the engine or goes to Zenodo. Cloning the folder wholesale would drag 92 GB into a paper repo.

## Three buckets

### A. Goes into `macroir-validity` (small, git-tracked)

- **Manuscript.** `papers/macroir-elife-2025/docs/manuscript-drafts/elife_paper.tex` (the live head) + `biblio.bib` + `elife.cls` and template assets needed to build. Drop the superseded drafts (`elife-macroir*.tex`) — they stay as history in `macro_dr`.
- **Program + paper-facing docs.** `From molecular mechanisms to data back and forth PROGRAM.md`, `analysis_figure_S1_score_mean.md`, `docs/notation_map.md`, `docs/corrected covariance justification.md`.
- **Figure scripts.** `projects/eLife_2025/figures/paper/*.Rmd` + `*_caption.md`, plus the `in_progress/figure_5|6|7_*.Rmd` that become final. NOT `to_classify/` or `archive/`.
- **Run configs (reproducibility set only).** From `ops/`: the FINAL `.macroir` (`figure_1.macroir`, `figure_2.macroir`, `figure_3.macroir`, `figure_3_mle_G.macroir`, the Gaussian-Fisher dispatch) + the dispatch/slurm/cluster scripts (`ops/slurm/`, `ops/clusters/`, `build_cluster.sh`). NOT the ~20 debug/old/test/IRT/micro variants (`figure_2_debug*`, `figure_2_old`, `figure_2_previous`, `figure_2_v0`, `*_IRT*`, `figure_2_micro*`, `figure_2_error`, ...).
- **Model inputs.** `scheme_CO_par.csv`, `scheme_CO_prior.csv`, `scheme_CO_model_description` (the 2-state model this paper uses). NOT the other schemes.
- **Figure-ready data (small).** Only the reduced CSVs the final `.Rmd` consume to draw figures. See the handoff convention below.
- **New files to author at carve time.** `README.md` (what the repo is, how to build), repro instructions, `CITATION.cff`, `.gitignore` (excludes raw data), and a note pinning `macro_dr @ <tag> (<hash>)`.

### B. Stays in `macro_dr` (the engine, shared across papers)

- The C++ engine (`src/`, `include/`, `CMakeLists.txt`, `tests/`, build system).
- `theory/macroir/` — cross-paper program knowledge. The paper's SI is a *frozen copy* of the specific `.tex` it cites, placed in the paper repo at freeze; the living theory stays here. **[Q]** confirm: SI = snapshot-copy into the paper repo, theory stays in `macro_dr`.
- The superseded manuscript drafts and `archive/` folders (history is preserved here).
- Everything not specific to this paper.

### C. Neither repo (heavy → Zenodo or regenerable)

- `projects/eLife_2025/runs/` (19 GB) — raw run outputs. Regenerable from engine + configs; archive to Zenodo if you want them citable.
- `projects/eLife_2025/figures/data/` (73 GB) — raw/intermediate. The figure-ready subset (small) is reduced into bucket A; the rest is regenerable/Zenodo.
- Root `figure_2_dlikelihood_*.csv` (94 MB each) — obsolete Taylor-variant (MRV/IRV) data from March, cut from the paper anyway.
- `logs/`, `diag.log`, `temp_script_file.txt` — scratch.

## Data handoff convention (Dirac → paper repo)

```
macro_dr @ tag  +  run configs   ──run on Dirac──►  raw outputs (large)
                                                          │  gitignored / Zenodo
                                                          ▼
                                                   reduction script
                                                          │
                                                          ▼
                                            data/figure-ready/*.csv  (small, git-tracked)
                                                          │
                                                          ▼
                                                    *.Rmd  ──►  figures/
```

- `data/raw/` in the paper repo is **gitignored** (populated from Dirac/Zenodo when someone reproduces).
- A reduction/aggregation step (part of the repo pipeline) turns raw → `data/figure-ready/` (small, tracked).
- `.Rmd` read only from `data/figure-ready/`, never from raw. **[Q]** does a reduction step already exist, or do the `.Rmd` currently read raw directly? If the latter, splitting the read into raw→reduced is the one real rewiring cost.

## Open boundary decisions [Q]

1. **Process docs.** Do the internal board docs (`00_master_plan_v2`, `01_workboard`, `02_decision_log`, `06_repro_pipeline`, `07_code_tasks`, `08_sources_audio_notes`) go into the paper repo as a `process/` folder (full transparency), or stay in `macro_dr` as dev-process? Default: paper-facing docs go, pure-process board stays.
2. **eLife-specific assets.** `Authors/eLife_LaTeX_template.zip`, `docs/elife-author-instructions.md` — venue-specific. Keep in a `submission/` folder, or drop for a venue-agnostic repo? Default: keep in `submission/`, since the repo name is already venue-agnostic.
3. **SI theory snapshot** (see bucket B): confirm the copy-at-freeze convention.
4. **Reduction step** (see handoff): confirm whether it exists.

## The 3 steps (at freeze)

1. **Consolidate** the bucket-A material into one folder (rewires figure/data paths — do this once, calmly).
2. **`git init`** a fresh repo from that folder. History starts clean; the prehistory stays in `macro_dr` (git never forgets), so nothing is lost. Tag `macro_dr` at the freeze commit (`macroir-validity-submitted-v1.0`) and record the hash in the repo.
3. **Delete** the bucket-A folder from `macro_dr`. The engine and `theory/` stay. `macro_dr` becomes just the engine + shared knowledge.

Subsequent papers start in their own folder/repo from day one and skip the carve entirely; only this paper needs it because it was born inside the monorepo.
