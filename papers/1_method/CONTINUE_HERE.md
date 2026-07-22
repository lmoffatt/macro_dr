# Continue here — figure work, state at 2026-07-22

> A pointer, not a second copy. Every fact below is owned by a document that already holds it; this
> file exists so a session starting cold knows **which** document, **what is in flight**, and **which
> traps cost hours to find**. If it starts duplicating those documents, delete it.

## Read in this order

1. `figures_build_plan.md` — build order, what each figure needs, the dirac dispatch commands, and
   §3b/§3c/§3d for what was built on 21–22 July and what it measured.
2. `results.md` §Fig 1–5 — the claims each figure now carries, with numbers.
3. `decisions.md` — the Fig-1 roster call, and the VR branch, now **resolved**.
4. `../_program/provenance.md` — which run made which data.

## Where the figures stand

- **Fig 1** built: the four-column ladder R → MR → VR → IR, over `figure_1_panels.R`, shared with
  `figure_1_all.Rmd` (six columns, orientation only, its own PDFs).
- **Fig 2** built: the A-strict roster on the **Gaussian anchor**, VR read from `0ffbda7` and the rest
  from `1c2ae6f` by a search path, so a new algorithm joins by itself when its files land.
- **Figs 3 and 4** built on regenerated dumps (engine `0ffbda7`, **seed now fixed** at 20260722; it
  was `0`, i.e. irreproducible). NR and NMR stay in the dumps and out of the panels, which is what
  lets Fig 4 keep its argument.
- **Fig 5** rebuilt with a **new subject**: where IR itself stops being faithful, in both moments. The
  old `figure_5_master*` are archived under `figures/archive/paper_superseded_20260722/`.
- **Figs 6 and 7** not built. The decomposition of the distortion into sample and correlation parts
  came out of Fig 5 and has no home yet.

## In flight, and what to check when it lands

A dirac run launched 2026-07-22:

```
RUN_DIR=1c2ae6f NCHS="10 20" N_SIMS="100000 100000" N_NOISE="0.05 0.05" \
GROUP_SIZE="100 250 1000" INTERVALS="1 0.05" N_ALGO="macro_IR" \
  projects/eLife_2025/ops/slurm/dispatch_figure_3_G.sh dirac
```

It exists to answer **one** question: with 1000 fits instead of 100, does the N_ch direction of
Fig 5's panel A cross below one, agreeing with the map, or is the disagreement structural? Look at
`figure_5A point 1` in the knit output. Second question, from the same files: does the skewness of the
MLE fall with group size as the CLT predicts, and where does it become small enough for an ellipse to
mean anything.

`dispatch_figure_3_G.sh` now takes **`INTERVALS`**. The seven-value grid was not a table, it was a
formula (n_samp = 500·Δ, n_step = 2/Δ and 4/Δ), so any multiple of 1/500 works and 0.002 is the floor.
The default reproduces the historical grid exactly; that was verified before use.

## Traps. Each of these cost real time; none is obvious from the code

- **The score-based maps and the MLE ellipses are different objects.** The map compares J against G
  unmarginalised and fits nothing; the ellipse compares an estimate's covariance against the sandwich
  marginalised over the other parameters. The ellipse is always milder, and by N_ch 200 it reads
  *calibrated* where the map says 1.25. Do not expect them to match, and do not "fix" one to the other.
- **group_size 10 at few channels reverses a sign.** Its finite-sample inflation is comparable to the
  true deflation of the N_ch direction (1.58 measured against 0.65 true). Group 100 is the floor.
- **The MLE is not normal at few channels** (skewness 4.13 at N_ch 10, group 10). This is a property
  of the estimator, **not a defect of the likelihood**, and it limits what any covariance or ellipse
  can mean. Keep the two statements apart in the text.
- **Estimating shape needs many fits, not big groups.** SE of a skewness is √(6/n): the exploratory
  `in_progress/figure_normality_shape.Rmd` looks structured but its within-panel variation is pure
  noise (sd/SE ≈ 1 in 14 of 16 blocks). Only N_ch 10 at group 10 carries real signal. The honest
  version of that figure is skewness against N_ch with a band, not a contour.
- **The noise label is not the noise.** `current_noise = label/1000`, and the CSV column
  `noise_in_conductance_tau` stores the **label**.
- **Do not write "dex".** It is astronomy jargon; use log10 units and give the percentage.
- **`Likelihood_Information_Distortion_Reconstituted` is unusable**: it implements an identity that is
  false (the exact relation is the congruence GIDM = K·CDM·Kᵀ) and is NaN in every Gaussian run. What
  *is* exact is log det GIDM = log det GSDM + log det CDM, to 1e-14.
- **`82b956f` is nsim 1000** while everything else is 10000. A glob written `nsim_*` pools across it
  and walks straight into the Jensen hazard. Hardcode `nsim_10000`.
- **Figure-1 CSVs are not in git** (`.gitignore` excludes `*.csv`) and `seed = 0` means random, so a
  re-run of a figure-1 script destroys the previous realization irrecoverably.

## What is settled that used to be open

- **VR's sign.** It came out over-confident and **worse than MR** (2.18 against 1.97 on the kinetic
  pair, 1.77 against 1.53 on the amplitude pair). The mechanism claim stands and the Results arc no
  longer needs two branches. `decisions.md`.
- **The two anchors agree** at the Fig-2 cell, so moving the main text onto the Gaussian anchor costs
  no result.
- **"MR → VR changes only the variance, VR → IR only the gain" is wrong** and was retired from four
  documents. The predictive variance divides the gain, and MR and IR predict the *same* observable
  variance from the same state, so the whole MR-to-IR difference is the gain.

## Still open, and they are decisions rather than work

- What Figures 6 and 7 are. Proposal on the table: 6 = the R-versus-IR decision map, 7 = the same map
  in bias. The decomposition needs a home.
- Whether the paper stays gated at six figures (`01_writing_plan.md` §0).
- The `Q-n` open questions in `00_plan.md` §8.
