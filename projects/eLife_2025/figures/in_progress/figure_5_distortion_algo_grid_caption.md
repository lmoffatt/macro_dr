# Figure 5 — captions

## Figure 5. Information distortion across the five likelihood approximations.

Filled contours show the diagonal information distortion (the ratio of the sampling
variance of the score to the model Fisher information, so that 1, shown white, is exact
self-consistency) for each algorithm (columns: NR, MNR, R, MR, IR) over the design plane of
measurement interval $\Delta/\tau_{\min}$ (horizontal) and channel number $N_{\mathrm{ch}}$
(vertical). Rows pair two parameters ($\log_{10} k_{on}$ and $\log_{10} N_{\mathrm{ch}}$)
with three instrument-noise levels. The colour scale is symmetric in $\log$ about 1 and
shared across all 30 panels. Solid lines mark distortion 1.15 and 1/1.15 (the $\pm 15\%$
self-consistency band, about 93% coverage of a nominal 95% interval); dashed lines mark 1.5
and 1/1.5 (about 89% coverage), which macro-R stays inside and macro-MR crosses. Distortion is
evaluated at the pooled estimate $\theta_{\mathrm{pool}}$, with
$N_{\mathrm{ch}} \ge 10$. IR is self-consistent across essentially the whole plane, and its
residual departures are toward under-confidence (the conservative direction) at few channels
and low noise, disappearing as noise rises. R and MR are systematically over-confident, and
the non-recursive NR and MNR are over-confident by one to three orders of magnitude.

## Figure 5B. Self-consistency versus instrument noise.

Fraction of the design plane inside the $\pm 15\%$ self-consistency band
($1/1.15 < \mathrm{distortion} < 1.15$) as a function of instrument noise, per algorithm, for
the two parameters. IR self-corrects (0.88 to 1.00) as noise rises; R and MR start well below
and remain systematically over-confident; NR and MNR stay near zero.

---

### Notes (not for publication, verify before submitting)

- Confirm the statistic definition matches the code: written here as "ratio of the sampling
  variance of the score to the model Fisher information". Adjust if the distortion is defined
  the other way up or with a different normalisation.
- The "$\pm 15\%$ band, about 93% coverage" assumes the distortion is a variance ratio
  (standard error scales as its square root). If the code already reports it on a
  standard-error scale, the coverage figure changes.
- Source files: `Figure_5_distortion_algo_grid.pdf`, `Figure_5_distortion_vs_noise_summary.pdf`
  (both from `figure_5_distortion_algo_grid.Rmd`, data `433ed13`, battery_pool).
