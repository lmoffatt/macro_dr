# D-2: units of the six `scheme_CO` parameters

> Prefilled for Luciano. Task B-2. His action: correct the wrong cells (~15 min).
>
> The units are written **nowhere** in the repo. Every "proposed unit" below is inferred from the
> arithmetic that consumes the parameter, not read from a declaration. The confidence column says how
> firm each inference is. Where it says GUESS, the row is phrased as a question at the bottom.
>
> Two facts that frame the whole table:
>
> 1. **The CSVs carry both the natural value and its log10.** `scheme_CO_par.csv` /
>    `scheme_CO_prior.csv` have a `parameter_value` column (the natural value) **and** a
>    `transformed_mean` column (its base-10 logarithm), with `parameter_transformation = "Log10"`
>    (`scheme_CO_par.csv:2-7`). The model's *working* parameter is the log10 (`create_parameters(...,
>    "Log10", value)` then `get_standard_parameter_transformed_values`), but both forms sit in the
>    file, so "natural or log10" is answered per-column, not per-file.
> 2. **The CSV values are not the figure values.** `scheme_CO_par.csv` holds `on = 0.1`,
>    `Current_Baseline = 0`; every production `.macroir` overrides these inline to `on = 10`,
>    `Current_Baseline = 1` and never loads the CSV (`figure_3_mle.macroir:45-50`;
>    `figure_provenance.md` §3). Units are identical either way. The table reports the CSV column
>    names as asked, and flags the value mismatch where it matters.
>
> **The time and concentration scales that fix every rate/current/PSD unit are inferred, not
> declared.** Time is in **seconds** (the per-sample step is `dt = number_of_samples / fs` with
> `fs = 50e3` Hz, `qmodel.h:3094`, `figure_2.macroir:75`; and `off = 100` gives
> `1/off = 10 ms = τ`, the acquisition axis `interval_in_tau` label "1" being 500 samples =
> 10 ms, consistent, `methods_plan.md:55`). Concentration is in **micromolar (µM)** (the only place
> it is written: the agonist step `10.0` is annotated `@ 10 µM`, `figure_3_mle.macroir:53`,
> `figure_3_mle_G.macroir:53`). Current has **no unit anywhere**; pA is a physical guess from
> single-channel scale.

| Parameter | Column name in the CSVs | Natural or log10 | Proposed unit | Confidence | Evidence (`file:line`) |
|---|---|---|---|---|---|
| forward binding rate `k_on` | `on` (`parameter_value` = 0.1 natural; `transformed_mean` = −1 log10) | both columns present; working value is log10 | **µM⁻¹·s⁻¹** (per-concentration per-second rate constant) | INFERRED (structure from arithmetic; concentration = µM stated in a comment; time = s inferred) | `models_simple.h:22` (`Qa` holds `on`); `qmodel.h:998` (`Qx = Q0 + Qa·x`, so `on·[agonist]` is a rate); `qmodel.h:3094` + `figure_2.macroir:75` (`dt = n_samples/fs`, `fs = 50e3` Hz → rates in s⁻¹); `figure_3_mle.macroir:53` (`@ 10 µM`) |
| unbinding / closing rate `k_off` | `off` (`parameter_value` = 100 natural; `transformed_mean` = 2 log10) | both columns present; working value is log10 | **s⁻¹** (first-order rate constant) | INFERRED (rate is firm; time = s inferred, corroborated by τ = 1/off = 10 ms) | `models_simple.h:19,54` (`Q0` holds `off` directly, no concentration factor); `qmodel.h:998` (enters the generator `Qx`); `qmodel.h:3094` (`dt` in seconds); `methods_plan.md:55` (`τ = 10 ms = 1/k_off`) |
| unitary current | `unitary_current` (`parameter_value` = 1 natural; `transformed_mean` = 0 log10) | both columns present; working value is log10 | **pA** (single-channel current; sign-flipped to inward/negative) | INFERRED that it is a *current*; **pA specifically is a GUESS** | `models_simple.h:24,40,59` (`g = unitary_current`, `× −1.0` sign flip); `qmodel.h:3735` (`y_mean = N·(P·g) + baseline`, so `g` is current per channel and `y` the measured current) |
| current noise (white) | `Current_Noise` (`parameter_value` = 0.001 natural; `transformed_mean` = −3 log10; swept in figures) | both columns present; working value is log10 | **pA²·s (= pA²/Hz)**, a white-noise power spectral density | INFERRED that it is `current²·time` (a PSD); **absolute scale rides on the pA GUESS** | `qmodel.h:3730` (`e = Current_Noise·fs / number_of_samples`, so `Current_Noise = e·n_samples/fs` = variance × time); `qmodel.h:3742-3744` (`e` is the current *variance* `y_var`, in current²) |
| current baseline | `Current_Baseline` (`parameter_value` = 0 natural, `transformed_mean` = −inf; figures use 1) | both columns present; working value is log10 | **pA** (additive current offset) | INFERRED that it is a *current*; **pA specifically is a GUESS** | `qmodel.h:3729,3735` (`y_baseline` added to `y_mean`, same dimension as the current) |
| mean number of channels | `Num_ch_mean` (`parameter_value` = 5000 natural; `transformed_mean` = 3.6989 log10; swept in figures) | both columns present; working value is log10 | **dimensionless (count of channels)** | INFERRED, near-STATED (name + arithmetic agree) | `models_simple.h:60` (`N_Ch_mean`); `qmodel.h:3735` (`N` multiplies `P·g` to scale current, and `N·gSg` to scale variance, i.e. a channel count) |

## The swept "noise" axis is a DIMENSIONLESS unit (corrected 2026-07-15, with Luciano)

B-2's first pass called `noise_in_conductance_tau` a misnomer. **That was wrong.** The name is
correct: the swept axis is a *dimensionless* noise, the natural (nondimensional) unit Luciano
designed, and the `/1000` is the nondimensionalization, not a meaningless rescale.

- The engine parameter is `Current_Noise`, a white-noise **power spectral density** in pA²·s. It
  enters `e = Current_Noise·fs/n_samples = Current_Noise/dt` (`qmodel.h:3730`), the white variance over
  an interval.
- The **figure axis label** is the dimensionless quantity, mapped as `Current_Noise = label/1000`
  (`dispatch_figure_3_G.sh:143-152`, comment `vnoise = label/1000`). So the production **"noise 0.1"**
  cell is `Current_Noise = 1e-4` (the label, not the value, is on every axis and filename).
- **What the label means.** The white variance over one relaxation time τ = 1/k_off is
  `Current_Noise·k_off`; against the single-channel signal power g² (g = unitary_current), the
  dimensionless noise is `ν = Current_Noise·k_off/g²` — precisely "noise in conductance·τ" (noise PSD /
  conductance² × 1/τ). With the real values (k_off = 100, g = 1) and `Current_Noise = label/1000`, the
  label equals **10·ν**: the noise-equals-single-channel-signal crossover (ν = 1, σ_white(τ) = g) sits
  at **label = 10** (`Current_Noise = 1e-2`), and `interval_in_tau = 1` really is dt = 1/k_off (the
  experiment uses n_samp = 500 at fs = 50 kHz, `dispatch_figure_3_G.sh:169-171`, confirming τ = 1/k_off,
  the "lazy k_off" instead of the true eigenvalue λ = k_on·[A]+k_off = 200).
- **Two loosenesses in the name, not errors:** it drops the square (it is g², a *variance* ratio), and
  the label carries an extra decade (label = 10ν; a "clean" ν = label would need `/100`, not `/1000`).
  The extra decade is a design choice (label = 1 as a low-noise reference, crossover at 10), **not** the
  k_off-vs-λ laziness (that is only a factor 2 here).

Consequence for the Methods and captions: the swept axis is **dimensionless** (`noise_in_conductance_tau`,
= 10·Current_Noise·k_off/g²), *not* a misnomer; report it as such. The underlying engine parameter is
`Current_Noise` in pA²·s (value e.g. `1e-4`); the figure shows the dimensionless label = 1000·Current_Noise.
Full derivation: memory `project_design_axes_nondimensional`.

## Questions for Luciano (only the GUESS / unconfirmed parts)

1. **Current unit (the one real GUESS, and it propagates).** The code carries no current unit. Is the
   unitary current **1 pA**? If yes, `unitary_current` and `Current_Baseline` are in **pA**, and
   `Current_Noise` is in **pA²/Hz** (its value `1e-4` meaning `1e-4 pA²/Hz`). If the intended unit is
   nA (or A), all three rescale together. This is the single cell most likely to be wrong.

2. **Time unit = seconds?** Everything treats rates as s⁻¹ because `fs = 50e3` is read as 50 kHz and
   `1/off = 10 ms` matches the `interval_in_tau` axis. Confirm rates are **s⁻¹** (so `off = 100` is
   100 s⁻¹ and `on = 10` is 10 µM⁻¹·s⁻¹). If `fs` were ever meant as kHz-in-ms or similar, the rates
   shift by a decade.

3. **Concentration unit = µM?** The only statement of it is the comment `@ 10 µM`
   (`figure_3_mle.macroir:53`). If the agonist `10.0` is really mM (or M), `on`'s unit changes from
   µM⁻¹·s⁻¹ accordingly. Confirm µM.

4. **(Not a unit, but blocks the Methods value.)** `scheme_CO_par.csv` says `on = 0.1` and
   `Current_Baseline = 0`; the figures actually ran `on = 10`, `Current_Baseline = 1`
   (`figure_3_mle.macroir:45-50`). Confirm the figure values are the ones to report, and whether the
   stale CSV should be deleted before the repo is carved (`09_carve_plan.md`).
