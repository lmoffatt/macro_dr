# Audit of `approach.md`: where each claim comes from

> Written 2026-08-04 at Luciano's request, after he flagged that `approach.md` was being quoted back
> to him as his settled position when parts of it might have been written by an agent and never
> reviewed. **This file makes no new claims about the science.** It classifies what is already in
> `approach.md` by provenance. It does not edit that file.
>
> **The rule, stated by Luciano 2026-08-04:** *"Lo que origina todo son los audios, en cuanto a la
> orientacion."* Orientation comes from him. Measurement comes from data.
>
> **Correction, same day.** A first version of this audit read the `session` label in `approach.md`
> section 9 as "the agent invented it" and listed five beats, plus sections 2 and 6, as unbacked. That
> was wrong, and the error is worth recording because it is the same failure in the other direction:
> declaring something unsourced without opening the source. Luciano then pointed at
> `~/.claude/projects/`, where his typed interventions are stored. Nearly every `session` item is his,
> verbatim. What follows is the corrected accounting.

## 1. The label `session` means "typed conversation with Luciano", not "agent invention"

Timestamps below are UTC as stored; local time is three hours earlier. All from
`71004b7d-881f-4c59-891e-7206fa48c1c6.jsonl`, the session that wrote `approach.md`, whose last
message lands minutes before the file's mtime.

| Item in `approach.md` | Table said | Actually |
|---|---|---|
| §1, the spine: the tool is the algorithm **and the executable code** | `10.11.53 + session` | **His, typed 14:18:42.** *"la herramienta no es solo el algoritmo sino el codigo ejecutable... Eso a un biofisico le soluciona un problema: confiar en la herramienta."* |
| §2 row 2, per-sample INR looks calibrated, the autocorrelation is what shows the problem | `session` | **His, typed 14:49:43.** *"hay un deception ahi tambien: que si miras por sample NIR parece estar todo bien!! tenes que ver por Autocorrelacion para ver las cosas"* |
| §3 **beat 6**, noise equalises everyone; the interval members always win on bias and MacroIR on distortion | `session` | **His, typed 16:54:30**, almost the beat verbatim. *"la figura 4 lo que hace es ilustrar que el ruido nos iguala a todos Y que los Intervals ganan siempre en bias y solo macroir gana siempre en distortion."* |
| §3 beat 2, the non-interval members also fail in bias | `session` | **His, typed 17:00:06.** *"en el beat 2 tambien vemos que los no-intervalos fallan en el bias"* |
| §10, Figure 2 observes the anisotropy and a Figure 4 supplement maps it | `session` | **His, typed 15:10:53.** *"Tenemos que señalar la anisotropia en la figura 2 y luego mostrarla en la 4 suplementaria."* Also 15:07:12 for the 15x / 10x / 8.8x reading. |
| The sample-versus-correlation distinction, which is axis 1 | not a row | **His, typed 17:06:50.** *"Hay que hacer explicita la diferncia entre sample y correlation (habria que evaluar si ordenar las filas de fig 3 para que sea mas explicito)."* |
| §6, the demarcation | no row; "Decided 2026-08-03", nobody named | **His, typed 21:42:14**, verbatim. *"O sea la discusion ambplia de cuando usar macroir no se puede dar por completo con los elementos del paper, hay que simplemente demarcarla."* |
| §3 beat 7, the limits are worth their price | none | **His, typed 21:50:54.** *"ojo, mostrar que macroir falla es importante, no le bajemos el precio"* |
| §7, the gift reduces to one measurement | none | **His, typed 20:15:16.** *"de ultima, todo se reduce a ver si tenes autocorrelacion y ahi no hay nada que estimar."* Plus 20:24:16 for the 1.15 threshold, and 21:12:15 for the caveat that τ_int must be taken at equilibrium or after a fit that tracks the mean well. |
| §8, M must be about 10⁴ | none | **His, typed 21:30:20.** *"con 500 samples no logras nada, tenes que hacer 10.000 samples para tener una estimacion decente del error"* |
| §8, do not commit to where the cost crossover sits | none | **His, typed 21:34:49 and 21:38:08.** *"el costo de macroir pareciera escalar como k^2... o sera k^3?... quizas recien con 10 estados ILSE seria competitivo?"* and *"no conviene casarse con numeros firmes"* |
| §8, where sweeps are plentiful least squares wins | none | **His, typed 21:38:56 and 21:40:46.** *"si LSE ganaria si podes promediar traces, ahi me parece que no hay mucha discusion"*, *"promediar en el sentido de fitearlos individualemente, o en conjuntos muestra boostrapeados como hice en mi paper"* |

He also commissioned the file itself, 19:38:45: *"ahora un documento md con toda nuestra argumentacion
acerca de la narrativa del paper, como implmentamos las ideas de los audios y que dudas quedarian."*
It was never claiming to be audio-only. It says "nuestra argumentación", and that is what it is.

## 2. What is genuinely the agent's, and unreviewed

Two things, and neither was found in any transcript or audio.

1. **§2's organising sentence, "everything local is right and everything non-local is wrong."** The
   four instances under it are his or measured; the sentence that unifies them is the file's. The file
   even records its own two drafts of it, which shows it was being composed rather than transcribed.
2. **§1's cast: "the villain is OPACITY", "the protagonist is the READER", "the paper is the mentor."**

The nine-beat enumeration with act boundaries, a MIDPOINT and a CLIMAX is a middle case. The
emotional shape is his, from the audios of 09.59.31 and 10.04.51, almost verbatim. The numbering is
the file's, but he argued inside it rather than inheriting it, 16:49:44: *"a mi no que no me parece
fuerte es acto II y climax, en cuanto a si la figura 5 es el climax... discutamos ese punto."*

## 3. Measurements checked against the repo

| Claim | Where | Verdict |
|---|---|---|
| `seed = 0` is the `random_device` sentinel | §11 | **VERIFIED.** `legacy/mcmc.h:38-39`: `if (initseed == 0) { std::random_device rd; ...` |
| Two states is the best case: the occupancy is exactly binomial and its propagation exact | §4 | **SUPPORTED.** The covariance update in `legacy/qmodel.h` (~4389) is `AT_B_A(t_P(), p_P_Cov() - diag(p_P_mean())) + diag(p_P_mean() * t_P())`, the exact multinomial propagation. |
| `max_lag` is fixed at 10, so τ_int has a ceiling near 18.9 | §7 | **PARTLY WRONG, two ways.** (i) The citation `legacy/moment_statistics.h:62-70` shows `max_lag` as a *function parameter*; the fixing happens in the configs, so the pointer is at the wrong place. (ii) "every production config sets it" is false. `max_lag = 10` is what the figure-3 scripts use (`1f7138b-dirty/.../script.macroir:203` and siblings), but a grep over `projects/eLife_2025/ops/` also returns `max_lag=100` and `max_lag=15`. The 18.9 ceiling follows from 10 and does not apply uniformly. **Resolve before the Figure 6B caption.** |

**And the ceiling is the wrong worry at the long end (Luciano, 2026-08-04).** The record is 10 τ at
every interval, so at Δt·k_off = 1 it holds ten observations, five on the rise and five on the decay.
A lag window of 10 then spans the whole record and the high-lag ρ̂ are estimated from one pair or
none. τ_int there is not a noisy estimator, it is an estimator without data, and the sample
autocorrelation's negative bias goes as 1/N, so with ten points τ_int is pulled down. §7 already
records that a short record deflates τ_int and that this is the direction that does not err safe; what
is new is that at the long-interval end the deflation is not a small correction.

**A defence, and the single number that decides whether it holds.** τ_int is large when there is
correlation, and there is correlation when the interval is short, which is where the record holds a
thousand points. Where the estimator starves, at long intervals, the true answer is near 1 and the
decision is easy. Required precision and available data move together. **What breaks that defence is
the K table of §7**, K = 3, 6, 13, 27, 42, 93, 160 across the seven intervals. If that list runs
short-to-long, then K = 160 sits at Δt·k_off = 1, meaning 160 sweeps are needed exactly where the
answer is trivially 1, and the table is measuring the starvation rather than a precision requirement.
If it runs long-to-short, the defence stands. **The ordering is not recoverable from `approach.md`
and comes from `wf_c7d9da53-c4b`. Check it before Figure 6B is built**, because it decides whether
the panel needs a masked or truncated long-interval end.

## 4. Measurements whose only source is a prior agent run

`approach.md` names `wf_cd2f528f-2c3` (§5) and `wf_c7d9da53-c4b` (§7). The numbers below rest on one
of those or on an unnamed session computation, and none was re-derived here. Not a claim that they are
wrong; a statement that the backing is an agent's output rather than an inspectable artifact.

- §5: the refutation of the five candidate breakdown criteria; the N_ch exponent of −0.25; the
  one-sided exclusion.
- §7: K = 3, 6, 13, 27, 42, 93, 160; the 67 per cent failure dropping to 3 at K = 25; the
  self-centring bias of 9 per cent at K = 2 and 1 per cent at K = 5.
- §8: the cost ratios (1.012 s vs 0.200; 0.0103 vs 0.0072); the nsim series 1.270 / 1.265 / 1.149 /
  1.127 / 1.064 against a floor of 1.043. Note that Luciano set the 10⁴ requirement by argument
  before these were run.
- §2 and §3: r² 0.9995 and accumulated 20.8 with score autocorrelation 0.868; per-sample 1.06 against
  a total of 1.60; IR's maxima 1.42, 1.077, 1.78 over 294 cells.
- §11: median |bias| 0.0019 / 0.0119 / 0.0119 / 0.0009.

Several have an independent trail in the manuscript's `% src:` comments (`04_results.tex`,
`figure_6.html`, `decisions/D-4_ranking_verdict.md`), produced by a different route. Check those first.

## 5. The one claim with no backing, and it was mine

**That the map shows least squares to be adequate over a large part of routine electrophysiology.**
Not in `approach.md`. Asserted in conversation on 2026-08-04 by sliding "calibrated" into "adequate".
Calibrated means the reported interval matches the delivered spread; it says nothing about how much
information is delivered, and in the high-noise regime where least squares calibrates, the information
matrix is rank deficient and only the product N·i is identified, so it is honest about strictly less.

The least-squares arm is also unfinished: standing blocker 3 has it at n_sims 1000 against 10⁴ for the
likelihood arm; the re-runs at 10⁴ ordered by the scope upgrade of 2026-07-23 are pending; §10 item 2
is not done and §10 item 6 is not started. No extensive study stands behind any statement of
least-squares adequacy.

## 6. Orientation stated by Luciano on 2026-08-04, not yet in `approach.md`

1. **Least-squares confidence intervals are not credible when the residuals are autocorrelated.**
2. **Empirical intervals for least squares are obtainable** from independent samples, by bootstrap.
   His own reservation in the same breath: resampling and taking global least-squares fits from the
   same population may not serve fully. It does not separate stochastic variability from variability
   of the rates, and it does give intervals. Not done in this work; done in the JGP 2007 paper. This
   extends what he already said on 2026-08-03 at 21:40:46.
3. **MacroIR will be equal or more precise, and that may offset the higher per-evaluation cost by
   needing fewer optimisation steps to reach the maximum**, plausibly in complex models with different
   apparent conductances. His own label: speculative, with a plausible construction. Do not report it
   as measured.
4. **MacroIR is a valid option; when each is preferable is not clear; it is one more tool.**
5. **Why the methods were not used:** there were no fast algorithms, in part because what was
   available before the 2025 work does not integrate the current over the acquisition interval. This
   is sharper than the two explanations in `00_abstract.tex` note 5 and is partly checkable against
   the prior art: the fast members (2007 recursive, Stepanyuk 2011) are wrong about the observable,
   and the member that is right about the observable (covariance fitting) cannot run on a whole
   record. It remains a hypothesis about other people's behaviour and should be written as one.

## 7. What to do with `approach.md`

Not done here, pending his word.

1. **Fix the `session` label**, which is what made this audit go wrong the first time. Split it into
   `luciano-typed` with a timestamp, and `agent`. The transcripts under `~/.claude/projects/` make
   this mechanical.
2. Mark the two genuine agent framings of section 2 above in place, without deleting them.
3. Amend the header: recording where the session corrected an audio is half the job; it also has to
   record what the session originated.
4. Fold section 6 of this audit into `approach.md`.
