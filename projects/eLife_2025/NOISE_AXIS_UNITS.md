# The noise axis: what the label means, what it was meant to mean, and the factor 10 between them

Written 2026-07-27, after the discrepancy was found while placing real recording configurations on
Figure 6. Nothing in the simulations is wrong. What is wrong is the stated meaning of the axis.

## One line

The sweep label is **ten times** the dimensionless noise it is named after. To get the intended
quantity, multiply the label by 0.1.

## Notation

- `S` is the instrumental noise power spectral density, in pA² per Hz. In the code it is the model
  parameter `Current_Noise`.
- `i` is the unitary current of one open channel, in pA. In the code, `unitary_current`.
- `k_off` is the channel closing rate in s⁻¹, and `τ = 1/k_off` is the channel time constant.
- `Δ` is the duration of one recorded sample, so `Δ = number_of_samples / fs`.
- `S̃` is the dimensionless noise defined below.
- "label" is the value that appears in `axis(name = "noise_in_conductance_tau", labels = [...])`, in
  every output file name, and on the y axis of every figure that sweeps noise.

## What was intended, and it is the right definition

Measure time in units of `τ` and current in units of `i`. Since `S` has units of current² × time,
there is exactly one dimensionless combination available, with no free constant:

    S̃ = S · k_off / i²

`S̃ = 1` means `S = i²·τ`, so the variance of a recorded point of duration `τ` is `S/τ = i²`. In words:
**the noise measured over one channel time constant equals the unitary current.** That is a good
reference point, it is what the axis name promises, and it is what the axis was meant to carry.

## What is implemented

The program never computes `S̃`. The axis label is a string, and the `Current_Noise` handed to the
model is written next to it by hand. In `ops/slurm/dispatch_figure_3_G.sh:209` the pairing is a shell
case statement, with its own comment saying so:

    case "$nnoise" in                    # label -> current_noise (vnoise = label / 1000)
        0.05) vnoise=0.00005;;
        0.1)  vnoise=0.0001;;
        ...

and in the `.macroir` scripts the same pairing appears literally, two lines apart:

    axis_noise    = axis(name= "noise_in_conductance_tau", labels= ["0.1"])
    current_noise = indexed_double_by(axis= axis_noise, values=[0.0001])

So `label = 1000 · S`. And `S` reaches the likelihood only through the variance of one recorded point,

    e = Current_Noise · fs / number_of_samples = S / Δ

at `legacy/qmodel.h:3740` (also 3920 and 4408, and the four sites in `legacy/micro_monoid.h`). Flat,
independent of the channel number, with no correlation time in it despite the axis name.

## Where the factor 10 comes from

Every run in this project uses `k_off = 100 s⁻¹` (from `{"off","Log10",100}`) and `i = 1 pA`. So

    S̃ = S · k_off / i² = (label/1000) · 100 / 1 = label / 10

The hardcoded 1000 is not `k_off/i² = 100`. It is `1/Δ` at the reference sampling interval
`Δ = 0.1·τ = 1 ms`. **The reference sampling interval got folded into the definition of the noise**,
where it does not belong: how much instrumental noise a rig has cannot depend on how often you choose
to sample it, and the sampling interval already has its own axis, `interval_in_tau`.

The arithmetic that shows it, at `label = 1`, so `S = 10⁻³ pA²/Hz`:

| averaging window | variance `S/Δ` | σ | σ / i |
|---|---|---|---|
| `Δ = 0.1·τ = 1 ms` | 1 pA² | 1 pA | **1** |
| `Δ = τ = 10 ms` | 0.1 pA² | 0.316 pA | 0.316 |

The equality `σ = i` holds one decade of averaging early, at `0.1·τ` instead of at `τ`. That factor
`√10` in σ is the factor 10 in the label.

## Why it went unnoticed

At `k_off = 100 s⁻¹` and `i = 1 pA`, the number `10·k_off/i²` equals 1000 exactly. Every run uses those
two values, so the implementation and the name agree numerically at the single operating point, and no
run in the project can tell them apart. The label is right by coincidence rather than by construction.
It also means the agreement breaks silently: change `k_off` or `i` and every label keeps its value while
its meaning moves, because the nondimensionalization is performed by a shell case statement instead of
by the program, which is the only place that knows the model's actual `k_off` and `i`.

## Conversion

    S̃ = 0.1 × label

| label | `S̃` | σ over one `τ`, in units of `i` | note |
|---|---|---|---|
| 0.05 | 0.005 | 0.071 | lowest simulated |
| 0.1 | 0.01 | 0.1 | the reference cell of most runs |
| 0.2 | 0.02 | 0.141 | |
| 0.3 | 0.03 | 0.173 | |
| 0.5 | 0.05 | 0.224 | |
| 1 | 0.1 | 0.316 | |
| 10 | 1 | **1** | this is the intended "noise = 1" |
| 100 | 10 | 3.16 | |
| 1000 | 100 | 10 | |
| 10000 | 1000 | 31.6 | |
| 100000 | 10⁴ | 100 | |
| 1000000 | 10⁵ | 316 | highest simulated |

The row worth remembering: **the sweep's `noise = 10` is the case the definition calls `noise = 1`**,
where the instrumental noise over one channel time constant equals the unitary current.

## The patch: convert on display, never in the data

The label is written into 2049 file names under this project, 1817 of them CSV, and appears in the
contents of 2292 files across `figures/data/` and every directory under `runs/`. Those strings are the
provenance: they are what ties a figure to the sweep cell that produced it. Rewriting them would break
that link and would buy nothing, because the simulations consumed `Current_Noise` and are unaffected by
what the label is called.

So: **leave every stored label alone.** Convert only where a number is shown to a reader.

Places that need the conversion:

- the y axis title and tick labels of every figure that shows the noise sweep
- any absolute noise value quoted in a caption or in prose. For example "0.05 is the lowest noise
  simulated" becomes `S̃ = 0.005`
- `figures/paper_both/figure_6.Rmd` already prints both columns, `label` and `S_tilde`, as of
  2026-07-27, so its render log carries the conversion without needing a rerun

Places that do **not** need anything:

- the recording-configuration boxes on Figure 6. They are computed in label units by
  `lbl(S, tau, i) = 10 * S / (tau * i^2)`, which is exactly `10·S̃`, the same units as the measured
  boundaries. Boxes and boundaries both divide by 10 together.
- consequently, every statement about orderings, region crossings, and the ratio
  `r = label / N_ch`. Those are unchanged, because numerator and denominator of every comparison move
  by the same factor.

## Verifying this in a minute

    grep -n -A16 'label -> current_noise' projects/eLife_2025/ops/slurm/dispatch_figure_3_G.sh
    grep -n -A2 'axis_noise' projects/eLife_2025/ops/local/figure_3.macroir
    grep -n 'Current_Noise>(m)' legacy/qmodel.h
    grep -n '"off","Log10"' projects/eLife_2025/ops/local/figure_3.macroir

The first two show the hand-written pairing, the third shows that `Current_Noise` enters only as
`S·fs/number_of_samples`, and the fourth shows `k_off = 100`.

## Still to decide

Whether the published axis becomes `S̃` with the labels divided by ten, or stays as the label with the
conversion stated once. Either is defensible. What is not defensible is the present state, where the
axis name says `S̃` and the numbers on it are `10·S̃`.

If it is worth closing properly rather than documenting, the fix is to compute `Current_Noise` inside
the program from `S̃` and the model's own `k_off` and `i`, so that the label cannot drift away from its
name again when either of those changes.
