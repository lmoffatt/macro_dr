## noise_units.R — the one definition of the dimensionless instrumental noise.
##
## Sourced by any notebook that puts a noise number in front of a reader. Nothing on disk changes:
## the stored label stays the label, in every file name and in the CSV column. This file exists so
## the factor 10 below is written once instead of being retyped per figure.
##
## THE THREE QUANTITIES, and why they are not the same number.
##
##   Current_Noise   the engine parameter, a white-noise power spectral density in pA^2*s. It reaches
##                   the likelihood only as the variance of one recorded point, e = Current_Noise/Delta
##                   (qmodel.h:3730). No channel count, no conductance, no correlation time in it.
##
##   label           what the dispatchers write. `current_noise = label/1000`, a shell case statement
##                   repeated in six scripts (dispatch_figure_{2,3,3_G,3_LSE,4}.sh,
##                   run_figure_3_G_local.sh). This is the number stored in the CSV column
##                   `noise_in_conductance_tau` and stamped into every file name as `noise_<label>`.
##
##   S_tilde         the dimensionless noise the axis is NAMED after: S_tilde = Current_Noise*k_off/i^2.
##                   S_tilde = 1 means the instrumental noise over one channel time constant equals the
##                   unitary current. THIS is the quantity to show a reader.
##
## The label is TEN TIMES S_tilde. The dispatcher divides by 1000 where the definition asks for
## k_off/i^2 = 100; the 1000 is 1/Delta at the reference sampling interval Delta = 1 ms = 0.1*tau, so
## the sampling interval got folded into a quantity that is supposed to be independent of it. The two
## agree at no operating point; they only look like one number because every run in the project uses
## k_off = 100 s^-1 and i = 1 pA, which is also why it went unnoticed. Full account, with the
## arithmetic and the decision to leave the stored data alone:  ../../NOISE_AXIS_UNITS.md
##
## RULES.
##  - Join, filter and match files on the RAW label. That is what the file names carry; converting
##    before a join silently fails to match.
##  - Convert only where a number is shown: axis titles, tick labels, printed tables, caption prose.
##  - Do not write a converted column back to any CSV under data/ or runs/.

# label -> dimensionless noise.  S_tilde = label/10.
S_TILDE_PER_LABEL <- 0.1

dimensionless_noise <- function(label) label * S_TILDE_PER_LABEL

# The inverse, for placing an externally computed S_tilde on an axis that is drawn in label units
# (figure_6 does this to put real recording configurations on the map).
noise_label <- function(S_tilde) S_tilde / S_TILDE_PER_LABEL

# Axis titles for scales drawn in dimensionless units. The tilde is the notation used in
# NOISE_AXIS_UNITS.md and in the Methods, where it has one meaning: the quantity in the natural
# units of the model, time in tau = 1/k_off and current in the unitary current i. So S_tilde for
# the noise and Delta_tilde = Delta*k_off for the acquisition interval, which the figures used to
# spell out as the product.
NOISE_AXIS_LAB <- expression(paste("instrumental noise  ", widetilde(S), " = S ", k[off], " / ", i^2))

INTERVAL_AXIS_LAB <- expression(widetilde(Delta) == Delta %.% k[off])
