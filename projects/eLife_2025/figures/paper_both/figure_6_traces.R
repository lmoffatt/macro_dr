## Figure 6 -- figure supplement 2 (VARIANT=N5 DHAT=0.01): what the record looks like on both sides of the resolvability limit (2026-08-29).
## Six cells of the sweep's scheme: N_ch = 10 (top) and 100 (bottom), S~ = 0.005, 0.1, 1 (columns), simulated with
## the project's exact simulator (ops/local/figure_6_traces*.macroir, seed 20260829, uniformization). Each panel: the
## true channel current (grey steps, open channels times i, from the substeps) and the sampled record (black, baseline
## removed), in units of the unitary current. Label: N_ch S~ (the noise on a conductance level, units of i^2) and the
## anisotropy a measured at that cell at the same interval. DHAT selects the sampling interval: 0.01 (one sample per
## level at N = 100, the sweep's shortest, 3 tau shown) or 0.1 (the reference interval, 20 tau shown).
## Style: Figures 5 / 6 (theme_bw 8/7 pt Helvetica, grey92 bold strips), pdf device; 5.6 in wide.
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
DHAT <- if (nzchar(Sys.getenv("DHAT"))) as.numeric(Sys.getenv("DHAT")) else 0.01
DAT <- "../data/figure_6_traces"; DIG <- "../in_progress/ir_reliability_20260828/digest"
VARIANT <- Sys.getenv("VARIANT")   # "" = the 2 x (3|4) grid; "N5" = rows N = 5, 10, 100 at S~ = 1e-4, 0.01, 0.1 (all measured but 1e-4)
if (VARIANT == "N5") { NCHS <- c(10, 100); STS <- c(0.0001, 0.01, 0.1) } else {   # the S2 grid: N = 5 dropped (it adds nothing to N = 10), 1 tau shown
  NCHS <- c(10, 100); STS <- if (DHAT < 0.05) c(0.0001, 0.005, 0.1, 1) else c(0.005, 0.1, 1) }   # 1e-4 is BEYOND the sweep (no measured a)
cells <- expand.grid(nch = NCHS, St = STS)
meas <- read.csv(file.path(DIG, "ir_pair_boot.csv")) |> mutate(St = 0.1 * z) |> filter(abs(interval - DHAT) < 1e-9) |> select(nch, St, s, wmax)
C50 <- 0.094; C95 <- 0.118   # median and 95% constants of a = exp(c x^-1/5), x = N S~ sqrt(Dhat), from figure_6_S1.R (printed there)
TAU <- 0.01; T0 <- 20 * DHAT * TAU; WIN <- if (VARIANT == "N5") 1 * TAU else if (DHAT < 0.05) 3 * TAU else 20 * TAU
pre <- if (DHAT < 0.05) "trace_D0.01_" else "trace_"
rd <- function(nch, St) {
  f <- file.path(DAT, sprintf("%snch_%d_S_%s.csv", pre, nch, format(St, scientific = FALSE, drop0trailing = TRUE)))
  d <- read.csv(f, check.names = FALSE)
  samp <- d |> filter(scope == "simulation") |> transmute(t = step_middle - T0, y = -(patch_current - 1))   # baseline 1 pA removed, open channels up
  sub  <- d |> filter(scope == "simulation_sub") |> transmute(t = step_middle - T0, y = -patch_current)
  list(samp = samp |> filter(t >= 0, t <= WIN), sub = sub |> filter(t >= 0, t <= WIN))
}
S <- list(); U <- list(); A <- list()
for (k in seq_len(nrow(cells))) {
  z <- rd(cells$nch[k], cells$St[k]); m <- meas |> filter(nch == cells$nch[k], abs(St - cells$St[k]) < 1e-9)
  S[[k]] <- z$samp |> mutate(nch = cells$nch[k], St = cells$St[k])
  U[[k]] <- z$sub  |> mutate(nch = cells$nch[k], St = cells$St[k])
  A[[k]] <- data.frame(nch = cells$nch[k], St = cells$St[k],
                       lab = if (nrow(m)) sprintf("N[ch]*tilde(S)==%s*','~~a==%.2f", format(cells$nch[k] * cells$St[k]), exp(m$s))
                             else sprintf("N[ch]*tilde(S)==%s*','~~a%%~~%%%.2f~'(law)'", format(cells$nch[k] * cells$St[k]),
                                          exp(C50 * (cells$nch[k] * cells$St[k] * sqrt(DHAT))^-0.2)))   # the median law of Figure 6 S1, extrapolated
}
S <- bind_rows(S); U <- bind_rows(U); A <- bind_rows(A)
# one strip per panel and its OWN y range (facet_wrap): a shared row range let the S~ = 1 noise crush the corner panel
fmtS <- function(x) ifelse(x < 1e-3, "10^-4", format(x, scientific = FALSE, drop0trailing = TRUE))
lv <- with(expand.grid(St = STS, nch = NCHS), sprintf("N[ch]==%d~~tilde(S)==%s", nch, fmtS(St)))
fac <- function(d) d |> mutate(pf = factor(sprintf("N[ch]==%d~~tilde(S)==%s", nch, fmtS(St)), levels = lv))
S <- fac(S); U <- fac(U); A <- fac(A)
ypos <- S |> group_by(pf) |> summarise(ytop = max(y) + 0.06 * diff(range(y)), .groups = "drop"); A <- A |> left_join(ypos, by = "pf")
thm <- theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 7), axis.title = element_text(size = 8),
        strip.background = element_rect(fill = "grey92", colour = NA), strip.text = element_text(size = 7, face = "bold"),
        panel.spacing = unit(0.15, "cm"), plot.margin = margin(2, 2, 2, 2))
xb <- if (VARIANT == "N5") c(0, 0.25, 0.5, 0.75, 1) else if (DHAT < 0.05) c(0, 1, 2, 3) else c(0, 5, 10, 15, 20)
g <- ggplot() +
  geom_line(data = S, aes(t / TAU, y), colour = "grey30", linewidth = 0.3) +
  geom_step(data = U, aes(t / TAU, y), colour = "#D55E00", linewidth = 0.45, alpha = 0.55) +   # the true open-channel count, on top, so the reader can see whether the record resolves its steps
  geom_text(data = A, aes(x = 0.02 * WIN / TAU, y = ytop, label = lab), parse = TRUE, hjust = 0, vjust = 1, size = 2.4) +
  facet_wrap(~ pf, ncol = length(STS), scales = "free_y", labeller = label_parsed) +
  scale_x_continuous(breaks = xb, expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  labs(x = expression("time after the agonist jump,  " * t / tau), y = "open channels  (current / i)") + thm
out <- if (VARIANT == "N5") "Figure_6_S2.pdf" else if (DHAT < 0.05) "Figure_6_traces_D0.01.pdf" else "Figure_6_traces.pdf"
ggsave(out, g, width = if (length(STS) == 4) 7.0 else 5.6, height = 1.6 * length(NCHS), device = "pdf", family = "Helvetica")
cat("written", out, "\n")
