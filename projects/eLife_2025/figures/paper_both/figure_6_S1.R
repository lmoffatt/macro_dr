## Figure 6 -- figure supplement 1: the boundary-conditioned filter's departure, as a law with its spread, in the
## paper's own coordinates (magnitude m and anisotropy a of Eq. magnitude-anisotropy; D_AI reported nowhere).
## One point per cell, all 560 IR cells at the simulation truth; colour = sampling interval (viridis, as Figure 5);
## grey45 reference at the null; the MEDIAN law (solid grey) and the 95% quantile line (dashed "22"), both with the
## slope fixed at -1/5 (-1 for the bias) and the constant set as the median / 95th percentile of v * x^(1/5) over
## the cells above 1.5x the resolution; dotted = the resolution of 10^4 recordings. Bars: 95% bootstrap over the
## recordings. theme_bw 8/7 pt Helvetica through the pdf device (as Figures 6 and 2 S1); 5.6 in.
## Why a median and a quantile and not a max envelope (2026-08-29, Luciano): a max over 560 noisy cells is set by
## one cell and its bootstrap noise, and the top of the cloud flattens toward the corner, so a fixed-slope envelope
## sat x1.7-2.6 above the data there. The quantile line, with the bars, is what a referee can check.
##   (A) a = e^s against x = N_ch S~ sqrt(Delta~).   (B) m = e^ebar, with the lines of |ebar| mirrored.
##   (C) the worst direction, exp(max|log lambda_i|), the factor by which that direction's information is
##       misreported (Figure 4's unit; right axis on the error bar).
##   (D) b'Hb per recording against N_ch S~^(3/4); right axis 1/(b'Hb), recordings to one joint standard error.
## Error bars: (B) the emitted bootstrap band of ebar; (A),(C) the emitted half-widths of sqrt(var log lambda) and
## of log lambda_max centred on the matrix value; (D) +-1.96 SE by the delta method (ir_bias_Hb.R).
## Data: digests of projects/eLife_2025/figures/in_progress/ir_reliability_20260828 (README = canonical).
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork)})
DIG <- "../in_progress/ir_reliability_20260828/digest"
lab7 <- function(v) factor(sprintf("%.2g", v), levels = sprintf("%.2g", sort(unique(v))))
r <- read.csv(file.path(DIG, "ir_pair_boot.csv")) |>
  mutate(St = 0.1 * z, x = nch * St * sqrt(interval), dlab = lab7(interval),
         a = exp(s), a_hw = (sqrt(s2_hi) - sqrt(s2_lo)) / 2, a_lo = exp(pmax(s - a_hw, 0)), a_hi = exp(s + a_hw),
         m = exp(ebar), m_lo = exp(ebar_lo), m_hi = exp(ebar_hi),
         w = exp(wmax), w_hw = (log(lmax_hi) - log(lmax_lo)) / 2, w_lo = exp(pmax(wmax - w_hw, 0)), w_hi = exp(wmax + w_hw))
b <- read.csv(file.path(DIG, "ir_bias_Hb_se.csv")) |>
  mutate(St = 0.1 * z, xb = nch * St^0.75, lo = pmax(q - 1.96 * se_q, 2e-5), hi = q + 1.96 * se_q, dlab = lab7(interval))
FL_s <- max(r$s[r$x > 1e4]); FL_e <- max(abs(r$ebar[r$x > 1e4])); FL_w <- max(r$wmax[r$x > 1e4]); FL_B <- max(b$q[b$xb > 300])
# law constants: v = c * x^slope; c50 = median, c95 = 95th percentile of v * x^-slope over cells above 1.5x the resolution
law <- function(v, x, fl, slope) { k <- (v * x^-slope)[v > 1.5 * fl]; c(c50 = unname(median(k)), c95 = unname(quantile(k, .95)), n = length(k)) }
L_s <- law(r$s, r$x, FL_s, -0.2); L_e <- law(abs(r$ebar), r$x, FL_e, -0.2); L_w <- law(r$wmax, r$x, FL_w, -0.2); L_b <- law(b$q, b$xb, FL_B, -1)
cat(sprintf("constants (median, q95, n): s %.3f %.3f %d | |ebar| %.3f %.3f %d | wmax %.3f %.3f %d | b'Hb %.4f %.4f %d\n",
            L_s[1], L_s[2], L_s[3], L_e[1], L_e[2], L_e[3], L_w[1], L_w[2], L_w[3], L_b[1], L_b[2], L_b[3]))
cat(sprintf("read-out: x=1  a %.2f (q95 %.2f) | worst info %.2f (q95 %.2f; bar %.2f, %.2f) || x=100  a %.2f (%.2f) | worst %.2f (%.2f) || worst cell info %.2f (bar %.2f), a %.2f\n",
            exp(L_s[1]), exp(L_s[2]), exp(L_w[1]), exp(L_w[2]), exp(L_w[1] / 2), exp(L_w[2] / 2),
            exp(L_s[1] * 100^-0.2), exp(L_s[2] * 100^-0.2), exp(L_w[1] * 100^-0.2), exp(L_w[2] * 100^-0.2), max(r$w), sqrt(max(r$w)), max(r$a)))
cat(sprintf("bias: n* = 1/(b'Hb) median %.0f, q95 %.0f  times N_ch S~^(3/4); share of above-floor cells above the q95 lines: s %.3f, w %.3f, b %.3f\n",
            1 / L_b[1], 1 / L_b[2], mean((r$s * r$x^0.2)[r$s > 1.5 * FL_s] > L_s[2]), mean((r$wmax * r$x^0.2)[r$wmax > 1.5 * FL_w] > L_w[2]), mean((b$q * b$xb)[b$q > 1.5 * FL_B] > L_b[2])))

thm <- theme_bw(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7), axis.title = element_text(size = 8),
        legend.position = "right", legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        legend.key.height = unit(9, "pt"), legend.spacing.y = unit(3, "pt"),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(face = "bold", size = 7, hjust = 0),
        plot.margin = margin(2, 2, 2, 2))
XB <- 10^seq(-2, 10, 3); XL <- parse(text = sprintf("10^%d", seq(-2, 10, 3)))
XLAB <- expression(N[ch] * widetilde(S) * sqrt(widetilde(Delta)))
xg <- 10^seq(log10(min(r$x)), log10(max(r$x)), length = 200)
ln <- function(y, lt, x = xg) geom_line(data = data.frame(x = x, y = y), aes(x, y), inherit.aes = FALSE, linetype = lt, linewidth = 0.35, colour = "grey35")
med <- function(y, x = xg) ln(y, "solid", x); q95 <- function(y, x = xg) ln(y, "22", x)
base <- function(p) p + facet_wrap(~ panel) + scale_x_log10(breaks = XB, labels = XL) +
  scale_colour_viridis_d(name = expression(widetilde(Delta)), option = "D", end = 0.92) + thm

pA <- base(ggplot(r |> mutate(panel = "A   anisotropy"), aes(x, a, colour = dlab)) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = "grey45") +
  geom_hline(yintercept = exp(FL_s), linewidth = 0.3, colour = "grey60", linetype = "dotted") +
  med(exp(L_s[1] * xg^-0.2)) + q95(exp(L_s[2] * xg^-0.2)) +
  geom_linerange(aes(ymin = a_lo, ymax = a_hi), linewidth = 0.25, alpha = 0.45) + geom_point(size = 0.55, alpha = 0.9) +
  scale_y_log10(breaks = c(1, 1.1, 1.2, 1.3, 1.4), labels = c("1", "1.1", "1.2", "1.3", "1.4"), limits = c(0.99, 1.45)) +
  labs(x = XLAB, y = expression(a == e^{s})))
pB <- base(ggplot(r |> mutate(panel = "B   magnitude"), aes(x, m, colour = dlab)) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = "grey45") +
  med(exp(L_e[1] * xg^-0.2)) + med(exp(-L_e[1] * xg^-0.2)) + q95(exp(L_e[2] * xg^-0.2)) + q95(exp(-L_e[2] * xg^-0.2)) +
  geom_linerange(aes(ymin = m_lo, ymax = m_hi), linewidth = 0.25, alpha = 0.45) + geom_point(size = 0.55, alpha = 0.9) +
  scale_y_log10(breaks = c(0.8, 0.9, 1, 1.1, 1.2), labels = c("0.8", "0.9", "1", "1.1", "1.2"), limits = c(0.78, 1.28)) +
  labs(x = XLAB, y = expression(m == e^{bar(e)})))
pC <- base(ggplot(r |> mutate(panel = "C   worst direction"), aes(x, w, colour = dlab)) +
  geom_hline(yintercept = 1, linewidth = 0.3, colour = "grey45") +
  geom_hline(yintercept = exp(FL_w), linewidth = 0.3, colour = "grey60", linetype = "dotted") +
  med(exp(L_w[1] * xg^-0.2)) + q95(exp(L_w[2] * xg^-0.2)) +
  geom_linerange(aes(ymin = w_lo, ymax = w_hi), linewidth = 0.25, alpha = 0.45) + geom_point(size = 0.55, alpha = 0.9) +
  scale_y_log10(breaks = c(1, 1.2, 1.5, 2), labels = c("1", "1.2", "1.5", "2"), limits = c(0.99, 2.3),
                sec.axis = sec_axis(~ sqrt(.), name = "on the error bar", breaks = c(1, 1.1, 1.2, 1.3, 1.5), labels = c("1", "1.1", "1.2", "1.3", "1.5"))) +
  labs(x = XLAB, y = "factor the information is off by"))
xgB <- 10^seq(log10(min(b$xb)), log10(max(b$xb)), length = 200)
pD <- ggplot(b |> mutate(panel = "D   bias in the Fisher metric"), aes(xb, q, colour = dlab)) +
  geom_hline(yintercept = FL_B, linewidth = 0.3, colour = "grey60", linetype = "dotted") +
  med(L_b[1] / xgB, xgB) + q95(L_b[2] / xgB, xgB) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.25, alpha = 0.45) + geom_point(size = 0.55, alpha = 0.9) +
  facet_wrap(~ panel) +
  scale_x_log10(breaks = 10^seq(-1, 9, 2), labels = parse(text = sprintf("10^%d", seq(-1, 9, 2)))) +
  scale_y_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1), labels = c("1e-4", "0.001", "0.01", "0.1"), limits = c(5e-5, 0.15),
                sec.axis = sec_axis(~ 1 / ., name = "recordings to one joint standard error", breaks = c(10, 100, 1000, 10000), labels = c("10", "100", "1000", "10000"))) +
  scale_colour_viridis_d(guide = "none", option = "D", end = 0.92) +
  labs(x = expression(N[ch] * widetilde(S)^{3/4}), y = expression(b^{T} * H * b)) + thm
g <- (pA + pB) / (pC + pD) + plot_layout(guides = "collect") & theme(legend.position = "right")
ggsave("Figure_6_S1.pdf", g, width = 5.6, height = 4.7, device = "pdf", family = "Helvetica")   # the pdf device, as Figures 6 and 2 S1; NOT cairo
cat("written Figure_6_S1.pdf (rasterise with pdftoppm to look at it)\n")
