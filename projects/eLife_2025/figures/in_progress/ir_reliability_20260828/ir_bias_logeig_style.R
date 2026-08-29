setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
# The bias treated EXACTLY like the log-eigenvalues of the distortion:
# per cell, 6 components z_p = bias_p/SE_p -> panel A: mean (signed), panel B: sd.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
b <- read.csv("digest/ir_bias_cis.csv") |> mutate(zr = m/se, S=0.1*z)
cell <- b |> group_by(nch, z, interval_in_tau) |>
  summarise(ebar_b = mean(zr), s_b = sd(zr), .groups="drop") |>
  rename(interval=interval_in_tau) |> mutate(S=0.1*z, x1=S*nch,
  dlab=factor(sprintf("%.2g",interval), levels=sprintf("%.2g",sort(unique(interval)))))
cat("ebar_b quantiles:\n"); print(round(quantile(cell$ebar_b, c(0,.05,.5,.95,1)),4))
cat("s_b quantiles:\n");    print(round(quantile(cell$s_b,   c(0,.05,.5,.95,1)),4))
# envelopes, exponent -1/3 (the bias pair-product law), floors from the high-x plateau
sA <- cell |> filter(abs(ebar_b) > 0.012)
aA <- max(abs(sA$ebar_b) * sA$x1^(1/3))
sB <- cell |> filter(s_b > 0.025)
aB <- max(sB$s_b * sB$x1^(1/3))
EAb <- function(x) pmax(aA*x^(-1/3), 0.012)
EBb <- function(x) pmax(aB*x^(-1/3), 0.025)
cat(sprintf("envelope constants: amplitude %.3f  anisotropy %.3f\n", aA, aB))
cat("violations A:", sum(abs(cell$ebar_b)>EAb(cell$x1)), " B:", sum(cell$s_b>EBb(cell$x1)), "of", nrow(cell), "\n")

COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(cell$dlab)), levels(cell$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") +
  theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.32,"cm"),
        panel.grid.minor=element_blank())
xg <- 10^seq(log10(min(cell$x1)), log10(max(cell$x1)), length=200)
XBRK <- 10^seq(-2,10,2); XLABS <- parse(text=sprintf("10^%d",seq(-2,10,2)))
XL <- expression(tilde(S)%.%N[ch])

pA <- ggplot(cell, aes(x1, ebar_b, colour=dlab)) +
  geom_hline(yintercept=0, colour="grey40", linewidth=.3) +
  geom_line(data=tibble(x1=xg,y= EAb(xg)),aes(x1,y),inherit.aes=FALSE,linetype="dashed",colour="grey30",linewidth=.4) +
  geom_line(data=tibble(x1=xg,y=-EAb(xg)),aes(x1,y),inherit.aes=FALSE,linetype="dashed",colour="grey30",linewidth=.4) +
  geom_point(size=.9, alpha=.75) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_colour_manual(values=COLd, guide="none") +
  labs(title="A   bias net amplitude", x=XL, y=expression(bar(e)[b]==mean~over~p~(bias[p]/SE[p]))) + thm

pB <- ggplot(cell, aes(x1, s_b, colour=dlab)) +
  geom_line(data=tibble(x1=xg,y=EBb(xg)),aes(x1,y),inherit.aes=FALSE,linetype="dashed",colour="grey30",linewidth=.4) +
  geom_hline(yintercept=0.025, colour="grey70", linetype="dotted", linewidth=.35) +
  geom_point(size=.9, alpha=.75) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_y_log10() +
  scale_colour_manual(values=COLd, name=expression(Delta%.%k[off])) +
  annotate("text", x=1e10, y=0.026, label="envelope floor (2.5x ensemble noise)", size=2.2, colour="grey45", hjust=1, vjust=0) +
  annotate("text", x=3e-2, y=max(cell$s_b), label="slope -1/3", size=2.4, colour="grey30", hjust=0) +
  labs(title="B   bias anisotropy", x=XL, y=expression(s[b]==sd~over~p~(bias[p]/SE[p]))) + thm

g <- pA + pB + plot_layout(widths=c(1,1.12))
ggsave("figures/Figure_IR_bias_logeig_style.pdf", g, width=6.2, height=2.7)
ggsave("figures/Figure_IR_bias_logeig_style.png", g, width=6.2, height=2.7, dpi=170)
cat("written figures/Figure_IR_bias_logeig_style.{pdf,png}\n")
