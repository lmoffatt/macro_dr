# Anatomy of ebar: residual witnesses (size, kurtosis, autocorrelation) and per-parameter J/G split into per-step shape vs temporal correlation.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
PN <- c("k_on","k_off","i","S","baseline","N_ch"); PLAB <- c(k_on="k[on]", k_off="k[off]", i="i", S="S", baseline="baseline", N_ch="N[ch]")
cw <- read.csv("digest/ir_residual_shape.csv"); pp <- read.csv("digest/ir_score_JG_perparam.csv")
d  <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |> select(nch,z,interval,ebar)
mk <- function(df) df |> mutate(St=0.1*z, x=nch*St*sqrt(interval), dlab=factor(sprintf("%.2g", interval), levels=sprintf("%.2g", sort(unique(interval)))))
cw <- mk(cw |> inner_join(d, by=c("nch","z","interval")))
pp <- mk(pp) |> mutate(shape=log(Jstep/G), tcorr=log(J/Jstep), param=factor(param, levels=PN)); levels(pp$param) <- PLAB[levels(pp$param)]
COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(cw$dlab)), levels(cw$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") + theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"),
  panel.grid.minor=element_blank(), strip.background=element_blank(), strip.text=element_text(size=8))
XBRK <- 10^seq(-2,10,4); XLABS <- parse(text=sprintf("10^%d", seq(-2,10,4))); XL <- expression(N[ch]%.%tilde(S)%.%hat(Delta)^0.5)
wit <- function(y, ttl, yl, lim) ggplot(cw, aes(x, .data[[y]], colour=dlab)) + geom_hline(yintercept=0, colour="grey40", linewidth=.3) +
  geom_point(size=.8, alpha=.75) + scale_x_log10(breaks=XBRK, labels=XLABS) + scale_colour_manual(values=COLd, name=expression(hat(Delta))) +
  coord_cartesian(ylim=lim) + labs(title=ttl, x=XL, y=yl) + thm
pA <- wit("r2bar", "A   size: residual variance misfit", expression(bar(r)[std]^2 - 1), c(-0.1,0.1))
pB <- wit("kurt",  "B   shape: residual excess kurtosis", expression(kappa[r]), c(-0.6,0.6))
pC <- wit("acf1",  "C   memory: residual autocorrelation term", expression(Var(sum(r))/sum(Var(r[t])) - 1), c(-0.1,0.25))
fac <- function(y, ttl, yl) ggplot(pp, aes(x, .data[[y]], colour=dlab)) + geom_hline(yintercept=0, colour="grey40", linewidth=.3) +
  geom_point(size=.6, alpha=.6) + facet_wrap(~param, nrow=1, labeller=label_parsed) + scale_x_log10(breaks=XBRK, labels=XLABS) +
  scale_colour_manual(values=COLd, guide="none") + coord_cartesian(ylim=c(-0.15,0.2)) + labs(title=ttl, x=XL, y=yl) + thm
pD <- fac("shape", "D   score variance vs Fisher, per step (shape of one sample):  log(J_step/G)", expression(log(J[step]/G)))
pE <- fac("tcorr", "E   score variance, total vs per step (temporal correlation):  log(J/J_step)", expression(log(J/J[step])))
g <- (pA + pB + pC + plot_layout(guides="collect")) / pD / pE + plot_layout(heights=c(1,1,1)) & theme(legend.position="right")
ggsave("figures/Figure_IR_ebar_anatomy.pdf", g, width=7.4, height=7.2); ggsave("figures/Figure_IR_ebar_anatomy.png", g, width=7.4, height=7.2, dpi=150)
cat("written figures/Figure_IR_ebar_anatomy.{pdf,png}\n")
