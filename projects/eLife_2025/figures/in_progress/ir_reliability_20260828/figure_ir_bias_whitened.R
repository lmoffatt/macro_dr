# Whitened bias v = G^{1/2} b given the log-lambda treatment: components per parameter, mean, sd. All vs x, colour = Dhat.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
PN <- c("k_on","k_off","i","S","baseline","N_ch"); PLAB <- c(k_on="k[on]", k_off="k[off]", i="i", S="S", baseline="baseline", N_ch="N[ch]")
w <- read.csv("digest/ir_bias_whitened.csv") |> mutate(St=0.1*z, x=nch*St*sqrt(interval),
       dlab=factor(sprintf("%.2g", interval), levels=sprintf("%.2g", sort(unique(interval)))),
       param=factor(param, levels=PN))
cell <- w |> group_by(nch,z,interval,x,dlab) |>
  summarise(vbar=mean(v), sv=sqrt(mean((v-mean(v))^2)), norm=sqrt(sum(v^2)), .groups="drop")
cat("check |v|^2 = 6(vbar^2+sv^2): max rel err", signif(max(abs(cell$norm^2 - 6*(cell$vbar^2+cell$sv^2))/cell$norm^2),3), "\n")
cat("vbar quantiles:\n"); print(round(quantile(cell$vbar, c(0,.05,.5,.95,1)),4))
cat("sv   quantiles:\n"); print(round(quantile(cell$sv, c(0,.05,.5,.95,1)),4))
cat("subordination |vbar| <= sv in", sum(abs(cell$vbar)<=cell$sv), "of", nrow(cell), "cells\n")
cat("\nsign of whitened component by parameter, cells with x<1 (big-bias corner): share positive\n")
print(w |> filter(x<1) |> group_by(param) |> summarise(pos=round(mean(v>0),2), med=round(median(v),3), q95abs=round(quantile(abs(v),.95),3), n=n()) |> as.data.frame(), row.names=FALSE)
cat("\nwhich parameter carries the largest |v| (all cells):\n"); print(table(w |> group_by(nch,z,interval) |> slice_max(abs(v), n=1, with_ties=FALSE) |> pull(param)))
s <- cell |> filter(sv > 0.03)
m1 <- lm(log10(sv) ~ log10(x), data=s); m2 <- lm(log10(sv) ~ log10(nch), data=s)
cat(sprintf("\nsv: vs x R2 %.3f slope %+.3f | vs N R2 %.3f slope %+.3f\n", summary(m1)$r.squared, coef(m1)[2], summary(m2)$r.squared, coef(m2)[2]))

COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(w$dlab)), levels(w$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") +
  theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank(),
        strip.background=element_blank(), strip.text=element_text(size=8))
XBRK <- 10^seq(-2,10,4); XLABS <- parse(text=sprintf("10^%d", seq(-2,10,4)))
XL <- expression(N[ch]%.%tilde(S)%.%hat(Delta)^0.5)
levels(w$param) <- PLAB[levels(w$param)]

pA <- ggplot(w, aes(x, v, colour=dlab)) +
  geom_hline(yintercept=0, colour="grey40", linewidth=.3) +
  geom_point(size=.6, alpha=.6) + facet_wrap(~param, nrow=1, labeller=label_parsed) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_colour_manual(values=COLd, name=expression(hat(Delta))) +
  coord_cartesian(ylim=c(-0.3,0.3)) +
  labs(title="A   whitened bias by parameter,  v = G^{1/2} b", x=XL, y=expression(v[p])) + thm
pB <- ggplot(cell, aes(x, vbar, colour=dlab)) +
  geom_hline(yintercept=0, colour="grey40", linewidth=.3) +
  geom_point(size=.9, alpha=.75) + scale_x_log10(breaks=XBRK, labels=XLABS) +
  scale_colour_manual(values=COLd, guide="none") + coord_cartesian(ylim=c(-0.16,0.16)) +
  labs(title="B   net (mean of components)", x=XL, y=expression(bar(v)==mean~v[p])) + thm
pC <- ggplot(cell, aes(x, sv, colour=dlab)) +
  geom_point(size=.9, alpha=.75) + scale_x_log10(breaks=XBRK, labels=XLABS) +
  scale_y_log10(breaks=c(0.003,0.01,0.03,0.1,0.3), labels=c("0.003","0.01","0.03","0.1","0.3")) +
  scale_colour_manual(values=COLd, guide="none") +
  labs(title="C   anisotropy (sd of components)", x=XL, y=expression(s[v]==sd~v[p])) + thm
g <- pA / (pB + pC) + plot_layout(heights=c(1,1.1), guides="collect") & theme(legend.position="right")
ggsave("figures/Figure_IR_bias_whitened.pdf", g, width=7.0, height=4.6)
ggsave("figures/Figure_IR_bias_whitened.png", g, width=7.0, height=4.6, dpi=170)
cat("written figures/Figure_IR_bias_whitened.{pdf,png}\n")
