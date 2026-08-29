suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
w <- read.csv("digest/ir_bias_whitened.csv")
cell <- w |> group_by(nch,z,interval) |> summarise(vsum=sum(v), vbar=mean(v), sv=sqrt(mean((v-mean(v))^2)), mahal=sqrt(sum(v^2)), .groups="drop")
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |> select(nch,z,interval,aniso,ebar,dAI)
cell <- inner_join(cell, d, by=c("nch","z","interval")) |> mutate(St=0.1*z, x=nch*St*sqrt(interval),
  dlab=factor(sprintf("%.2g", interval), levels=sprintf("%.2g", sort(unique(interval)))))
cat("cor(ebar, sum v_p) =", round(cor(cell$ebar, cell$vsum),3), "  cor(ebar, vbar) same; Spearman", round(cor(cell$ebar, cell$vsum, method="spearman"),3), "\n")
cat("cor(s, mahal) =", round(cor(cell$aniso, cell$mahal),3), "  on log scale", round(cor(log(cell$aniso), log(cell$mahal)),3), "\n")
cat("cor(dAI, mahal) log =", round(cor(log(cell$dAI), log(cell$mahal)),3), "\n")
big <- cell |> filter(abs(ebar) > 0.03)
cat("cells |ebar|>0.03: n=", nrow(big), " cor(ebar, sum v)=", round(cor(big$ebar, big$vsum),3), " ratio sum v / ebar quantiles:\n"); print(round(quantile(big$vsum/big$ebar, c(.1,.25,.5,.75,.9)),2))
big2 <- cell |> filter(aniso>0.06, mahal>0.06)
cat("ratio mahal/s (both>0.06) quantiles:"); print(round(quantile(big2$mahal/big2$aniso, c(.1,.25,.5,.75,.9)),2))
cat("sign agreement ebar vs sum v (|ebar|>0.03):", round(mean(sign(big$ebar)==sign(big$vsum)),2), "\n")

COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(cell$dlab)), levels(cell$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") + theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank())
XBRK <- 10^seq(-2,10,4); XLABS <- parse(text=sprintf("10^%d", seq(-2,10,4))); XL <- expression(N[ch]%.%tilde(S)%.%hat(Delta)^0.5)
YB <- c(0.01,0.03,0.1,0.3); YL <- c(0.008, 0.35)
p1 <- ggplot(cell, aes(x, aniso, colour=dlab)) + geom_point(size=.9, alpha=.75) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_y_log10(breaks=YB, limits=YL) + scale_colour_manual(values=COLd, name=expression(hat(Delta))) +
  labs(title="A   distortion:  s = sd log lambda", x=XL, y="s") + thm
p2 <- ggplot(cell, aes(x, mahal, colour=dlab)) + geom_point(size=.9, alpha=.75) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_y_log10(breaks=YB, limits=YL) + scale_colour_manual(values=COLd, guide="none") +
  labs(title="B   bias:  |v| = sqrt(b'Gb)", x=XL, y="|v|") + thm
p3 <- ggplot(cell, aes(ebar, vsum, colour=dlab)) + geom_abline(slope=1, intercept=0, colour="grey40", linewidth=.3) +
  geom_hline(yintercept=0, colour="grey80", linewidth=.3) + geom_vline(xintercept=0, colour="grey80", linewidth=.3) +
  geom_point(size=.9, alpha=.75) + scale_colour_manual(values=COLd, guide="none") + coord_equal(xlim=c(-0.1,0.16), ylim=c(-0.1,0.16)) +
  labs(title=sprintf("C   ebar vs sum v_p   (r = %.2f)", cor(cell$ebar, cell$vsum)), x=expression(bar(e)==mean~log~lambda), y=expression(sum(v[p]))) + thm
p4 <- ggplot(cell, aes(aniso, mahal, colour=dlab)) + geom_abline(slope=1, intercept=0, colour="grey40", linewidth=.3) +
  geom_point(size=.9, alpha=.75) + scale_colour_manual(values=COLd, guide="none") +
  scale_x_log10(breaks=YB, limits=YL) + scale_y_log10(breaks=YB, limits=YL) + coord_equal() +
  labs(title=sprintf("D   s vs |v|   (r = %.2f, log)", cor(log(cell$aniso), log(cell$mahal))), x="s", y="|v|") + thm
g <- (p1 + p2) / (p3 + p4) + plot_layout(guides="collect") & theme(legend.position="right")
ggsave("figures/Figure_IR_bias_vs_distortion.pdf", g, width=6.6, height=5.6)
ggsave("figures/Figure_IR_bias_vs_distortion.png", g, width=6.6, height=5.6, dpi=170)
cat("written figures/Figure_IR_bias_vs_distortion.{pdf,png}\n")
