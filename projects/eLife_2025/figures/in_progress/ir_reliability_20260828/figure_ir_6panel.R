# A ebar vs x | B s vs x | C sum v vs x | D |v| vs x | E sum v vs N*S | F |v| vs N*S     (x = N*S*sqrt(Dhat), colour = Dhat)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
w <- read.csv("digest/ir_bias_whitened.csv")
cell <- w |> group_by(nch,z,interval) |> summarise(vsum=sum(v), mahal=sqrt(sum(v^2)), .groups="drop")
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |> select(nch,z,interval,aniso,ebar)
cell <- inner_join(cell, d, by=c("nch","z","interval")) |> mutate(St=0.1*z, x=nch*St*sqrt(interval), x0=nch*St,
  dlab=factor(sprintf("%.2g", interval), levels=sprintf("%.2g", sort(unique(interval)))))
FLOOR_V <- sqrt(6/10000)   # whitened-mean-score noise with 1e4 sims: E|chi_6|/100 ~ 0.0235
# trends above floor, for the caption
fx <- function(y, xs, lab){ m <- lm(log10(y) ~ log10(xs)); r <- 10^abs(resid(m))
  cat(sprintf("  %-22s slope %+.3f  R2 %.3f  spread x%.2f  n=%d\n", lab, coef(m)[2], summary(m)$r.squared, median(r), length(y))) }
s <- cell |> filter(mahal > 0.06)
cat("|v| above floor (>0.06):\n"); fx(s$mahal, s$x, "vs x"); fx(s$mahal, s$x0, "vs N*S"); fx(s$mahal, s$nch, "vs N")
sa <- cell |> filter(aniso > 0.06)
cat("s above floor (>0.06):\n"); fx(sa$aniso, sa$x, "vs x"); fx(sa$aniso, sa$x0, "vs N*S")
# envelope for F: |v| <= max(c (N S)^-b, floor) with b from the fit, c = min with zero violations above floor
mF <- lm(log10(mahal) ~ log10(x0), data=s); bF <- -coef(mF)[2]; cF <- ceiling(100*max(s$mahal*s$x0^bF))/100
cat(sprintf("F envelope: |v| <= max(%.2f (N S)^-%.2f, 0.055): violations %d\n", cF, bF, sum(cell$mahal > pmax(cF*cell$x0^-bF, 0.055))))
mD <- lm(log10(mahal) ~ log10(x), data=s); bD <- -coef(mD)[2]; cD <- ceiling(100*max(s$mahal*s$x^bD))/100
cat(sprintf("D envelope: |v| <= max(%.2f x^-%.2f, 0.055): violations %d\n", cD, bD, sum(cell$mahal > pmax(cD*cell$x^-bD, 0.055))))

COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(cell$dlab)), levels(cell$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") + theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank())
XBRK <- 10^seq(-2,10,2); XLABS <- parse(text=sprintf("10^%d", seq(-2,10,2)))
XL <- expression(N[ch]%.%tilde(S)%.%hat(Delta)^0.5); XL0 <- expression(N[ch]%.%tilde(S))
xg <- 10^seq(log10(min(cell$x)), log10(max(cell$x)), length=200); xg0 <- 10^seq(log10(min(cell$x0)), log10(max(cell$x0)), length=200)
EA <- function(x) pmax(0.05*x^-0.2, 0.02); EB <- function(x) pmax(0.21*x^-0.2, 0.045)
YB <- c(0.01,0.03,0.1,0.3); YL <- c(0.008,0.35)
env <- function(df, xv) geom_line(data=df, aes(x=.data[[xv]], y=y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4)

pA <- ggplot(cell, aes(x, ebar, colour=dlab)) + geom_hline(yintercept=0, colour="grey40", linewidth=.3) +
  env(tibble(x=xg, y=EA(xg)), "x") + env(tibble(x=xg, y=-EA(xg)), "x") + geom_point(size=.9, alpha=.75) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_colour_manual(values=COLd, name=expression(hat(Delta))) + coord_cartesian(ylim=c(-0.16,0.16)) +
  labs(title="A   distortion: mean log lambda", x=XL, y=expression(bar(e))) + thm
pB <- ggplot(cell, aes(x, aniso, colour=dlab)) + env(tibble(x=xg, y=EB(xg)), "x") +
  geom_hline(yintercept=0.045, colour="grey70", linewidth=.3, linetype="dotted") + geom_point(size=.9, alpha=.75) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_y_log10(breaks=YB, limits=YL) + scale_colour_manual(values=COLd, guide="none") +
  labs(title="B   distortion: sd log lambda", x=XL, y="s") + thm
pC <- ggplot(cell, aes(x, vsum, colour=dlab)) + geom_hline(yintercept=0, colour="grey40", linewidth=.3) + geom_point(size=.9, alpha=.75) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_colour_manual(values=COLd, guide="none") + coord_cartesian(ylim=c(-0.25,0.25)) +
  labs(title="C   bias: sum of whitened components", x=XL, y=expression(sum(v[p]))) + thm
pD <- ggplot(cell, aes(x, mahal, colour=dlab)) + env(tibble(x=xg, y=pmax(cD*xg^-bD, 0.055)), "x") +
  geom_hline(yintercept=FLOOR_V, colour="grey70", linewidth=.3, linetype="dotted") + geom_point(size=.9, alpha=.75) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_y_log10(breaks=YB, limits=YL) + scale_colour_manual(values=COLd, guide="none") +
  labs(title="D   bias: Fisher metric", x=XL, y=expression(sqrt(b^T*G*b))) + thm
pE <- ggplot(cell, aes(x0, vsum, colour=dlab)) + geom_hline(yintercept=0, colour="grey40", linewidth=.3) + geom_point(size=.9, alpha=.75) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_colour_manual(values=COLd, guide="none") + coord_cartesian(ylim=c(-0.25,0.25)) +
  labs(title="E   bias: sum of whitened components", x=XL0, y=expression(sum(v[p]))) + thm
pF <- ggplot(cell, aes(x0, mahal, colour=dlab)) + env(tibble(x0=xg0, y=pmax(cF*xg0^-bF, 0.055)), "x0") +
  geom_hline(yintercept=FLOOR_V, colour="grey70", linewidth=.3, linetype="dotted") + geom_point(size=.9, alpha=.75) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_y_log10(breaks=YB, limits=YL) + scale_colour_manual(values=COLd, guide="none") +
  annotate("text", x=min(cell$x0)*3, y=0.3, label=sprintf("%.2f*x^{-%.2f}", cF, bF), parse=TRUE, size=2.4, colour="grey30", hjust=0) +
  labs(title="F   bias: Fisher metric", x=XL0, y=expression(sqrt(b^T*G*b))) + thm
g <- (pA + pB) / (pC + pD) / (pE + pF) + plot_layout(guides="collect") & theme(legend.position="right")
ggsave("figures/Figure_IR_6panel.pdf", g, width=6.6, height=7.6)
ggsave("figures/Figure_IR_6panel.png", g, width=6.6, height=7.6, dpi=160)
cat("written figures/Figure_IR_6panel.{pdf,png}\n")
