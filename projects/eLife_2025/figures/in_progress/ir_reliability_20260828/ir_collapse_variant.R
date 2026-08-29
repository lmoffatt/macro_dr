setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
# Generic collapse variant: x = Nch^pN * S~^pS * Dhat^pD. Usage: Rscript ir_collapse_variant.R pN pS pD tag
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
a <- commandArgs(TRUE); pN <- as.numeric(a[1]); pS <- as.numeric(a[2]); pD <- as.numeric(a[3]); tag <- a[4]
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |>
  mutate(S=0.1*z, x=nch^pN*S^pS*interval^pD,
         dlab = factor(sprintf("%.2g", interval), levels=sprintf("%.2g", sort(unique(interval)))))
s <- d |> filter(aniso > 0.055)
m <- lm(log10(aniso) ~ log10(x), data=s); r <- 10^abs(resid(m))
cat(sprintf("x = N^%g S^%g D^%g : slope %+.3f  R2 %.3f  spread median x%.2f  q90 x%.2f  (n=%d above floor)\n",
    pN,pS,pD, coef(m)[2], summary(m)$r.squared, median(r), quantile(r,.9), nrow(s)))
# envelope constants at slope -1/5 (min a, zero violations above floor)
cells <- d |> filter(aniso > 0.045); aB <- max(cells$aniso*cells$x^0.2)
cellsA <- d |> filter(abs(ebar) > 0.02); aA <- max(abs(cellsA$ebar)*cellsA$x^0.2)
aB <- ceiling(aB*100)/100; aA <- ceiling(aA*100)/100
cat(sprintf("envelopes at -1/5: s <= max(%.2f x^-0.2, 0.045)   |ebar| <= max(%.2f x^-0.2, 0.02)\n", aB, aA))
# also: interval stratification test = R2 gain from adding interval factor to the 1-D fit
m2 <- lm(log10(aniso) ~ log10(x) + dlab, data=s)
cat(sprintf("residual interval stratification: R2 %.3f -> %.3f with Delta as factor (gain %.3f)\n",
    summary(m)$r.squared, summary(m2)$r.squared, summary(m2)$r.squared-summary(m)$r.squared))

EA <- function(x) pmax(aA*x^-0.2, 0.02); EB <- function(x) pmax(aB*x^-0.2, 0.045)
COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(d$dlab)), levels(d$dlab))
xg <- 10^seq(log10(min(d$x)), log10(max(d$x)), length=200)
thm <- theme_bw(base_size=8, base_family="Helvetica") +
  theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.32,"cm"), panel.grid.minor=element_blank())
lo <- floor(log10(min(d$x))); hi <- ceiling(log10(max(d$x))); XBRK <- 10^seq(lo,hi,2); XLABS <- parse(text=sprintf("10^%d", seq(lo,hi,2)))
fx <- function(p, sym) if(p==1) sym else if(p==0) "" else sprintf("%s^%g", sym, p)
XL <- parse(text=paste(Filter(nzchar, c(fx(pN,"N[ch]"), fx(pS,"tilde(S)"), fx(pD,"hat(Delta)"))), collapse="%.%"))

pA <- ggplot(d, aes(x, ebar, colour=dlab)) +
  geom_hline(yintercept=0, colour="grey40", linewidth=.3) +
  geom_line(data=tibble(x=xg, y= EA(xg)), aes(x,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
  geom_line(data=tibble(x=xg, y=-EA(xg)), aes(x,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
  geom_point(size=.9, alpha=.75) + scale_x_log10(breaks=XBRK, labels=XLABS) +
  scale_colour_manual(values=COLd, guide="none") + coord_cartesian(ylim=c(-0.16,0.16)) +
  labs(title="A   net amplitude", x=XL, y=expression(bar(e)==mean~log~lambda)) + thm
pB <- ggplot(d, aes(x, aniso, colour=dlab)) +
  geom_line(data=tibble(x=xg, y=EB(xg)), aes(x,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
  geom_hline(yintercept=0.045, colour="grey70", linewidth=.3, linetype="dotted") +
  geom_point(size=.9, alpha=.75) + scale_x_log10(breaks=XBRK, labels=XLABS) + scale_y_log10(breaks=c(0.02,0.05,0.1,0.2,0.3)) +
  scale_colour_manual(values=COLd, name=expression(hat(Delta)==Delta%.%k[off])) +
  annotate("text", x=min(d$x)*3, y=0.30, label=sprintf("%.2f*x^{-1/5}", aB), parse=TRUE, size=2.4, colour="grey30", hjust=0) +
  labs(title=sprintf("B   anisotropy   (R2 %.2f, spread x%.2f)", summary(m)$r.squared, median(r)), x=XL, y=expression(s==sd~log~lambda)) + thm
g <- pA + pB + plot_layout(widths=c(1,1.2))
ggsave(sprintf("figures/Figure_IR_collapse_%s.pdf", tag), g, width=5.2, height=2.7)
ggsave(sprintf("figures/Figure_IR_collapse_%s.png", tag), g, width=5.2, height=2.7, dpi=170)
cat(sprintf("written figures/Figure_IR_collapse_%s.{pdf,png}\n", tag))
