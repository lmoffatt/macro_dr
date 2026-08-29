setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
# Variant: A/B against the TRIPLE product x2 = S~ * Nch * (Delta*koff); C unchanged.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |>
  mutate(S=0.1*z, x=S*nch*interval,
         dlab = factor(sprintf("%.2g", interval), levels=sprintf("%.2g", sort(unique(interval)))))
ci <- read.csv("digest/ir_bias_cis.csv") |>
  mutate(zp = abs(m)/se, sig = (lo>0 | hi<0)) |>
  group_by(nch, z, interval_in_tau) |> slice_max(zp, n=1, with_ties=FALSE) |>
  ungroup() |> transmute(nch, z, interval=interval_in_tau, zmax_ci=zp, sig)
dc <- d |> inner_join(ci, by=c("nch","z","interval")) |> mutate(Slab = log10(S))

EA <- function(x) pmax(0.04*x^-0.2, 0.02)
EB <- function(x) pmax(0.21*x^-0.2, 0.045)
EC <- function(n) pmax(0.8/sqrt(n), 0.04)
cat("A violations:", sum(abs(d$ebar) > EA(d$x)), " B violations:", sum(d$aniso > EB(d$x)),
    " C violations (sig):", nrow(dc |> filter(sig, zmax_ci > EC(nch))), "\n")

COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(d$dlab)), levels(d$dlab))
xg <- 10^seq(log10(min(d$x)), log10(max(d$x)), length=200)
thm <- theme_bw(base_size=8, base_family="Helvetica") +
  theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.32,"cm"),
        panel.grid.minor=element_blank())
XBRK <- 10^seq(-4,10,2); XLABS <- parse(text=sprintf("10^%d", seq(-4,10,2)))
XL <- expression(tilde(S)%.%N[ch]%.%Delta%.%k[off])

pA <- ggplot(d, aes(x, ebar, colour=dlab)) +
  geom_hline(yintercept=0, colour="grey40", linewidth=.3) +
  geom_line(data=tibble(x=xg, y= EA(xg)), aes(x,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
  geom_line(data=tibble(x=xg, y=-EA(xg)), aes(x,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
  geom_point(size=.9, alpha=.75) +
  scale_x_log10(breaks=XBRK, labels=XLABS) +
  scale_colour_manual(values=COLd, guide="none") +
  coord_cartesian(ylim=c(-0.16,0.16)) +
  labs(title="A   net amplitude: subordinate", x=XL, y=expression(bar(e)==mean~log~lambda)) + thm

pB <- ggplot(d, aes(x, aniso, colour=dlab)) +
  geom_line(data=tibble(x=xg, y=EB(xg)), aes(x,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
  geom_hline(yintercept=0.045, colour="grey70", linewidth=.3, linetype="dotted") +
  geom_point(size=.9, alpha=.75) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_y_log10(breaks=c(0.02,0.05,0.1,0.2,0.3)) +
  scale_colour_manual(values=COLd, name=expression(Delta%.%k[off])) +
  annotate("text", x=1e10, y=0.049, label="resolution floor", size=2.2, colour="grey45", hjust=1, vjust=0) +
  annotate("text", x=1e-2, y=0.30, label="0.21*x^{-1/5}", parse=TRUE, size=2.4, colour="grey30", hjust=0) +
  labs(title="B   anisotropy vs the triple product", x=XL, y=expression(s==sd~log~lambda)) + thm

ng <- 10^seq(log10(5), 4, length=100)
pC <- ggplot(dc, aes(nch, zmax_ci)) +
  geom_line(data=tibble(nch=ng, y=EC(ng)), aes(nch,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
  geom_point(aes(fill=Slab, shape=sig), size=1.2, alpha=.8, colour="grey25", stroke=.25) +
  scale_shape_manual(values=c(`TRUE`=21, `FALSE`=1), name="bias CI\nexcludes 0",
                     labels=c(`TRUE`="yes",`FALSE`="no"),
                     guide=guide_legend(override.aes=list(fill=c("grey55","white"), size=2))) +
  scale_fill_gradient(low="#FEE8C8", high="#7F0000", name=expression(log[10]~tilde(S))) +
  scale_x_log10(breaks=c(10,100,1000,10000), labels=c("10","100","1000","10000")) +
  scale_y_log10(breaks=c(0.003,0.01,0.03,0.1,0.3), labels=c("0.003","0.01","0.03","0.1","0.3")) +
  annotate("text", x=25, y=0.24, label="0.8/sqrt(N[ch])", parse=TRUE, size=2.4, colour="grey30", hjust=0) +
  annotate("text", x=9000, y=0.047, label="0.04", size=2.4, colour="grey30", hjust=1, vjust=0) +
  labs(title="C   bias, channels alone", x=expression(N[ch]), y=expression(max[p]~"|bias|"/SE)) + thm

g <- pA + pB + pC + plot_layout(widths=c(1,1.15,1.15))
ggsave("figures/Figure_IR_reliability_SNd.pdf", g, width=7.0, height=2.7)
ggsave("figures/Figure_IR_reliability_SNd.png", g, width=7.0, height=2.7, dpi=170)
cat("written figures/Figure_IR_reliability_SNd.{pdf,png}\n")
