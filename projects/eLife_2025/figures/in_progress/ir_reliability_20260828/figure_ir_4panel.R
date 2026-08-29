# Four panels, all vs x = Nch * S~ * sqrt(Dhat): A ebar, B aniso, C Mahalanobis bias sqrt(b'Gb), D max_p |b_p|/SE_p
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |> select(-zmax, -znrm) |>
  mutate(S=0.1*z, x=nch*S*sqrt(interval),
         dlab=factor(sprintf("%.2g", interval), levels=sprintf("%.2g", sort(unique(interval)))))
ci <- read.csv("digest/ir_bias_cis.csv") |> mutate(zp=abs(m)/se, sig=(lo>0|hi<0)) |>
  group_by(nch,z,interval_in_tau) |> summarise(zmax=max(zp), sig=sig[which.max(zp)], anysig=any(sig), .groups="drop") |>
  rename(interval=interval_in_tau)
mh <- read.csv("digest/ir_bias_mahal.csv") |> select(nch,z,interval,mahal)
dd <- d |> inner_join(ci, by=c("nch","z","interval")) |> inner_join(mh, by=c("nch","z","interval")) |> mutate(lN=log10(nch))
cat("cells:", nrow(dd), "  mahal >= zmax in", sum(dd$mahal >= dd$zmax), "\n")

fit <- function(y, xs, lab){ m <- lm(log10(y) ~ log10(xs)); r <- 10^abs(resid(m))
  sprintf("%-14s slope %+.3f R2 %.3f spread x%.2f", lab, coef(m)[2], summary(m)$r.squared, median(r)) }
sB <- dd |> filter(aniso > 0.055); sC <- dd |> filter(anysig); sD <- dd |> filter(sig)
cat("B aniso   :", fit(sB$aniso, sB$x, "vs x"), "|", fit(sB$aniso, sB$nch, "vs N"), "\n")
cat("C mahal   :", fit(sC$mahal, sC$x, "vs x"), "|", fit(sC$mahal, sC$nch, "vs N"), "  (n sig", nrow(sC), ")\n")
cat("D zmax    :", fit(sD$zmax,  sD$x, "vs x"), "|", fit(sD$zmax,  sD$nch, "vs N"), "  (n sig", nrow(sD), ")\n")
mC <- lm(log10(mahal) ~ log10(nch)+log10(S)+log10(interval), data=sC)
cat(sprintf("C mahal free exponents: N %+.3f S %+.3f D %+.3f R2 %.3f\n", coef(mC)[2],coef(mC)[3],coef(mC)[4],summary(mC)$r.squared))
# envelopes
EA <- function(x) pmax(0.05*x^-0.2, 0.02); EB <- function(x) pmax(0.21*x^-0.2, 0.045)
aC <- max(sC$mahal*sqrt(sC$nch)); cat(sprintf("C: significant mahal*sqrt(N) max = %.2f ; overall mahal max = %.3f\n", aC, max(dd$mahal)))
cat("A viol:", sum(abs(dd$ebar)>EA(dd$x)), " B viol:", sum(dd$aniso>EB(dd$x)), "\n")

COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(dd$dlab)), levels(dd$dlab))
xg <- 10^seq(log10(min(dd$x)), log10(max(dd$x)), length=200)
thm <- theme_bw(base_size=8, base_family="Helvetica") +
  theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank())
XBRK <- 10^seq(-2,10,2); XLABS <- parse(text=sprintf("10^%d", seq(-2,10,2)))
XL <- expression(N[ch]%.%tilde(S)%.%hat(Delta)^0.5)

pA <- ggplot(dd, aes(x, ebar, colour=dlab)) +
  geom_hline(yintercept=0, colour="grey40", linewidth=.3) +
  geom_line(data=tibble(x=xg,y=EA(xg)), aes(x,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
  geom_line(data=tibble(x=xg,y=-EA(xg)), aes(x,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
  geom_point(size=.9, alpha=.75) + scale_x_log10(breaks=XBRK, labels=XLABS) +
  scale_colour_manual(values=COLd, name=expression(hat(Delta))) + coord_cartesian(ylim=c(-0.16,0.16)) +
  labs(title="A   distortion: net amplitude", x=XL, y=expression(bar(e)==mean~log~lambda)) + thm
pB <- ggplot(dd, aes(x, aniso, colour=dlab)) +
  geom_line(data=tibble(x=xg,y=EB(xg)), aes(x,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
  geom_hline(yintercept=0.045, colour="grey70", linewidth=.3, linetype="dotted") +
  geom_point(size=.9, alpha=.75) + scale_x_log10(breaks=XBRK, labels=XLABS) + scale_y_log10(breaks=c(0.02,0.05,0.1,0.2,0.3)) +
  scale_colour_manual(values=COLd, guide="none") +
  annotate("text", x=min(dd$x)*3, y=0.30, label="0.21*x^{-1/5}", parse=TRUE, size=2.4, colour="grey30", hjust=0) +
  labs(title="B   distortion: anisotropy", x=XL, y=expression(s==sd~log~lambda)) + thm
pC <- ggplot(dd, aes(x, mahal, colour=lN)) +
  geom_hline(yintercept=0.3, linetype="dashed", colour="grey30", linewidth=.4) +
  geom_point(aes(shape=anysig), size=1.1, alpha=.8) +
  scale_shape_manual(values=c(`TRUE`=16, `FALSE`=1), guide="none") +
  scale_colour_viridis_c(name=expression(log[10]~N[ch]), option="C", end=.9) +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_y_log10(breaks=c(0.01,0.03,0.1,0.3), labels=c("0.01","0.03","0.1","0.3")) +
  labs(title="C   bias in the Fisher metric", x=XL, y=expression(sqrt(b^T*G*b))) + thm
pD <- ggplot(dd, aes(x, zmax, colour=lN)) +
  geom_hline(yintercept=0.3, linetype="dashed", colour="grey30", linewidth=.4) +
  geom_point(aes(shape=sig), size=1.1, alpha=.8) +
  scale_shape_manual(values=c(`TRUE`=16, `FALSE`=1), name="CI excludes 0", labels=c(`TRUE`="yes",`FALSE`="no")) +
  scale_colour_viridis_c(option="C", end=.9, guide="none") +
  scale_x_log10(breaks=XBRK, labels=XLABS) + scale_y_log10(breaks=c(0.003,0.01,0.03,0.1,0.3), labels=c("0.003","0.01","0.03","0.1","0.3")) +
  labs(title="D   bias per parameter, worst", x=XL, y=expression(max[p]~"|"*b[p]*"|"/SE[p])) + thm
g <- (pA + pB) / (pC + pD) + plot_layout(guides="collect") & theme(legend.position="right")
ggsave("figures/Figure_IR_4panel.pdf", g, width=6.4, height=5.0)
ggsave("figures/Figure_IR_4panel.png", g, width=6.4, height=5.0, dpi=170)
cat("written figures/Figure_IR_4panel.{pdf,png}\n")
