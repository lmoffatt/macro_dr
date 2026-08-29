suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
d <- read.csv("digest/ir_components_cells.csv") |> filter(comp=="sample") |> mutate(q=6*(ebar^2+s^2), St=0.1*z, D=interval, xref=nch*St*sqrt(D),
  dlab=factor(sprintf("%.2g", D), levels=sprintf("%.2g", sort(unique(D)))))
floor <- max(d$q[d$xref>1e4]); s <- d |> filter(q > 1.5*floor)
D0 <- 0.1; gpk <- function(D) log10(D0/D + D/D0)
# saturating noise fraction: q ~ N^a * (1 + c*S/(N*D))^b * peak(D)   -> grid over c
cat("SATURATING form  log q = a logN + b log(1 + c*S/(N*D)) + k*g(D):\n")
for(cc in c(0.03,0.1,0.3,1,3,10,30,100)){ m <- lm(log10(q) ~ log10(nch) + log10(1 + cc*St/(nch*D)) + gpk(D), data=s); co <- coef(m)
  m2 <- lm(log10(q) ~ log10(nch) + log10(1 + cc*St/(nch*D)) + gpk(D) + factor(D), data=s)
  cat(sprintf("  c=%6.2f  R2 %.3f gainD %.3f   a %+.3f  b %+.3f  k %+.3f\n", cc, summary(m)$r.squared, summary(m2)$r.squared-summary(m)$r.squared, co[2], co[3], co[4])) }
cat("\nfor scale, product form with peak:  log q = a logN + b logS + k g(D)\n")
m <- lm(log10(q) ~ log10(nch) + log10(St) + gpk(D), data=s); co <- coef(m); cat(sprintf("  R2 %.3f   a %+.3f b %+.3f k %+.3f\n", summary(m)$r.squared, co[2], co[3], co[4]))
# residual structure of the product+peak form: by N and by S
s$res <- resid(m)
cat("residual (dex) by N:\n"); print(s |> group_by(nch) |> summarise(med=round(median(res),2), n=n()) |> as.data.frame(), row.names=FALSE)
cat("residual (dex) by S~:\n"); print(s |> group_by(St) |> summarise(med=round(median(res),2), n=n()) |> as.data.frame(), row.names=FALSE)

# figure: sample sum log^2 lambda vs x_peak = N * S^1.7 * (D0/D + D/D0)^1.3   (rounded: S^1.5, exponent 1.3 -> keep fitted)
eS <- 1.7; eG <- 1.3
d <- d |> mutate(xpk = nch * St^eS * (D0/D + D/D0)^eG, xpow = nch*St^2)
sfit <- d |> filter(q > 1.5*floor)
for(v in c("xpow","xpk")){ m1 <- lm(log10(sfit$q) ~ log10(sfit[[v]])); m2 <- lm(log10(sfit$q) ~ log10(sfit[[v]]) + sfit$dlab)
  cat(sprintf("1-D fit vs %-5s slope %+.3f R2 %.3f gainD %.3f spread x%.2f\n", v, coef(m1)[2], summary(m1)$r.squared, summary(m2)$r.squared-summary(m1)$r.squared, median(10^abs(resid(m1))))) }
COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(d$dlab)), levels(d$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") + theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank())
YSC <- scale_y_log10(breaks=c(1e-5,1e-4,1e-3,1e-2,1e-1,1), labels=c("1e-5","1e-4","0.001","0.01","0.1","1"), limits=c(2e-6,1.5))
pA <- ggplot(d, aes(xpow, q, colour=dlab)) + geom_hline(yintercept=floor, colour="grey70", linewidth=.3, linetype="dotted") + geom_point(size=.8, alpha=.7) + scale_x_log10() + YSC +
  scale_colour_manual(values=COLd, name=expression(hat(Delta))) + labs(title="A   sample distortion vs power law", x=expression(N[ch]%.%tilde(S)^2), y=expression(sum(log^2*lambda))) + thm
pB <- ggplot(d, aes(xpk, q, colour=dlab)) + geom_hline(yintercept=floor, colour="grey70", linewidth=.3, linetype="dotted") + geom_point(size=.8, alpha=.7) + scale_x_log10() + YSC +
  scale_colour_manual(values=COLd, guide="none") + labs(title="B   sample distortion vs peaked variable", x=expression(N[ch]%.%tilde(S)^1.7%.%(0.1/hat(Delta)+hat(Delta)/0.1)^1.3), y=expression(sum(log^2*lambda))) + thm
prof <- d |> group_by(nch,z) |> filter(n()==7, max(q) > 3*floor) |> mutate(qn=q/max(q)) |> ungroup()
pC <- ggplot(prof, aes(D, qn, group=interaction(nch,z))) + geom_line(colour="grey60", linewidth=.3, alpha=.6) +
  stat_summary(aes(group=1), fun=median, geom="line", colour="#B2182B", linewidth=.8) + scale_x_log10(breaks=c(0.01,0.1,1)) +
  labs(title="C   profile in the interval, 46 (N, S) cells", x=expression(hat(Delta)), y=expression(q/q[max])) + thm
g <- pA + pB + pC + plot_layout(widths=c(1,1,0.8), guides="collect") & theme(legend.position="right")
ggsave("figures/Figure_IR_sample_peak.pdf", g, width=8.6, height=2.7); ggsave("figures/Figure_IR_sample_peak.png", g, width=8.6, height=2.7, dpi=150)
cat("written figures/Figure_IR_sample_peak.{pdf,png}\n")
