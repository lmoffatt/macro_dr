suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
d <- read.csv("digest/ir_components_cells.csv") |> filter(comp=="sample") |> mutate(q=6*(ebar^2+s^2), St=0.1*z, D=interval, xref=nch*St*sqrt(D),
  dlab=factor(sprintf("%.2g", D), levels=sprintf("%.2g", sort(unique(D)))))
# 1. floor by interval in the calibrated region
flD <- d |> filter(xref > 1e4) |> group_by(D) |> summarise(n=n(), med=median(q), q90=quantile(q,.9), .groups="drop")
cat("floor of the SAMPLE component by interval (calibrated region xref>1e4):\n"); print(flD |> mutate(med=signif(med,2), q90=signif(q90,2), med_over_D=signif(med/D,2)) |> as.data.frame(), row.names=FALSE)
m <- lm(log10(med) ~ log10(D), data=flD); cat(sprintf("  floor ~ D^%.2f  (R2 %.3f);  floor = %.2e * D\n", coef(m)[2], summary(m)$r.squared, 10^coef(m)[1]))
cf <- median(flD$med/flD$D); cat("  using floor(D) = median(med/D) * D =", signif(cf,3), "* D\n")
# 2. subtract the floor (in quadrature-additive sense: q_obs = q_true + floor)
d <- d |> mutate(floorD = cf*D, qc = ifelse(q > floorD, q - floorD, NA_real_), above = q > 3*floorD)
cat("cells clearly above their own floor (q > 3 floor):", sum(d$above), "of", nrow(d), "\n")
s <- d |> filter(above)
D0 <- 0.1; gpk <- function(D) log10(D0/D + D/D0)
fit3 <- function(y, lab){ m <- lm(log10(y) ~ log10(nch) + log10(St) + gpk(D), data=s); co <- coef(m); m2 <- lm(log10(y) ~ log10(nch) + log10(St) + gpk(D) + factor(D), data=s)
  cat(sprintf("  %-18s R2 %.3f gainD %.3f   N %+.3f  S %+.3f  peak %+.3f   -> x = N * S^%.2f * (D0/D+D/D0)^%.2f\n", lab, summary(m)$r.squared, summary(m2)$r.squared-summary(m)$r.squared, co[2], co[3], co[4], co[3]/co[2], co[4]/co[2])); m }
cat("\nproduct+peak fits on the SAME cells:\n"); m_raw <- fit3(s$q, "raw q"); m_cor <- fit3(s$qc, "floor-subtracted")
cat("\nresiduals of the corrected fit by N and by S:\n"); s$res <- resid(m_cor)
print(s |> group_by(nch) |> summarise(med=round(median(res),2), n=n()) |> as.data.frame(), row.names=FALSE)
print(s |> group_by(St) |> summarise(med=round(median(res),2), n=n()) |> as.data.frame(), row.names=FALSE)
# also: is D0 still 0.1 after subtraction?
cat("\nD0 scan on corrected q:\n"); for(D0x in c(0.05,0.07,0.1,0.14,0.2)){ g <- log10(D0x/s$D + s$D/D0x); mm <- lm(log10(qc) ~ log10(nch)+log10(St)+g, data=s); cat(sprintf("  D0=%.2f R2 %.3f\n", D0x, summary(mm)$r.squared)) }
# 3. figure: raw vs corrected against the fitted composite
co <- coef(m_cor); eS <- co[3]/co[2]; eG <- co[4]/co[2]
d <- d |> mutate(xpk = nch*St^eS*(D0/D+D/D0)^eG)
sf <- d |> filter(above); m1 <- lm(log10(qc) ~ log10(xpk), data=sf); cat(sprintf("\n1-D corrected: slope %+.3f R2 %.3f spread x%.2f (n=%d)\n", coef(m1)[2], summary(m1)$r.squared, median(10^abs(resid(m1))), nrow(sf)))
COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(d$dlab)), levels(d$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") + theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank())
YSC <- scale_y_log10(breaks=c(1e-5,1e-4,1e-3,1e-2,1e-1,1), labels=c("1e-5","1e-4","0.001","0.01","0.1","1"), limits=c(2e-6,1.5))
XL <- parse(text=sprintf("N[ch]%%.%%tilde(S)^%.2g%%.%%(0.1/hat(Delta)+hat(Delta)/0.1)^%.2g", eS, eG))
pA <- ggplot(d, aes(xpk, q, colour=dlab)) + geom_point(size=.8, alpha=.7) + scale_x_log10() + YSC + scale_colour_manual(values=COLd, name=expression(hat(Delta))) +
  labs(title="A   raw: floor is proportional to the interval", x=XL, y=expression(sum(log^2*lambda))) + thm
xg <- 10^seq(log10(min(sf$xpk)), log10(max(sf$xpk)), length=100)
pB <- ggplot(sf, aes(xpk, qc, colour=dlab)) + geom_line(data=tibble(x=xg, y=10^predict(m1, newdata=tibble(xpk=xg))), aes(x,y), inherit.aes=FALSE, colour="grey30", linewidth=.4, linetype="dashed") +
  geom_point(size=.8, alpha=.7) + scale_x_log10() + YSC + scale_colour_manual(values=COLd, guide="none") +
  labs(title=sprintf("B   floor-subtracted, cells above 3x their floor   (R2 %.2f)", summary(m1)$r.squared), x=XL, y=expression(sum(log^2*lambda) - floor(hat(Delta)))) + thm
pC <- ggplot(flD, aes(D, med)) + geom_line(data=tibble(D=c(0.01,1), med=cf*c(0.01,1)), colour="grey40", linetype="dashed") + geom_point(colour="#B2182B", size=1.5) +
  scale_x_log10(breaks=c(0.01,0.1,1)) + scale_y_log10() + labs(title="C   floor vs interval (calibrated region)", x=expression(hat(Delta)), y="median q") + thm
g <- pA + pB + pC + plot_layout(widths=c(1,1,0.7), guides="collect") & theme(legend.position="right")
ggsave("figures/Figure_IR_sample_floor.pdf", g, width=9, height=2.8); ggsave("figures/Figure_IR_sample_floor.png", g, width=9, height=2.8, dpi=150)
cat("written figures/Figure_IR_sample_floor.{pdf,png}\n")
