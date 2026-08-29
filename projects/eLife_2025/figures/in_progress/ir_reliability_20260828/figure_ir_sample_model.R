suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
d <- read.csv("digest/ir_components_cells.csv") |> filter(comp=="sample") |> mutate(q=6*(ebar^2+s^2), St=0.1*z, D=interval,
  dlab=factor(sprintf("%.2g", D), levels=sprintf("%.2g", sort(unique(D)))))
cf <- 0.00139; fit <- readRDS("digest/ir_sample_nls.rds")$p
# ROUNDED model: q = A / (N^(1/4) * S * (0.1/D + D/0.1)) + cf*D, only A free
xr <- with(d, nch^0.25 * St * (0.1/D + D/0.1)); TSS <- sum((log(d$q)-mean(log(d$q)))^2)
oA <- optimize(function(lA) sum((log(d$q) - log(exp(lA)/xr + cf*d$D))^2), c(-12, 0))
A <- exp(oA$minimum); R2r <- 1 - oA$objective/TSS
cat(sprintf("ROUNDED model  q = %.4f / (N^(1/4) * S * (0.1/D + D/0.1)) + %.2e * D :  R2(log) %.3f  (free-exponent model 0.847; no floor 0.610)\n", A, cf, R2r))
d <- d |> mutate(xr = xr, pred = A/xr + cf*D, signal = A/xr, qc = q - cf*D, dom = signal > 3*cf*D)
cat("signal-dominated cells:", sum(d$dom), "; 1:1 spread (all cells) x", round(exp(median(abs(log(d$q/d$pred)))),2), "; signal-dominated only x", round(exp(median(abs(log(d$q[d$dom]/d$pred[d$dom])))),2), "\n")
COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(d$dlab)), levels(d$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") + theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank())
BR <- c(1e-5,1e-4,1e-3,1e-2,1e-1,1); LB <- c("1e-5","1e-4","0.001","0.01","0.1","1")
pA <- ggplot(d, aes(pred, q, colour=dlab)) + geom_abline(slope=1, intercept=0, colour="grey40", linewidth=.3) + geom_point(size=.8, alpha=.7) +
  scale_x_log10(breaks=BR, labels=LB) + scale_y_log10(breaks=BR, labels=LB) + coord_equal() + scale_colour_manual(values=COLd, name=expression(hat(Delta))) +
  labs(title=sprintf("A   sample distortion: model vs observed  (R2 %.2f)", R2r), x=expression(0.0028/(N[ch]^{1/4}~tilde(S)~(0.1/hat(Delta)+hat(Delta)/0.1)) + 0.0014~hat(Delta)), y=expression(sum(log^2*lambda))) + thm
s <- d |> filter(dom, qc > 0) |> mutate(ix = 1/xr); xg <- 10^seq(log10(min(s$ix)), log10(max(s$ix)), length=50)
pB <- ggplot(s, aes(ix, qc, colour=dlab)) + geom_line(data=tibble(x=xg, y=A*xg), aes(x,y), inherit.aes=FALSE, colour="grey30", linewidth=.4, linetype="dashed") +
  geom_point(size=.8, alpha=.7) + scale_x_log10() + scale_y_log10(breaks=BR, labels=LB) + scale_colour_manual(values=COLd, guide="none") +
  labs(title="B   signal part: floor subtracted, slope fixed at -1", x=expression(1/(N[ch]^{1/4}~tilde(S)~(0.1/hat(Delta)+hat(Delta)/0.1))), y=expression(sum(log^2*lambda) - 0.0014~hat(Delta))) + thm
g <- pA + pB + plot_layout(guides="collect") & theme(legend.position="right")
ggsave("figures/Figure_IR_sample_model.pdf", g, width=7.6, height=3.4); ggsave("figures/Figure_IR_sample_model.png", g, width=7.6, height=3.4, dpi=150)
cat("written figures/Figure_IR_sample_model.{pdf,png}\n")
