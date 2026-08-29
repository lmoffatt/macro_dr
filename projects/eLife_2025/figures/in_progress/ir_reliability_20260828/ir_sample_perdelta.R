# Sample distortion fitted SEPARATELY per interval: at fixed Dhat the floor is a constant, no subtraction needed.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
d <- read.csv("digest/ir_components_cells.csv") |> filter(comp=="sample") |> mutate(q=6*(ebar^2+s^2), St=0.1*z, D=interval, xref=nch*St*sqrt(D),
  dlab=factor(sprintf("%.2g", D), levels=sprintf("%.2g", sort(unique(D)))))
fl <- d |> filter(xref > 1e4) |> group_by(D) |> summarise(floor=max(q), .groups="drop")   # per-interval floor = max in the calibrated region
d <- d |> inner_join(fl, by="D")
per <- d |> filter(q > 3*floor) |> group_by(D, dlab, floor) |> group_modify(~{
  m <- lm(log10(q) ~ log10(nch) + log10(St), data=.x); co <- coef(m); r <- 10^abs(resid(m))
  pr <- predict(m, newdata=data.frame(nch=100, St=0.01), se.fit=TRUE)
  tibble(n=nrow(.x), a=co[2], b=co[3], A=10^co[1], R2=summary(m)$r.squared, spread=median(r), se_a=summary(m)$coef[2,2], se_b=summary(m)$coef[3,2],
         lq_ref=pr$fit, se_ref=pr$se.fit) }) |> ungroup() |> mutate(q_lo=10^(lq_ref-se_ref), q_hi=10^(lq_ref+se_ref))
cat("SAMPLE distortion, per-interval fits  log10 q = log10 A + a log10 N + b log10 S~   (cells above 3x that interval's floor)\n")
print(per |> transmute(D, n, floor=signif(floor,2), a=round(a,2), se_a=round(se_a,2), b=round(b,2), se_b=round(se_b,2), b_over_a=round(b/a,2), A=signif(A,2), R2=round(R2,2), spread=round(spread,2)) |> as.data.frame(), row.names=FALSE)
# amplitude at a reference point (N=100, S=0.01): the profile in Dhat, free of the floor
per <- per |> mutate(q_ref = A * 100^a * 0.01^b)
cat("\namplitude at N=100, S~=0.01 by interval (the peak, no subtraction):\n"); print(per |> transmute(D, q_ref=signif(q_ref,2)) |> as.data.frame(), row.names=FALSE)
# for comparison, the same per-interval fit for TOTAL and CORR (are their exponents stable across intervals?)
for(cn in c("total","corr")){ dd <- read.csv("digest/ir_components_cells.csv") |> filter(comp==cn) |> mutate(q=6*(ebar^2+s^2), St=0.1*z, D=interval, xref=nch*St*sqrt(D))
  f2 <- dd |> filter(xref>1e4) |> group_by(D) |> summarise(floor=max(q), .groups="drop"); dd <- dd |> inner_join(f2, by="D") |> filter(q > 3*floor)
  pp <- dd |> group_by(D) |> group_modify(~{ m <- lm(log10(q) ~ log10(nch)+log10(St), data=.x); co <- coef(m); tibble(n=nrow(.x), a=round(co[2],2), b=round(co[3],2), A=signif(10^co[1],2), R2=round(summary(m)$r.squared,2)) }) |> ungroup()
  cat(sprintf("\n%s, per interval:\n", cn)); print(as.data.frame(pp), row.names=FALSE) }

# figure: 7 facets, q vs N*S^(b/a) with the fit and the floor; plus amplitude profile and exponents
d <- d |> inner_join(per |> select(D, a, b, A), by="D") |> mutate(xd = nch * St^(b/a))
COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(d$dlab)), levels(d$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") + theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank(), strip.background=element_blank())
d <- d |> mutate(pred = A * nch^a * St^b)
pA <- ggplot(d, aes(pred, q, colour=dlab)) + geom_abline(slope=1, intercept=0, colour="grey30", linewidth=.4, linetype="dashed") +
  geom_hline(data=fl |> mutate(dlab=factor(sprintf("%.2g", D), levels=levels(d$dlab))), aes(yintercept=floor), colour="grey60", linetype="dotted", linewidth=.3) +
  geom_point(size=.7, alpha=.75) + facet_wrap(~dlab, nrow=1, labeller=label_bquote(hat(Delta)==.(as.character(dlab))), scales="free_x") +
  scale_x_log10() + scale_y_log10(breaks=c(1e-5,1e-3,1e-1), labels=c("1e-5","0.001","0.1")) + scale_colour_manual(values=COLd, guide="none") +
  labs(title="A   sample distortion, one power-law fit per interval: observed vs fitted A*N^a*S^b (dotted = that interval's floor)", x=expression(A%.%N[ch]^a%.%tilde(S)^b), y=expression(sum(log^2*lambda))) + thm
pB <- ggplot(per, aes(D, q_ref)) + geom_line(colour="grey50") + geom_errorbar(aes(ymin=q_lo, ymax=q_hi), width=.12, linewidth=.35, colour="#B2182B") +
  geom_point(colour="#B2182B", size=1.6) + scale_x_log10(breaks=c(0.01,0.1,1)) + scale_y_log10() +
  labs(title="B   amplitude at N=100, S=0.01  (+-1 SE)", x=expression(hat(Delta)), y="fitted q") + thm
pC <- ggplot(per |> pivot_longer(c(a,b), names_to="exp", values_to="v") |> mutate(se=ifelse(exp=="a", se_a, se_b)), aes(D, v, colour=exp)) +
  geom_hline(yintercept=0, colour="grey70", linewidth=.3) + geom_errorbar(aes(ymin=v-se, ymax=v+se), width=.1, linewidth=.3) + geom_point(size=1.4) + geom_line(linewidth=.3) +
  scale_x_log10(breaks=c(0.01,0.1,1)) + scale_colour_manual(values=c(a="#0072B2", b="#E69F00"), labels=c(a=expression(a~(N[ch])), b=expression(b~(tilde(S)))), name=NULL) +
  labs(title="C   exponents by interval", x=expression(hat(Delta)), y="exponent") + thm
g <- pA / (pB + pC + plot_spacer() + plot_layout(widths=c(1,1,1.2))) + plot_layout(heights=c(1,0.9))
ggsave("figures/Figure_IR_sample_perdelta.pdf", g, width=9.0, height=4.6); ggsave("figures/Figure_IR_sample_perdelta.png", g, width=9.0, height=4.6, dpi=150)
cat("written figures/Figure_IR_sample_perdelta.{pdf,png}\n")
