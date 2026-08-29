# Per-interval power-law fits (log q = log A + a log N + b log S) for bias (b'Gb), total, sample, correlation (sum log^2 lambda).
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
comp <- read.csv("digest/ir_components_cells.csv") |> mutate(q=6*(ebar^2+s^2)) |> select(nch,z,interval,comp,q)
bias <- read.csv("digest/ir_bias_whitened.csv") |> group_by(nch,z,interval) |> summarise(q=sum(v^2), .groups="drop") |> mutate(comp="bias")
all <- bind_rows(bias, comp) |> mutate(St=0.1*z, D=interval, xref=nch*St*sqrt(D), comp=factor(comp, levels=c("bias","total","sample","corr")))
fl <- all |> filter(xref > 1e4) |> group_by(comp, D) |> summarise(floor=max(q), .groups="drop")
all <- all |> inner_join(fl, by=c("comp","D"))
per <- all |> filter(q > 3*floor) |> group_by(comp, D, floor) |> group_modify(~{
  if(nrow(.x) < 6) return(tibble(n=nrow(.x), a=NA, b=NA, A=NA, R2=NA, spread=NA, se_a=NA, se_b=NA))
  m <- lm(log10(q) ~ log10(nch) + log10(St), data=.x); co <- coef(m); sm <- summary(m)$coef
  pr <- predict(m, newdata=data.frame(nch=100, St=0.01), se.fit=TRUE)   # log10 q at the reference point, with its SE
  tibble(n=nrow(.x), a=co[2], b=co[3], A=10^co[1], R2=summary(m)$r.squared, spread=median(10^abs(resid(m))), se_a=sm[2,2], se_b=sm[3,2],
         lq_ref=pr$fit, se_ref=pr$se.fit) }) |> ungroup() |>
  mutate(q_ref = 10^lq_ref, q_lo = 10^(lq_ref - se_ref), q_hi = 10^(lq_ref + se_ref))
for(cn in levels(all$comp)){ cat(sprintf("\n=== %s ===\n", cn))
  print(per |> filter(comp==cn) |> transmute(D, n, floor=signif(floor,2), a=round(a,2), se_a=round(se_a,2), b=round(b,2), se_b=round(se_b,2), q_ref=signif(q_ref,2), q_lo=signif(q_lo,2), q_hi=signif(q_hi,2), R2=round(R2,2), spread=round(spread,2)) |> as.data.frame(), row.names=FALSE) }
write.csv(per, "digest/ir_perdelta_all.csv", row.names=FALSE)

thm <- theme_bw(base_size=8, base_family="Helvetica") + theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank())
TT <- c(bias="bias  b'Gb", total="total distortion", sample="sample distortion", corr="correlation distortion")
amp <- function(cn, lab){ p <- per |> filter(comp==cn, !is.na(A))
  ggplot(p, aes(D, q_ref)) + geom_line(colour="grey50") + geom_errorbar(aes(ymin=q_lo, ymax=q_hi), width=.12, linewidth=.35, colour="#B2182B") +
    geom_point(colour="#B2182B", size=1.6) + scale_x_log10(breaks=c(0.01,0.1,1)) + scale_y_log10() +
    labs(title=paste0(lab, "   ", TT[[cn]]), x=expression(hat(Delta)), y=expression(q~at~N==100*","~tilde(S)==0.01~~(""%+-%""*1~SE))) + thm }
expo <- function(cn, lab){ p <- per |> filter(comp==cn, !is.na(a)) |> pivot_longer(c(a,b), names_to="exp", values_to="v") |> mutate(se=ifelse(exp=="a", se_a, se_b))
  ggplot(p, aes(D, v, colour=exp)) + geom_hline(yintercept=0, colour="grey70", linewidth=.3) + geom_errorbar(aes(ymin=v-se, ymax=v+se), width=.1, linewidth=.3) +
    geom_point(size=1.4) + geom_line(linewidth=.3) + scale_x_log10(breaks=c(0.01,0.1,1)) + coord_cartesian(ylim=c(-1.5,0.3)) +
    scale_colour_manual(values=c(a="#0072B2", b="#E69F00"), labels=c(a=expression(a~(N[ch])), b=expression(b~(tilde(S)))), name=NULL) +
    labs(title=lab, x=expression(hat(Delta)), y="exponent") + thm }
g <- (amp("bias","A") + amp("total","B") + amp("sample","C") + amp("corr","D") + plot_layout(nrow=1)) /
     (expo("bias","E") + expo("total","F") + expo("sample","G") + expo("corr","H") + plot_layout(nrow=1, guides="collect")) & theme(legend.position="right")
ggsave("figures/Figure_IR_perdelta_all.pdf", g, width=9.4, height=4.6); ggsave("figures/Figure_IR_perdelta_all.png", g, width=9.4, height=4.6, dpi=150)
cat("written figures/Figure_IR_perdelta_all.{pdf,png}\n")
