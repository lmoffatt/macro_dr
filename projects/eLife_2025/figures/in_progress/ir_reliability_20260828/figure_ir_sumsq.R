# One row: sum of squares  (bias: b'Gb = sum v_p^2 ; distortions: sum log^2 lambda = d_AI^2 = 6(ebar^2+s^2)) vs the rounded composites.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
comp <- read.csv("digest/ir_components_cells.csv") |> mutate(q = 6*(ebar^2 + s^2))   # = sum log^2 lambda (dAI column is NA: tibble masking bug)
bias <- read.csv("digest/ir_bias_whitened.csv") |> group_by(nch,z,interval) |> summarise(q=sum(v^2), .groups="drop") |> mutate(comp="bias")
all <- bind_rows(bias, comp |> select(nch,z,interval,comp,q)) |>
  mutate(St=0.1*z, D=interval, xref=nch*St*sqrt(D), dlab=factor(sprintf("%.2g", D), levels=sprintf("%.2g", sort(unique(D)))),
         comp=factor(comp, levels=c("bias","total","sample","corr")))
fits <- tibble(comp=factor(c("bias","total","sample","corr"), levels=levels(all$comp)), eS=c(0.75,1,2,1), eD=c(0,0.5,0,1), slope=c(-1, -2/5, -1/4, -2/5))
all <- all |> inner_join(fits, by="comp") |> mutate(xb = nch*St^eS*D^eD)
fl <- all |> filter(xref > 1e4) |> group_by(comp) |> summarise(floor=max(q), .groups="drop"); all <- all |> inner_join(fl, by="comp")
chk <- all |> filter(q > 1.5*floor) |> group_by(comp) |> group_modify(~{
  m1 <- lm(log10(q) ~ log10(xb), data=.x); m2 <- lm(log10(q) ~ log10(xb)+dlab, data=.x); r <- 10^abs(resid(m1))
  tibble(slope_free=coef(m1)[2], R2=summary(m1)$r.squared, gainD=summary(m2)$r.squared-summary(m1)$r.squared, spread=median(r), n=nrow(.x),
         a=ceiling(1000*max(.x$q*.x$xb^-.x$slope[1]))/1000) }) |> ungroup()
tab <- fits |> inner_join(chk, by="comp") |> inner_join(fl, by="comp")
viol <- all |> inner_join(tab |> select(comp,a), by="comp") |> mutate(v = q > pmax(a*xb^slope, 1.5*floor)) |> group_by(comp) |> summarise(viol=sum(v))
cat("SUM OF SQUARES vs rounded composites (slopes doubled):\n")
for(i in seq_len(nrow(tab))) with(tab[i,], cat(sprintf("  %-7s x=N*S^%.2g*D^%.2g  slope fixed %+.2f free %+.3f  R2 %.3f  gainD %.3f  spread x%.2f  a=%.3f floor=%.4f  n=%d\n", comp, eS,eD, slope, slope_free, R2, gainD, spread, a, floor, n)))
print(as.data.frame(viol), row.names=FALSE)
# for comparison: same fits on s (sd) and on sqrt(q)
cat("\n(for reference: R2 on sd log lambda from the previous figure: bias 0.889 total 0.819 sample 0.520 corr 0.741)\n")

COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(all$dlab)), levels(all$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") + theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank())
TT <- c(bias="bias:  b'Gb", total="total:  sum log^2 lambda", sample="sample:  sum log^2 lambda", corr="correlation:  sum log^2 lambda")
YL <- c(bias="b^T*G*b", total="sum(log^2*lambda)", sample="sum(log^2*lambda)", corr="sum(log^2*lambda)")
xlab_of <- function(eS, eD){ f <- function(e, sym) if(abs(e) < 0.05) "" else if(abs(e-1) < 0.05) sym else sprintf("%s^%.2g", sym, e)
  parse(text=paste(Filter(nzchar, c("N[ch]", f(eS,"tilde(S)"), f(eD,"hat(Delta)"))), collapse="%.%")) }
pan <- function(cn, lab){ d <- all |> filter(comp==cn); t <- tab |> filter(comp==cn); xg <- 10^seq(log10(min(d$xb)), log10(max(d$xb)), length=200)
  ggplot(d, aes(xb, q, colour=dlab)) + geom_hline(yintercept=t$floor, colour="grey70", linewidth=.3, linetype="dotted") +
    geom_line(data=tibble(x=xg, y=pmax(t$a*xg^t$slope, 1.5*t$floor)), aes(x,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
    geom_point(size=.8, alpha=.7) + scale_x_log10() +
    scale_y_log10(breaks=c(1e-5,1e-4,1e-3,1e-2,1e-1,1), labels=c("1e-5","1e-4","0.001","0.01","0.1","1"), limits=c(2e-6,1.5)) +
    scale_colour_manual(values=COLd, name=expression(hat(Delta))) +
    annotate("text", x=min(d$xb)*2, y=1.1, label=sprintf("%.3f*x^{%s}~~(R^2==%.2f)", t$a, c("-1","-2/5","-1/4","-2/5")[as.integer(t$comp)], t$R2), parse=TRUE, size=2.3, colour="grey30", hjust=0) +
    labs(title=paste0(lab, "   ", TT[[cn]]), x=xlab_of(t$eS, t$eD), y=parse(text=YL[[cn]])) + thm }
g <- wrap_plots(pan("bias","A"), pan("total","B"), pan("sample","C"), pan("corr","D"), nrow=1) + plot_layout(guides="collect") & theme(legend.position="right")
ggsave("figures/Figure_IR_sumsq.pdf", g, width=9.0, height=2.6); ggsave("figures/Figure_IR_sumsq.png", g, width=9.0, height=2.6, dpi=150)
cat("written figures/Figure_IR_sumsq.{pdf,png}\n")
