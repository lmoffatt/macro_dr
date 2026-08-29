# Each quantity against its own best-fitting composite  x_best = N * S^(bS/bN) * D^(bD/bN)  (fit on cells above the resolution floor)
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
comp <- read.csv("digest/ir_components_cells.csv")
bias <- read.csv("digest/ir_bias_whitened.csv") |> group_by(nch,z,interval) |> summarise(ebar=mean(v), s=sqrt(mean((v-mean(v))^2)), .groups="drop") |> mutate(comp="bias")
all <- bind_rows(bias, comp |> select(nch,z,interval,comp,ebar,s)) |>
  mutate(St=0.1*z, D=interval, xref=nch*St*sqrt(D), dlab=factor(sprintf("%.2g", D), levels=sprintf("%.2g", sort(unique(D)))),
         comp=factor(comp, levels=c("bias","total","sample","corr")))
fl <- all |> filter(xref > 1e4) |> group_by(comp) |> summarise(floor=max(s), .groups="drop")
all <- all |> inner_join(fl, by="comp")
# ROUNDED exponents (grid-searched, tmp/ir_round_exponents.R) and rounded slopes
fits <- tibble(comp=factor(c("bias","total","sample","corr"), levels=levels(all$comp)),
               eS=c(0.75, 1, 2, 1), eD=c(0, 0.5, 0, 1), slope_fix=c(-1/2, -1/5, -1/8, -1/5),
               n=NA_integer_, bN=NA_real_, bS=NA_real_, bD=NA_real_, R2=NA_real_)
all <- all |> inner_join(fits, by="comp") |> mutate(xb = nch * St^eS * D^eD)
# 1-D check vs x_best and interval-residual, envelope constant
chk <- all |> filter(s > 1.5*floor) |> group_by(comp) |> group_modify(~{
  m1 <- lm(log10(s) ~ log10(xb), data=.x); m2 <- lm(log10(s) ~ log10(xb) + dlab, data=.x); r <- 10^abs(resid(m1))
  tibble(slope_free=coef(m1)[2], R2_1d=summary(m1)$r.squared, gainD=summary(m2)$r.squared-summary(m1)$r.squared, spread=median(r), n=nrow(.x)) }) |> ungroup() |>
  inner_join(fits |> select(comp, slope=slope_fix), by="comp")
envc <- all |> filter(s > 1.5*floor) |> inner_join(chk, by="comp") |> group_by(comp) |> summarise(a=ceiling(100*max(s*xb^-slope[1]))/100, .groups="drop")
tab <- fits |> select(comp,eS,eD) |> inner_join(chk, by="comp") |> inner_join(envc, by="comp") |> inner_join(fl, by="comp")
cat("component  x = N*S^eS*D^eD   slope(fixed) slope(free)  R2(1d)  gainD  spread  envelope a  floor\n")
for(i in seq_len(nrow(tab))) with(tab[i,], cat(sprintf("%-8s  N*S^%.2f*D^%.1f   %+.3f   %+.3f   %.3f   %.3f   x%.2f   %.2f   %.3f  n=%d\n", comp, eS,eD, slope, slope_free, R2_1d, gainD, spread, a, floor, n)))
viol <- all |> inner_join(tab |> select(comp,a,slope), by="comp") |> mutate(v = s > pmax(a*xb^slope, 1.5*floor)) |> group_by(comp) |> summarise(viol=sum(v)); print(as.data.frame(viol), row.names=FALSE)

COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(all$dlab)), levels(all$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") + theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank())
TT <- c(bias="bias  (v = G^{1/2} b)", total="total distortion", sample="sample distortion", corr="correlation distortion")
YM <- c(bias="mean~v[p]", total="mean~log~lambda", sample="mean~log~lambda", corr="mean~log~lambda")
YS <- c(bias="sd~v[p]", total="sd~log~lambda", sample="sd~log~lambda", corr="sd~log~lambda")
xlab_of <- function(eS, eD){ f <- function(e, sym) if(abs(e) < 0.05) "" else if(abs(e-1) < 0.05) sym else sprintf("%s^%.2g", sym, e)
  parse(text=paste(Filter(nzchar, c("N[ch]", f(eS,"tilde(S)"), f(eD,"hat(Delta)"))), collapse="%.%")) }
top <- function(cn, lab){ d <- all |> filter(comp==cn); t <- tab |> filter(comp==cn)
  ggplot(d, aes(xb, ebar, colour=dlab)) + geom_hline(yintercept=0, colour="grey40", linewidth=.3) + geom_point(size=.8, alpha=.7) +
    scale_x_log10() + scale_colour_manual(values=COLd, name=expression(hat(Delta))) + coord_cartesian(ylim=c(-0.15,0.15)) +
    labs(title=paste0(lab, "   ", TT[[cn]]), x=xlab_of(t$eS, t$eD), y=parse(text=YM[[cn]])) + thm }
bot <- function(cn, lab){ d <- all |> filter(comp==cn); t <- tab |> filter(comp==cn); xg <- 10^seq(log10(min(d$xb)), log10(max(d$xb)), length=200)
  ggplot(d, aes(xb, s, colour=dlab)) + geom_hline(yintercept=t$floor, colour="grey70", linewidth=.3, linetype="dotted") +
    geom_line(data=tibble(x=xg, y=pmax(t$a*xg^t$slope, 1.5*t$floor)), aes(x,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
    geom_point(size=.8, alpha=.7) + scale_x_log10() +
    scale_y_log10(breaks=c(0.003,0.01,0.03,0.1,0.3), labels=c("0.003","0.01","0.03","0.1","0.3"), limits=c(0.0015,0.4)) +
    scale_colour_manual(values=COLd, guide="none") +
    annotate("text", x=min(d$xb)*2, y=0.33, label=sprintf("%.2f*x^{-1/%d}~~(R^2==%.2f)", t$a, round(-1/t$slope), t$R2_1d), parse=TRUE, size=2.3, colour="grey30", hjust=0) +
    labs(title=lab, x=xlab_of(t$eS, t$eD), y=parse(text=YS[[cn]])) + thm }
g <- wrap_plots(top("bias","A"), top("total","B"), top("sample","C"), top("corr","D"), nrow=1) /
     wrap_plots(bot("bias","E"), bot("total","F"), bot("sample","G"), bot("corr","H"), nrow=1) + plot_layout(guides="collect") & theme(legend.position="right")
ggsave("figures/Figure_IR_8panel_round.pdf", g, width=9.0, height=4.8); ggsave("figures/Figure_IR_8panel_round.png", g, width=9.0, height=4.8, dpi=150)
cat("written figures/Figure_IR_8panel_round.{pdf,png}\n")
