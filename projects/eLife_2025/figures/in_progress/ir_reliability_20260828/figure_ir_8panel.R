# 8 panels: columns bias | total distortion | sample distortion | correlation distortion; rows: mean (vbar / ebar) and sd (s_v / s). All vs x, colour Dhat.
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
comp <- read.csv("digest/ir_components_cells.csv")
bias <- read.csv("digest/ir_bias_whitened.csv") |> group_by(nch,z,interval) |>
  summarise(ebar=mean(v), s=sqrt(mean((v-mean(v))^2)), .groups="drop") |> mutate(comp="bias")
all <- bind_rows(bias, comp |> select(nch,z,interval,comp,ebar,s)) |>
  mutate(St=0.1*z, x=nch*St*sqrt(interval), dlab=factor(sprintf("%.2g", interval), levels=sprintf("%.2g", sort(unique(interval)))),
         comp=factor(comp, levels=c("bias","total","sample","corr")))
# floors from the flat region (nch>=5000 | z>=100): 90th percentile of s there
fl <- all |> filter(x > 1e4) |> group_by(comp) |> summarise(floor=max(s), .groups="drop")
print(fl |> mutate(floor=round(floor,3)) |> as.data.frame(), row.names=FALSE)
# envelopes at slope -1/5 vs x, zero violations above floor, per component (for reference lines)
env <- all |> inner_join(fl, by="comp") |> filter(s > floor) |> group_by(comp) |> summarise(a=ceiling(100*max(s*x^0.2))/100, .groups="drop") |> inner_join(fl, by="comp")
print(env |> mutate(floor=round(floor,3)) |> as.data.frame(), row.names=FALSE)
# also: interval stratification residual per component (R2 gain from adding Dhat to the 1-D fit vs x)
for(cn in levels(all$comp)){ s2 <- all |> inner_join(fl, by="comp") |> filter(comp==cn, s > 1.3*floor)
  m1 <- lm(log10(s) ~ log10(x), data=s2); m2 <- lm(log10(s) ~ log10(x)+dlab, data=s2)
  cat(sprintf("  %-6s vs x: slope %+.3f R2 %.3f ; with Dhat factor R2 %.3f (gain %.3f)  n=%d\n", cn, coef(m1)[2], summary(m1)$r.squared, summary(m2)$r.squared, summary(m2)$r.squared-summary(m1)$r.squared, nrow(s2))) }

COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(all$dlab)), levels(all$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") + theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank())
XBRK <- 10^seq(-2,10,4); XLABS <- parse(text=sprintf("10^%d", seq(-2,10,4))); XL <- expression(N[ch]%.%tilde(S)%.%hat(Delta)^0.5)
TT <- c(bias="bias  (v = G^{1/2} b)", total="total distortion", sample="sample distortion", corr="correlation distortion")
YM <- c(bias="mean~v[p]", total="mean~log~lambda", sample="mean~log~lambda", corr="mean~log~lambda")
YS <- c(bias="sd~v[p]", total="sd~log~lambda", sample="sd~log~lambda", corr="sd~log~lambda")
xg <- 10^seq(log10(min(all$x)), log10(max(all$x)), length=200)
top <- function(cn, lab){ d <- all |> filter(comp==cn)
  ggplot(d, aes(x, ebar, colour=dlab)) + geom_hline(yintercept=0, colour="grey40", linewidth=.3) + geom_point(size=.8, alpha=.7) +
    scale_x_log10(breaks=XBRK, labels=XLABS) + scale_colour_manual(values=COLd, name=expression(hat(Delta))) + coord_cartesian(ylim=c(-0.15,0.15)) +
    labs(title=paste0(lab, "   ", TT[[cn]]), x=NULL, y=parse(text=YM[[cn]])) + thm }
bot <- function(cn, lab){ d <- all |> filter(comp==cn); e <- env |> filter(comp==cn)
  p <- ggplot(d, aes(x, s, colour=dlab)) + geom_hline(yintercept=e$floor, colour="grey70", linewidth=.3, linetype="dotted") + geom_point(size=.8, alpha=.7)
  if(cn %in% c("bias","total")) p <- p + geom_line(data=tibble(x=xg, y=pmax(e$a*xg^-0.2, e$floor)), aes(x,y), inherit.aes=FALSE, linetype="dashed", colour="grey30", linewidth=.4) +
    annotate("text", x=min(d$x)*3, y=0.33, label=sprintf("%.2f*x^{-1/5}", e$a), parse=TRUE, size=2.3, colour="grey30", hjust=0)
  p +
    scale_x_log10(breaks=XBRK, labels=XLABS) + scale_y_log10(breaks=c(0.003,0.01,0.03,0.1,0.3), labels=c("0.003","0.01","0.03","0.1","0.3"), limits=c(0.0015,0.4)) +
    scale_colour_manual(values=COLd, guide="none") +
    labs(title=lab, x=XL, y=parse(text=YS[[cn]])) + thm }
g <- wrap_plots(top("bias","A"), top("total","B"), top("sample","C"), top("corr","D"), nrow=1) /
     wrap_plots(bot("bias","E"), bot("total","F"), bot("sample","G"), bot("corr","H"), nrow=1) + plot_layout(guides="collect") & theme(legend.position="right")
ggsave("figures/Figure_IR_8panel.pdf", g, width=9.0, height=4.8); ggsave("figures/Figure_IR_8panel.png", g, width=9.0, height=4.8, dpi=150)
cat("written figures/Figure_IR_8panel.{pdf,png}\n")
