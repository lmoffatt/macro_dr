# Third witness: standardized-residual second moment r2bar = Sum_r2_std / n, n = 10/interval (null 1).
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(readr);library(ggplot2);library(patchwork)})
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/data")
DIRS <- c("1c2ae6f","0ffbda7","87889e6","82b956f")
files <- unlist(lapply(DIRS, function(d) list.files(d, pattern="macro_IR.*battery_sim_G\\.csv$", full.names=TRUE)))
meta <- tibble(f=files) |> mutate(nch=as.numeric(sub(".*nch_([0-9]+)_.*","\\1",f)), z=as.numeric(sub(".*noise_([0-9.]+)_battery.*","\\1",f))) |>
  group_by(nch,z) |> slice(1) |> ungroup()
one <- function(f,nch,z){
  d <- suppressWarnings(read_csv(f, skip=1, show_col_types=FALSE, progress=FALSE)) |>
    filter(component_path %in% c("Probit_statistics_Moment_statistics_Sum_r2_std","Probit_statistics_Moment_statistics_Sum_r_std"),
           statistic=="mean", probit %in% c("mean","quantile")) |>
    mutate(key=ifelse(probit=="mean","m", ifelse(quantile_level<0.5,"lo","hi"))) |> filter(key=="m" | quantile_level %in% c(0.025,0.975)) |>
    transmute(interval=interval_in_tau, comp=ifelse(grepl("r2",component_path),"r2","r1"), key, value) |>
    pivot_wider(names_from=c(comp,key), values_from=value) |> mutate(nch=nch, z=z)
  d
}
r <- pmap_dfr(meta, function(f,nch,z) one(f,nch,z)) |> mutate(n=10/interval, r2=r2_m/n-1, r2lo=r2_lo/n-1, r2hi=r2_hi/n-1, r1=r1_m/n)
setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828"); write_csv(r, "digest/ir_r2std_cells.csv")
w <- read.csv("digest/ir_bias_whitened.csv") |> group_by(nch,z,interval) |> summarise(vsum=sum(v), .groups="drop")
d <- read.csv("digest/ir_collapse_cells.csv") |> filter(anchor=="sim") |> select(nch,z,interval,ebar,aniso)
cell <- r |> inner_join(w, by=c("nch","z","interval")) |> inner_join(d, by=c("nch","z","interval")) |>
  mutate(St=0.1*z, x=nch*St*sqrt(interval), Dcl=cut(interval, c(0,0.05,0.2,1.01), labels=c("short","mid","long")),
         dlab=factor(sprintf("%.2g", interval), levels=sprintf("%.2g", sort(unique(interval)))))
cat("cells:", nrow(cell), "\n")
cat("r2bar-1 quantiles:"); print(round(quantile(cell$r2, c(0,.05,.5,.95,1)),4))
cat("\ncorrelations (all cells):  cor(ebar, r2-1) =", round(cor(cell$ebar, cell$r2),3), "  cor(vsum, r2-1) =", round(cor(cell$vsum, cell$r2),3), "  cor(ebar, vsum) =", round(cor(cell$ebar, cell$vsum),3), "\n")
big <- cell |> filter(x < 1)
cat("corner x<1 (n=", nrow(big), "): cor(ebar, r2-1) =", round(cor(big$ebar, big$r2),3), "  cor(vsum, r2-1) =", round(cor(big$vsum, big$r2),3), "\n")
cat("\nCORNER x<1, by interval class: sign shares and medians\n")
print(big |> group_by(Dcl) |> summarise(n=n(), ebar_pos=round(mean(ebar>0),2), vsum_pos=round(mean(vsum>0),2), r2_pos=round(mean(r2>0),2),
  ebar_med=round(median(ebar),3), vsum_med=round(median(vsum),3), r2_med=round(median(r2),4), r2_sig=round(mean(r2lo>0 | r2hi<0),2)) |> as.data.frame(), row.names=FALSE)
cat("\nsign agreement in corner: ebar~r2", round(mean(sign(big$ebar)==sign(big$r2)),2), " vsum~r2", round(mean(sign(big$vsum)==sign(big$r2)),2), " ebar~vsum", round(mean(sign(big$ebar)==sign(big$vsum)),2), "\n")
cat("restricted to |ebar|>0.03:\n"); b2 <- big |> filter(abs(ebar)>0.03)
cat("  n=", nrow(b2), " ebar~r2", round(mean(sign(b2$ebar)==sign(b2$r2)),2), " vsum~r2", round(mean(sign(b2$vsum)==sign(b2$r2)),2), "\n")

COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(cell$dlab)), levels(cell$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") + theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.3,"cm"), panel.grid.minor=element_blank())
XBRK <- 10^seq(-2,10,4); XLABS <- parse(text=sprintf("10^%d", seq(-2,10,4))); XL <- expression(N[ch]%.%tilde(S)%.%hat(Delta)^0.5)
pp <- function(y, ttl, yl, lim) ggplot(cell, aes(x, .data[[y]], colour=dlab)) + geom_hline(yintercept=0, colour="grey40", linewidth=.3) +
  geom_point(size=.9, alpha=.75) + scale_x_log10(breaks=XBRK, labels=XLABS) + scale_colour_manual(values=COLd, name=expression(hat(Delta))) +
  coord_cartesian(ylim=lim) + labs(title=ttl, x=XL, y=yl) + thm
g <- pp("ebar", "A   distortion: mean log lambda", expression(bar(e)), c(-0.16,0.16)) +
     pp("vsum", "B   bias: sum of whitened components", expression(sum(v[p])), c(-0.25,0.25)) +
     pp("r2", "C   residual variance misfit", expression(bar(r)[std]^2 - 1), c(-0.12,0.12)) +
     plot_layout(guides="collect") & theme(legend.position="right")
ggsave("figures/Figure_IR_sign_witness.pdf", g, width=7.4, height=2.6); ggsave("figures/Figure_IR_sign_witness.png", g, width=7.4, height=2.6, dpi=170)
cat("written figures/Figure_IR_sign_witness.{pdf,png}\n")
