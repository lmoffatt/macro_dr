setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
# Bias magnitude + anisotropy, mirroring the distortion panels.
# A: Mahalanobis |bias| (joint sigmas) vs S~*Nch (pair product; Delta-independent)
# B: share of |bias|^2 by covariance eigen-rank, real cells vs null reference
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork)})
d <- read.csv("digest/ir_bias_direction.csv") |> mutate(S=0.1*z, x1=S*nch,
   dlab=factor(sprintf("%.2g",interval), levels=sprintf("%.2g",sort(unique(interval)))),
   real = mag > 2*floor)
r <- read.csv("digest/ir_bias_rank.csv") |>
  inner_join(d |> select(nch,z,interval,floor,real), by=c("nch","z","interval"))

cat("median share by eigen-rank (1=softest), real vs null (isotropy 0.167):\n")
tb <- r |> group_by(real,rank) |> summarise(med=round(median(share),3), .groups="drop") |>
      pivot_wider(names_from=rank, values_from=med); print(as.data.frame(tb), row.names=FALSE)

s <- d |> filter(real)
aENV <- max(s$mag * s$x1^(1/3)); cat(sprintf("envelope constant a (slope -1/3): %.3f\n", aENV))
EM <- function(x) pmax(aENV*x^(-1/3), 0.055)
cat("violations mag > envelope:", sum(d$mag > EM(d$x1)), "of", nrow(d), "\n")

COLd <- setNames(colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#F4A582","#D6604D","#B2182B","#67001F"))(nlevels(d$dlab)), levels(d$dlab))
thm <- theme_bw(base_size=8, base_family="Helvetica") +
  theme(plot.title=element_text(size=8, face="bold"), legend.key.height=unit(0.32,"cm"),
        panel.grid.minor=element_blank())
xg <- 10^seq(log10(min(d$x1)), log10(max(d$x1)), length=200)

pA <- ggplot(d, aes(x1, mag, colour=dlab, alpha=real)) +
  geom_line(data=tibble(x1=xg, y=EM(xg)), aes(x1,y), inherit.aes=FALSE,
            linetype="dashed", colour="grey30", linewidth=.4) +
  geom_hline(yintercept=median(d$floor), colour="grey70", linetype="dotted", linewidth=.35) +
  geom_point(size=.9) +
  scale_alpha_manual(values=c(`TRUE`=.85, `FALSE`=.25), guide="none") +
  scale_x_log10(breaks=10^seq(-2,10,2), labels=parse(text=sprintf("10^%d",seq(-2,10,2)))) +
  scale_y_log10(breaks=c(0.02,0.05,0.1,0.2,0.4), labels=c("0.02","0.05","0.1","0.2","0.4")) +
  scale_colour_manual(values=COLd, name=expression(Delta%.%k[off])) +
  annotate("text", x=1e8, y=0.0235, label="ensemble noise floor", size=2.2, colour="grey45", hjust=1, vjust=-0.4) +
  annotate("text", x=3e-2, y=0.42, label="slope -1/3", size=2.4, colour="grey30", hjust=0) +
  labs(title="A   bias magnitude, pair product, no interval dependence",
       x=expression(tilde(S)%.%N[ch]), y=expression("|bias|"[M]~~"(joint sigmas)")) + thm

rr <- r |> mutate(grp=ifelse(real,"real bias (114 cells)","at the noise floor (446)"))
pB <- ggplot(rr, aes(factor(rank), share, fill=grp)) +
  geom_hline(yintercept=1/6, colour="grey40", linetype="dashed", linewidth=.35) +
  geom_boxplot(outlier.size=.4, outlier.alpha=.3, linewidth=.25, width=.7,
               position=position_dodge(.8)) +
  scale_fill_manual(values=c(`real bias (114 cells)`="#B2182B", `at the noise floor (446)`="grey75"), name=NULL) +
  scale_x_discrete(labels=c("1\nsoftest","2","3","4","5","6\nstiffest")) +
  coord_cartesian(ylim=c(0,0.85)) +
  annotate("text", x=6.2, y=1/6, label="isotropy 1/6", size=2.2, colour="grey40", hjust=1, vjust=-0.5) +
  labs(title="B   bias direction: avoids soft and stiff",
       x="covariance eigen-rank", y=expression(share~of~"|bias|"^2)) +
  thm + theme(legend.position="top", legend.margin=margin(0,0,-4,0))

g <- pA + pB + plot_layout(widths=c(1.15,1))
ggsave("figures/Figure_IR_bias_maganiso.pdf", g, width=7.0, height=2.9)
ggsave("figures/Figure_IR_bias_maganiso.png", g, width=7.0, height=2.9, dpi=170)
cat("written figures/Figure_IR_bias_maganiso.{pdf,png}\n")
