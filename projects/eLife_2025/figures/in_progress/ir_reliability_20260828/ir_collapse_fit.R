setwd("/home/lmoffatt/Code/macro_dr/macro_dr/projects/eLife_2025/figures/in_progress/ir_reliability_20260828")
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
d <- read.csv("digest/ir_collapse_cells.csv")
for(a in c("sim","pool")){
 s <- d |> filter(anchor==a, is.finite(dAI))
 cat("\n=== anchor:",a," cells:",nrow(s),"  Nch:",paste(sort(unique(s$nch)),collapse=","),
     "\n    noise:",paste(sort(unique(s$z)),collapse=","),"\n")
 for(y in c("dAI","aniso","zmax")){
   s2 <- s |> filter(is.finite(.data[[y]]), .data[[y]]>0)
   m <- lm(log10(s2[[y]]) ~ log10(z)+log10(nch)+log10(interval), data=s2)
   co <- coef(m)
   cat(sprintf("  %-6s  exp(noise)=%+.3f exp(Nch)=%+.3f exp(int)=%+.3f  R2=%.3f  n=%d\n",
       y, co[2], co[3], co[4], summary(m)$r.squared, nrow(s2)))
   # collapse quality on 1 predictor, controlling interval as factor
   for(nm in c("prod","ratio","nch","noise")){
     x <- switch(nm, prod=log10(s2$z*s2$nch), ratio=log10(s2$z/s2$nch), nch=log10(s2$nch), noise=log10(s2$z))
     r2 <- summary(lm(log10(s2[[y]]) ~ x + factor(s2$interval)))$r.squared
     cat(sprintf("      x=%-5s  R2(with interval factor)=%.3f\n", nm, r2))
   }
 }
}
