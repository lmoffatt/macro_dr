# fig3_claims.R -- checks on the panel-by-panel reading of Figure 3, 2026-08-12. Feeds the numbers
# cited in the residual and correlation paragraphs of 04_results.tex: the mean squared standardized
# residual by segment (0.010 before the agonist, 1.91 over the pulse, 0.69 through the washout for
# either least-squares arm), r^2 at the stepping interval, and the residual autocorrelation at lag
# one.
#
# Run from the repository root:  Rscript papers/1_method/decisions/recompute/fig3_claims.R
# PROMOTED 2026-08-26 from tmp/, which is gitignored; only the paths changed.
# Checks on Luciano's panel-by-panel description of Figure 3, 2026-08-12.
suppressPackageStartupMessages({library(data.table)})
DTOK <- c(LSE="LSE_av0", ILSE="LSE", NR="NR", INR="INR", R="R", MR="MR", VR="VR", IR="IR")
PSEL <- c(k_off=1L, N_ch=5L); FMIN <- 1e-6
wide <- function(d,v) as.matrix(dcast(d, simulation_index ~ sample_index, value.var=v)[,-1])
colVar <- function(M){n<-nrow(M);m<-colMeans(M);(colMeans(M*M)-m^2)*n/(n-1)}

cat("=== B: time course of r^2. pulse = samples 21-60, washout = 61-100 ===\n")
cat(sprintf("%-5s %8s %8s %8s %8s %8s\n","key","pre(1-20)","early(21-30)","late(31-60)","wash","overall"))
for (k in names(DTOK)) {
  s <- as.data.table(readRDS(sprintf("projects/eLife_2025/figures/data/digest/figure_3_digest_%s.rds", DTOK[[k]]))$scalar)
  s <- s[is.finite(r_std)][, r2 := r_std^2]
  f <- function(a,b) s[sample_index %in% a:b, mean(r2)]
  cat(sprintf("%-5s %8.3f %8.3f %8.3f %8.3f %8.3f\n", k, f(1,20), f(21,30), f(31,60), f(61,100), mean(s$r2)))
}

cat("\n=== C: after the pulse, is F_t(k_off) proportional to the mean current? ===\n")
for (k in c("NR","INR","R","IR")) {
  dg <- readRDS(sprintf("projects/eLife_2025/figures/data/digest/figure_3_digest_%s.rds", DTOK[[k]]))
  P <- as.data.table(dg$param)[param_index %in% PSEL]
  S <- as.data.table(dg$scalar)
  cur <- S[, .(y = mean(y_mean)), by = sample_index]
  Fk <- P[param_index==1L, .(F = mean(I)), by = sample_index]
  Fn <- P[param_index==5L, .(F = mean(I)), by = sample_index]
  d <- merge(merge(cur, Fk, by="sample_index"), Fn, by="sample_index", suffixes=c("_koff","_Nch"))
  w <- d[sample_index %in% 61:100 & F_koff > FMIN]
  cat(sprintf("%-4s washout: cor(log F_koff, log |y|) = %+.3f over %d steps | slope = %+.2f | F_Nch drops %.1e -> %.1e (steps 61 -> 70)\n",
      k, cor(log(w$F_koff), log(abs(w$y))), nrow(w),
      coef(lm(log(F_koff) ~ log(abs(y)), data = w))[2],
      d[sample_index==61, F_Nch], d[sample_index==70, F_Nch]))
}

cat("\n=== E: does the per-interval ratio Var(s_t)/F_t rise during the decay? (k_off) ===\n")
cat(sprintf("%-5s %10s %10s %10s\n","key","pulse","early wash","late wash"))
for (k in names(DTOK)) {
  P <- as.data.table(readRDS(sprintf("projects/eLife_2025/figures/data/digest/figure_3_digest_%s.rds", DTOK[[k]]))$param)[param_index==1L]
  Sm <- wide(P,"s"); Im <- wide(P,"I"); Fb <- colMeans(Im); r <- colVar(Sm)/Fb
  smp <- as.integer(colnames(Sm)); ok <- is.finite(r) & Fb > FMIN
  m <- function(a,b) mean(r[ok & smp>=a & smp<=b])
  cat(sprintf("%-5s %10.3f %10.3f %10.3f\n", k, m(21,60), m(61,80), m(81,100)))
}

cat("\n=== G: shape of the score ACF (k_off): lag-1, and the lag where it first falls below 0.1 ===\n")
acf_one <- function(x) as.numeric(acf(x, lag.max=25, plot=FALSE)$acf)
for (k in names(DTOK)) {
  P <- as.data.table(readRDS(sprintf("projects/eLife_2025/figures/data/digest/figure_3_digest_%s.rds", DTOK[[k]]))$param)[param_index==1L]
  A <- colMeans(do.call(rbind, P[, .(a=list(acf_one(s))), by=simulation_index]$a), na.rm=TRUE)
  below <- which(A[-1] < 0.1)[1]
  cat(sprintf("%-5s acf1 %+.3f  acf2 %+.3f  acf5 %+.3f  acf10 %+.3f  | first lag < 0.1: %s | sum_{1..25} %+.2f\n",
      k, A[2], A[3], A[6], A[11], ifelse(is.na(below),">25",below), sum(A[2:26])))
}
