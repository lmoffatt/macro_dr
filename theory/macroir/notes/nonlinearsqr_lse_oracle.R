# nonlinearsqr LSE — R validation oracle for scheme_CO (2-state C<->O)
# =====================================================================
# Companion to nonlinearsqr_lse_plan.md + nonlinearsqr_cpp_spec.md.
#
# Purpose: an INDEPENDENT reference for the LSE driver's mean trajectory,
# logL, score and Fisher, so the C++ can be validated on ONE recording at a
# fixed theta BEFORE running any figure. Two layers (see the plan's test ladder):
#
#   (A) finalize check  -- logL/score/Fisher recomputed from residuals; needs
#                          only mu_i, dmu_i, y_i. Fully independent of Q/timing.
#   (B) mean check      -- mu_i recomputed analytically from theta + the
#                          per-interval (dt, agonist) table. Cross-checks the
#                          conductance-averaging.
#   (C) AD-trap check   -- the Fisher must be the GN block (n/SSE) J^T J, NOT
#                          the true Hessian. They DIFFER off the optimum; this
#                          file reports both so a driver that used the 2nd-order
#                          AD Hessian is caught. Evaluate at a NON-optimal theta.
#
# scheme_CO (legacy/models_simple.h:10): state 0 = C (closed), 1 = O (open).
#   Q[0,1] = agonist * on ;  Q[1,0] = off ;  g_open = -unitary_current .
#   mu(t) = N * g_open * P_open(t) + baseline .
#   P_open(t) = Pinf + (P0 - Pinf) exp(-t/tau),  Pinf = kon/(kon+koff),
#     tau = 1/(kon+koff),  kon = agonist*on_nat,  koff = off_nat .
#   P_initial = [C=1, O=0]  ->  P_open(0) = 0 .
#   All 6 params are Log10-transformed: natural = 10^theta.
#   Order: on, off, unitary, curr_noise, baseline, Nch .  curr_noise does NOT
#   enter the mean. Per the identifiability fixes the FREE set for the fit is
#   {on, off, baseline, Nch}; unitary + curr_noise are Fixed.
#
# Averaging (verified qmodel.h:4367-4395):
#   av=1  interval-AVERAGED conductance -> mu = N*g_open*<P_open>_interval + base
#   av=0  INSTANTANEOUS at interval midpoint (P_half) -> mu = N*g_open*P_mid + base
#   av=2  gives the SAME mean as av=1  ->  excluded from the LSE.
#
# LSE finalize (plan §1; Jeffreys-marginalized sigma^2):
#   SSE = sum (y-mu)^2 over finite y ;  n = count finite y ;  sigma2_hat = SSE/n
#   logL   = lgamma(n/2) - (n/2) log(pi) - (n/2) log(SSE)
#   score_a  = (n/SSE) * sum_i (y_i-mu_i) dmu_i/dtheta_a        [ = d logL/dtheta_a ]
#   Fisher_ab= (n/SSE) * sum_i (dmu_i/dtheta_a)(dmu_i/dtheta_b) [ GN block ]
#
# Usage:
#   Rscript nonlinearsqr_lse_oracle.R          # runs self-test + a demo
#   source("nonlinearsqr_lse_oracle.R")        # loads the functions only
#   compare_driver("figure_..._replicates.csv", theta, av=1, ...)  # vs a dump

# ---- model ----------------------------------------------------------------

# Natural parameter value from the Log10-transformed theta entry.
nat <- function(theta, name) 10^unname(theta[[name]])

# P_open trajectory over a table of intervals (each row: dt, agonist).
# Returns, per interval, P_open at the start, the midpoint, the time-average
# over the interval, and the end (= next interval's start). P_open(0)=0.
p_open_trajectory <- function(theta, intervals) {
  on  <- nat(theta, "on")
  off <- nat(theta, "off")
  n   <- nrow(intervals)
  Pstart <- Pmid <- Pavg <- Pend <- numeric(n)
  cur <- 0.0                              # P_initial = [C=1, O=0] -> P_open=0
  for (i in seq_len(n)) {
    A  <- intervals$agonist[i]
    dt <- intervals$dt[i]
    kon <- A * on; koff <- off; ksum <- kon + koff
    if (ksum > 0) { Pinf <- kon / ksum; tau <- 1 / ksum } else { Pinf <- cur; tau <- Inf }
    Pstart[i] <- cur
    ed   <- if (is.finite(tau)) exp(-dt / tau)       else 1
    emid <- if (is.finite(tau)) exp(-(dt / 2) / tau) else 1
    Pend[i] <- Pinf + (cur - Pinf) * ed
    Pmid[i] <- Pinf + (cur - Pinf) * emid
    Pavg[i] <- if (is.finite(tau)) Pinf + (cur - Pinf) * (tau / dt) * (1 - ed) else cur
    cur <- Pend[i]
  }
  data.frame(Pstart, Pmid, Pavg, Pend)
}

# Analytic LSE mean mu_i for av in {0,1}.
mu_analytic <- function(theta, intervals, av = 1) {
  g_open <- -nat(theta, "unitary")
  N      <-  nat(theta, "Nch")
  base   <-  nat(theta, "baseline")
  tr     <- p_open_trajectory(theta, intervals)
  Popen  <- if (av == 1) tr$Pavg else if (av == 0) tr$Pmid else
            stop("av must be 0 or 1 (av=2 has the same mean as av=1)")
  as.numeric(N * g_open * Popen + base)
}

DEFAULT_FREE <- c("on", "off", "baseline", "Nch")

# Central-difference Jacobian J[i,a] = dmu_i/dtheta_a over the free params.
dmu_dtheta <- function(theta, intervals, av = 1, free = DEFAULT_FREE, h = 1e-6) {
  mu0 <- mu_analytic(theta, intervals, av)
  J <- matrix(0, length(mu0), length(free), dimnames = list(NULL, free))
  for (a in free) {
    tp <- theta; tp[[a]] <- tp[[a]] + h
    tm <- theta; tm[[a]] <- tm[[a]] - h
    J[, a] <- (mu_analytic(tp, intervals, av) - mu_analytic(tm, intervals, av)) / (2 * h)
  }
  J
}

# ---- finalize (independent of the mean) -----------------------------------

# logL/score/Fisher from mu, its Jacobian J, and y. NaN y are dropped (n counts
# only finite y). This is the pure §1 finalize -- the highest-value check.
lse_stats <- function(mu, J, y) {
  ok  <- is.finite(y)
  r   <- (y - mu)[ok]
  Jok <- J[ok, , drop = FALSE]
  n   <- length(r)
  SSE <- sum(r^2)
  list(
    n      = n,
    SSE    = SSE,
    sigma2 = SSE / n,
    logL   = lgamma(n / 2) - (n / 2) * log(pi) - (n / 2) * log(SSE),
    score  = as.numeric((n / SSE) * crossprod(Jok, r)),  # (n/SSE) J^T r
    fisher = (n / SSE) * crossprod(Jok)                  # (n/SSE) J^T J  (GN block)
  )
}

# Scalar logL at theta (for the FD-of-value score check).
lse_loglik_value <- function(theta, y, intervals, av = 1) {
  mu  <- mu_analytic(theta, intervals, av)
  ok  <- is.finite(y); r <- (y - mu)[ok]; n <- length(r); SSE <- sum(r^2)
  lgamma(n / 2) - (n / 2) * log(pi) - (n / 2) * log(SSE)
}

# Score by finite-differencing the logL VALUE (independent of the analytic
# score assembly). Must match lse_stats()$score.
fd_score <- function(theta, y, intervals, av = 1, free = DEFAULT_FREE, h = 1e-6) {
  vapply(free, function(a) {
    tp <- theta; tp[[a]] <- tp[[a]] + h
    tm <- theta; tm[[a]] <- tm[[a]] - h
    (lse_loglik_value(tp, y, intervals, av) - lse_loglik_value(tm, y, intervals, av)) / (2 * h)
  }, numeric(1))
}

# Observed information -d^2 logL/dtheta^2 by FD of the analytic score. This is
# the TRUE Hessian; it EQUALS the GN Fisher only at the optimum. Used to show
# the AD-trap check discriminates (GN Fisher vs Hessian differ off-optimum).
fd_neg_hessian <- function(theta, y, intervals, av = 1, free = DEFAULT_FREE, h = 1e-5) {
  score_at <- function(th) {
    mu <- mu_analytic(th, intervals, av)
    J  <- dmu_dtheta(th, intervals, av, free)
    lse_stats(mu, J, y)$score
  }
  p <- length(free); H <- matrix(0, p, p, dimnames = list(free, free))
  for (k in seq_along(free)) {
    tp <- theta; tp[[free[k]]] <- tp[[free[k]]] + h
    tm <- theta; tm[[free[k]]] <- tm[[free[k]]] - h
    H[, k] <- (score_at(tp) - score_at(tm)) / (2 * h)
  }
  -0.5 * (H + t(H))   # symmetrize; observed information = -d score/dtheta
}

# One-shot: everything the driver reports, from theta + y + intervals.
oracle <- function(theta, y, intervals, av = 1, free = DEFAULT_FREE) {
  mu <- mu_analytic(theta, intervals, av)
  J  <- dmu_dtheta(theta, intervals, av, free)
  s  <- lse_stats(mu, J, y)
  s$mu <- mu; s$J <- J; s
}

# ---- self-test (validates THIS oracle, no driver needed) ------------------

# A small scheme_CO experiment: pre 0uM, agonist 10uM, post 0uM. dt chosen so
# tau (=5 ms at 10uM with on=10,off=100) is well resolved. Timing here is only
# for the self-test; the driver comparison reads (dt,agonist) from the dump.
demo_intervals <- function(n_per_seg = 40, dt = 5e-4) {
  seg <- function(A) data.frame(dt = rep(dt, n_per_seg), agonist = rep(A, n_per_seg))
  rbind(seg(0.0), seg(10.0), seg(0.0))
}

demo_theta <- function() c(on = 1, off = 2, unitary = 0, curr_noise = -4, baseline = 0, Nch = log10(5000))
# natural: on=10, off=100, unitary=1, curr_noise=1e-4, baseline=1, Nch=5000.

self_test <- function(tol_grad = 1e-4, tol_disc = 1e-6, seed = 1) {
  set.seed(seed)
  theta <- demo_theta()
  iv    <- demo_intervals()
  free  <- DEFAULT_FREE

  mu_true <- mu_analytic(theta, iv, av = 1)
  # emission sd ~ current-noise scale; small relative to the signal.
  y <- mu_true + rnorm(length(mu_true), 0, 0.02 * sd(mu_true))
  y[1:5] <- NaN                                  # pre-agonist NaN prefix

  # Evaluate at a NON-optimal theta (perturb away from the truth).
  theta_off <- theta; theta_off["on"] <- theta_off["on"] + 0.15
  theta_off["off"] <- theta_off["off"] - 0.10

  ok <- TRUE
  say <- function(name, pass, detail = "") {
    cat(sprintf("  [%s] %-34s %s\n", if (pass) "PASS" else "FAIL", name, detail))
    ok <<- ok && pass
  }

  # 1. NaN prefix dropped, n = count of finite y.
  s <- oracle(theta_off, y, iv, av = 1, free = free)
  say("n = count(finite y)", s$n == sum(is.finite(y)), sprintf("n=%d", s$n))

  # 2. logL matches the closed form recomputed here (sanity of lse_stats).
  ok_y <- is.finite(y); r <- (y - s$mu)[ok_y]; n <- length(r); SSE <- sum(r^2)
  logL_ref <- lgamma(n/2) - (n/2)*log(pi) - (n/2)*log(SSE)
  say("logL == closed form", abs(s$logL - logL_ref) < 1e-10)

  # 3. analytic score == FD of the logL value (gradient consistency), off-optimum.
  g_fd <- fd_score(theta_off, y, iv, av = 1, free = free)
  rel  <- max(abs(s$score - g_fd) / (abs(g_fd) + 1e-8))
  say("score == FD(logL value)", rel < tol_grad, sprintf("max rel %.2e", rel))

  # 4. AD-trap discrimination: GN Fisher != observed information OFF the optimum.
  H <- fd_neg_hessian(theta_off, y, iv, av = 1, free = free)
  gap <- max(abs(s$fisher - H)) / max(abs(s$fisher))
  say("GN Fisher != Hessian (off-opt)", gap > tol_disc,
      sprintf("rel gap %.2e -> the Fisher test is discriminating", gap))

  # 5. av=0 and av=1 give DIFFERENT means (both are real variants).
  d <- max(abs(mu_analytic(theta, iv, 0) - mu_analytic(theta, iv, 1)))
  say("av=0 mean != av=1 mean", d > 0, sprintf("max |dmu| %.3g", d))

  cat(if (ok) "self_test: ALL PASS\n" else "self_test: FAILURES ABOVE\n")
  invisible(ok)
}

# ---- driver comparison ----------------------------------------------------

# Compare the driver's per-interval dump against the oracle.
#   csv          : path to the driver dump (one row per interval).
#   theta        : named Log10 theta the driver was evaluated at.
#   av           : 0 or 1, matching the driver run.
#   cols         : column names in the dump. Adjust to the real schema once the
#                  driver emits it. Required: y (obs), dt, agonist. Optional:
#                  y_mean (driver mu, to cross-check the analytic mean).
#   driver_logL, driver_score, driver_fisher : the driver's reported scalars/
#                  vectors/matrix, if available, to diff against the oracle.
compare_driver <- function(csv, theta, av = 1, free = DEFAULT_FREE,
                           cols = list(y = "y", dt = "dt", agonist = "agonist",
                                       y_mean = "y_mean"),
                           driver_logL = NULL, driver_score = NULL,
                           driver_fisher = NULL, skip = 0) {
  d  <- utils::read.csv(csv, skip = skip, check.names = FALSE)
  y  <- d[[cols$y]]
  iv <- data.frame(dt = d[[cols$dt]], agonist = d[[cols$agonist]])
  s  <- oracle(theta, y, iv, av = av, free = free)

  cat(sprintf("intervals=%d  n(finite y)=%d  SSE=%.6g  sigma2=%.6g  logL=%.8g\n",
              nrow(d), s$n, s$SSE, s$sigma2, s$logL))

  if (!is.null(cols$y_mean) && cols$y_mean %in% names(d)) {
    dm <- max(abs(s$mu - d[[cols$y_mean]]))
    cat(sprintf("mean check   : max |mu_oracle - y_mean_driver| = %.3e\n", dm))
  }
  if (!is.null(driver_logL))
    cat(sprintf("logL   diff  : %.3e\n", abs(s$logL - driver_logL)))
  if (!is.null(driver_score))
    cat(sprintf("score  diff  : max |.| = %.3e (rel %.2e)\n",
                max(abs(s$score - driver_score)),
                max(abs(s$score - driver_score) / (abs(driver_score) + 1e-8))))
  if (!is.null(driver_fisher))
    cat(sprintf("Fisher diff  : max |.| = %.3e (rel %.2e)\n",
                max(abs(s$fisher - driver_fisher)),
                max(abs(s$fisher - driver_fisher) / (max(abs(driver_fisher)) + 1e-30))))
  invisible(s)
}

# ---- main -----------------------------------------------------------------

if (identical(environment(), globalenv()) && !interactive()) {
  cat("== nonlinearsqr LSE oracle: self-test ==\n")
  self_test()
  cat("\n== demo (scheme_CO, theta_sim, av=1) ==\n")
  th <- demo_theta(); iv <- demo_intervals()
  mu <- mu_analytic(th, iv, av = 1)
  cat(sprintf("n_intervals=%d  mu range [%.3f, %.3f]  (N=%.0f, g_open=%.0f, base=%.0f)\n",
              length(mu), min(mu), max(mu),
              nat(th, "Nch"), -nat(th, "unitary"), nat(th, "baseline")))
}
