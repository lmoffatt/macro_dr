#ifndef MACRODR_OPTIMIZATION_GAUSS_NEWTON_H
#define MACRODR_OPTIMIZATION_GAUSS_NEWTON_H

// Generic Gauss-Newton / Levenberg-Marquardt optimizer.
//
// Maximizes a smooth scalar objective whose evaluation returns a structured
// type with three named accessors:
//
//   get<logL>(state)()                 -> double                value
//   get<Grad>(state)().value()         -> Matrix<double>        gradient
//   get<CurvTag>(state)().value()      -> SymPosDefMatrix       curvature (PSD)
//
// `CurvTag` is the only template parameter the caller specifies; `logL` and
// `Grad` are hard-coded because the codebase uses a single canonical name for
// each. The curvature is templated because there are three mathematically
// distinct PSD curvature estimators in macro-dr:
//
//   - Gaussian_Fisher_Information           (= G_lik, the moment-matched FIM)
//   - Likelihood_Numerical_Fisher_Information (= F, the observed Hessian; can
//                                              be indefinite, requires strong
//                                              LM damping or projection)
//   - Score_Covariance_Matrix               (= J, empirical score covariance)
//
// For MLE with the Gaussian-moment-matched likelihood, instantiate with
// `<Gaussian_Fisher_Information>` (the slot populated by update_macro_state
// in the AD pipeline, accumulating G_lik per timestep). For MAP, a wrapper
// objective composes G_post = G_lik + H_prior before passing it back here.
//
// Step: standard LM with adaptive damping
//   (curv + lambda * diag(curv)) * delta = grad
//   theta_new = theta + delta
// Accept if value increased; otherwise grow lambda and retry from same theta.
// Shrink lambda on accept, grow on reject. Convergence on ||grad|| < tol or
// |dvalue| < tol.

#include "distributions.h"    // logL, Grad, Gaussian_Fisher_Information, Vector_Space, get<>
#include "lapack_headers.h"   // Lapack_SymmPosDef_inv
#include "matrix.h"           // Matrix, SymPosDefMatrix
#include "maybe_error.h"      // Maybe_error, error_message

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>
#include <utility>

namespace macrodr::optimization {





struct gauss_newton_options {
    // Damping policy: lambda walks along the geometric sequence
    //   0 ↔ kickoff ↔ kickoff·factor ↔ kickoff·factor² ↔ … ↔ lambda_max
    //
    // Loop starts at lambda = 0 (pure Newton).
    //   - On reject/failure: step UP in the sequence (0 → kickoff → kickoff·factor → …).
    //   - On accept:         step DOWN in the sequence (… → kickoff·factor → kickoff → 0).
    //
    // Symmetric — no "snap" asymmetry between accept and reject. Gradual —
    // factor=3 ≈ √10 gives ~2× finer control than factor=10 at similar
    // convergence speed; the sequence 0,1,3,9,27,… has roughly log₂ steps to
    // reach a target damping versus log₁₀ for factor=10.
    //
    // Intermediate values between 0 and kickoff (e.g. 0.001, 0.01, 0.1) are
    // skipped on purpose — they are "almost Newton" without adding meaningful
    // regularisation.
    double lambda_kickoff = 1.0;   // smallest non-zero value in the sequence;
                                   // configurable — start higher if you suspect
                                   // θ_init is far from the optimum
    double lambda_factor  = 3.0;   // geometric ratio between successive damped values
    double lambda_max     = 1e10;  // give-up bound

    std::size_t max_iter  = 100;

    // Convergence (PRIMARY): Newton decrement  ½·gᵀH⁻¹g < newton_dec2_tol.
    // ½·gᵀH⁻¹g is the quadratic-model estimate of the objective sub-optimality
    // f(θ) − f*, in objective (log-likelihood) units. It is SCALE-INVARIANT:
    // unlike ‖grad‖ — which scales ∝ N for an MLE, so a fixed grad_rtol becomes
    // effectively N× tighter at large N and can fall BELOW the gradient's
    // numerical floor (→ the optimiser churns near the optimum and escalates λ
    // to lambda_max instead of converging) — the Newton decrement's numerical
    // floor shrinks ∝ 1/N, so one fixed tolerance converges robustly across data
    // sizes. H is the (undamped) CurvTag curvature.
    double newton_dec2_tol = 1e-8;

    // Convergence (FALLBACK): ||grad|| < grad_rtol * max(1, ||theta||).
    // Kept for back-compat; at large N the Newton decrement above is preferred.
    double grad_rtol      = 1e-6;

    // Convergence: |dvalue| < dvalue_tol (absolute on the objective value)
    double dvalue_tol     = 1e-10;

    bool verbose          = false;
};

template <class Parameters, class Evaluted_Object>
struct gauss_newton_result {
    Parameters argmax;
    Evaluted_Object initial_value;
    Evaluted_Object result_value;
    double max_value = 0.0;
    std::size_t n_iter = 0;
    std::string status;                    // "converged_grad" / "converged_value" /
                                           // "max_iter_reached" / "stalled"
};

namespace detail {

template <class Mat>
inline double frobenius_norm(Mat const& m) {
    double s = 0.0;
    for (std::size_t i = 0; i < m.nrows(); ++i)
        for (std::size_t j = 0; j < m.ncols(); ++j)
            s += m(i, j) * m(i, j);
    return std::sqrt(s);
}

}  // namespace detail


// ----------------------------------------------------------------------------
// gauss_newton_maximize
// ----------------------------------------------------------------------------
//
// Generic Gauss-Newton-with-LM-damping maximizer.
//
// Template parameters:
//   CurvTag      : the tag used to extract the PSD curvature matrix from the
//                  objective's state (e.g. FIM, Gaussian_Fisher_Information,
//                  Likelihood_Numerical_Fisher_Information, Score_Covariance_Matrix).
//   Objective    : callable with signature
//                    Maybe_error<State> operator()(Parameters const&) const;
//                  where State satisfies:
//                    - get<logL>(state)()                 -> double
//                    - get<Grad>(state)().value()         -> Matrix<double>
//                    - get<CurvTag>(state)().value()      -> SymPosDefMatrix<double>
//                    - get<Grad>(state)().parameters_ptr()    -> ptr (optional, for names)
//                    - get<CurvTag>(state)().parameters_ptr() -> ptr (optional, for names)
//   Parameters   : the parameter type (typically var::Parameters_transformed),
//                  supporting operator() to extract values as Matrix and
//                  create(Matrix) to produce a new instance with updated values.
//
// Behaviour:
//   - Iterates LM-damped Gauss-Newton steps until convergence or max_iter.
//   - Adapts lambda: shrink on accept, grow on reject.
//   - Returns Maybe_error<gauss_newton_result>: error if lambda exceeds
//     lambda_max while still rejecting / failing.
//
template <class CurvTag, class Objective, class Parameters>
auto gauss_newton_maximize(Objective const& objective, Parameters initial,
                           gauss_newton_options const& opts = {}) {

    auto theta  = std::move(initial);
    double lambda = 0.0;  // start in Newton mode; jumps to opts.lambda_kickoff on first failure

    // Initial objective evaluation.
    auto Maybe_state = objective(theta);

    // Deduce the evaluated-object type and the result type from the objective's
    // return. EvalT is whatever the objective puts inside its Maybe_error.
    using EvalT   = std::decay_t<decltype(Maybe_state.value())>;
    using ResultT = gauss_newton_result<Parameters, EvalT>;
    using ReturnT = Maybe_error<ResultT>;

    if (!Maybe_state)
        return ReturnT(error_message(
            std::string("gauss_newton: initial objective failed: ") + Maybe_state.error()()));
    EvalT initial_state = Maybe_state.value();
    double value = get<logL>(Maybe_state.value())();

    // Stable reference scale for the grad-norm fallback: the WARMSTART norm, NOT
    // the current theta(). Using the current theta is self-defeating — a non-
    // identifiable parameter that runs away inflates ‖θ‖, which inflates the
    // tolerance (grad_rtol·‖θ‖), so a diverged, NON-converged point (large
    // gradient) gets falsely accepted as "converged_grad". The warmstart is stable.
    const double theta_ref_norm = std::max(1.0, detail::frobenius_norm(theta()));

    std::size_t iter = 0;
    for (; iter < opts.max_iter; ++iter) {

        // Extract gradient and curvature at current theta (for the LM step).
        auto const& grad = get<Grad>(Maybe_state.value())().value();
        auto const& curv = get<CurvTag>(Maybe_state.value())().value();

        // -------- Convergence checks --------
        const double grad_norm  = detail::frobenius_norm(grad);
        const double theta_norm = std::max(1.0, detail::frobenius_norm(theta()));

        // Newton decrement λ² = gᵀ H⁻¹ g on the UNDAMPED curvature (scale-
        // invariant; see newton_dec2_tol). curv is PSD by construction; if it is
        // singular the decrement is left at +inf and we fall back to the
        // grad-norm criterion below. Cheap (p×p, p = #params); the cholesky is
        // recomputed for the damped solve below but that is negligible for small p.
        double newton_dec2 = std::numeric_limits<double>::infinity();
        {
            auto Maybe_cholH = cholesky(curv);
            if (Maybe_cholH) {
                auto Maybe_LHinv = inv(Maybe_cholH.value());
                if (Maybe_LHinv) {
                    auto Hinv   = XXT(tr(Maybe_LHinv.value()));
                    auto Hinv_g = Hinv * grad;
                    double acc  = 0.0;
                    for (std::size_t i = 0; i < grad.nrows(); ++i)
                        acc += grad(i, 0ul) * Hinv_g(i, 0ul);
                    newton_dec2 = acc;
                }
            }
        }

        if (opts.verbose) {
            std::cerr << "[gn] iter=" << iter
                      << " lambda=" << lambda
                      << " value=" << value
                      << " grad_norm=" << grad_norm
                      << " newton_dec2=" << newton_dec2
                      << " theta_norm=" << theta_norm
                      << "\n";
        }

        // PRIMARY: Newton-decrement convergence (objective sub-optimality ½·λ²).
        if (0.5 * newton_dec2 < opts.newton_dec2_tol) {
            return ReturnT(ResultT{
                std::move(theta), std::move(initial_state),
                std::move(Maybe_state.value()),
                value, iter, "converged_newton_dec"
            });
        }
        // FALLBACK: raw gradient norm, relative to the STABLE warmstart scale
        // (theta_ref_norm), NEVER the current ‖θ‖ — a diverging parameter would
        // inflate the current norm and falsely satisfy this test (see above).
        if (grad_norm < opts.grad_rtol * theta_ref_norm) {
            return ReturnT(ResultT{
                std::move(theta), std::move(initial_state),
                std::move(Maybe_state.value()),
                value, iter, "converged_grad"
            });
        }

        // -------- Hybrid Marquardt + Levenberg damping --------
        // damped(i,i) = curv(i,i) * (1 + lambda) + lambda * max_diag
        //                ^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^
        //                  Marquardt (scale-       Levenberg shift
        //                   invariant)              (lifts near-zero eigs)
        //
        // Pure Marquardt (the multiplicative term alone) preserves the
        // condition number of curv — it can't rescue near-singular FIMs
        // (e.g. figure_2 ~ scheme_CO with σ_n = 1e-4 gives κ ~ 500 from
        // the Current_Noise vs Num_ch_mean diagonal-ratio alone, and worse
        // off-diagonal near-colinearities). Without the Levenberg term, the
        // LM step in the near-null direction blows up regardless of lambda,
        // producing overflow in the AD parameter transform (Log10⁻¹ → inf).
        //
        // The Levenberg term lambda·max_diag floods small eigenvalues with a
        // floor proportional to lambda. Effect at moderate lambda (~0.1·max_diag):
        // all damped diagonals become comparable, conditioning improves to O(10).
        // At small lambda, the Marquardt term dominates (good near-Newton step);
        // at large lambda, the Levenberg term dominates (well-conditioned
        // gradient-descent-like step). PSD preserved (curv PSD + nonnegative
        // shift to diagonals).
        //
        // Off-diagonals untouched in either term.
        SymPosDefMatrix<double> damped = curv;
        double max_diag = 0.0;
        for (std::size_t i = 0; i < damped.nrows(); ++i) {
            max_diag = std::max(max_diag, damped(i, i));
        }
        const double lev_shift = lambda * max_diag;
        if (opts.verbose) {
            std::cerr << "[gn]   pre-damping: curv(0,0)=" << curv(0,0)
                      << " curv(3,3)=" << curv(3,3)
                      << " max_diag=" << max_diag
                      << " lev_shift=" << lev_shift
                      << "\n";
        }
        for (std::size_t i = 0; i < damped.nrows(); ++i) {
            damped.set(i, i, damped(i, i) * (1.0 + lambda) + lev_shift);
        }
        if (opts.verbose) {
            std::cerr << "[gn]   post-damping: damped(0,0)=" << damped(0,0)
                      << " damped(3,3)=" << damped(3,3)
                      << "\n";
        }

        // Solve damped * delta = grad via Cholesky factorisation. Mirrors the
        // pattern used by parallel_levenberg_tempering (the only place in the
        // codebase that actually inverts a damped SPD matrix in production).
        // The direct `inv(damped)` route via Lapack_SymmPosDef_inv (dpotrf+dpotri)
        // is declared but unexercised — it returned garbage in practice, leaving
        // step magnitude essentially insensitive to lambda.
        //
        // Steps:
        //   1. damped = Lᵀ·L   (cholesky → upper-triangular L)
        //   2. L_inv = L⁻¹      (inv of triangular, tested)
        //   3. damped⁻¹ = L⁻ᵀ·L⁻¹ = XXT(tr(L_inv))
        // Helpers for damping walk along the sequence
        //   0 ↔ kickoff ↔ kickoff·factor ↔ kickoff·factor² ↔ …
        // step_up = reject/failure direction; step_down = accept direction.
        // Intermediate values (0, kickoff) are skipped on both sides — they
        // are "almost Newton" with no regularisation gain.
        auto step_up = [&]() {
            return (lambda == 0.0) ? opts.lambda_kickoff : lambda * opts.lambda_factor;
        };
        auto step_down = [&]() {
            return (lambda <= opts.lambda_kickoff) ? 0.0 : (lambda / opts.lambda_factor);
        };

        auto Maybe_chol = cholesky(damped);
        if (!Maybe_chol) {
            // Cholesky failed -> escalate damping, retry.
            if (opts.verbose) {
                std::cerr << "[gn]   cholesky(damped) FAILED — escalate lambda " << lambda
                          << " -> " << step_up() << "\n";
            }
            lambda = step_up();
            if (lambda > opts.lambda_max) {
                // lambda_max reached: return the current best theta as "stalled"
                // (a numerical stationary point) so EVERY group yields an MLE
                // result instead of being dropped. Quality is read downstream
                // from the gradient (Grad ~ 0 => genuine optimum).
                return ReturnT(ResultT{std::move(theta), std::move(initial_state),
                                       std::move(Maybe_state.value()), value, iter, "stalled"});
            }
            continue;
        }
        auto Maybe_L_inv = inv(Maybe_chol.value());
        if (!Maybe_L_inv) {
            if (opts.verbose) {
                std::cerr << "[gn]   inv(L) FAILED — escalate lambda " << lambda
                          << " -> " << step_up() << "\n";
            }
            lambda = step_up();
            if (lambda > opts.lambda_max) {
                // lambda_max reached: return the current best theta as "stalled"
                // (see the Cholesky branch) so the group still yields an MLE.
                return ReturnT(ResultT{std::move(theta), std::move(initial_state),
                                       std::move(Maybe_state.value()), value, iter, "stalled"});
            }
            continue;
        }
        auto damped_inv = XXT(tr(Maybe_L_inv.value()));

        // delta_theta = damped^{-1} * grad
        Matrix<double> delta_theta = damped_inv * grad;
        if (opts.verbose) {
            const double step_norm = detail::frobenius_norm(delta_theta);
            std::cerr << "[gn]   step_norm=" << step_norm
                      << " damped_inv(0,0)=" << damped_inv(0,0)
                      << " damped_inv(3,3)=" << damped_inv(3,3)
                      << "\n";
            std::cerr << "[gn]   step =";
            for (std::size_t i = 0; i < delta_theta.nrows(); ++i)
                std::cerr << " " << delta_theta(i, 0);
            std::cerr << "\n";
        }

        // Propose new theta. Parameters_transformed::create preserves the
        // parameter-name infrastructure.
        auto theta_proposed = theta.create(theta() + delta_theta);

        // -------- Evaluate objective at proposed theta --------
        auto Maybe_state_new = objective(theta_proposed);
        if (!Maybe_state_new) {
            // Objective failure (e.g. AD pipeline NaN, out-of-domain step):
            // escalate damping, retry from same theta.
            if (opts.verbose) {
                std::cerr << "[gn]   objective FAILED: " << Maybe_state_new.error()()
                          << " — escalate lambda " << lambda
                          << " -> " << step_up() << "\n";
            }
            lambda = step_up();
            if (lambda > opts.lambda_max) {
                // lambda_max reached with every proposal failing: return the last
                // good theta as "stalled" (see the Cholesky branch) so the group
                // still yields an MLE.
                return ReturnT(ResultT{std::move(theta), std::move(initial_state),
                                       std::move(Maybe_state.value()), value, iter, "stalled"});
            }
            continue;
        }

        const double value_new = get<logL>(Maybe_state_new.value())();
        const double dvalue    = value_new - value;

        // -------- Armijo-trivial accept/reject --------
        if (dvalue > 0.0) {
            // Accept: step DOWN one notch in the damping sequence. At kickoff
            // (smallest damped value) the next step is 0 (Newton); higher up
            // we divide by factor. Intermediate values between 0 and kickoff
            // are skipped — they would be "almost Newton" without regularisation.
            if (opts.verbose) {
                std::cerr << "[gn]   ACCEPT value " << value << " -> " << value_new
                          << " (dvalue=" << dvalue
                          << ") — step_down lambda " << lambda
                          << " -> " << step_down() << "\n";
            }
            theta       = std::move(theta_proposed);
            Maybe_state = std::move(Maybe_state_new);
            value       = value_new;
            lambda      = step_down();

            // Convergence check on value change.
            if (std::abs(dvalue) < opts.dvalue_tol) {
                return ReturnT(ResultT{
                    std::move(theta), std::move(initial_state),
                    std::move(Maybe_state.value()),
                    value, iter + 1, "converged_value"
                });
            }
        } else {
            // Reject: escalate damping and retry from same theta.
            if (opts.verbose) {
                std::cerr << "[gn]   REJECT value " << value << " -> " << value_new
                          << " (dvalue=" << dvalue
                          << ") — escalate lambda " << lambda
                          << " -> " << step_up() << "\n";
            }
            lambda = step_up();
            if (lambda > opts.lambda_max) {
                // No step improves even at maximal damping: the proposed step is
                // ~0 yet value_new <= value (sub-FP-resolution progress on a very
                // sharp/flat surface) — the current theta IS the numerical
                // optimum. Return it as "stalled" (see the Cholesky branch) so
                // every group yields an MLE instead of failing the whole refit.
                return ReturnT(ResultT{std::move(theta), std::move(initial_state),
                                       std::move(Maybe_state.value()), value, iter, "stalled"});
            }
        }
    }

    // Max iterations reached without convergence.
    return ReturnT(ResultT{
        std::move(theta), std::move(initial_state),
        std::move(Maybe_state.value()),
        value, iter, "max_iter_reached"
    });
}

}  // namespace macrodr::optimization

#endif  // MACRODR_OPTIMIZATION_GAUSS_NEWTON_H
