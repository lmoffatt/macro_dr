#pragma once

// acquisition_filter.h — analog acquisition-filter machinery for the Bessel
// member (MacroIR-Bessel plan, theory/macroir/notes/macroir_bessel_plan.md,
// milestone M0 "pole tables and pole-shifted divided differences").
//
// Self-contained: standard library only, no macro_dr dependencies, so it is
// testable in isolation (tests/math/test_acquisition_filter.cpp) and usable
// from both the likelihood side (calc_Qdtf_eig) and the simulator twin.
//
// Scope of this header (phase 0 of the implementation plan):
//   * reverse Bessel polynomials by exact integer recurrence,
//     θ_n(s) = (2n−1)·θ_{n−1}(s) + s²·θ_{n−2}(s);
//   * pole finding (Durand–Kerner) + −3 dB normalization by bisection —
//     no transcribed pole constants anywhere, everything is computed and
//     canaried, which kills the delay-norm vs −3 dB ×2.114 trap at the root;
//   * partial-fraction residues, impulse response, H(0)=1 canary;
//   * complex-capable divided differences of exp: dd1 (uniformly stable
//     exp((x+y)/2)·sinhc((x−y)/2) form) and dd2 (Opitz: (0,2) entry of the
//     exponential of a 3×3 upper-triangular divided-difference matrix, which
//     is uniformly accurate through every coincidence pattern — this is the
//     kernel the "near-coincident shifted arguments" M0 gate stresses);
//   * dd1_pair: the same dd1 written in explicit real-pair arithmetic with
//     no std::complex — the template seed for derivative mode (plan
//     decision 2: real 2×2 pole-pair arithmetic under Der<double>);
//   * exp_rot: the real 2×2 decay–rotation block e^{−σΔ}·R(ωΔ) that
//     propagates one conjugate pole pair;
//   * closed forms for the boxcar-read variance and lag-1 covariance of
//     filtered white instrument noise — validated in the tests against the
//     manuscript's appendix bound table (61.5% deficit / 0.641 lag-1
//     correlation at f_c·Δ = 0.2).
//
// Deliberately NOT here (later phases): the Qdtf window object, the
// pole-shifted gtotal tables (they live next to calc_Qdtm_eig in qmodel.h
// and consume dd1/dd2 from here), the augmented recursion, DSL plumbing.

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <numbers>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace macrodr::acqf {

using cdouble = std::complex<double>;

//--------------------------------------------------------------------------
// Reverse Bessel polynomial coefficients, ascending: a[0] + a[1]s + ... + s^n.
// Exact integers for n ≤ 8 (max coefficient 2027025 ≪ 2^53).
//--------------------------------------------------------------------------
inline std::vector<double> reverse_bessel_coeffs(int n) {
    std::vector<double> th_prev{1.0};        // θ0 = 1
    std::vector<double> th{1.0, 1.0};        // θ1 = 1 + s
    if (n == 0)
        return th_prev;
    for (int k = 2; k <= n; ++k) {
        std::vector<double> next(k + 1, 0.0);
        const double c = 2.0 * k - 1.0;
        for (std::size_t i = 0; i < th.size(); ++i) next[i] += c * th[i];
        for (std::size_t i = 0; i < th_prev.size(); ++i) next[i + 2] += th_prev[i];
        th_prev = std::move(th);
        th = std::move(next);
    }
    return th;
}

inline cdouble poly_eval(const std::vector<double>& a, cdouble z) {
    cdouble r = 0.0;
    for (std::size_t i = a.size(); i-- > 0;) r = r * z + a[i];
    return r;
}

//--------------------------------------------------------------------------
// Durand–Kerner root finder for a monic polynomial with simple roots.
// Bessel polynomials have well-separated simple roots; convergence is fast.
//--------------------------------------------------------------------------
inline std::vector<cdouble> poly_roots(const std::vector<double>& a_in) {
    std::vector<double> a = a_in;
    const std::size_t n = a.size() - 1;
    for (auto& c : a) c /= a_in[n];  // make monic

    double radius = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        radius = std::max(radius, std::pow(std::abs(a[i]), 1.0 / double(n - i)));
    radius = 1.0 + radius;

    std::vector<cdouble> z(n);
    for (std::size_t k = 0; k < n; ++k) {
        double ang = 2.0 * std::numbers::pi * double(k) / double(n) + 0.4;
        z[k] = radius * cdouble(std::cos(ang), std::sin(ang));
    }
    for (int iter = 0; iter < 500; ++iter) {
        double step = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            cdouble denom = 1.0;
            for (std::size_t j = 0; j < n; ++j)
                if (j != k)
                    denom *= (z[k] - z[j]);
            cdouble dz = poly_eval(a, z[k]) / denom;
            z[k] -= dz;
            step = std::max(step, std::abs(dz) / (1.0 + std::abs(z[k])));
        }
        if (step < 1e-15)
            break;
    }
    return z;
}

//--------------------------------------------------------------------------
// Divided differences of exp, complex-capable.
//
// dd1[x,y] = (e^x − e^y)/(x − y), continuous value e^x at x = y.
// Stable form: e^{(x+y)/2} · sinhc((x−y)/2),  sinhc(d) = sinh(d)/d.
//--------------------------------------------------------------------------
inline cdouble sinhc(cdouble d) {
    if (std::norm(d) < 0.25) {  // |d| < 0.5 → series, error < 1e-18 at 8 terms
        cdouble d2 = d * d, term = 1.0, sum = 1.0;
        double fact = 1.0;
        for (int k = 1; k <= 8; ++k) {
            fact = fact * (2.0 * k) * (2.0 * k + 1.0);
            term *= d2;
            sum += term / fact;
        }
        return sum;
    }
    return std::sinh(d) / d;
}

// Wide real separation: the sinhc form multiplies e^{(x+y)/2} (underflowing
// to 0) by sinh((x−y)/2) (overflowing to inf once |Re(x−y)/2| ≳ 710) and
// yields 0·inf = NaN, although the value itself, ≈ e^{max(x,y)}/(x−y), is
// representable. Beyond |Re d| = 50 the two exponentials differ by e^{100}, so
// the direct quotient has no cancellation and is used instead. Reached by
// the filter poles over long windows: |p|Δ ≈ 10³ for a 10 kHz Bessel over a
// 10 ms window (tests/macroir/test_qdtf_member.cpp, 2026-09-06).
inline constexpr double dd1_direct_threshold = 50.0;

inline cdouble dd1(cdouble x, cdouble y) {
    const cdouble d = 0.5 * (x - y);
    if (std::abs(d.real()) > dd1_direct_threshold)
        return (std::exp(x) - std::exp(y)) / (x - y);
    return std::exp(0.5 * (x + y)) * sinhc(d);
}

// dd2[x,y,z]: second divided difference of exp. Opitz: it is the (0,2)
// entry of expm(T) with T = [[x,1,0],[0,y,1],[0,0,z]] upper triangular.
// Computed by scaling-and-squaring with a Taylor series on the scaled
// matrix; uniformly accurate through every coincidence pattern, which is
// what the near-coincident-arguments gate needs. Cost is irrelevant here
// (3×3, called per mode pair per distinct window length, memoized upstream).
namespace detail {
struct UT3 {  // upper-triangular 3×3 complex
    cdouble a00, a01, a02, a11, a12, a22;
};
inline UT3 mul(const UT3& A, const UT3& B) {
    UT3 C;
    C.a00 = A.a00 * B.a00;
    C.a11 = A.a11 * B.a11;
    C.a22 = A.a22 * B.a22;
    C.a01 = A.a00 * B.a01 + A.a01 * B.a11;
    C.a12 = A.a11 * B.a12 + A.a12 * B.a22;
    C.a02 = A.a00 * B.a02 + A.a01 * B.a12 + A.a02 * B.a22;
    return C;
}
}  // namespace detail

inline cdouble dd2(cdouble x, cdouble y, cdouble z) {
    using detail::UT3;
    double nrm = std::max({std::abs(x), std::abs(y), std::abs(z), 1.0});
    int s = std::max(0, int(std::ceil(std::log2(nrm))) + 1);  // scaled norm ≤ 0.5
    double sc = std::ldexp(1.0, -s);
    UT3 T{x * sc, sc, 0.0, y * sc, sc, z * sc};
    // expm(T) by Taylor: I + T + T²/2! + ... ; scaled ‖T‖ ≤ ~0.5 → 20 terms ≫ enough
    UT3 E{1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
    UT3 P = T;
    double fact = 1.0;
    for (int k = 1; k <= 20; ++k) {
        fact *= k;
        E.a00 += P.a00 / fact;
        E.a01 += P.a01 / fact;
        E.a02 += P.a02 / fact;
        E.a11 += P.a11 / fact;
        E.a12 += P.a12 / fact;
        E.a22 += P.a22 / fact;
        P = detail::mul(P, T);
    }
    for (int k = 0; k < s; ++k) E = detail::mul(E, E);
    return E.a02;
}

//--------------------------------------------------------------------------
// dd1_pair: dd1 with real x and complex y = (y_re, y_im), in explicit
// real-pair arithmetic (no std::complex). This is the algebra that the
// derivative-mode overloads templated on C_double will reuse verbatim:
// every operation below is +,−,×,÷,exp,sin,cos,sinh,cosh on reals.
//--------------------------------------------------------------------------
inline void dd1_pair(double x, double y_re, double y_im, double& out_re, double& out_im) {
    const double m_re = 0.5 * (x + y_re), m_im = 0.5 * y_im;
    const double d_re = 0.5 * (x - y_re), d_im = -0.5 * y_im;

    if (std::abs(d_re) > dd1_direct_threshold) {
        // wide separation: (e^x − e^y)/(x − y) in real-pair arithmetic (see dd1)
        const double ey = std::exp(y_re);
        const double n_re = std::exp(x) - ey * std::cos(y_im), n_im = -ey * std::sin(y_im);
        const double q_re = x - y_re, q_im = -y_im;
        const double q2 = q_re * q_re + q_im * q_im;
        out_re = (n_re * q_re + n_im * q_im) / q2;
        out_im = (n_im * q_re - n_re * q_im) / q2;
        return;
    }

    // e^m
    const double em = std::exp(m_re);
    const double em_re = em * std::cos(m_im), em_im = em * std::sin(m_im);

    // sinhc(d)
    double s_re, s_im;
    const double d2 = d_re * d_re + d_im * d_im;
    if (d2 < 0.25) {
        // series in d²: d² = (d_re²−d_im², 2 d_re d_im)
        const double d2_re = d_re * d_re - d_im * d_im;
        const double d2_im = 2.0 * d_re * d_im;
        double t_re = 1.0, t_im = 0.0;  // running d^{2k}
        s_re = 1.0;
        s_im = 0.0;
        double fact = 1.0;
        for (int k = 1; k <= 8; ++k) {
            fact = fact * (2.0 * k) * (2.0 * k + 1.0);
            const double n_re = t_re * d2_re - t_im * d2_im;
            const double n_im = t_re * d2_im + t_im * d2_re;
            t_re = n_re;
            t_im = n_im;
            s_re += t_re / fact;
            s_im += t_im / fact;
        }
    } else {
        // sinh(d) = (sinh(d_re)cos(d_im), cosh(d_re)sin(d_im)); then /d
        const double sh_re = std::sinh(d_re) * std::cos(d_im);
        const double sh_im = std::cosh(d_re) * std::sin(d_im);
        s_re = (sh_re * d_re + sh_im * d_im) / d2;
        s_im = (sh_im * d_re - sh_re * d_im) / d2;
    }
    out_re = em_re * s_re - em_im * s_im;
    out_im = em_re * s_im + em_im * s_re;
}

//--------------------------------------------------------------------------
// exp_rot: the real 2×2 propagator of one conjugate pole pair over dt.
// State (u,v) = (2·Re x, 2·Im x) of the complex modal state x with pole
// p = −σ + iω:  d/dt (u,v) = [[−σ, −ω],[ω, −σ]](u,v) + input terms.
// exp of that block = e^{−σ dt}·[[cos ωdt, −sin ωdt],[sin ωdt, cos ωdt]].
// Returned row-major {m00, m01, m10, m11}.
//--------------------------------------------------------------------------
inline std::array<double, 4> exp_rot(double sigma_dt, double omega_dt) {
    const double e = std::exp(-sigma_dt);
    const double c = std::cos(omega_dt), s = std::sin(omega_dt);
    return {e * c, -e * s, e * s, e * c};
}

//--------------------------------------------------------------------------
// Filter realization: poles (as conjugate pairs) + residues, scaled so the
// −3 dB point of |H| sits at f_c (Hz). H(s) = θ(0)/θ(s/ω_scale) has unit DC
// gain by construction; the H(0)=1 canary re-checks it from the residues.
//--------------------------------------------------------------------------
struct ModePair {
    double sigma;  // > 0; pole p = −sigma + i·omega (and its conjugate)
    double omega;  // ≥ 0
    double r_re;   // residue at p = −sigma + i·omega
    double r_im;
};

struct FilterRealization {
    int order = 0;
    double f_c_hz = 0.0;
    double omega3_norm = 0.0;             // −3 dB frequency of the normalized θ_n
    std::vector<ModePair> pairs;          // order/2 pairs (even order)
    std::vector<cdouble> poles;           // all `order` poles, physical scale (rad/s)
    std::vector<cdouble> residues;        // matching residues, h(t) = Σ r_k e^{p_k t}
};

inline double impulse_area(const FilterRealization& f) {  // = H(0)
    cdouble s = 0.0;
    for (std::size_t k = 0; k < f.poles.size(); ++k) s += f.residues[k] / (-f.poles[k]);
    return s.real();
}

inline double h_impulse(const FilterRealization& f, double t) {
    double h = 0.0;
    for (const auto& m : f.pairs)
        h += 2.0 * std::exp(-m.sigma * t) *
             (m.r_re * std::cos(m.omega * t) - m.r_im * std::sin(m.omega * t));
    return h;
}

inline std::optional<std::string> validate(const FilterRealization& f, double tol = 1e-8) {
    for (const auto& p : f.poles)
        if (!(p.real() < 0.0))
            return "unstable pole (Re ≥ 0)";
    if (std::abs(impulse_area(f) - 1.0) > tol)
        return "H(0) != 1: " + std::to_string(impulse_area(f));
    return std::nullopt;
}

// Bessel low-pass of even order (4 or 8 in production), −3 dB cutoff f_c_hz.
inline FilterRealization make_bessel_filter(int order, double f_c_hz) {
    FilterRealization f;
    f.order = order;
    f.f_c_hz = f_c_hz;

    const auto coeffs = reverse_bessel_coeffs(order);
    const double th0 = coeffs[0];

    // −3 dB of the normalized filter: |θ(iω)|² = 2·θ(0)², monotone in ω.
    double lo = 0.0, hi = 4.0 * order;
    for (int it = 0; it < 200; ++it) {
        double mid = 0.5 * (lo + hi);
        double m2 = std::norm(poly_eval(coeffs, cdouble(0.0, mid)));
        (m2 < 2.0 * th0 * th0 ? lo : hi) = mid;
    }
    f.omega3_norm = 0.5 * (lo + hi);

    const double scale = 2.0 * std::numbers::pi * f_c_hz / f.omega3_norm;

    auto roots = poly_roots(coeffs);
    for (auto& z : roots) z *= scale;
    f.poles = roots;

    // Residues of H(s) = A/Π(s − p_j), A = Π(−p_j) (unit DC gain):
    // r_k = A / Π_{j≠k}(p_k − p_j).
    cdouble A = 1.0;
    for (const auto& p : f.poles) A *= (-p);
    f.residues.resize(f.poles.size());
    for (std::size_t k = 0; k < f.poles.size(); ++k) {
        cdouble denom = 1.0;
        for (std::size_t j = 0; j < f.poles.size(); ++j)
            if (j != k)
                denom *= (f.poles[k] - f.poles[j]);
        f.residues[k] = A / denom;
    }

    for (std::size_t k = 0; k < f.poles.size(); ++k)
        if (f.poles[k].imag() > 0.0)
            f.pairs.push_back({-f.poles[k].real(), f.poles[k].imag(), f.residues[k].real(),
                               f.residues[k].imag()});
    return f;
}

//--------------------------------------------------------------------------
// White instrument noise of spectral density S0 (E[ξ(t)ξ(t')] = S0·δ(t−t'))
// through the filter, read as the boxcar average over a window of length dt:
//   z_j = (1/dt)·∫_{j·dt}^{(j+1)·dt} (h*ξ)(t) dt.
// Closed forms from the stationary autocovariance
//   R(τ) = S0 · Σ_{k,l} r_k r_l · e^{p_l τ} / (−(p_k+p_l)),  τ ≥ 0.
// These are the numbers behind the manuscript appendix's bound table
// (variance deficit and lag-one correlation vs f_c·Δ); the box model
// assigns Var = S0/dt and zero lag covariance.
//--------------------------------------------------------------------------
inline double noise_autocov(const FilterRealization& f, double S0, double tau) {
    cdouble s = 0.0;
    for (std::size_t k = 0; k < f.poles.size(); ++k)
        for (std::size_t l = 0; l < f.poles.size(); ++l)
            s += f.residues[k] * f.residues[l] * std::exp(f.poles[l] * tau) /
                 (-(f.poles[k] + f.poles[l]));
    return S0 * s.real();
}

inline double noise_read_variance(const FilterRealization& f, double S0, double dt) {
    cdouble s = 0.0;
    for (std::size_t k = 0; k < f.poles.size(); ++k)
        for (std::size_t l = 0; l < f.poles.size(); ++l) {
            const cdouble p = f.poles[l];
            s += f.residues[k] * f.residues[l] / (-(f.poles[k] + p)) *
                 (std::exp(p * dt) - 1.0 - p * dt) / (p * p);
        }
    return S0 * 2.0 / (dt * dt) * s.real();
}

// Covariance of reads j and j+m, m ≥ 1.
inline double noise_read_cov(const FilterRealization& f, double S0, double dt, int m) {
    cdouble s = 0.0;
    for (std::size_t k = 0; k < f.poles.size(); ++k)
        for (std::size_t l = 0; l < f.poles.size(); ++l) {
            const cdouble p = f.poles[l];
            const cdouble e1 = std::exp(p * dt) - 1.0;
            s += f.residues[k] * f.residues[l] / (-(f.poles[k] + p)) *
                 std::exp(p * dt * double(m - 1)) * e1 * e1 / (p * p);
        }
    return S0 / (dt * dt) * s.real();
}

}  // namespace macrodr::acqf
