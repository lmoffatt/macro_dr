#pragma once

// qdtf_engine.h — the mathematical core of the Bessel member: the spectral
// window producer (what calc_Qdtf_eig computes) and the augmented recursion
// (predict / update), in plain double precision.
//
// Architecture note. This is deliberately a self-contained header (stdlib +
// acquisition_filter.h only) so it compiles in seconds and is gated in
// isolation (tests/math/test_qdtf_engine.cpp) against quadrature of the
// defining integrals and against the boundary-state (Σ^bnd) reference. The
// qmodel.h calc_Qdtf_eig is then a thin adapter: extract (V, λ, U, γ, Δ)
// from Qx_eig/Patch_Model, call calc_qdtf_window, wrap the result in the
// Qdtf Vector_Space. The recursion sibling adapts predict/update the same
// way. Equations: theory/macroir/notes/bessel_reference/
// augmented_recursion_blocks.md, validated by Monte Carlo in
// bessel_oracle.py (gates O1-O5, 2026-08-31).
//
// Representation. The observation system is carried in PHYSICAL real states:
// one (u,v) block per conjugate pole pair ((u,v) = (2Re x, 2Im x) of the
// complex modal state), plus one accumulator state u_acc for the grouped and
// box reads (reset at each window start — the recorder's own reset). Window
// integrals are computed SPECTRALLY over the full complex mode list
// {filter poles (both conjugates)} ∪ {0 (the accumulator's input mode)} via
// the pole-shifted divided differences dd1/dd2, then assembled to physical
// states through the fixed complex map T:
//     pair block:  u = x + x̄,   v = −i(x − x̄)
//     accumulator: u_acc = Σ_j x_j / p_j + w₀      (w₀ = ∫ i, the zero mode;
//                  Σ_j r_j/p_j = −1 because H(0)=1 makes this exact)
// Free evolution of the accumulator has the closed form
//     u_acc(Δ) = u_acc(0) + Σ_j x_j(0)·(e^{p_jΔ}−1)/p_j,
// which is the Phi row of the accumulator; the zero mode has no free growth.
//
// Channel-side inputs: the eigendecomposition of the generator,
// P(t) = V·diag(e^{λt})·U with U = V⁻¹ and real eigenvalues λ, conductance
// vector γ, window length Δ. The three spectral formulas (derivation:
// simplex/Hermite–Genocchi, see the notes):
//   F1c[mode][i][j] = ρΔ Σ_{a,b} V_ia Mg_ab U_bj · dd1(λ_aΔ, (λ_b+q)Δ)
//   S2c[m1][m2][i]  = ρ₁ρ₂Δ² Σ_{a,b} V_ia Mg_ab g̃_b ·
//                     [dd2(λ_aΔ, (λ_b+q₁)Δ, (q₁+q₂)Δ) + (q₁ ↔ q₂)]
//   Nc[m1][m2]      = S0·ρ₁ρ₂·Δ·E1((q₁+q₂)Δ)
// with Mg = U·diag(γ)·V, g̃ = U·γ. Setting every q = 0 reproduces the box
// member's Ee/E3 tables identically (the regression anchor).

#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

#include "acquisition_filter.h"

namespace macrodr::qdtf {

using acqf::cdouble;
using acqf::dd1;
using acqf::dd2;

//--------------------------------------------------------------------------
// Minimal dense real matrix (row-major); K and m are single digits here.
//--------------------------------------------------------------------------
struct Mat {
    std::size_t nr = 0, nc = 0;
    std::vector<double> d;
    Mat() = default;
    Mat(std::size_t r, std::size_t c) : nr(r), nc(c), d(r * c, 0.0) {}
    double& operator()(std::size_t r, std::size_t c) {
        return d[r * nc + c];
    }
    double operator()(std::size_t r, std::size_t c) const {
        return d[r * nc + c];
    }
};

inline Mat mul(const Mat& A, const Mat& B) {
    Mat C(A.nr, B.nc);
    for (std::size_t i = 0; i < A.nr; ++i)
        for (std::size_t k = 0; k < A.nc; ++k) {
            const double a = A(i, k);
            if (a == 0.0)
                continue;
            for (std::size_t j = 0; j < B.nc; ++j) C(i, j) += a * B(k, j);
        }
    return C;
}

inline Mat trans(const Mat& A) {
    Mat T(A.nc, A.nr);
    for (std::size_t i = 0; i < A.nr; ++i)
        for (std::size_t j = 0; j < A.nc; ++j) T(j, i) = A(i, j);
    return T;
}

inline Mat diag_of(const std::vector<double>& v) {
    Mat D(v.size(), v.size());
    for (std::size_t i = 0; i < v.size(); ++i) D(i, i) = v[i];
    return D;
}

//--------------------------------------------------------------------------
// Mode set: full complex mode list + the physical-state map T + read info.
//--------------------------------------------------------------------------
struct QdtfModeSet {
    std::vector<cdouble> q;    // poles, conjugates explicit; zero mode last if present
    std::vector<cdouble> rho;  // input loadings (residues; 1.0 for the zero mode)
    std::size_t n_pairs = 0;   // physical (u,v) blocks; pair k ↔ modes 2k, 2k+1
    bool has_acc = false;      // accumulator state present (grouped / box reads)

    std::size_t m() const {
        return 2 * n_pairs + (has_acc ? 1 : 0);
    }
    std::size_t acc_row() const {
        return m() - 1;
    }  // valid when has_acc
    // T[a][j]: physical state a as Σ_j T[a][j]·(modal j); reality holds by
    // construction, the assembled tables take the real part.
    std::vector<std::vector<cdouble>> T;
};

namespace detail {
inline void build_T(QdtfModeSet& ms) {
    const std::size_t nfull = ms.q.size();
    ms.T.assign(ms.m(), std::vector<cdouble>(nfull, 0.0));
    for (std::size_t k = 0; k < ms.n_pairs; ++k) {
        ms.T[2 * k][2 * k] = 1.0;  // u = x + x̄
        ms.T[2 * k][2 * k + 1] = 1.0;
        ms.T[2 * k + 1][2 * k] = cdouble(0.0, -1.0);  // v = −i(x − x̄)
        ms.T[2 * k + 1][2 * k + 1] = cdouble(0.0, 1.0);
    }
    if (ms.has_acc) {
        auto& row = ms.T[ms.acc_row()];
        for (std::size_t j = 0; j < 2 * ms.n_pairs; ++j) row[j] = 1.0 / ms.q[j];
        row[nfull - 1] = 1.0;  // the zero mode w₀
    }
}
}  // namespace detail

// Native read (one raw sample): states = filter pairs, read = filter output.
inline QdtfModeSet make_native_modes(const acqf::FilterRealization& f) {
    QdtfModeSet ms;
    for (const auto& pr : f.pairs) {
        ms.q.push_back(cdouble(-pr.sigma, pr.omega));
        ms.q.push_back(cdouble(-pr.sigma, -pr.omega));
        ms.rho.push_back(cdouble(pr.r_re, pr.r_im));
        ms.rho.push_back(cdouble(pr.r_re, -pr.r_im));
    }
    ms.n_pairs = f.pairs.size();
    ms.has_acc = false;
    detail::build_T(ms);
    return ms;
}

// Grouped read: filter pairs + zero mode + accumulator of the filter OUTPUT.
inline QdtfModeSet make_grouped_modes(const acqf::FilterRealization& f) {
    QdtfModeSet ms = make_native_modes(f);
    ms.q.push_back(0.0);
    ms.rho.push_back(1.0);
    ms.has_acc = true;
    detail::build_T(ms);
    return ms;
}

// Box read: the current member — accumulator of the raw input, no filter.
inline QdtfModeSet make_box_modes() {
    QdtfModeSet ms;
    ms.q.push_back(0.0);
    ms.rho.push_back(1.0);
    ms.n_pairs = 0;
    ms.has_acc = true;
    detail::build_T(ms);
    return ms;
}

// Read vector on physical states. Native: the filter output Σ u_k.
// Grouped/box: the accumulator average u_acc/Δ.
inline std::vector<double> read_vector(const QdtfModeSet& ms, double dt) {
    std::vector<double> r(ms.m(), 0.0);
    if (ms.has_acc) {
        r[ms.acc_row()] = 1.0 / dt;
    } else {
        for (std::size_t k = 0; k < ms.n_pairs; ++k) r[2 * k] = 1.0;
    }
    return r;
}

// Native (single-raw-sample) read on a layout that carries the accumulator:
// sample the filter output, ignore the accumulator. Lets one state layout
// serve records that mix native-rate and grouped windows.
inline std::vector<double> native_read_vector(const QdtfModeSet& ms) {
    std::vector<double> r(ms.m(), 0.0);
    for (std::size_t k = 0; k < ms.n_pairs; ++k) r[2 * k] = 1.0;
    return r;
}

//--------------------------------------------------------------------------
// The window object: what calc_Qdtf_eig produces per (Δ, agonist, mode set).
//--------------------------------------------------------------------------
struct QdtfWindow {
    Mat P;                // K×K transition matrix
    std::vector<Mat> F1;  // K entries; F1[i](j, a) = E[1{end=j}·w_a | start=i]
    Mat f1;               // K×m, start-conditioned firsts (row sums of F1)
    std::vector<Mat> S2;  // K entries m×m; E[w wᵀ | start=i] (end-marginalized)
    Mat Phi;              // m×m free-evolution propagator over Δ
    Mat Sig;              // m×m instrument-noise Gramian
};

inline QdtfWindow calc_qdtf_window(const Mat& V, const std::vector<double>& lam, const Mat& U,
                                   const std::vector<double>& gamma, double dt,
                                   const QdtfModeSet& ms, double S0) {
    const std::size_t K = lam.size();
    const std::size_t nfull = ms.q.size();
    const std::size_t m = ms.m();
    QdtfWindow w;

    // P = V diag(e^{λΔ}) U
    std::vector<double> elam(K);
    for (std::size_t a = 0; a < K; ++a) elam[a] = std::exp(lam[a] * dt);
    w.P = mul(mul(V, diag_of(elam)), U);

    // Mg = U diag(γ) V,  g̃ = U γ
    Mat Mg(K, K);
    std::vector<double> gt(K, 0.0);
    for (std::size_t a = 0; a < K; ++a)
        for (std::size_t b = 0; b < K; ++b) {
            double s = 0.0;
            for (std::size_t k = 0; k < K; ++k) s += U(a, k) * gamma[k] * V(k, b);
            Mg(a, b) = s;
        }
    for (std::size_t a = 0; a < K; ++a) {
        double s = 0.0;
        for (std::size_t k = 0; k < K; ++k) s += U(a, k) * gamma[k];
        gt[a] = s;
    }

    // modal firsts F1c[mode](i,j) and start-marginalized seconds S2c[m1][m2][i]
    std::vector<std::vector<cdouble>> F1c(nfull, std::vector<cdouble>(K * K, 0.0));
    for (std::size_t jm = 0; jm < nfull; ++jm) {
        const cdouble q = ms.q[jm], rho = ms.rho[jm];
        for (std::size_t a = 0; a < K; ++a)
            for (std::size_t b = 0; b < K; ++b) {
                const cdouble ker =
                    dd1(cdouble(lam[a] * dt, 0.0), (cdouble(lam[b], 0.0) + q) * dt);
                const cdouble w_ab = rho * dt * Mg(a, b) * ker;
                for (std::size_t i = 0; i < K; ++i)
                    for (std::size_t j = 0; j < K; ++j)
                        F1c[jm][i * K + j] += V(i, a) * w_ab * U(b, j);
            }
    }
    std::vector<std::vector<std::vector<cdouble>>> S2c(
        nfull, std::vector<std::vector<cdouble>>(nfull, std::vector<cdouble>(K, 0.0)));
    for (std::size_t j1 = 0; j1 < nfull; ++j1)
        for (std::size_t j2 = j1; j2 < nfull; ++j2) {
            const cdouble q1 = ms.q[j1], q2 = ms.q[j2];
            for (std::size_t a = 0; a < K; ++a)
                for (std::size_t b = 0; b < K; ++b) {
                    const cdouble k12 = dd2(cdouble(lam[a] * dt, 0.0),
                                            (cdouble(lam[b], 0.0) + q1) * dt, (q1 + q2) * dt);
                    const cdouble k21 = dd2(cdouble(lam[a] * dt, 0.0),
                                            (cdouble(lam[b], 0.0) + q2) * dt, (q1 + q2) * dt);
                    const cdouble w_ab =
                        ms.rho[j1] * ms.rho[j2] * dt * dt * Mg(a, b) * gt[b] * (k12 + k21);
                    for (std::size_t i = 0; i < K; ++i) S2c[j1][j2][i] += V(i, a) * w_ab;
                }
            if (j2 != j1)
                S2c[j2][j1] = S2c[j1][j2];
        }

    // assemble physical tables through T
    w.F1.assign(K, Mat(K, m));
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j)
            for (std::size_t a = 0; a < m; ++a) {
                cdouble s = 0.0;
                for (std::size_t jm = 0; jm < nfull; ++jm)
                    s += ms.T[a][jm] * F1c[jm][i * K + j];
                w.F1[i](j, a) = s.real();
            }
    w.f1 = Mat(K, m);
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t a = 0; a < m; ++a) {
            double s = 0.0;
            for (std::size_t j = 0; j < K; ++j) s += w.F1[i](j, a);
            w.f1(i, a) = s;
        }
    w.S2.assign(K, Mat(m, m));
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t a = 0; a < m; ++a)
            for (std::size_t b = 0; b < m; ++b) {
                cdouble s = 0.0;
                for (std::size_t j1 = 0; j1 < nfull; ++j1)
                    for (std::size_t j2 = 0; j2 < nfull; ++j2)
                        s += ms.T[a][j1] * ms.T[b][j2] * S2c[j1][j2][i];
                w.S2[i](a, b) = s.real();
            }

    // noise Gramian
    w.Sig = Mat(m, m);
    for (std::size_t a = 0; a < m; ++a)
        for (std::size_t b = 0; b < m; ++b) {
            cdouble s = 0.0;
            for (std::size_t j1 = 0; j1 < nfull; ++j1)
                for (std::size_t j2 = 0; j2 < nfull; ++j2) {
                    const cdouble z = (ms.q[j1] + ms.q[j2]) * dt;
                    s += ms.T[a][j1] * ms.T[b][j2] * ms.rho[j1] * ms.rho[j2] * dt *
                         dd1(z, cdouble(0.0, 0.0));
                }
            w.Sig(a, b) = S0 * s.real();
        }

    // free-evolution propagator: exp_rot pair blocks; accumulator closed form
    w.Phi = Mat(m, m);
    for (std::size_t k = 0; k < ms.n_pairs; ++k) {
        const double sg = -ms.q[2 * k].real(), om = ms.q[2 * k].imag();
        auto R = acqf::exp_rot(sg * dt, om * dt);
        w.Phi(2 * k, 2 * k) = R[0];
        w.Phi(2 * k, 2 * k + 1) = R[1];
        w.Phi(2 * k + 1, 2 * k) = R[2];
        w.Phi(2 * k + 1, 2 * k + 1) = R[3];
    }
    if (ms.has_acc) {
        const std::size_t ua = ms.acc_row();
        w.Phi(ua, ua) = 1.0;
        for (std::size_t k = 0; k < ms.n_pairs; ++k) {
            const cdouble p = ms.q[2 * k];
            const cdouble c = (std::exp(p * dt) - 1.0) / p;
            // Σ_j x_j(0)(e^{p_jΔ}−1)/p_j over the pair = 2Re(c·x) = Re(c)u − Im(c)v
            w.Phi(ua, 2 * k) = c.real();
            w.Phi(ua, 2 * k + 1) = -c.imag();
        }
    }
    return w;
}

//--------------------------------------------------------------------------
// The augmented recursion, count scale (equations and validation:
// augmented_recursion_blocks.md; the box subcase reproduces MacroIR).
//--------------------------------------------------------------------------
struct QdtfState {
    Mat muN;  // K×1
    Mat SNN;  // K×K
    Mat mx;   // m×1
    Mat CNx;  // K×m  = Cov(N, x)
    Mat Sxx;  // m×m
};

inline QdtfState make_initial_state(const std::vector<double>& p_start, double Nch,
                                    std::size_t m) {
    const std::size_t K = p_start.size();
    QdtfState s;
    s.muN = Mat(K, 1);
    s.SNN = Mat(K, K);
    for (std::size_t i = 0; i < K; ++i) {
        s.muN(i, 0) = Nch * p_start[i];
        for (std::size_t j = 0; j < K; ++j)
            s.SNN(i, j) = Nch * ((i == j ? p_start[i] : 0.0) - p_start[i] * p_start[j]);
    }
    s.mx = Mat(m, 1);
    s.CNx = Mat(K, m);
    s.Sxx = Mat(m, m);
    return s;
}

inline QdtfState predict(const QdtfState& st, const QdtfWindow& w, int reset_row = -1) {
    const std::size_t K = st.muN.nr, m = st.mx.nr;
    QdtfState p = st;
    if (reset_row >= 0) {  // the recorder's own reset of the accumulator
        const std::size_t r = std::size_t(reset_row);
        p.mx(r, 0) = 0.0;
        for (std::size_t i = 0; i < K; ++i) p.CNx(i, r) = 0.0;
        for (std::size_t a = 0; a < m; ++a) {
            p.Sxx(r, a) = 0.0;
            p.Sxx(a, r) = 0.0;
        }
    }
    QdtfState out;
    // μ_N⁺ = Pᵀ μ_N ;   S_NN⁺ = Pᵀ(S_NN − diag μ_N)P + diag(μ_N⁺)   [σ_pre]
    Mat Pt = trans(w.P);
    out.muN = mul(Pt, p.muN);
    Mat SmD = p.SNN;
    for (std::size_t i = 0; i < K; ++i) SmD(i, i) -= p.muN(i, 0);
    out.SNN = mul(mul(Pt, SmD), w.P);
    for (std::size_t j = 0; j < K; ++j) out.SNN(j, j) += out.muN(j, 0);

    // m_x⁺ = Φ m_x + Σ_i μ_i f1_i
    out.mx = mul(w.Phi, p.mx);
    for (std::size_t a = 0; a < m; ++a)
        for (std::size_t i = 0; i < K; ++i) out.mx(a, 0) += p.muN(i, 0) * w.f1(i, a);

    // C_Nx⁺[j] = (Pᵀ C_Nx Φᵀ)[j] + Σ_i μ_i (F1[i](j,·) − P_ij f1_i)
    //          + Σ_{i,k} P_ij S_NN[i,k] f1_k                    [the gS shape]
    out.CNx = mul(Pt, mul(p.CNx, trans(w.Phi)));
    for (std::size_t j = 0; j < K; ++j)
        for (std::size_t a = 0; a < m; ++a) {
            double s = 0.0;
            for (std::size_t i = 0; i < K; ++i) {
                s += p.muN(i, 0) * (w.F1[i](j, a) - w.P(i, j) * w.f1(i, a));
                for (std::size_t k = 0; k < K; ++k)
                    s += w.P(i, j) * p.SNN(i, k) * w.f1(k, a);
            }
            out.CNx(j, a) += s;
        }

    // S_xx⁺ = Φ S_xx Φᵀ + Φ M1 + (Φ M1)ᵀ + V_w + Σ_ξ
    Mat M1(m, m);  // Cov(x0, Σw) = Σ_i outer(C_Nx[i], f1_i)
    for (std::size_t a = 0; a < m; ++a)
        for (std::size_t b = 0; b < m; ++b) {
            double s = 0.0;
            for (std::size_t i = 0; i < K; ++i) s += p.CNx(i, a) * w.f1(i, b);
            M1(a, b) = s;
        }
    Mat PhiM1 = mul(w.Phi, M1);
    out.Sxx = mul(mul(w.Phi, p.Sxx), trans(w.Phi));
    for (std::size_t a = 0; a < m; ++a)
        for (std::size_t b = 0; b < m; ++b) {
            double s = PhiM1(a, b) + PhiM1(b, a) + w.Sig(a, b);
            for (std::size_t i = 0; i < K; ++i) {
                s += p.muN(i, 0) * (w.S2[i](a, b) - w.f1(i, a) * w.f1(i, b));
                for (std::size_t k = 0; k < K; ++k)
                    s += p.SNN(i, k) * w.f1(i, a) * w.f1(k, b);
            }
            out.Sxx(a, b) += s;
        }
    return out;
}

struct QdtfUpdateResult {
    double zhat = 0.0;
    double v = 0.0;
    double logl = 0.0;
};

inline QdtfUpdateResult update(QdtfState& st, const std::vector<double>& read, double z,
                               double floor_var = 0.0) {
    const std::size_t K = st.muN.nr, m = st.mx.nr;
    QdtfUpdateResult r;
    for (std::size_t a = 0; a < m; ++a) r.zhat += read[a] * st.mx(a, 0);
    std::vector<double> kx(m, 0.0), kN(K, 0.0);
    for (std::size_t a = 0; a < m; ++a)
        for (std::size_t b = 0; b < m; ++b) kx[a] += st.Sxx(a, b) * read[b];
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t a = 0; a < m; ++a) kN[i] += st.CNx(i, a) * read[a];
    r.v = floor_var;
    for (std::size_t a = 0; a < m; ++a) r.v += read[a] * kx[a];
    const double d = z - r.zhat;
    for (std::size_t i = 0; i < K; ++i) st.muN(i, 0) += kN[i] * d / r.v;
    for (std::size_t a = 0; a < m; ++a) st.mx(a, 0) += kx[a] * d / r.v;
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t j = 0; j < K; ++j) st.SNN(i, j) -= kN[i] * kN[j] / r.v;
    for (std::size_t i = 0; i < K; ++i)
        for (std::size_t a = 0; a < m; ++a) st.CNx(i, a) -= kN[i] * kx[a] / r.v;
    for (std::size_t a = 0; a < m; ++a)
        for (std::size_t b = 0; b < m; ++b) st.Sxx(a, b) -= kx[a] * kx[b] / r.v;
    r.logl = -0.5 * (std::log(2.0 * std::numbers::pi * r.v) + d * d / r.v);
    return r;
}

}  // namespace macrodr::qdtf
