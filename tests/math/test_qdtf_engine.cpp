// Gates for legacy/qdtf_engine.h — the spectral window producer (the math of
// calc_Qdtf_eig) and the augmented recursion of the Bessel member.
//
// Gate map:
//   [tables]   spectral window tables (F1, S2, Phi, Sig, P) against direct
//              quadrature of their defining integrals on the real coupled
//              observation system, for the native (K+4), grouped (K+5) and
//              box (K+1) mode sets, symmetric and asymmetric rates.
//   [anchor]   box mode set reproduces the classic Ee/E3 route of
//              calc_Qdt_eig (gtotal_ij and the collapsed marginal second
//              moment) to machine precision — the ν = 0 regression anchor.
//   [assembly] predict/update against an independently coded boundary-state
//              (Σ^bnd) assembly (mr_vs_ir_from_macror.md, count scale),
//              two windows deep with an update in between.
//   [run]      multi-window recursion smoke: finite logL, positive v,
//              symmetric covariances, simplex-consistent μ.
//
// The Monte Carlo validation of the same equations lives in
// theory/macroir/notes/bessel_reference/bessel_oracle.py (gates O1-O5).

#include <cmath>
#include <complex>
#include <vector>

#include "catch_amalgamated.hpp"
#include "qdtf_engine.h"

using namespace macrodr::qdtf;
using macrodr::acqf::FilterRealization;
using macrodr::acqf::make_bessel_filter;
using std::abs;
using std::size_t;

namespace {

//------------------------------------------------------------------ helpers
Mat expm_mat(const Mat& M) {  // scaling-and-squaring Taylor, small dense
    double nrm = 0.0;
    for (size_t i = 0; i < M.nr; ++i) {
        double r = 0.0;
        for (size_t j = 0; j < M.nc; ++j) r += std::abs(M(i, j));
        nrm = std::max(nrm, r);
    }
    int s = std::max(0, int(std::ceil(std::log2(std::max(1.0, nrm)))) + 1);
    const double sc = std::ldexp(1.0, -s);
    Mat A(M.nr, M.nc);
    for (size_t i = 0; i < M.d.size(); ++i) A.d[i] = M.d[i] * sc;
    Mat E(M.nr, M.nc), P = A;
    for (size_t i = 0; i < M.nr; ++i) E(i, i) = 1.0;
    double fact = 1.0;
    for (int k = 1; k <= 18; ++k) {
        fact *= k;
        for (size_t i = 0; i < E.d.size(); ++i) E.d[i] += P.d[i] / fact;
        P = mul(P, A);
    }
    for (int k = 0; k < s; ++k) E = mul(E, E);
    return E;
}

struct ChannelEig {  // K = 2 two-state channel, P(t) = V diag(e^{λt}) U
    Mat V, U;
    std::vector<double> lam;
    Mat Q;
};

ChannelEig make_channel(double k1, double k2) {
    ChannelEig c;
    c.Q = Mat(2, 2);
    c.Q(0, 0) = -k1;
    c.Q(0, 1) = k1;
    c.Q(1, 0) = k2;
    c.Q(1, 1) = -k2;
    c.lam = {0.0, -(k1 + k2)};
    c.V = Mat(2, 2);
    c.V(0, 0) = 1.0;
    c.V(0, 1) = k1;
    c.V(1, 0) = 1.0;
    c.V(1, 1) = -k2;
    c.U = Mat(2, 2);
    const double s = 1.0 / (k1 + k2);
    c.U(0, 0) = k2 * s;
    c.U(0, 1) = k1 * s;
    c.U(1, 0) = s;
    c.U(1, 1) = -s;
    return c;
}

// real coupled observation system (A, b) for a mode set: pair blocks plus,
// when the accumulator is present, the row du/dt = (filter output) = Σ u_k
// (grouped) or du/dt = input (box: no pairs, b_u = 1).
void real_obs_system(const QdtfModeSet& ms, Mat& A, std::vector<double>& b) {
    const size_t m = ms.m();
    A = Mat(m, m);
    b.assign(m, 0.0);
    for (size_t k = 0; k < ms.n_pairs; ++k) {
        const double sg = -ms.q[2 * k].real(), om = ms.q[2 * k].imag();
        A(2 * k, 2 * k) = -sg;
        A(2 * k, 2 * k + 1) = -om;
        A(2 * k + 1, 2 * k) = om;
        A(2 * k + 1, 2 * k + 1) = -sg;
        b[2 * k] = 2.0 * ms.rho[2 * k].real();
        b[2 * k + 1] = 2.0 * ms.rho[2 * k].imag();
    }
    if (ms.has_acc) {
        const size_t u = ms.acc_row();
        if (ms.n_pairs > 0) {
            for (size_t k = 0; k < ms.n_pairs; ++k) A(u, 2 * k) = 1.0;  // output row
        } else {
            b[u] = 1.0;  // box: accumulate the raw input
        }
    }
}

Mat Pmat(const ChannelEig& c, double t) {
    Mat Qt(2, 2);
    for (size_t i = 0; i < 4; ++i) Qt.d[i] = c.Q.d[i] * t;
    return expm_mat(Qt);
}

// G(s) = e^{A(Δ−s)} b : the loading of input at time s onto the state at Δ
std::vector<double> Gvec(const Mat& A, const std::vector<double>& b, double dt, double s) {
    Mat As(A.nr, A.nc);
    for (size_t i = 0; i < A.d.size(); ++i) As.d[i] = A.d[i] * (dt - s);
    Mat E = expm_mat(As);
    std::vector<double> g(A.nr, 0.0);
    for (size_t i = 0; i < A.nr; ++i)
        for (size_t j = 0; j < A.nc; ++j) g[i] += E(i, j) * b[j];
    return g;
}

}  // namespace

TEST_CASE("spectral window tables match quadrature", "[qdtf][tables]") {
    const double S0 = 2e-5;
    auto filt = make_bessel_filter(4, 10e3);
    struct Case {
        QdtfModeSet ms;
        double dt;
        const char* name;
    };
    std::vector<Case> cases;
    cases.push_back({make_native_modes(filt), 20e-6, "native"});
    cases.push_back({make_grouped_modes(filt), 100e-6, "grouped"});
    cases.push_back({make_box_modes(), 20e-6, "box"});

    for (auto rates : {std::pair<double, double>{100.0, 100.0}, {37.0, 211.0}}) {
        auto ch = make_channel(rates.first, rates.second);
        const std::vector<double> gamma{0.0, 1.0};
        for (auto& cs : cases) {
            INFO(cs.name << " k1=" << rates.first);
            const double dt = cs.dt;
            auto w = calc_qdtf_window(ch.V, ch.lam, ch.U, gamma, dt, cs.ms, S0);
            const size_t m = cs.ms.m();

            // P
            Mat Pref = Pmat(ch, dt);
            for (size_t i = 0; i < 4; ++i) CHECK(abs(w.P.d[i] - Pref.d[i]) < 1e-12);

            Mat A;
            std::vector<double> b;
            real_obs_system(cs.ms, A, b);

            // Phi = expm(A dt), closed form vs expm
            {
                Mat Adt(m, m);
                for (size_t i = 0; i < A.d.size(); ++i) Adt.d[i] = A.d[i] * dt;
                Mat E = expm_mat(Adt);
                for (size_t i = 0; i < m * m; ++i) CHECK(abs(w.Phi.d[i] - E.d[i]) < 1e-11);
            }

            // Sig vs Van Loan
            {
                Mat M(2 * m, 2 * m);
                for (size_t i = 0; i < m; ++i)
                    for (size_t j = 0; j < m; ++j) {
                        M(i, j) = -A(i, j) * dt;
                        M(i, m + j) = b[i] * b[j] * dt;
                        M(m + i, m + j) = A(j, i) * dt;
                    }
                Mat E = expm_mat(M);
                double scale = 1e-300;
                for (size_t i = 0; i < m * m; ++i) scale = std::max(scale, abs(w.Sig.d[i]));
                for (size_t a = 0; a < m; ++a)
                    for (size_t bb = 0; bb < m; ++bb) {
                        double s = 0.0;
                        for (size_t k = 0; k < m; ++k) s += E(m + k, m + a) * E(k, m + bb);
                        // 1e-7·scale: the modal partial-fraction assembly of the
                        // Gramian cancels ~4 digits at Bessel residue magnitudes
                        CHECK(abs(w.Sig(a, bb) - S0 * s) < 1e-7 * scale);
                    }
            }

            // F1[i](j,·) = ∫ [P(s) diag(γ) P(Δ−s)]_{ij} G(s) ds   (Simpson)
            {
                const int M = 800;
                const double h = dt / M;
                std::vector<Mat> Ps(M + 1);
                std::vector<std::vector<double>> Gs(M + 1);
                for (int a = 0; a <= M; ++a) {
                    Ps[a] = Pmat(ch, a * h);
                    Gs[a] = Gvec(A, b, dt, a * h);
                }
                double scale = 1e-300;
                for (size_t i = 0; i < 2; ++i)
                    for (size_t z = 0; z < 2 * m; ++z) scale = std::max(scale, abs(w.F1[i].d[z]));
                for (size_t i = 0; i < 2; ++i)
                    for (size_t j = 0; j < 2; ++j)
                        for (size_t a = 0; a < m; ++a) {
                            double acc = 0.0;
                            for (int n = 0; n <= M; ++n) {
                                double corr = 0.0;
                                for (size_t k = 0; k < 2; ++k)
                                    corr += Ps[n](i, k) * gamma[k] * Ps[M - n](k, j);
                                const double wgt = (n == 0 || n == M) ? 1.0
                                                   : (n % 2 ? 4.0 : 2.0);
                                acc += wgt * corr * Gs[n][a];
                            }
                            acc *= h / 3.0;
                            CHECK(abs(w.F1[i](j, a) - acc) < 1e-8 * scale);
                        }
            }

            // S2[i] = ∫∫_{s<t} (G(s)G(t)ᵀ + G(t)G(s)ᵀ)·E[γ_s γ_t | i] ds dt  (trapezoid)
            {
                const int M = 240;
                const double h = dt / M;
                std::vector<Mat> Ps(M + 1);
                std::vector<std::vector<double>> Gs(M + 1);
                for (int a = 0; a <= M; ++a) {
                    Ps[a] = Pmat(ch, a * h);
                    Gs[a] = Gvec(A, b, dt, a * h);
                }
                double scale = 1e-300;
                for (size_t i = 0; i < 2; ++i)
                    for (size_t z = 0; z < m * m; ++z) scale = std::max(scale, abs(w.S2[i].d[z]));
                for (size_t i = 0; i < 2; ++i) {
                    Mat acc(m, m);
                    for (int na = 0; na <= M; ++na) {
                        const double wa = (na == 0 || na == M) ? 0.5 : 1.0;
                        for (int nb = na; nb <= M; ++nb) {
                            const double wb = (nb == 0 || nb == M) ? 0.5 : 1.0;
                            double corr = 0.0;  // E[γ(s_na) γ(s_nb) | i], ordered
                            for (size_t k = 0; k < 2; ++k) {
                                if (gamma[k] == 0.0)
                                    continue;
                                double inner = 0.0;
                                for (size_t l = 0; l < 2; ++l)
                                    inner += Ps[nb - na](k, l) * gamma[l];
                                corr += Ps[na](i, k) * gamma[k] * inner;
                            }
                            const double f = wa * wb * h * h * corr * (na == nb ? 0.5 : 1.0);
                            for (size_t a = 0; a < m; ++a)
                                for (size_t bb = 0; bb < m; ++bb)
                                    acc(a, bb) += f * (Gs[na][a] * Gs[nb][bb] +
                                                       Gs[nb][a] * Gs[na][bb]);
                        }
                    }
                    for (size_t z = 0; z < m * m; ++z)
                        CHECK(abs(w.S2[i].d[z] - acc.d[z]) < 5e-4 * scale);
                }
            }
        }
    }
}

TEST_CASE("box mode set reproduces the classic Ee/E3 route", "[qdtf][anchor]") {
    auto ch = make_channel(100.0, 100.0);
    const std::vector<double> gamma{0.0, 1.0};
    const double dt = 20e-6;
    auto ms = make_box_modes();
    auto w = calc_qdtf_window(ch.V, ch.lam, ch.U, gamma, dt, ms, 0.0);

    // Mg = U diag(γ) V, g̃ = U γ (as in calc_Qdt_eig's WgV)
    Mat Mg(2, 2);
    std::vector<double> gt(2, 0.0);
    for (size_t a = 0; a < 2; ++a) {
        for (size_t b = 0; b < 2; ++b) {
            double s = 0.0;
            for (size_t k = 0; k < 2; ++k) s += ch.U(a, k) * gamma[k] * ch.V(k, b);
            Mg(a, b) = s;
        }
        for (size_t k = 0; k < 2; ++k) gt[a] += ch.U(a, k) * gamma[k];
    }
    auto Ee = [](double x, double y) {
        return std::abs(x - y) < 1e-12 ? std::exp(x) : (std::exp(x) - std::exp(y)) / (x - y);
    };
    // gtotal_ij = V (Mg ∘ Ee) U : interval-AVERAGED first moment; engine F1/Δ
    for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 2; ++j) {
            double g = 0.0;
            for (size_t a = 0; a < 2; ++a)
                for (size_t b = 0; b < 2; ++b)
                    g += ch.V(i, a) * Mg(a, b) * Ee(ch.lam[a] * dt, ch.lam[b] * dt) * ch.U(b, j);
            CHECK(abs(w.F1[i](j, 0) / dt - g) < 1e-13 * (1.0 + abs(g)));
        }
    // marginal second moment: classic full contraction 2·Σ V Mg Mg (U𝟙) E3
    // versus the engine's collapsed g̃·dd2(·,·,0) — the collapse identity live
    auto E3c = [&](double x, double y, double z) {
        return dd2(cdouble(x, 0), cdouble(y, 0), cdouble(z, 0)).real();
    };
    std::vector<double> U1(2, 0.0);
    for (size_t a = 0; a < 2; ++a)
        for (size_t k = 0; k < 2; ++k) U1[a] += ch.U(a, k);
    for (size_t i = 0; i < 2; ++i) {
        double classic = 0.0;
        for (size_t a = 0; a < 2; ++a)
            for (size_t b = 0; b < 2; ++b)
                for (size_t c = 0; c < 2; ++c)
                    classic += 2.0 * ch.V(i, a) * Mg(a, b) * Mg(b, c) * U1[c] *
                               E3c(ch.lam[a] * dt, ch.lam[b] * dt, ch.lam[c] * dt);
        classic *= dt * dt;
        CHECK(abs(w.S2[i](0, 0) - classic) < 1e-12 * (1.0 + abs(classic)));
    }
}

TEST_CASE("predict/update match the boundary-state assembly", "[qdtf][assembly]") {
    auto ch = make_channel(100.0, 100.0);
    const std::vector<double> gamma{0.0, 1.0};
    const double dt = 20e-6, S0 = 2e-5, Nch = 10.0;
    auto ms = make_box_modes();
    auto w = calc_qdtf_window(ch.V, ch.lam, ch.U, gamma, dt, ms, S0);
    auto read = read_vector(ms, dt);

    auto st = make_initial_state({0.5, 0.5}, Nch, ms.m());
    for (int win = 0; win < 2; ++win) {
        // independently coded Σ^bnd reference from the same spectral tables
        // (Γ̄ = F1/(PΔ); V̄ enters only through the S2 marginal, so the whole
        // reference is expressible in engine tables — the four-boxes algebra)
        Mat Gbar(2, 2);
        for (size_t i = 0; i < 2; ++i)
            for (size_t j = 0; j < 2; ++j) Gbar(i, j) = w.F1[i](j, 0) / (w.P(i, j) * dt);
        double yhat_ref = 0.0, v_ref = S0 / dt;
        std::vector<double> gbar0(2, 0.0);
        for (size_t i = 0; i < 2; ++i)
            for (size_t j = 0; j < 2; ++j) {
                yhat_ref += st.muN(i, 0) * w.P(i, j) * Gbar(i, j);
                gbar0[i] += w.P(i, j) * Gbar(i, j);
            }
        // Γ̄ᵀ Cov(B) Γ̄  +  Σ μP·V̄, via the contracted route. The +Σ μ_ij Γ̄²
        // diagonal of Cov(B) cancels exactly against the −Σ μ_ij Γ̄² inside
        // Σ μP·V̄ (the four-boxes bookkeeping), so neither is written out.
        for (size_t i = 0; i < 2; ++i)
            for (size_t k = 0; k < 2; ++k) {
                const double SmD = st.SNN(i, k) - (i == k ? st.muN(i, 0) : 0.0);
                v_ref += gbar0[i] * SmD * gbar0[k];
            }
        for (size_t i = 0; i < 2; ++i) v_ref += st.muN(i, 0) * w.S2[i](0, 0) / (dt * dt);
        // κ_j = Σ_i P_ij (SmD ḡ⁰)_i + Σ_i μ_i P_ij Γ̄_ij
        std::vector<double> kap_ref(2, 0.0);
        for (size_t j = 0; j < 2; ++j) {
            double s = 0.0;
            for (size_t i = 0; i < 2; ++i) {
                double smdg = 0.0;
                for (size_t k = 0; k < 2; ++k)
                    smdg += (st.SNN(i, k) - (i == k ? st.muN(i, 0) : 0.0)) * gbar0[k];
                s += w.P(i, j) * smdg + st.muN(i, 0) * w.P(i, j) * Gbar(i, j);
            }
            kap_ref[j] = s;
        }

        auto stp = predict(st, w, int(ms.acc_row()));
        QdtfState probe = stp;
        auto res = update(probe, read, yhat_ref);
        std::vector<double> kN(2, 0.0);
        for (size_t i = 0; i < 2; ++i) kN[i] = stp.CNx(i, 0) * read[0];

        INFO("window " << win);
        CHECK(abs(res.zhat - yhat_ref) < 1e-11 * (1.0 + abs(yhat_ref)));
        CHECK(abs(res.v - v_ref) < 1e-11 * (1.0 + abs(v_ref)));
        for (size_t j = 0; j < 2; ++j)
            CHECK(abs(kN[j] - kap_ref[j]) < 1e-11 * (1.0 + abs(kap_ref[j])));

        st = stp;
        update(st, read, yhat_ref + 0.3 * std::sqrt(res.v));
    }
}

TEST_CASE("multi-window recursion runs clean for all three reads", "[qdtf][run]") {
    auto filt = make_bessel_filter(4, 10e3);
    auto ch = make_channel(100.0, 100.0);
    const std::vector<double> gamma{0.0, 1.0};
    const double S0 = 2e-5, Nch = 10.0;

    struct Cfg {
        QdtfModeSet ms;
        double dt;
        const char* name;
    };
    std::vector<Cfg> cfgs;
    cfgs.push_back({make_native_modes(filt), 20e-6, "native"});
    cfgs.push_back({make_grouped_modes(filt), 100e-6, "grouped"});
    cfgs.push_back({make_box_modes(), 20e-6, "box"});

    for (auto& cf : cfgs) {
        INFO(cf.name);
        auto w = calc_qdtf_window(ch.V, ch.lam, ch.U, gamma, cf.dt, cf.ms, S0);
        auto read = read_vector(cf.ms, cf.dt);
        auto st = make_initial_state({0.5, 0.5}, Nch, cf.ms.m());
        const int reset = cf.ms.has_acc ? int(cf.ms.acc_row()) : -1;
        double logl = 0.0;
        for (int j = 0; j < 50; ++j) {
            st = predict(st, w, reset);
            auto r = update(st, read, 4.9 + 0.01 * j);
            CHECK(r.v > 0.0);
            logl += r.logl;
        }
        CHECK(std::isfinite(logl));
        // μ stays a (count-scaled) probability vector
        double tot = 0.0;
        for (size_t i = 0; i < 2; ++i) {
            CHECK(st.muN(i, 0) > -1e-9);
            tot += st.muN(i, 0);
        }
        CHECK(abs(tot - Nch) < 1e-8);
        // covariances stay symmetric; variances nonnegative
        for (size_t a = 0; a < cf.ms.m(); ++a) {
            CHECK(st.Sxx(a, a) > -1e-12);
            for (size_t b = 0; b < cf.ms.m(); ++b)
                CHECK(abs(st.Sxx(a, b) - st.Sxx(b, a)) < 1e-9);
        }
    }
}

TEST_CASE("native read on the grouped layout equals the native layout", "[qdtf][run]") {
    // The accumulator does not feed back into the filter pairs, so sampling
    // the filter output must give the same predictive on either layout.
    auto filt = make_bessel_filter(4, 10e3);
    auto ch = make_channel(100.0, 100.0);
    const std::vector<double> gamma{0.0, 1.0};
    const double S0 = 2e-5, Nch = 10.0, dt = 20e-6;

    auto msn = make_native_modes(filt);
    auto msg = make_grouped_modes(filt);
    auto wn = calc_qdtf_window(ch.V, ch.lam, ch.U, gamma, dt, msn, S0);
    auto wg = calc_qdtf_window(ch.V, ch.lam, ch.U, gamma, dt, msg, S0);
    auto stn = make_initial_state({0.5, 0.5}, Nch, msn.m());
    auto stg = make_initial_state({0.5, 0.5}, Nch, msg.m());
    stn = predict(stn, wn);
    stg = predict(stg, wg, int(msg.acc_row()));
    auto rn = read_vector(msn, dt);
    auto rg = native_read_vector(msg);
    QdtfState pn = stn, pg = stg;
    auto an = update(pn, rn, 5.0);
    auto ag = update(pg, rg, 5.0);
    CHECK(abs(an.zhat - ag.zhat) < 1e-11 * (1.0 + abs(an.zhat)));
    CHECK(abs(an.v - ag.v) < 1e-11 * (1.0 + abs(an.v)));
}
