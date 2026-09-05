#pragma once

// qdtf_member.h — the Bessel member's likelihood loop over a recording,
// wiring the verified engine (legacy/qdtf_engine.h) into macro_dr's model /
// experiment / recording types.
//
// What this member is: the MacroIR recursion with the amplifier's state
// carried inside it (theory/macroir/notes/macroir_bessel_plan.md; block
// equations and their Monte Carlo validation in theory/macroir/notes/
// bessel_reference/). n_poles = 0 selects the box configuration, which is
// mathematically the existing av=2 member and serves as the regression
// anchor (tests/macroir/test_qdtf_member.cpp); n_poles = 4 or 8 selects the
// Bessel filter with the given −3 dB cutoff.
//
// Design choices, and their reasons:
//   * One state layout for the whole record (filter pairs + accumulator).
//     Native-rate windows read the filter output; grouped windows read the
//     accumulator average. The accumulator is reset at every window start
//     (the recorder's own reset), so mixing window kinds is exact.
//   * Windows containing agonist sub-steps are handled by running the
//     prediction per sub-step and updating once at the window end. Carrying
//     the state composes the sub-intervals exactly; no table-level monoid
//     is needed.
//   * The filter has been stationary since before the record, so the state
//     is warmed up by predict-only windows at the first agonist before the
//     first sample (filter memory ≈ 0.69/f_c; the warm-up covers ≈ 20/f_c).
//     The baseline current passes the filter with H(0) = 1 exactly, so it
//     enters as an offset on the observed value, never as filter input.
//   * Current_Noise is read as the pre-filter power spectral density S0
//     (same number as the box member's convention, where the per-sample
//     variance is S0/Δ). Pink_Noise stays an additive per-observation floor.
//     Proportional_Noise is NOT inherited: its box scaling does not
//     transfer under a filter (plan, section 8).
//   * Plain (non-derivative) path only for now; the engine is double
//     precision. Fitting via gradients and thermodynamic evidence hook in
//     later; direct evaluation, MCMC-style samplers and quadrature over the
//     R bindings need only this entry point.
//   * Time-varying N_Ch_mean is rejected explicitly rather than silently
//     mishandled: the engine state is at count scale for one fixed N.

#include <cmath>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "qdtf_engine.h"
#include "qmodel.h"

namespace macrodr {

template <class Model, class C_Parameters>
auto qdtf_log_Likelihood(const Model& model, const C_Parameters& par, const Recording& y,
                         const Experiment& e, int n_poles, double cutoff_hz)
    -> Maybe_error<logLs> {
    auto Maybe_m = model(par);
    if (!is_valid(Maybe_m)) {
        return get_error(Maybe_m);
    }
    auto m = std::move(get_value(Maybe_m));
    auto md = Macro_DMR{};
    const double fs = get<Frequency_of_Sampling>(e).value();

    auto Maybe_ini = md.init(m);
    if (!Maybe_ini)
        return Maybe_ini.error();
    auto& t_ini = Maybe_ini.value();
    auto& p0 = get<P_mean>(t_ini())();
    auto& C0 = get<P_Cov>(t_ini())();
    const std::size_t K = p0.size();

    qdtf::QdtfModeSet ms;
    if (n_poles > 0) {
        auto filt = acqf::make_bessel_filter(n_poles, cutoff_hz);
        if (auto err = acqf::validate(filt); err.has_value())
            return error_message("qdtf filter construction: " + err.value());
        ms = qdtf::make_grouped_modes(filt);
    } else {
        ms = qdtf::make_box_modes();
    }
    const std::size_t mst = ms.m();
    const int reset = int(ms.acc_row());

    const double S0 = get<Current_Noise>(m).value();
    const double pink = get<Pink_Noise>(m).value();
    const double baseline = get<Current_Baseline>(m)();

    auto& t_g = get<g>(m)();
    std::vector<double> gam(K);
    for (std::size_t k = 0; k < K; ++k) gam[k] = t_g[k];

    auto Nchs = get<N_Ch_mean>(m)();
    const double Nch = Nchs[0];
    for (std::size_t i = 0; i < Nchs.size(); ++i)
        if (std::abs(Nchs[i] - Nch) > 1e-9 * (1.0 + std::abs(Nch)))
            return error_message("qdtf member: time-varying N_Ch_mean is not supported yet");

    // per-agonist channel eigendecomposition, copied into engine matrices
    struct EigR {
        qdtf::Mat V, U;
        std::vector<double> lam;
    };
    std::map<double, EigR> eig_cache;
    auto get_eig = [&](const Agonist_concentration& x) -> Maybe_error<const EigR*> {
        auto it = eig_cache.find(x());
        if (it != eig_cache.end())
            return &it->second;
        auto Me = md.calc_eigen(m, x);
        if (!Me)
            return Me.error();
        auto& tV = get<V>(Me.value())();
        auto& tW = get<W>(Me.value())();  // W = V⁻¹; the engine calls it U
        auto& tl = get<lambda>(Me.value())();
        EigR r;
        r.V = qdtf::Mat(K, K);
        r.U = qdtf::Mat(K, K);
        r.lam.resize(K);
        for (std::size_t i = 0; i < K; ++i) {
            r.lam[i] = tl[i];
            for (std::size_t j = 0; j < K; ++j) {
                r.V(i, j) = tV(i, j);
                r.U(i, j) = tW(i, j);
            }
        }
        auto pos = eig_cache.emplace(x(), std::move(r)).first;
        return &pos->second;
    };

    // window tables memoized on (Δ, agonist); the pole set is fixed per run
    std::map<std::pair<double, double>, qdtf::QdtfWindow> win_cache;
    auto get_window = [&](double dt, const Agonist_concentration& x)
        -> Maybe_error<const qdtf::QdtfWindow*> {
        const auto key = std::make_pair(dt, x());
        auto it = win_cache.find(key);
        if (it != win_cache.end())
            return &it->second;
        auto Me = get_eig(x);
        if (!Me)
            return Me.error();
        auto pos = win_cache
                       .emplace(key, qdtf::calc_qdtf_window(Me.value()->V, Me.value()->lam,
                                                            Me.value()->U, gam, dt, ms, S0))
                       .first;
        return &pos->second;
    };

    // initial state at count scale, from the model's own initial law
    qdtf::QdtfState st;
    st.muN = qdtf::Mat(K, 1);
    st.SNN = qdtf::Mat(K, K);
    st.mx = qdtf::Mat(mst, 1);
    st.CNx = qdtf::Mat(K, mst);
    st.Sxx = qdtf::Mat(mst, mst);
    for (std::size_t i = 0; i < K; ++i) {
        st.muN(i, 0) = Nch * p0[i];
        for (std::size_t j = 0; j < K; ++j) st.SNN(i, j) = Nch * C0(i, j);
    }

    auto& steps = get<Recording_conditions>(e)();
    const std::size_t n_steps = std::min(y().size(), steps.size());

    // warm-up: the amplifier has been at the pre-record condition since long
    // before the first sample, so bring the filter blocks to stationarity
    if (n_poles > 0 && n_steps > 0) {
        const Agonist_concentration* x0 = nullptr;
        const Agonist_evolution& first_ev = get<Agonist_evolution>(steps[0]);
        for (std::size_t k = 0; k < first_ev().size() && x0 == nullptr; ++k)
            if (get<number_of_samples>(first_ev()[k])() > 0.0)
                x0 = &get<Agonist_concentration>(first_ev()[k]);
        if (x0 != nullptr) {
            auto Mw = get_window(1.0 / cutoff_hz, *x0);
            if (!Mw)
                return Mw.error();
            for (int w = 0; w < 20; ++w) st = qdtf::predict(st, *Mw.value(), reset);
        }
    }

    double sum_logl = 0.0;
    for (std::size_t i_step = 0; i_step < n_steps; ++i_step) {
        const Agonist_evolution& t_step = get<Agonist_evolution>(steps[i_step]);
        double dt_total = 0.0, nsamp_total = 0.0;
        bool first_sub = true;
        for (std::size_t k = 0; k < t_step().size(); ++k) {
            auto& sub = t_step()[k];
            const double ns = get<number_of_samples>(sub)();
            if (ns <= 0.0)
                continue;
            const double dt = ns / fs;
            auto Mw = get_window(dt, get<Agonist_concentration>(sub));
            if (!Mw)
                return Mw.error();
            st = qdtf::predict(st, *Mw.value(), first_sub ? reset : -1);
            first_sub = false;
            dt_total += dt;
            nsamp_total += ns;
        }
        const double yv = y()[i_step].value();
        if (!std::isnan(yv) && dt_total > 0.0) {
            std::vector<double> rd;
            if (n_poles > 0 && std::abs(nsamp_total - 1.0) < 1e-9)
                rd = qdtf::native_read_vector(ms);
            else {
                rd.assign(mst, 0.0);
                rd[std::size_t(reset)] = 1.0 / dt_total;
            }
            auto res = qdtf::update(st, rd, yv - baseline, pink);
            if (!std::isfinite(res.v) || res.v <= 0.0)
                return error_message("qdtf member: invalid predictive variance at step " +
                                     std::to_string(i_step));
            sum_logl += res.logl;
        }
    }
    return logLs(logL(sum_logl), elogL(0.0), vlogL(0.0));
}

}  // namespace macrodr
