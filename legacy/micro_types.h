#pragma once
// NOTE: <micro_full.h> is included at the very END of this file, NOT here.
// micro_full.h includes qmodel.h, which (after Macro_DMR is defined) re-includes
// micro_types.h to pick up the unified dMacro_State_Ev_gradient_* aliases. If
// we included micro_full.h up here, the chain
//   micro_types.h:2 → micro_full.h → qmodel.h:5719 → micro_types.h (skipped, in guard)
// would arrive at qmodel.h:6295's use of dMacro_State_Ev_gradient_all *before*
// micro_types.h had finished defining it. Pushing the micro_full.h include to
// the bottom guarantees the unified aliases are visible by the time qmodel.h
// reaches line 6295.
#include <variables.h>
#include <cstddef>
#include <vector>

#include "general_algorithm_on_containers.h"
// Derivative<F, X> specializations must be visible BEFORE qmodel_types.h is
// parsed, because qmodel_types.h's dMacro_State constructor body (around
// line 778) uses var::get_dx_of_dfdx(Derivative<Patch_State, …>), which
// requires Derivative<Patch_State, …> to be complete. The old include chain
// brought these in transitively via micro_full.h → qmodel.h →
// parallel_levenberg_tempering.h → variables_derivative.h; now that
// micro_full.h lives at the bottom of this file, the specialization headers
// are included explicitly here.
//   - variables_derivative.h: Derivative<Var<Id, T>, X>, Derivative<Vector_Space<…>, X>
//   - matrix_derivative.h:    Derivative<Matrix<double>, double>
//   - parameters_derivative.h: Derivative<double, Parameters_transformed>,
//                               Derivative<Matrix<double>, Parameters_transformed>
#include "variables_derivative.h"
#include "matrix_derivative.h"
#include "parameters_derivative.h"
#include "qmodel_types.h"

namespace macrodr {

class micro_Nchannels: public var::Constant<micro_Nchannels, std::size_t> {
   public:
    using var::Constant<micro_Nchannels, std::size_t>::Constant;
    friend std::string className(micro_Nchannels) { return "micro_Nchannels"; }
};




// Probability distribution over microstates.
// Shape: (1, num_full_states). One probability per row of Micro_state_Num_ch.
// Parallel to P_mean (which at the macro level has shape (1, k_states)).
class micro_P_state: public var::Var<micro_P_state, Matrix<double>> {
   public:
    using var::Var<micro_P_state, Matrix<double>>::Var;
    friend std::string className(micro_P_state) { return "micro_P_state"; }
};




class micro_State_p : public var::Var<micro_State_p, Matrix<double>> {};



class micro_g : public var::Var<micro_g, Matrix<double>> {};

class micro_g_var : public var::Var<micro_g_var, Matrix<double>> {};


class micro_P_state_t2_y0 : public var::Var<micro_P_state_t2_y0, Matrix<double>> {
   public:
    friend std::string className(micro_P_state_t2_y0) { return "micro_P_state_t2_y0"; }
};

class micro_P_state_t2_y1 : public var::Var<micro_P_state_t2_y1, Matrix<double>> {
   public:
    friend std::string className(micro_P_state_t2_y1) { return "micro_P_state_t2_y1"; }
};

class micro_P_state_t15_y0 : public var::Var<micro_P_state_t15_y0, Matrix<double>> {
   public:
    friend std::string className(micro_P_state_t15_y0) { return "micro_P_state_t15_y0"; }
};

class micro_P_state_t15_y1 : public var::Var<micro_P_state_t15_y1, Matrix<double>> {
   public:
    friend std::string className(micro_P_state_t15_y1) { return "micro_P_state_t15_y1"; }
};


class micro_P_state_t1_y1 : public var::Var<micro_P_state_t1_y1, Matrix<double>> {
   public:
    friend std::string className(micro_P_state_t1_y1) { return "micro_P_state_t1_y1"; }
};

class micro_P_state_t20_y1 : public var::Var<micro_P_state_t20_y1, Matrix<double>> {
   public:
    friend std::string className(micro_P_state_t20_y1) { return "micro_P_state_t20_y1"; }
};

class micro_P_state_t11_y0 : public var::Var<micro_P_state_t11_y0, Matrix<double>> {
   public:
    friend std::string className(micro_P_state_t11_y0) { return "micro_P_state_t11_y0"; }
};

class micro_P_state_t10_y1 : public var::Var<micro_P_state_t10_y1, Matrix<double>> {
   public:
    friend std::string className(micro_P_state_t10_y1) { return "micro_P_state_t10_y1"; }
};

class micro_P_state_0t_y0 : public var::Var<micro_P_state_0t_y0, Matrix<double>> {
   public:
    // Boundary-state prior mean over the interval: (i0,it) -> P(X0=i0, Xt=it).
    friend std::string className(micro_P_state_0t_y0) { return "micro_P_state_0t_y0"; }
};

class micro_P_state_0t_y1 : public var::Var<micro_P_state_0t_y1, Matrix<double>> {
   public:
    // Boundary-state posterior mean over the interval after conditioning on y_{0->t}.
    friend std::string className(micro_P_state_0t_y1) { return "micro_P_state_0t_y1"; }
};

class micro_P : public Var<micro_P, Matrix<double>> {
    friend std::string className(const micro_P&) { return "micro_P_ij"; }
};

class micro_P_half : public Var<micro_P_half, Matrix<double>> {
    friend std::string className(const micro_P_half&) { return "micro_Ph_ij"; }
};

class micro_P_state_to_P_mean : public Var<micro_P_state_to_P_mean, Matrix<double>> {
    friend std::string className(micro_P_state_to_P_mean) { return "micro_P_state_to_P_mean"; }
};


class micro_gmean_i : public Var<micro_gmean_i, Matrix<double>> {
    friend std::string className(micro_gmean_i) { return "micro_gmean_i"; }
};
class micro_gtotal_ij : public Var<micro_gtotal_ij, Matrix<double>> {
    friend std::string className(micro_gtotal_ij) { return "micro_gtotal_ij"; }
};
class micro_gmean_ij : public Var<micro_gmean_ij, Matrix<double>> {
    friend std::string className(micro_gmean_ij) { return "micro_gmean_ij"; }
};
class micro_gtotal_sqr_ij : public Var<micro_gtotal_sqr_ij, Matrix<double>> {
    friend std::string className(micro_gtotal_sqr_ij) { return "micro_gtotal_sqr_ij"; }
};
class micro_gsqr_i : public Var<micro_gsqr_i, Matrix<double>> {
    friend std::string className(micro_gsqr_i) { return "micro_gsqr_i"; }
};
class micro_gvar_i : public Var<micro_gvar_i, Matrix<double>> {
    friend std::string className(micro_gvar_i) { return "micro_gvar_i"; }
};
class micro_gtotal_var_ij : public Var<micro_gtotal_var_ij, Matrix<double>> {
    friend std::string className(micro_gtotal_var_ij) { return "micro_gtotal_var_ij"; }
};
class micro_gvar_ij : public Var<micro_gvar_ij, Matrix<double>> {
    friend std::string className(micro_gvar_ij) { return "micro_gvar_ij"; }
};

class micro_y_mean : public var::Var<micro_y_mean, Matrix<double>> {
    friend std::string className(micro_y_mean) { return "micro_y_mean"; }
};
class micro_y_var : public var::Var<micro_y_var, Matrix<double>> {
    friend std::string className(micro_y_var) { return "micro_y_var"; }
};
    


class micro_r_std : public var::Var<micro_r_std, double> {
    friend std::string className(micro_r_std) { return "micro_r_std"; }
};

class micro_Chi2 : public var::Var<micro_Chi2, double> {
    friend std::string className(micro_Chi2) { return "micro_Chi2"; }
};


// gaussian_logL: future diagnostic comparator — log φ(y; μ, σ²) evaluated against
// the moment-matched marginal of the micro state.  Reserved; not currently wired.
// The per-step logL that Micror folds into the running likelihood is just `logL`
// (computed in mixture form by `calculate_micro_logL`), mirroring how macror
// uses `logL` for its Gaussian-moment-match form via `calculate_logL`.
class gaussian_logL : public var::Var<gaussian_logL, double> {
    friend std::string className(gaussian_logL) { return "gaussian_logL"; }
};

using micro_Qdtg = Vector_Space<number_of_samples, min_P, micro_Nchannels, micro_P_half, micro_P_state_to_P_mean, micro_g>;

using micro_Qdtm =
    Vector_Space<number_of_samples, min_P, micro_Nchannels,micro_P, micro_P_state_to_P_mean,micro_gmean_i, micro_gvar_i>;

using micro_Qdt = Vector_Space<number_of_samples, min_P, micro_Nchannels, micro_P, micro_P_state_to_P_mean, micro_gmean_ij, micro_gvar_ij>;


struct Calc_micro_Qdt {
    friend std::string ToString(Calc_micro_Qdt) { return "Calc_micro_Qdt"; }
};

struct Calc_micro_Qdt_step {
    friend std::string ToString(Calc_micro_Qdt_step) { return "Calc_micro_Qdt_step"; }
};

struct Calc_micro_Qdtm_step {
    friend std::string ToString(Calc_micro_Qdtm_step) { return "Calc_micro_Qdtm_step"; }
};

struct Calc_micro_Qdtg_step {
    friend std::string ToString(Calc_micro_Qdtg_step) { return "Calc_micro_Qdtg_step"; }
};

struct Calc_micro_Qx {
    friend std::string ToString(Calc_micro_Qx) { return "Calc_micro_Qx"; }
};





using micro_Algo_State_Dynamic_Space =
var::push_back_var_t<Algo_State_Dynamic_Space,
logL, elogL, micro_r_std, micro_Chi2, gaussian_logL,
          micro_P,micro_P_half,micro_gmean_i,micro_gvar_i,micro_gmean_ij,
          micro_gtotal_ij,micro_P_state_t2_y0, micro_P_state_t2_y1,
          micro_P_state_t15_y0, micro_P_state_t15_y1,
                       micro_P_state_t1_y1, micro_P_state_t20_y1, micro_P_state_t11_y0, micro_P_state_t10_y1, 
                       micro_P_state_0t_y0,micro_P_state_0t_y1>; 
 
 
class micro_Algo_State_Dynamic
    : public var::Var<
          Algo_State_Dynamic, micro_Algo_State_Dynamic_Space> {
   public:
    Matrix<double> const& get_micro_P_state() const {
        if (get<micro_P_state_t2_y1>((*this)())().size() > 0) {
            return get<micro_P_state_t2_y1>((*this)())();
        }
        if (get<micro_P_state_t2_y0>((*this)())().size() > 0) {
            return get<micro_P_state_t2_y0>((*this)())();
        }
        
        return get<micro_P_state_t20_y1>((*this)())();
    }
    Matrix<double> const& get_P_mean() const {
        if (get<P_mean_t2_y1>((*this)())().size() > 0) {
            return get<P_mean_t2_y1>((*this)())();
        }
        if (get<P_mean_t2_y0>((*this)())().size() > 0) {
            return get<P_mean_t2_y0>((*this)())();
        }
        
        return get<P_mean_t20_y1>((*this)())();
    }
    SymmetricMatrix<double> const& get_P_Cov() const {
        if (get<P_Cov_t2_y1>((*this)())().size() > 0) {
            return get<P_Cov_t2_y1>((*this)())();
        }if (get<P_Cov_t2_y0>((*this)())().size() > 0) {
            return get<P_Cov_t2_y0>((*this)())();
        }

        
        return get<P_Cov_t20_y1>((*this)())();
    }
};


using micro_Algo_State_space=
var::push_back_var_t<Algo_State_space,
 logL, elogL, micro_r_std, micro_Chi2, gaussian_logL, micro_P_state>;


class micro_Algo_State
    : public var::Var<micro_Algo_State, micro_Algo_State_space> {
   public:
    using base_type =
        var::Var<micro_Algo_State, micro_Algo_State_space>;
    micro_Algo_State(const micro_Algo_State_Dynamic& p)
        : base_type{Vector_Space(get<y_mean>(p()), get<y_var>(p()), get<trust_coefficient>(p()),
                                 get<taylor_trust_coefficient>(p()), get<taylor_vSv>(p()),
                                 get<taylor_strength>(p()),
                                 get<r_std>(p()), get<Chi2>(p()),
                                 P_mean(p.get_P_mean()), P_Cov(p.get_P_Cov()),
                                 trust_mean_coefficient{}, trust_psd_coefficient{},
                                 get<logL>(p()), get<elogL>(p()),
                                 get<micro_r_std>(p()), get<micro_Chi2>(p()),
                                 get<gaussian_logL>(p()),
                                 micro_P_state(p.get_micro_P_state()))} {}

    using base_type::Var;
};


//using micro_Patch_Model = push_back_var_t<Patch_Model, micro_P_state, micro_g, micro_g_var>;

struct micro_Patch_State : public var::Var<micro_Patch_State, Vector_Space<micro_P_state>> {};

template <typename... Vars>
struct MMicro_State : public Vector_Space<logL, micro_Patch_State, Vars...> {
    MMicro_State() = default;
    MMicro_State(micro_Patch_State&& ps) { get<micro_Patch_State>((*this)) = std::move(ps); }
    MMicro_State(logL&& l, micro_Patch_State&& ps, Vars&&... vars)
        : Vector_Space<logL, micro_Patch_State, Vars...>(std::move(l), std::move(ps),
                                                         std::forward<Vars>(vars)...) {}
};

template <typename... Vars>struct dMicro_State
    : public Vector_Space<var::Derivative<logL, var::Parameters_transformed>,
                          var::Derivative<micro_Patch_State, var::Parameters_transformed>, 
                          Vars...> {
    dMicro_State(var::Derivative<micro_Patch_State, var::Parameters_transformed>&& dps) {
        auto const& dx = var::get_dx_of_dfdx(dps);

        // Seed the prior patch state (carries dx).
        get<var::Derivative<micro_Patch_State, var::Parameters_transformed>>(*this) = std::move(dps);
        // Accumulators start at zero/empty but share the same dx.
        auto seed_with_dx = [&](auto& component) {
            using Comp = std::decay_t<decltype(component)>;
            if constexpr (std::constructible_from<Comp, decltype(dx) const&>) {
                component = Comp(dx);
            } else {
                component = Comp{};
                if constexpr (requires { component.derivative().set_dx(dx); }) {
                    component.derivative().set_dx(dx);
                }
            }
        };
        seed_with_dx(get<var::Derivative<logL, var::Parameters_transformed>>(*this));
        if constexpr (sizeof...(Vars) > 0)
            ((seed_with_dx(get<Vars>(*this))), ...);
    }
    dMicro_State(var::Derivative<logL, var::Parameters_transformed>&& dl,
                 var::Derivative<micro_Patch_State, var::Parameters_transformed>&& dps,
                 Vars&&... vars)
        : Vector_Space<var::Derivative<logL, var::Parameters_transformed>,
                       var::Derivative<micro_Patch_State, var::Parameters_transformed>,
                          Vars...>(
              std::move(dl), std::move(dps),
              std::forward<Vars>(vars)...) {}
    dMicro_State() = default;       
                    
};

template <typename... Vars>
struct ddMicro_State
    : public var::Derivative<Vector_Space<logL, micro_Patch_State, Vars...>, var::Parameters_transformed> {
    ddMicro_State(var::Derivative<micro_Patch_State, var::Parameters_transformed>&& dps) {
        auto const& dx = var::get_dx_of_dfdx(dps);
        get<micro_Patch_State>(*this) = std::move(dps);

        auto seed_with_dx = [&](auto& component) {
            using Comp = std::decay_t<decltype(component)>;
            if constexpr (std::constructible_from<Comp, decltype(dx) const&>) {
                component = Comp(dx);
            } else {
                component = Comp{};
                if constexpr (requires { component.derivative().set_dx(dx); }) {
                    component.derivative().set_dx(dx);
                }
            }
        };
        seed_with_dx(get<logL>(*this));
        if constexpr (sizeof...(Vars) > 0)
            ((seed_with_dx(get<Vars>(*this))), ...);
    }
    ddMicro_State(var::Derivative<logL, var::Parameters_transformed>&& dl,
                  var::Derivative<micro_Patch_State, var::Parameters_transformed>&& dps,
                  var::Derivative_t<Vars, var::Parameters_transformed>&&... vars)
        : var::Derivative<Vector_Space<logL, micro_Patch_State, Vars...>, var::Parameters_transformed>(
              std::move(dl), std::move(dps), std::forward<Vars>(vars)...) {}
};
}
namespace var {
template <class... Vars>
struct transformation_type<macrodr::ddMicro_State<Vars...>> {
    using type = Derivative_Op<Parameters_transformed>;
};

template <class... Vars>
struct is_derivative<macrodr::ddMicro_State<Vars...>> : std::true_type {};

template <class... Vars, class... Ds>
struct dx_of_dfdx<macrodr::ddMicro_State<Vars...>, Ds...> {
    using type = Parameters_transformed;
};

template <class G, class... Vars, class... Ds>
    requires(!is_derivative_v<G>)
struct dx_of_dfdx<G, macrodr::ddMicro_State<Vars...>, Ds...> {
    using type = Parameters_transformed;
};

template <class F, class... Vars, class... Ds>
struct dx_of_dfdx<Derivative<F, Parameters_transformed>, macrodr::ddMicro_State<Vars...>, Ds...> {
    using type = Parameters_transformed;
};
}  // namespace var

namespace macrodr {

using micro_predictions_element =add_t<predictions_element,var::please_include<micro_r_std,micro_P>>;

using micro_diagnostic_element = var::please_include<logL, elogL, vlogL, micro_Algo_State_Dynamic>;

// Structured as supersets of the macro per-step elements: a macro filter writes
// only the macro slots (r_std / Derivative<r_std>) and leaves the micro-only
// slots (micro_r_std / Derivative<micro_r_std>) as zero-sized defaults; a micro
// filter populates both. This is the same empty-matrix-as-null convention used
// by micro_Algo_State_Dynamic::get_micro_P_state() / get_P_mean() at
// micro_types.h:199-228 — pick the first slot that has been populated.
using micro_gradient_minimal_element =
    add_t<gradient_minimal_element, var::please_include<micro_r_std>>;

using micro_gradient_all_element =
    add_t<gradient_all_element,
          var::please_include<var::Derivative<micro_r_std, var::Parameters_transformed>>>;

// Removed: Micro_State_minimal / Micro_State_reg / Micro_State_Ev_predictions /
// Micro_State_Ev_diagnostic. They referenced Micro_State<> from micro_full.h
// (the deprecated capital-M state). They had no users outside their own
// definitions and pulled an include cycle when qmodel.h re-includes
// micro_types.h after Macro_DMR is defined: test_micro_full.cpp →
// micro_full.h:14 → qmodel.h → micro_types.h reaches these aliases before
// micro_full.h has finished expanding Micro_State.

using dMicro_State_Ev_gradient_minimal =
    add_t<dMicro_State<>,
          var::please_include<Evolution_of<add_t<Vector_Space<>, micro_gradient_minimal_element>>,
                              evaluation_time>>;

using dMicro_State_Ev_gradient_all =
    add_t<dMicro_State<>,
          var::please_include<Evolution_of<add_t<Vector_Space<>, micro_gradient_all_element>>,
                              evaluation_time>>;

// Unified output types: same Evolution element shape as the dMicro variants
// (the per-step element is now a superset of the macro one), but carry
// Derivative<Patch_State> so the macro filter writes its native shape.
// micro algorithms convert their dMicro_State_Ev_gradient_all output to one of
// these via widen_to_dMacro_Ev_gradient_all (defined below).
using dMacro_State_Ev_gradient_minimal =
    add_t<dMacro_State<>,
          var::please_include<Evolution_of<add_t<Vector_Space<>, micro_gradient_minimal_element>>,
                              evaluation_time>>;

using dMacro_State_Ev_gradient_all =
    add_t<dMacro_State<>,
          var::please_include<Evolution_of<add_t<Vector_Space<>, micro_gradient_all_element>>,
                              evaluation_time, Gaussian_Fisher_Information>>;
// Note: Gaussian_Fisher_Information slot added at the top level so this state
// can be used directly as the gauss_newton_result evaluator type for MLE

// Per-replica per-sample numerical Fisher Information variant. Each timestep
// of the Evolution_of<> slot carries the per-step (NOT cumulative) FD-computed
// F contribution F_t = -∂²(logL_t)/∂θ². By additivity of logL across timesteps,
// Σ_t F_t equals the global F returned by
// calculate_n_simulation_mnumerical_fisher_information.
// Used by calculate_per_sample_n_simulation_mnumerical_fisher_information for
// diagnostic localization of FD instabilities in score derivatives — the bug
// hunt looking for which timestep triggers the discontinuity on specific
// parameter directions.
using dMacro_State_Ev_per_sample_F =
    add_t<dMacro_State<>,
          var::please_include<Evolution_of<Vector_Space<Likelihood_Numerical_Fisher_Information>>>>;

// Detailed per-sample diagnostic state: per-step Evolution of the significant
// recursion variables as Derivatives (logL, y_mean, y_var, P_mean, P_Cov,
// trust_coefficient). Used by calc_per_sample_numerical_fisher_information
// _detailed to dump, per sample, the value X AND the regular AD gradient ∂X/∂θ
// (over all parameters) at θ, so an intermediate variable whose gradient
// misbehaves can be spotted by comparing sick vs sane replicas.
using dMacro_State_Ev_detailed =
    add_t<dMacro_State<>,
          var::please_include<Evolution_of<add_t<Vector_Space<>, detailed_element>>,
                              evaluation_time>>;
// per-replicate analysis (Path B / "complete"). The slot is accumulated by
// `update_macro_state` via has_var_c<>, no extra integration code needed.

// Adapter: micro_*_element is a superset of gradient_*_element, so the
// per-step Evolution slot has the same shape on both sides — copy it through.
// The carried prior is dropped (Derivative<micro_Patch_State> → default
// Derivative<Patch_State>); downstream consumers read only the Evolution slot
// and the carried logL gradient. This is the empty-matrix-as-null convention:
// "Patch_State unpopulated" signals "consult the Evolution slot instead."
inline auto widen_to_dMacro_Ev_gradient_all(dMicro_State_Ev_gradient_all&& src)
    -> dMacro_State_Ev_gradient_all {
    dMacro_State_Ev_gradient_all out;
    using DLogL = var::Derivative<logL, var::Parameters_transformed>;
    using EvSlot = Evolution_of<add_t<Vector_Space<>, micro_gradient_all_element>>;
    get<DLogL>(out) = std::move(get<DLogL>(src));
    get<EvSlot>(out) = std::move(get<EvSlot>(src));
    get<evaluation_time>(out) = std::move(get<evaluation_time>(src));
    return out;
}



}  // namespace macrodr

// micro_full.h depends on the types defined above (and on qmodel.h's Macro_DMR
// once qmodel.h re-includes micro_types.h via its end-of-file include). By
// pulling micro_full.h here, all types that micro_full.h or its transitively
// included qmodel.h needs from micro_types.h are already defined.
#include <micro_full.h>
