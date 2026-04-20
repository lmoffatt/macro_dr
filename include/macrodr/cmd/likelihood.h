#pragma once

#include <derivative_fwd.h>
#include <distributions.h>
#include <macrodr/interface/IModel.h>
#include <macrodr/cmd/detail/write_csv_common.h>

#include <type_traits>
#include <utility>

#include "patch_model.h"
#include "qmodel.h"

namespace macrodr::cmd {

inline auto build_likelihood_function(const ModelPtr& model0,
                                      bool adaptive_approximation, bool recursive_approximation,
                                      int averaging_approximation, bool variance_approximation,
                                      bool taylor_variance_correction_approximation) {
    auto nsub = Simulation_n_sub_dt(100);
    const interface::IModel<var::Parameters_values>& model_ref = *model0;

    return merge_Maybe_variant(
        Likelihood_Model_regular<
               var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
               var::constexpr_Var_domain<bool, uses_recursive_aproximation, false,true>,
               var::constexpr_Var_domain<int, uses_averaging_aproximation,0,1, 2>,
               var::constexpr_Var_domain<bool, uses_variance_aproximation, false,true>,
               var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, false>,
               decltype(model_ref)>(model_ref, nsub,
                                 uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
        .get_variant(),
        Likelihood_Model_regular<
               var::constexpr_Var_domain<bool, uses_adaptive_aproximation, false>,
               var::constexpr_Var_domain<bool, uses_recursive_aproximation,true>,
               var::constexpr_Var_domain<int, uses_averaging_aproximation,1, 2>,
               var::constexpr_Var_domain<bool, uses_variance_aproximation, true>,
               var::constexpr_Var_domain<bool, uses_taylor_variance_correction_aproximation, true>,
               decltype(model_ref)>(model_ref, nsub,
                                 uses_adaptive_aproximation_value(adaptive_approximation),
                                 uses_recursive_aproximation_value(recursive_approximation),
                                 uses_averaging_aproximation_value(averaging_approximation),
                                 uses_variance_aproximation_value(variance_approximation),
                                 uses_taylor_variance_correction_aproximation_value(
                                     taylor_variance_correction_approximation))
        .get_variant());
;
}

using likelihood_algorithm_type = var::untransformed_type_t<decltype(build_likelihood_function(
    std::declval<const ModelPtr&>(), false, false, 2, true, false))>;

auto calculate_mlikelihood(const likelihood_algorithm_type& likelihood_algorithm,
                           const var::Parameters_transformed& par, const Experiment& e,
                           const Recording& r) -> Maybe_error<Vector_Space<logL, elogL, vlogL>>;

auto calculate_mdlikelihood(const likelihood_algorithm_type& likelihood_algorithm,
                            const var::Parameters_transformed& par, const Experiment& e,
                            const Recording& r) -> Maybe_error<dMacro_State_Hessian_minimal>;

auto calculate_mdiff_likelihood(const likelihood_algorithm_type& likelihood_algorithm,
                                const var::Parameters_transformed& par, const Experiment& e,
                                const Recording& r, double delta_param)
    -> Maybe_error<diff_Macro_State_Gradient_Hessian>;

auto calculate_mlikelihood_predictions(const likelihood_algorithm_type& likelihood_algorithm,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Recording& r)
    -> Maybe_error<Macro_State_Ev_predictions>;

inline auto calculate_simulation_mlikelihood_predictions(const likelihood_algorithm_type& likelihood_algorithm,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Simulated_Recording<var::please_include<>>& r)
    -> Maybe_error<Macro_State_Ev_predictions>
    {
        return calculate_mlikelihood_predictions   (likelihood_algorithm, par, e,
                                             get<Recording>(r()));

    }



auto calculate_mlikelihood_diagnostics(const likelihood_algorithm_type& likelihood_algorithm,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Recording& r)
    -> Maybe_error<Macro_State_Ev_diagnostic>;

inline auto calculate_simulated_mlikelihood_diagnostics(const likelihood_algorithm_type& likelihood_algorithm,
                                       const var::Parameters_transformed& par,const Experiment& e,
                                       const Simulated_Recording<var::please_include<>>& simulation)
    -> Maybe_error<Macro_State_Ev_diagnostic> {
    return calculate_mlikelihood_diagnostics(likelihood_algorithm, par, e,
                                             get<Recording>(simulation()));
}

auto calculate_mdlikelihood_predictions(const likelihood_algorithm_type& likelihood_algorithm,
                                        const var::Parameters_transformed& par, const Experiment& e,
                                        const Recording& r)
    -> Maybe_error<dMacro_State_Ev_gradient_all>;

inline auto calculate_simulation_mdlikelihood_predictions(
    const likelihood_algorithm_type& likelihood_algorithm, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_mdlikelihood_predictions(likelihood_algorithm, par, e,
                                              get<Recording>(simulation()));
}

inline auto calculate_n_simulation_mdlikelihood_predictions(
    const likelihood_algorithm_type& likelihood_algorithm, const var::Parameters_transformed& par,
    const Experiment& e, const std::vector<Simulated_Recording<var::please_include<>>>& simulation)
    -> Maybe_error<std::vector<dMacro_State_Ev_gradient_all>> {
    
    std::vector<dMacro_State_Ev_gradient_all> results;
    results.reserve(simulation.size());
    for (const  auto& sim : simulation){
        auto res = calculate_mdlikelihood_predictions(likelihood_algorithm, par, e,
                                              get<Recording>(sim()));
        if (!res){
            return res.error();
        }
        results.push_back(std::move(res.value()));        
    }
    return results;
}


auto calculate_likelihood(const ModelPtr& model0,
                          const var::Parameters_transformed& par, const Experiment& e,
                          const Recording& r, bool adaptive_approximation,
                          bool recursive_approximation, int averaging_approximation,
                          bool variance_approximation,
                          bool taylor_variance_correction_approximation)
    -> Maybe_error<Vector_Space<logL, elogL, vlogL>>;

inline auto calculate_simulation_likelihood(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<Vector_Space<logL, elogL, vlogL>> {
    return calculate_likelihood(model0, par, e, get<Recording>(simulation()),
                                adaptive_approximation, recursive_approximation,
                                averaging_approximation, variance_approximation,
                                taylor_variance_correction_approximation);
}

auto calculate_dlikelihood(const ModelPtr& model0,
                           const var::Parameters_transformed& par, const Experiment& e,
                           const Recording& r, bool adaptive_approximation,
                           bool recursive_approximation, int averaging_approximation,
                           bool variance_approximation,
                           bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Hessian_minimal>;

inline auto calculate_simulation_dlikelihood(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Hessian_minimal> {
    return calculate_dlikelihood(model0, par, e, get<Recording>(simulation()),
                                 adaptive_approximation, recursive_approximation,
                                 averaging_approximation, variance_approximation,
                                 taylor_variance_correction_approximation);
}
auto calculate_diff_likelihood(const ModelPtr& model0,
                               const var::Parameters_transformed& par, const Experiment& e,
                               const Recording& r, bool adaptive_approximation,
                               bool recursive_approximation, int averaging_approximation,
                               bool variance_approximation,
                               bool taylor_variance_correction_approximation, double delta_param)
    -> Maybe_error<diff_Macro_State_Gradient_Hessian>;

inline auto calculate_simulation_diff_likelihood(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation, double delta_param)
    -> Maybe_error<diff_Macro_State_Gradient_Hessian> {
    return calculate_diff_likelihood(model0, par, e, get<Recording>(simulation()),
                                     adaptive_approximation, recursive_approximation,
                                     averaging_approximation, variance_approximation,
                                     taylor_variance_correction_approximation, delta_param);
}

auto calculate_likelihood_predictions(const ModelPtr& model0,
                                      const var::Parameters_transformed& par, const Experiment& e,
                                      const Recording& r, bool adaptive_approximation,
                                      bool recursive_approximation, int averaging_approximation,
                                      bool variance_approximation,
                                      bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_predictions>;

inline auto calculate_simulation_likelihood_predictions(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_predictions> {
    return calculate_likelihood_predictions(model0, par, e, get<Recording>(simulation()),
                                            adaptive_approximation, recursive_approximation,
                                            averaging_approximation, variance_approximation,
                                            taylor_variance_correction_approximation);
}

auto calculate_likelihood_diagnostics(const ModelPtr& model0,
                                      const var::Parameters_transformed& par, const Experiment& e,
                                      const Recording& r, bool adaptive_approximation,
                                      bool recursive_approximation, int averaging_approximation,
                                      bool variance_approximation,
                                      bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_diagnostic>;

inline auto calculate_simulation_likelihood_diagnostics(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<Macro_State_Ev_diagnostic> {
    return calculate_likelihood_diagnostics(model0, par, e, get<Recording>(simulation()),
                                            adaptive_approximation, recursive_approximation,
                                            averaging_approximation, variance_approximation,
                                            taylor_variance_correction_approximation);
}

auto calculate_dlikelihood_predictions(const ModelPtr& model0,
                                       const var::Parameters_transformed& par, const Experiment& e,
                                       const Recording& r, bool adaptive_approximation,
                                       bool recursive_approximation, int averaging_approximation,
                                       bool variance_approximation,
                                       bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all>;

inline auto calculate_simulation_dlikelihood_predictions(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e, const Simulated_Recording<var::please_include<>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_dlikelihood_predictions(model0, par, e, get<Recording>(simulation()),
                                             adaptive_approximation, recursive_approximation,
                                             averaging_approximation, variance_approximation,
                                             taylor_variance_correction_approximation);
}

inline auto calculate_simulation_sub_dlikelihood_predictions(
    const ModelPtr& model0, const var::Parameters_transformed& par,
    const Experiment& e,
    const Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>& simulation,
    bool adaptive_approximation, bool recursive_approximation, int averaging_approximation,
    bool variance_approximation, bool taylor_variance_correction_approximation)
    -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_dlikelihood_predictions(model0, par, e, get<Recording>(simulation()),
                                             adaptive_approximation, recursive_approximation,
                                             averaging_approximation, variance_approximation,
                                             taylor_variance_correction_approximation);
}

auto calculate_dlikelihood_predictions_model(
    const std::string& model0, const var::Parameters_transformed& par, const Experiment& e,
    const Recording& r, bool adaptive_approximation, bool recursive_approximation,
    int averaging_approximation, bool variance_approximation,
    bool taylor_variance_correction_approximation) -> Maybe_error<dMacro_State_Ev_gradient_all>;

inline auto calculate_simulation_dlikelihood_predictions_model(
    const std::string& model0, const var::Parameters_transformed& par, const Experiment& e,
    const Simulated_Recording<var::please_include<>>& simulation, bool adaptive_approximation,
    bool recursive_approximation, int averaging_approximation, bool variance_approximation,
    bool taylor_variance_correction_approximation) -> Maybe_error<dMacro_State_Ev_gradient_all> {
    return calculate_dlikelihood_predictions_model(model0, par, e, get<Recording>(simulation()),
                                                   adaptive_approximation, recursive_approximation,
                                                   averaging_approximation, variance_approximation,
                                                   taylor_variance_correction_approximation);
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   Simulated_Recording<SimTag> const& simulation,
                                   TMacro_State<vVars...> const& lik, std::string path);

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   var::Indexed<Simulated_Recording<SimTag>> const& simulation,
                                   TMacro_State<vVars...> const& lik, std::string path);

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   Simulated_Recording<SimTag> const& simulation,
                                   var::Indexed<TMacro_State<vVars...>> const& lik,
                                   std::string path);

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e,
                                   var::Indexed<Simulated_Recording<SimTag>> const& simulation,
                                   var::Indexed<TMacro_State<vVars...>> const& lik,
                                   std::string path);

template <class... Vs>
inline Maybe_error<std::string> write_csv(var::Vector_Space<Vs...> const& lik, std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}

template <class T>
inline Maybe_error<std::string> write_csv(var::Indexed<T> const& indexed, std::string path) {
    return detail::write_summary_csv(indexed, std::move(path), "summary");
}

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(Experiment const& e, std::vector<Simulated_Recording<SimTag>> const& simulation,
                                   std::vector<TMacro_State<vVars...>> const& liks, std::string path);

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    Experiment const& e, var::Indexed<std::vector<Simulated_Recording<SimTag>>> const& simulation,
    std::vector<TMacro_State<vVars...>> const& liks, std::string path);

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    Experiment const& e, std::vector<Simulated_Recording<SimTag>> const& simulation,
    var::Indexed<std::vector<TMacro_State<vVars...>>> const& liks, std::string path);

template <typename SimTag, template <typename...> class TMacro_State, typename... vVars>
    requires(macrodr::has_var_c<TMacro_State<vVars...> const&, Evolution>)
Maybe_error<std::string> write_csv(
    Experiment const& e, var::Indexed<std::vector<Simulated_Recording<SimTag>>> const& simulation,
    var::Indexed<std::vector<TMacro_State<vVars...>>> const& liks, std::string path);

template <template <typename...> class TMacro_State, typename... vVars>
Maybe_error<std::string> write_csv(TMacro_State<vVars...> const& lik, std::string path);


auto calculate_boot_Likelihood_diagnostics(const std::vector<dMacro_State_Ev_gradient_all>& dy, 
               const std::vector<Simulated_Recording<var::please_include<>>>& simulation)
    -> Maybe_error<std::vector<Macro_State_Ev_diagnostic>>;


using Analisis_derivative_diagnostic_base = var::Vector_Space
<
        Probit_statistics<Moment_statistics<Sum<logL>, false>>,
        Probit_statistics<Moment_statistics<Sum<elogL>, false>>,
        Probit_statistics<Moment_statistics<Sum<macrodr::r_std>, false>>,
        Probit_statistics<Moment_statistics<Sum<macrodr::r2_std>, false>>,
        Probit_statistics<Moment_statistics<Sum<macrodr::trust_coefficient>, false>>,
        Probit_statistics<Moment_statistics<Sum<dlogL>, true>>,
        Probit_statistics<Moment_statistics<Sum<Gaussian_Fisher_Information>, false>>,
        Probit_statistics<Sum<Moment_statistics<macrodr::r_std, false>>>,
        Probit_statistics<Sum<Moment_statistics<dlogL, true>>>,
        Probit_statistics<Sum<Moment_statistics<Gaussian_Fisher_Information, false>>>,
        Probit_statistics<Information_Distortion_Matrix>,
        Probit_statistics<log_Det<Information_Distortion_Matrix>>,
        Probit_statistics<Information_Distortion_Reconstituted>,
        Probit_statistics<Sample_Distortion_Matrix>,
        Probit_statistics<log_Det<Sample_Distortion_Matrix>>,
        Probit_statistics<Correlation_Distortion_Matrix>,
        Probit_statistics<log_Det<Correlation_Distortion_Matrix>>,
        Probit_statistics<Fisher_Covariance>,
        Probit_statistics<log_Det<Fisher_Covariance>>,
        Probit_statistics<Correlation_Of<Fisher_Covariance>>,
        Probit_statistics<Eigenvalue_Spectrum<Fisher_Covariance>>,
        Probit_statistics<Effective_Rank<Fisher_Covariance>>,
        Probit_statistics<Spectrum_Condition_Number<Fisher_Covariance>>,
        Probit_statistics<Null_Space_Projector<Fisher_Covariance>>,
        Probit_statistics<Worst_Subspace_Projector<Fisher_Covariance>>,
        Probit_statistics<Distortion_Corrected_Covariance>,
        Probit_statistics<log_Det<Distortion_Corrected_Covariance>>,
        Probit_statistics<Correlation_Of<Distortion_Corrected_Covariance>>,
        Probit_statistics<Eigenvalue_Spectrum<Distortion_Corrected_Covariance>>,
        Probit_statistics<Effective_Rank<Distortion_Corrected_Covariance>>,
        Probit_statistics<Spectrum_Condition_Number<Distortion_Corrected_Covariance>>,
        Probit_statistics<Null_Space_Projector<Distortion_Corrected_Covariance>>,
        Probit_statistics<Worst_Subspace_Projector<Distortion_Corrected_Covariance>>,
        Probit_statistics<Distortion_Induced_Bias>>;

// Per-sample derived diagnostics: Sample_Distortion_Matrix and
// Distortion_Induced_Bias evaluated at each sample. Included by presets that
// ask for trace-local detail (series_var, series_cov, series_kernel_full).
using Per_sample_derived_diagnostics =
    Probit_statistics<macrodr::Evolution_of<var::Vector_Space<
        Sample_Distortion_Matrix, Distortion_Induced_Bias>>>;

// Preset 1: basic — one integral correlation lag per headline observable.
// Used for lattice-scale runs; keeps output at megabyte scale.
using Analisis_derivative_diagnostic_basic = var::concatenate_t<
    Analisis_derivative_diagnostic_base,
    var::Vector_Space<
        Probit_statistics<Report_integral<logL>>,
        Probit_statistics<Report_integral<macrodr::y_mean>>,
        Probit_statistics<Report_integral<macrodr::y_var>>,
        Probit_statistics<Report_integral<macrodr::r_std>>,
        Probit_statistics<Report_integral<dlogL>>>>;

// Preset 2: series_var — per-sample moments (diagonal variance only), plus
// forward and integral lags, for 7 leaf observables + Per_sample_derived.
using Analisis_derivative_diagnostic_series_var = var::concatenate_t<
    var::concatenate_t<
        Analisis_derivative_diagnostic_base,
        var::Vector_Space<
            Probit_statistics<Report_local_var<logL>>,
            Probit_statistics<Report_local_var<elogL>>,
            Probit_statistics<Report_local_var<macrodr::y_mean>>,
            Probit_statistics<Report_local_var<macrodr::y_var>>,
            Probit_statistics<Report_local_var<macrodr::r_std>>,
            Probit_statistics<Report_local_var<macrodr::trust_coefficient>>,
            Probit_statistics<Report_local_var<dlogL>>>>,
    var::Vector_Space<Per_sample_derived_diagnostics>>;

// Preset 3: series_cov — same as series_var but with full per-sample
// covariance blocks for vector V's; adds GFI at local_var level
// (matrix-of-matrix covariance is impractical for GFI).
using Analisis_derivative_diagnostic_series_cov = var::concatenate_t<
    var::concatenate_t<
        Analisis_derivative_diagnostic_base,
        var::Vector_Space<
            Probit_statistics<Report_local_cov<logL>>,
            Probit_statistics<Report_local_cov<elogL>>,
            Probit_statistics<Report_local_cov<macrodr::y_mean>>,
            Probit_statistics<Report_local_cov<macrodr::y_var>>,
            Probit_statistics<Report_local_cov<macrodr::r_std>>,
            Probit_statistics<Report_local_cov<macrodr::trust_coefficient>>,
            Probit_statistics<Report_local_cov<dlogL>>,
            Probit_statistics<Report_local_var<Gaussian_Fisher_Information>>>>,
    var::Vector_Space<Per_sample_derived_diagnostics>>;

// Preset 4: series_kernel — full cross-correlation kernel for the 5 core
// observables. No auxiliaries (elogL, trust, GFI), no Per_sample_derived.
// Focused single-condition mechanism view.
using Analisis_derivative_diagnostic_series_kernel = var::concatenate_t<
    Analisis_derivative_diagnostic_base,
    var::Vector_Space<
        Probit_statistics<Report_cross<logL>>,
        Probit_statistics<Report_cross<macrodr::y_mean>>,
        Probit_statistics<Report_cross<macrodr::y_var>>,
        Probit_statistics<Report_cross<macrodr::r_std>>,
        Probit_statistics<Report_cross<dlogL>>>>;

// Preset 5: series_kernel_full — exhaustive single-condition: kernel for 7
// leaf observables, GFI at local_var, and Per_sample_derived included.
using Analisis_derivative_diagnostic_series_kernel_full = var::concatenate_t<
    var::concatenate_t<
        Analisis_derivative_diagnostic_base,
        var::Vector_Space<
            Probit_statistics<Report_cross<logL>>,
            Probit_statistics<Report_cross<elogL>>,
            Probit_statistics<Report_cross<macrodr::y_mean>>,
            Probit_statistics<Report_cross<macrodr::y_var>>,
            Probit_statistics<Report_cross<macrodr::r_std>>,
            Probit_statistics<Report_cross<macrodr::trust_coefficient>>,
            Probit_statistics<Report_cross<dlogL>>,
            Probit_statistics<Report_local_var<Gaussian_Fisher_Information>>>>,
    var::Vector_Space<Per_sample_derived_diagnostics>>;


auto calculate_Likelihood_derivative_basic_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_basic;

auto calculate_Likelihood_derivative_series_var_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_series_var;

auto calculate_Likelihood_derivative_series_cov_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_series_cov;

auto calculate_Likelihood_derivative_series_kernel_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_series_kernel;

auto calculate_Likelihood_derivative_series_kernel_full_diagnostics(
    const std::vector<dMacro_State_Ev_gradient_all>& dy, std::size_t n_boostrap_samples,
    const std::set<double>& cis, std::size_t seed, std::size_t max_lag)
    -> Analisis_derivative_diagnostic_series_kernel_full;


inline Maybe_error<std::string> write_csv(Analisis_derivative_diagnostic_base const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}
inline Maybe_error<std::string> write_csv(Analisis_derivative_diagnostic_basic const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}
inline Maybe_error<std::string> write_csv(Analisis_derivative_diagnostic_series_var const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}
inline Maybe_error<std::string> write_csv(Analisis_derivative_diagnostic_series_cov const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}
inline Maybe_error<std::string> write_csv(Analisis_derivative_diagnostic_series_kernel const& lik,
                                          std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}
inline Maybe_error<std::string> write_csv(
    Analisis_derivative_diagnostic_series_kernel_full const& lik, std::string path) {
    return detail::write_summary_csv(lik, std::move(path), "summary");
}


}  // namespace macrodr::cmd
