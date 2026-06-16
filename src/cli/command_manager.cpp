#include <experiment.h>
#include <macrodr/cli/command_manager.h>
#include <macrodr/cmd/cli_meta.h>
#include <macrodr/cmd/indexed_construction.h>
#include <macrodr/cmd/likelihood.h>
#include <macrodr/cmd/load_experiment.h>
#include <macrodr/cmd/load_model.h>
#include <macrodr/cmd/load_parameters.h>
#include <macrodr/cmd/patch_model.h>
#include <macrodr/cmd/simulate.h>
#include <maybe_error.h>
#include <parameters.h>
#include <qmodel.h>
#include <variables.h>

#include <cstddef>
#include <set>

// Legacy registry builders used until migrated
#include "CLI_function_table.h"
#include "CLI_likelihood.h"
#include "CLI_macro_dr_base.h"
#include "CLI_thermo_evidence_dts.h"
#include "function_builder.h"
#include "lapack_headers.h"

// ---------------------------------------------------------------------------
// JSON (de)serialization stub for the DEFERRED Probit_Samples_at_Group_Size
// slot of MLE_Group_Cloud. That slot is a std::vector<Probit_Sample_Record>,
// and Probit_Sample_Record has no JSON conversion. The DSL registers every
// command return type for JSON round-trip (Compiler::ensure_type_registered),
// which would otherwise try to (de)serialize each Probit_Sample_Record and fail
// (the generic vector path does `T element{}`, ill-formed for the Matrix member,
// and there is no element from_json). In the driver-first cut this slot is
// ALWAYS EMPTY, so trivial (de)serialization is both sufficient and correct.
// Found by ADL on Probit_Sample_Record's namespace (macrodr); more specialized
// than the generic std::vector<T> overloads, so it intercepts before recursion.
// Replace with real serialization when the representative-group snapshotting
// (the rest of the figure_2 battery) is implemented.
namespace macrodr {
template <class State>
inline io::json::Json to_json(const std::vector<Probit_Sample_Record<State>>&,
                              io::json::conv::TagPolicy) {
    return io::json::Json::array();
}
template <class State>
inline Maybe_error<void> from_json(const io::json::Json&,
                                   std::vector<Probit_Sample_Record<State>>& out,
                                   const std::string&, io::json::conv::TagPolicy) {
    out.clear();
    return {};
}
}  // namespace macrodr

namespace macrodr::cli {
using std::size_t;
using ModelPtr = macrodr::cmd::ModelPtr;

namespace {

// DSL wrappers for the five diagnostic presets. Each takes max_lag since
// all presets now depend on the kernel-based integral_correlation_lag. F is
// the per-recording numerical Fisher information vector (parallel to dy);
// it must be computed upstream via calculate_numerical_fisher_information.
auto calculate_likelihood_derivative_basic_diagnostics_dsl(
    const std::vector<dMacro_State_Ev_gradient_all>& dlikelihood_predictions,
    const std::vector<parameter_spd_payload>& F_per_recording,
    std::size_t n_boostrap_samples, std::set<double> probits, std::size_t seed,
    std::size_t max_lag) {
    return cmd::calculate_Likelihood_derivative_basic_diagnostics(
        dlikelihood_predictions, F_per_recording, n_boostrap_samples, probits, seed, max_lag);
}

auto calculate_likelihood_derivative_basic_diagnostics_paired_dsl(
    const std::vector<dMacro_State_Ev_gradient_all>& dlikelihood_predictions,
    const std::vector<parameter_spd_payload>& F_per_recording,
    std::size_t n_boostrap_samples, std::set<double> probits, std::size_t seed,
    std::size_t max_lag, std::string samples_dir) {
    return cmd::calculate_Likelihood_derivative_basic_diagnostics_paired(
        dlikelihood_predictions, F_per_recording, n_boostrap_samples, probits, seed, max_lag,
        std::move(samples_dir));
}

auto calculate_likelihood_derivative_series_var_diagnostics_dsl(
    const std::vector<dMacro_State_Ev_gradient_all>& dlikelihood_predictions,
    const std::vector<parameter_spd_payload>& F_per_recording,
    std::size_t n_boostrap_samples, std::set<double> probits, std::size_t seed,
    std::size_t max_lag) {
    return cmd::calculate_Likelihood_derivative_series_var_diagnostics(
        dlikelihood_predictions, F_per_recording, n_boostrap_samples, probits, seed, max_lag);
}

auto calculate_likelihood_derivative_series_cov_diagnostics_dsl(
    const std::vector<dMacro_State_Ev_gradient_all>& dlikelihood_predictions,
    const std::vector<parameter_spd_payload>& F_per_recording,
    std::size_t n_boostrap_samples, std::set<double> probits, std::size_t seed,
    std::size_t max_lag) {
    return cmd::calculate_Likelihood_derivative_series_cov_diagnostics(
        dlikelihood_predictions, F_per_recording, n_boostrap_samples, probits, seed, max_lag);
}

auto calculate_likelihood_derivative_series_kernel_diagnostics_dsl(
    const std::vector<dMacro_State_Ev_gradient_all>& dlikelihood_predictions,
    const std::vector<parameter_spd_payload>& F_per_recording,
    std::size_t n_boostrap_samples, std::set<double> probits, std::size_t seed,
    std::size_t max_lag) {
    return cmd::calculate_Likelihood_derivative_series_kernel_diagnostics(
        dlikelihood_predictions, F_per_recording, n_boostrap_samples, probits, seed, max_lag);
}

auto calculate_likelihood_derivative_series_kernel_full_diagnostics_dsl(
    const std::vector<dMacro_State_Ev_gradient_all>& dlikelihood_predictions,
    const std::vector<parameter_spd_payload>& F_per_recording,
    std::size_t n_boostrap_samples, std::set<double> probits, std::size_t seed,
    std::size_t max_lag) {
    return cmd::calculate_Likelihood_derivative_series_kernel_full_diagnostics(
        dlikelihood_predictions, F_per_recording, n_boostrap_samples, probits, seed, max_lag);
}

// Build a gauss_newton_options struct by value (mirrors build_likelihood_function).
// Exposes the 7 tunable fields so a DSL script can configure the per-group MLE
// inner loop without recompiling.
auto build_gauss_newton_options_dsl(double lambda_kickoff, double lambda_factor,
                                     double lambda_max, std::size_t max_iter, double grad_rtol,
                                     double dvalue_tol, bool verbose)
    -> macrodr::optimization::gauss_newton_options {
    macrodr::optimization::gauss_newton_options o;
    o.lambda_kickoff = lambda_kickoff;
    o.lambda_factor = lambda_factor;
    o.lambda_max = lambda_max;
    o.max_iter = max_iter;
    o.grad_rtol = grad_rtol;
    o.dvalue_tol = dvalue_tol;
    o.verbose = verbose;
    return o;
}

// Extract the plain Recording from each loaded Simulated_recording. The
// calc_MLE_per_group_of_replicates input is std::vector<Recording>; load_simulations
// yields std::vector<Simulated_Recording<...>>, whose Vector_Space carries a
// Recording slot (get<Recording>(sim()), same access as calculate_simulation_*).
auto get_recordings_from_simulations_dsl(
    const std::vector<Simulated_Recording<var::please_include<>>>& simulations)
    -> std::vector<Recording> {
    std::vector<Recording> out;
    out.reserve(simulations.size());
    for (auto const& s : simulations) out.push_back(get<Recording>(s()));
    return out;
}

// Driver wrapper pinning State = dMacro_State_Hessian_minimal_param (Path A).
// F_h_relative is bound to its header default (1e-5) here rather than exposed.
auto calc_MLE_per_group_of_replicates_dsl(
    const cmd::likelihood_algorithm_type& likelihood_algorithm,
    const var::Parameters_transformed& theta_warmstart,
    const var::Parameters_transformed& theta_reference, const Experiment& experiment,
    const std::vector<Recording>& recordings, std::size_t group_size,
    std::size_t n_bootstrap_samples, std::size_t min_groups_for_bootstrap,
    std::set<double> probit_cis, std::set<double> probit_sample_heights,
    std::vector<std::string> ranking_variables, std::size_t seed,
    const macrodr::optimization::gauss_newton_options& gn_opts)
    -> Maybe_error<cmd::MLE_Group_Cloud<dMacro_State_Hessian_minimal_param>> {
    return cmd::calc_MLE_per_group_of_replicates<dMacro_State_Hessian_minimal_param>(
        likelihood_algorithm, theta_warmstart, theta_reference, experiment, recordings, group_size,
        n_bootstrap_samples, min_groups_for_bootstrap, probit_cis, probit_sample_heights,
        ranking_variables, seed, gn_opts, 1e-5);
}

// θ̄ handoff: rebuild the full-sample mean of θ̂ as a Parameters_transformed from
// the cloud's BARE Model_Parameters_Hat slot. State pinned to Path A.
auto get_parameters_mean_dsl(
    const cmd::MLE_Group_Cloud<dMacro_State_Hessian_minimal_param>& cloud)
    -> Maybe_error<var::Parameters_transformed> {
    return cmd::get_parameters_mean<dMacro_State_Hessian_minimal_param>(cloud);
}

// figure_3 Fase-2 empirical-vs-theoretical capstone. State pinned to Path A;
// rtol / atol bound to their header defaults (1e-10 / 0.0).
auto calc_empirical_distortion_dsl(
    const cmd::MLE_Group_Cloud<dMacro_State_Hessian_minimal_param>& cloud,
    const std::vector<parameter_spd_payload>& fim_sim,
    const std::vector<parameter_spd_payload>& fim_bar,
    const std::vector<dMacro_State_Ev_gradient_all>& dlik_bar,
    std::size_t n_bootstrap, std::size_t seed, std::set<double> probit_cis)
    -> Maybe_error<cmd::Empirical_Distortion_Bootstrap> {
    return cmd::calc_empirical_distortion<dMacro_State_Hessian_minimal_param>(
        cloud, fim_sim, fim_bar, dlik_bar, n_bootstrap, seed, probit_cis, 1e-10, 0.0);
}

}

inline macrodr::dsl::Compiler<macrodr::dsl::Lexer> make_simulations_compiler() {
    dsl::Compiler<dsl::Lexer> cm;
    cm.push_function("get_num_parameters",
                     dsl::to_typed_function<std::string>(&get_num_parameters, "model"));

    cm.push_function("simulate",
                     dsl::to_typed_function<std::string, cmd::recording_type, cmd::experiment_type,
                                            std::size_t, const std::string&,
                                            cmd::parameters_value_type, cmd::simulation_algo_type>(
                         &cmd::runsimulation, "filename_prefix", "recording", "experiment",
                         "init_seed", "modelName", "parameter_values", "simulation_algorithm"));

    {
        using SimFromValues = Maybe_error<Simulated_Recording<var::please_include<>>> (*)(
            const ModelPtr&, const var::Parameters_values&,
            const Experiment&, const Recording&, std::size_t, std::size_t);
        cm.push_function(
            "simulate",
            dsl::to_typed_function<const ModelPtr&,
                                   const var::Parameters_values&, const Experiment&,
                                   const Recording&, std::size_t, std::size_t>(
                static_cast<SimFromValues>(&cmd::run_simulations), "model", "parameter_values",
                "experiment", "observations", "number_of_substeps", "seed"));
    }

    {
        using SimUniformFromValues = Maybe_error<Simulated_Recording<var::please_include<>>> (*)(
            const ModelPtr&, const var::Parameters_values&, const Experiment&, const Recording&,
            std::string, std::size_t, std::size_t);
        cm.push_function(
            "simulate",
            dsl::to_typed_function<const ModelPtr&, const var::Parameters_values&,
                                   const Experiment&, const Recording&, std::string,std::size_t, std::size_t>
                                   (
                static_cast<SimUniformFromValues>(&cmd::run_simulations), "model",
                "parameter_values", "experiment", "observations", "simulation_algorithm","number_of_substeps",
                "seed"));
    }

    {
        using SimFromTransformed = Maybe_error<Simulated_Recording<var::please_include<>>> (*)(
            const ModelPtr&, const var::Parameters_transformed&,
            const Experiment&, const Recording&, std::size_t, std::size_t);
        cm.push_function(
            "simulate",
            dsl::to_typed_function<const ModelPtr&,
                                   const var::Parameters_transformed&, const Experiment&,
                                   const Recording&, std::size_t, std::size_t>(
                static_cast<SimFromTransformed>(&cmd::run_simulations), "model",
                "parameter_values", "experiment", "observations", "number_of_substeps", "seed"));
    }

    {
        using SimUniformFromTransformed = Maybe_error<Simulated_Recording<var::please_include<>>> (*)(
            const ModelPtr&, const var::Parameters_transformed&, const Experiment&,
            const Recording&, std::string, std::size_t, std::size_t);
        cm.push_function(
            "simulate",
            dsl::to_typed_function<const ModelPtr&, const var::Parameters_transformed&,
                                   const Experiment&, const Recording&, std::string,std::size_t, std::size_t>(
                static_cast<SimUniformFromTransformed>(&cmd::run_simulations), "model",
                "parameter_transformed", "experiment", "observations",
                "simulation_algorithm", "number_of_substeps", "seed"));
    }
    cm.push_function(
        "simulate",
        dsl::to_typed_return_function<
            Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>>,
            const ModelPtr&,
            const var::Parameters_transformed&, std::size_t,const Experiment&,
            const Recording&, std::size_t, std::size_t>(
                static_cast<Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> (*)(
                    const ModelPtr&, const var::Parameters_transformed&, std::size_t,
                    const Experiment&, const Recording&, std::size_t, std::size_t)>(
                    &cmd::run_n_simulations),
                "model",
                "parameter_transformed", "n_simulations", "experiment", "observations", "number_of_substeps", "seed"));

    cm.push_function(
        "simulate",
        dsl::to_typed_return_function<
            Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>>,
            const ModelPtr&, const var::Parameters_transformed&, std::size_t,
            const Experiment&, const Recording&, std::string, std::size_t, std::size_t>(
            static_cast<Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> (*)(
                const ModelPtr&, const var::Parameters_transformed&, std::size_t,
                const Experiment&, const Recording&, std::string, std::size_t, std::size_t)>(
                &cmd::run_n_simulations),
            "model", "parameter_transformed", "n_simulations", "experiment", "observations",
            "simulation_algorithm", "number_of_substeps", "seed"));
    
    cm.push_function(
        "simulate",
        dsl::to_typed_return_function<
            Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>>,
            const ModelPtr&,
            const var::Parameters_values&, std::size_t,const Experiment&,
            const Recording&, std::size_t, std::size_t>(
                static_cast<Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> (*)(
                    const ModelPtr&, const var::Parameters_values&, std::size_t,
                    const Experiment&, const Recording&, std::size_t, std::size_t)>(
                    &cmd::run_n_simulations),
                "model",
                "parameter_values", "n_simulations", "experiment", "observations", "number_of_substeps", "seed"));

    cm.push_function(
        "simulate",
        dsl::to_typed_return_function<
            Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>>,
            const ModelPtr&, const var::Parameters_values&, std::size_t,
            const Experiment&, const Recording&, std::string, std::size_t, std::size_t>(
            static_cast<Maybe_error<std::vector<Simulated_Recording<var::please_include<>>>> (*)(
                const ModelPtr&, const var::Parameters_values&, std::size_t,
                const Experiment&, const Recording&, std::string, std::size_t, std::size_t)>(
                &cmd::run_n_simulations),
            "model", "parameter_values", "n_simulations", "experiment", "observations",
            "simulation_algorithm", "number_of_substeps", "seed"));
    

    {
        using SimSubFromValues = Maybe_error<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> (*)(
            const ModelPtr&, const var::Parameters_values&,
            const Experiment&, const Recording&, std::size_t, std::size_t);
        cm.push_function(
            "simulate_with_sub_intervals",
            dsl::to_typed_function<const ModelPtr&,
                                   const var::Parameters_values&, const Experiment&,
                                   const Recording&, std::size_t, std::size_t>(
                static_cast<SimSubFromValues>(&cmd::run_simulations_with_sub_intervals), "model",
                "parameter_values", "experiment", "observations", "number_of_substeps", "seed"));
    }

    {
        using SimSubUniformFromValues = Maybe_error<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> (*)(
            const ModelPtr&, const var::Parameters_values&, const Experiment&, const Recording&,
            std::string, std::size_t);
        cm.push_function(
            "simulate_with_sub_intervals",
            dsl::to_typed_function<const ModelPtr&, const var::Parameters_values&,
                                   const Experiment&, const Recording&, std::string,
                                   std::size_t>(
                static_cast<SimSubUniformFromValues>(&cmd::run_simulations_with_sub_intervals),
                "model", "parameter_values", "experiment", "observations",
                "simulation_algorithm", "seed"));
    }

    {
        using SimSubFromTransformed = Maybe_error<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> (*)(
            const ModelPtr&, const var::Parameters_transformed&,
            const Experiment&, const Recording&, std::size_t, std::size_t);
        cm.push_function(
            "simulate_with_sub_intervals",
            dsl::to_typed_function<const ModelPtr&,
                                   const var::Parameters_transformed&, const Experiment&,
                                   const Recording&, std::size_t, std::size_t>(
                static_cast<SimSubFromTransformed>(&cmd::run_simulations_with_sub_intervals),
                "model", "parameters", "experiment", "observations", "number_of_substeps",
                "seed"));
    }

    {
        using SimSubUniformFromTransformed = Maybe_error<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>> (*)(
            const ModelPtr&, const var::Parameters_transformed&, const Experiment&,
            const Recording&, std::string, std::size_t);
        cm.push_function(
            "simulate_with_sub_intervals",
            dsl::to_typed_function<const ModelPtr&, const var::Parameters_transformed&,
                                   const Experiment&, const Recording&, std::string,
                                   std::size_t>(
                static_cast<SimSubUniformFromTransformed>(
                    &cmd::run_simulations_with_sub_intervals),
                "model", "parameter_transformed", "experiment", "observations",
                "simulation_algorithm", "seed"));
    }

    cm.push_function(
        "remove_intervals",
        dsl::to_typed_function<
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const&>(
            &cmd::remove_intervals, "simulation_with_sub_intervals"));  
  
    return cm;
}

dsl::Compiler<dsl::Lexer> make_compiler_new() {
    dsl::Compiler<dsl::Lexer> cm;
    // Meta commands (help/version) available everywhere
    cm.merge(macrodr::cmd::make_cli_meta_compiler());
    cm.merge(macrodr::cmd::make_utilities_compiler());
    cm.merge(macrodr::cmd::make_io_compiler());
    cm.merge(macrodr::cmd::make_experiment_compiler());
    cm.merge(macrodr::cmd::make_model_compiler());
    cm.merge(make_simulations_compiler());
    cm.merge(macrodr::cmd::make_likelihood_compiler());
    cm.merge(macrodr::cmd::make_dts_compiler());
    cm.push_function("load_model",
                     dsl::to_typed_function<std::string>(&cmd::load_model, "model_name"));
    cm.push_function("load_experiment", dsl::to_typed_function<std::string, double, double>(
                                            &macrodr::cmd::load_experiment, "filename",
                                            "frequency_of_sampling", "initial_agonist"));

    cm.push_function(
        "create_experiment",
        dsl::to_typed_function<std::vector<std::tuple<std::size_t, std::size_t, double>>, double,
                               double, double>(&macrodr::cmd::create_experiment,
                                               "experiment_structure", "frequency_of_sampling",
                                               "initial_agonist", "initial_time"));

    cm.push_function("load_observations", dsl::to_typed_function<std::string>(
                                              &macrodr::cmd::load_recording, "filename"));

    // Load simulations from a CSV (the _simulation.csv written by
    // write_csv(experiment, simulations, path)). Returns var::Indexed (native
    // Indexed return, like indexed_by) reconstructed from the CSV axis columns,
    // so the axis chain survives into downstream Indexed pipelines. Two overloads:
    //   load_simulations(filename)                  -> all simulations
    //   load_simulations(filename, replica_indices) -> only listed indices
    // Args by-value so the DSL accepts literal arguments inline.
    using LoadSimAll =
        Maybe_error<var::Indexed<std::vector<Simulated_Recording<var::please_include<>>>>> (*)(
            std::string);
    using LoadSimFiltered =
        Maybe_error<var::Indexed<std::vector<Simulated_Recording<var::please_include<>>>>> (*)(
            std::string, std::vector<std::size_t>);
    cm.push_function(
        "load_simulations",
        dsl::to_typed_function<std::string>(
            static_cast<LoadSimAll>(&macrodr::cmd::load_simulations), "filename"));
    cm.push_function(
        "load_simulations",
        dsl::to_typed_function<std::string, std::vector<std::size_t>>(
            static_cast<LoadSimFiltered>(&macrodr::cmd::load_simulations),
            "filename", "replica_indices"));

    cm.push_function("set_observations", dsl::to_typed_function<std::vector<double>>(
                                              &macrodr::cmd::define_recording, "values_set"));

    // Build a Recording with explicit missing-measurement flags. Indices in
    // `missing_indices` are set to NaN; the likelihood path treats NaN samples
    // as "predict, don't update" — the posterior propagates forward without a
    // Bayes update and the log-likelihood contribution at that step is zero.
    cm.push_function("set_observations_with_gaps",
                     dsl::to_typed_function<std::size_t, double, std::vector<std::size_t>>(
                         &macrodr::cmd::define_recording_with_gaps, "n_samples", "fill_value",
                         "missing_indices"));

    // Convenience for a contiguous gap (e.g. a long pre-equilibration lead-in):
    // marks [missing_start, missing_end) as NaN.
    cm.push_function(
        "set_observations_with_missing_range",
        dsl::to_typed_function<std::size_t, double, std::size_t, std::size_t>(
            &macrodr::cmd::define_recording_with_missing_range, "n_samples", "fill_value",
            "missing_start", "missing_end"));

    cm.push_function("axis",
                     dsl::to_typed_function<std::string, std::vector<std::string>>(
                         &cmd::axis, "name", "labels"));
    cm.push_function("indexed_bool_by",
                     dsl::to_typed_function<var::Axis, std::vector<bool>>(
                         &cmd::indexed_by<bool>, "axis", "values"));
    cm.push_function("indexed_int_by",
                     dsl::to_typed_function<var::Axis, std::vector<int>>(
                         &cmd::indexed_by<int>, "axis", "values"));
    cm.push_function("indexed_size_by",
                     dsl::to_typed_function<var::Axis, std::vector<std::size_t>>(
                         &cmd::indexed_by<std::size_t>, "axis", "values"));
    cm.push_function("indexed_double_by",
                     dsl::to_typed_function<var::Axis, std::vector<double>>(
                         &cmd::indexed_by<double>, "axis", "values"));
    cm.push_function("indexed_string_by",
                     dsl::to_typed_function<var::Axis, std::vector<std::string>>(
                         &cmd::indexed_by<std::string>, "axis", "values"));
    cm.push_function("indexed_by",
                     dsl::to_typed_function<var::Axis, std::vector<std::size_t>>(
                         &cmd::indexed_by<std::size_t>, "axis", "values"));
    cm.push_function("indexed_by",
                     dsl::to_typed_function<var::Axis, std::vector<int>>(
                         &cmd::indexed_by<int>, "axis", "values"));
    cm.push_function("indexed_by",
                     dsl::to_typed_function<var::Axis, std::vector<double>>(
                         &cmd::indexed_by<double>, "axis", "values"));
    cm.push_function("indexed_by",
                     dsl::to_typed_function<var::Axis, std::vector<std::string>>(
                         &cmd::indexed_by<std::string>, "axis", "values"));

    using SimulationBatch = std::vector<Simulated_recording>;
    using IndexedSimulation = var::Indexed<Simulated_recording>;
    using IndexedSimulationBatch = var::Indexed<SimulationBatch>;
    using SubSimulation =
        Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>;
    using IndexedSubSimulation = var::Indexed<SubSimulation>;
    using GradientBatch = std::vector<dMacro_State_Ev_gradient_all>;
    using IndexedPredictions = var::Indexed<Macro_State_Ev_predictions>;
    using IndexedDiagnostics = var::Indexed<Macro_State_Ev_diagnostic>;
    using IndexedGradients = var::Indexed<dMacro_State_Ev_gradient_all>;
    using IndexedGradientBatch = var::Indexed<GradientBatch>;
    // Per-sample numerical Fisher Information State (diagnostic): same shape
    // family as GradientBatch but the Evolution slot carries per-step F_t
    // contributions instead of the gradient_all element. Used by the bug-hunt
    // pipeline (calc_per_sample_numerical_fisher_information).
    using PerSampleFBatch = std::vector<dMacro_State_Ev_per_sample_F>;
    using IndexedPerSampleFBatch = var::Indexed<PerSampleFBatch>;
    // Detailed per-sample diagnostic State: Evolution slot carries the rich
    // detailed_element (P_mean/P_Cov/y_mean/y_var/trust_coefficient/logL as
    // Derivatives). Used by calc_per_sample_numerical_fisher_information_detailed.
    using DetailedBatch = std::vector<dMacro_State_Ev_detailed>;
    using IndexedDetailedBatch = var::Indexed<DetailedBatch>;

    cm.push_function("write_csv", dsl::to_typed_return_function<
        Maybe_error<std::string>,Experiment const&,Simulated_recording const&, std::string >(&macrodr::cmd::write_csv,
            "experiment","simulation", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "path"));
    cm.push_function("write_csv", dsl::to_typed_return_function<
        Maybe_error<std::string>,Experiment const&,std::vector<Simulated_recording> const&, std::string >(&macrodr::cmd::write_csv,
            "experiment","simulations", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulationBatch const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulations", "path"));
    cm.push_function("write_csv", dsl::to_typed_return_function<
        Maybe_error<std::string>,Experiment const&,Recording const&, std::string >(&macrodr::cmd::write_csv,
            "experiment","observations", "path"));

    cm.push_function("write_csv", 
        dsl::to_typed_return_function<
        Maybe_error<std::string>,
        Experiment const&,
        Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const&, 
        std::size_t,
        std::string >(
            &macrodr::cmd::write_csv,
            "experiment",
            "simulation", 
            "number_of_substates",
            "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSubSimulation const&, std::size_t, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "number_of_substates",
            "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_Recording<var::please_include<>> const&,
                                      Macro_State_Ev_predictions const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&,
                                      Macro_State_Ev_predictions const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_recording const&, IndexedPredictions const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&, IndexedPredictions const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_Recording<var::please_include<>> const&,
                                      Macro_State_Ev_diagnostic const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&,
                                      Macro_State_Ev_diagnostic const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_recording const&, IndexedDiagnostics const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&, IndexedDiagnostics const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_Recording<var::please_include<>> const&,
                                      dMacro_State_Ev_gradient_all const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&,
                                      dMacro_State_Ev_gradient_all const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      Simulated_recording const&, IndexedGradients const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulation const&, IndexedGradients const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      SimulationBatch const&, GradientBatch const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulations", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulationBatch const&, GradientBatch const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulations", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      SimulationBatch const&, IndexedGradientBatch const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulations", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulationBatch const&, IndexedGradientBatch const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulations", "likelihood", "path"));

    // write_csv for the per-sample numerical Fisher Information State. Same
    // 4-variant set as GradientBatch above (plain vs Indexed × simulations vs
    // states) — needed so DSL dispatch picks the right overload regardless of
    // how the inputs got wrapped through the axis-combo loop.
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      SimulationBatch const&, PerSampleFBatch const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulations", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulationBatch const&, PerSampleFBatch const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulations", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      SimulationBatch const&, IndexedPerSampleFBatch const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulations", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulationBatch const&, IndexedPerSampleFBatch const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulations", "likelihood", "path"));
    // Dedicated 2-arg overload (Indexed states + path, no Experiment/Simulation):
    // the DSL matches the Indexed param directly without lifting, so the axis
    // columns survive into the per_sample_F CSV.
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>,
                                      IndexedPerSampleFBatch const&, std::string>(
            &macrodr::cmd::write_csv, "likelihood", "path"));

    // DETAILED diagnostic write_csv: windowed (sample_min/sample_max) variant is
    // the primary one — both simulations and likelihood arrive Indexed over the
    // scenario/data axes. The window keeps the CSV small while preserving
    // absolute sample indices / times / agonist / patch_current.
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulationBatch const&, IndexedDetailedBatch const&,
                                      std::size_t, std::size_t, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulations", "likelihood", "sample_min",
            "sample_max", "path"));
    // Non-windowed sibling (full Evolution) for convenience.
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSimulationBatch const&, IndexedDetailedBatch const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulations", "likelihood", "path"));

    // Per-replica state + numerical Fisher dump (no Experiment, no Simulation).
    // Emits each state without its Evolution_of slot and the matched F as
    // Likelihood_Numerical_Fisher_Information. Designed for bug-hunting on
    // numeric_Fisher_information variability across replicates.
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>,
                                      GradientBatch const&,
                                      std::vector<parameter_spd_payload> const&,
                                      std::string>(
            &macrodr::cmd::write_csv,
            "dlikelihood_predictions", "numerical_fisher_information", "path"));

    // Indexed (multi-axis) variant: when called from a DSL context with
    // axes in scope (algorithm × interval × Num_ch × h_rel × ...), the DSL
    // hands these as Indexed values. This overload iterates the axis space
    // and emits the CSV with axis_values columns populated per coord (matching
    // the analysis CSV schema).
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>,
                                      var::Indexed<GradientBatch> const&,
                                      var::Indexed<std::vector<parameter_spd_payload>> const&,
                                      std::string>(
            &macrodr::cmd::write_csv,
            "dlikelihood_predictions", "numerical_fisher_information", "path"));


    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>, Experiment const&,
            Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>> const&,
            dMacro_State_Ev_gradient_all const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSubSimulation const&,
                                      dMacro_State_Ev_gradient_all const&, std::string>(
            &macrodr::cmd::write_csv, "experiment", "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      SubSimulation const&, IndexedGradients const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, Experiment const&,
                                      IndexedSubSimulation const&, IndexedGradients const&,
                                      std::string>(&macrodr::cmd::write_csv, "experiment",
                                                   "simulation", "likelihood", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>, macrodr::cmd::Analisis_derivative_diagnostic_base const&,
            std::string>(&macrodr::cmd::write_csv, "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>, macrodr::cmd::Analisis_derivative_diagnostic_basic const&,
            std::string>(&macrodr::cmd::write_csv, "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>, macrodr::cmd::Analisis_derivative_diagnostic_series_var const&,
            std::string>(&macrodr::cmd::write_csv, "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>, macrodr::cmd::Analisis_derivative_diagnostic_series_cov const&,
            std::string>(&macrodr::cmd::write_csv, "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>,
            macrodr::cmd::Analisis_derivative_diagnostic_series_kernel const&,
            std::string>(&macrodr::cmd::write_csv, "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>,
            macrodr::cmd::Analisis_derivative_diagnostic_series_kernel_full const&,
            std::string>(&macrodr::cmd::write_csv, "analysis", "path"));

    // Indexed variants for each preset type. The DSL needs a concrete function
    // pointer at each type — the generic write_csv(Indexed<T>, …) template at
    // include/macrodr/cmd/likelihood.h:314 supplies the instantiation.
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>,
                                      const var::Indexed<
                                          macrodr::cmd::Analisis_derivative_diagnostic_base>&,
                                      std::string>(
            &macrodr::cmd::write_csv<macrodr::cmd::Analisis_derivative_diagnostic_base>,
            "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>,
                                      const var::Indexed<
                                          macrodr::cmd::Analisis_derivative_diagnostic_basic>&,
                                      std::string>(
            &macrodr::cmd::write_csv<macrodr::cmd::Analisis_derivative_diagnostic_basic>,
            "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>,
                                      const var::Indexed<
                                          macrodr::cmd::Analisis_derivative_diagnostic_series_var>&,
                                      std::string>(
            &macrodr::cmd::write_csv<macrodr::cmd::Analisis_derivative_diagnostic_series_var>,
            "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>,
                                      const var::Indexed<
                                          macrodr::cmd::Analisis_derivative_diagnostic_series_cov>&,
                                      std::string>(
            &macrodr::cmd::write_csv<macrodr::cmd::Analisis_derivative_diagnostic_series_cov>,
            "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>,
            const var::Indexed<macrodr::cmd::Analisis_derivative_diagnostic_series_kernel>&,
            std::string>(
            &macrodr::cmd::write_csv<macrodr::cmd::Analisis_derivative_diagnostic_series_kernel>,
            "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<
            Maybe_error<std::string>,
            const var::Indexed<macrodr::cmd::Analisis_derivative_diagnostic_series_kernel_full>&,
            std::string>(&macrodr::cmd::write_csv<
                             macrodr::cmd::Analisis_derivative_diagnostic_series_kernel_full>,
                         "analysis", "path"));

    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, const var::Indexed<std::size_t>&,
                                      std::string>(&macrodr::cmd::write_csv<std::size_t>,
                                                   "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, const var::Indexed<int>&,
                                      std::string>(&macrodr::cmd::write_csv<int>, "analysis",
                                                   "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, const var::Indexed<double>&,
                                      std::string>(&macrodr::cmd::write_csv<double>, "analysis",
                                                   "path"));

            

                                              // Expose low-level model/qmodel helpers
    // Load parameters file (filename, sep) using model’s schema
    cm.push_function("load_parameters",
                     dsl::to_typed_function<ModelPtr, std::string>(&macrodr::cmd::load_parameters,
                                                                   "model", "parameter_file"));

    cm.push_function("create_parameters",
                     dsl::to_typed_function<ModelPtr,  std::vector<std::tuple<std::string, std::string, double>>>
                     (&macrodr::cmd::create_parameters,"model", "parameters_info"));

    // In-line Gaussian prior construction (counterpart to load_prior); each
    // tuple is (name, transformation, value, transformed_variance).
    cm.push_function("create_prior",
                     dsl::to_typed_function<ModelPtr, std::vector<std::tuple<std::string, std::string, double, double>>>
                     (&macrodr::cmd::create_prior, "model", "prior_info"));


                                                                   // load_parameter_values returns transformations (safer to persist)
    cm.push_function("get_standard_parameter_values",
                     dsl::to_typed_function<var::Parameters_Transformations const&>(
                         &macrodr::cmd::get_standard_parameter_values, "parameters"));

    cm.push_function("get_standard_parameter_transformed_values",
                     dsl::to_typed_function<var::Parameters_Transformations const&>(
                         &macrodr::cmd::get_standard_parameter_transformed_values, "parameters"));

    // Two-step perturbation for the per-sample FD diagnostic, built outside the
    // (agnostic) command. Step 1: by_parameter_coordinate → Indexed raw Δ (h·eᵢ)
    // over `parameter_coordinate` × `h_rel` (the `by_` prefix = returns Indexed
    // keyed by a derived axis; plain args + Indexed return → native indexed, no
    // lift). Step 2: apply_relative_perturbation applies θ′ᵢ = θᵢ + Δᵢ·max(|θᵢ|,1)
    // and returns PLAIN, so the DSL lifts it over Δ's axes and combines them.
    // h_rels is a plain vector ([+h,−h], or several scales).
    cm.push_function(
        "by_parameter_coordinate",
        dsl::to_typed_function<var::Parameters_transformed const&, std::vector<double>>(
            &macrodr::cmd::by_parameter_coordinate, "parameters", "h_rels"));
    cm.push_function("apply_relative_perturbation",
                     dsl::to_typed_function<var::Parameters_transformed const&,
                                            var::Parameters_transformed const&>(
                         &macrodr::cmd::apply_relative_perturbation, "parameters", "perturbation"));


    cm.push_function(
        "patch_model",
        dsl::to_typed_return_function<Maybe_error<macrodr::cmd::PatchModel>, ModelPtr const&,
                                      std::pair<std::string, std::string> const&>(
            &macrodr::cmd::patch_model, "model", "parameter_values"));
    // Overload: patch from in-memory values
    {
        using Return = Maybe_error<macrodr::cmd::PatchModel>;
        using FnPatchFromValues = Return (*)(const ModelPtr&, const var::Parameters_values&);
        cm.push_function(
            "patch_model",
            dsl::to_typed_function<ModelPtr, var::Parameters_values>(
                static_cast<FnPatchFromValues>(&macrodr::cmd::patch_model), "model", "values"));
    }
    {
        using Return = Maybe_error<macrodr::cmd::PatchModel>;
        using FnPatchFromTr = Return (*)(const ModelPtr&, const var::Parameters_Transformations&);
        cm.push_function("patch_model",
                         dsl::to_typed_function<ModelPtr, var::Parameters_Transformations>(
                             static_cast<FnPatchFromTr>(&macrodr::cmd::patch_model), "model",
                             "parameter_values"));
    }
    // Back-compat: accept 'parameter_values' as the argument name for in-memory values as well
    {
        using Return = Maybe_error<macrodr::cmd::PatchModel>;
        using FnPatchFromValues = Return (*)(const ModelPtr&, const var::Parameters_values&);
        cm.push_function("patch_model",
                         dsl::to_typed_function<ModelPtr, var::Parameters_values>(
                             static_cast<FnPatchFromValues>(&macrodr::cmd::patch_model), "model",
                             "parameter_values"));
    }

    using PatchModel = macrodr::Transfer_Op_to<var::Parameters_values, macrodr::Patch_Model>;
    cm.push_function("calc_eigen", dsl::to_typed_function<PatchModel const&, double>(
                                       &macrodr::cmd::calc_eigen, "patch_model", "agonist"));

    using QxEig = macrodr::Transfer_Op_to<PatchModel, macrodr::Eigs>;
    cm.push_function("calc_peq", dsl::to_typed_function<QxEig const&, PatchModel const&>(
                                     &macrodr::cmd::calc_peq, "qx_eig", "patch_model"));

    // Patch state initial condition helpers
    {
        using RetPS = Maybe_error<macrodr::Patch_State>;
        using FnPSfromPM = RetPS (*)(const macrodr::cmd::PatchModel&, double);
        cm.push_function("path_state", dsl::to_typed_function<PatchModel const&, double>(
                                           static_cast<FnPSfromPM>(&macrodr::cmd::path_state),
                                           "patch_model", "initial_agonist"));
    }
    {
        using RetPS = Maybe_error<macrodr::Patch_State>;
        using FnPSfromVals = RetPS (*)(const ModelPtr&, const var::Parameters_values&, double);
        cm.push_function("path_state",
                         dsl::to_typed_function<ModelPtr, var::Parameters_values const&, double>(
                             static_cast<FnPSfromVals>(&macrodr::cmd::path_state), "model",
                             "values", "initial_agonist"));
    }

    cm.push_function("build_likelihood_function",
                     dsl::to_typed_function<
                         ModelPtr const&, bool, bool, int, bool, bool, bool, bool>(
                         &macrodr::cmd::build_likelihood_function, "model",
                                "adaptive_approximation", "recursive_approximation",
                                "averaging_approximation", "variance_approximation",
                                "taylor_variance_correction", "micro_approximation",
                                "taylor_qdt_approximation")
                            );

    cm.push_function(
        "calc_likelihood",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&>(
            &cmd::calculate_mlikelihood, "likelihood_algorithm", "parameters", "experiment", "data"));

    cm.push_function(
        "calc_likelihood",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&>(
            &cmd::calculate_simulation_mlikelihood, "likelihood_algorithm", "parameters",
            "experiment", "data"));

    cm.push_function(
        "calc_dlikelihood",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&>(
            &cmd::calculate_mdlikelihood, "likelihood_algorithm", "parameters", "experiment",
            "data"));

    cm.push_function(
        "calc_diff_likelihood",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&, const var::Parameters_transformed&,
                               const Experiment&, const Recording&, double>(
            &cmd::calculate_mdiff_likelihood, "likelihood_algorithm", "parameters", "experiment",
            "data", "delta_param"));

    cm.push_function(
        "calc_likelihood_predictions",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&>(
            &cmd::calculate_mlikelihood_predictions, "likelihood_algorithm", "parameters",
            "experiment", "data"));

        cm.push_function(
        "calc_likelihood_predictions",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Simulated_Recording<var::please_include<>>&>(
            &cmd::calculate_simulation_mlikelihood_predictions, "likelihood_algorithm", "parameters",
            "experiment", "data")); 


    cm.push_function(
        "calc_likelihood_diagnostic",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&>(
            &cmd::calculate_mlikelihood_diagnostics, "likelihood_algorithm", "parameters",
            "experiment", "data"));
    cm.push_function(
        "calc_likelihood_diagnostic",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Simulated_Recording<var::please_include<>>&&>(
            &cmd::calculate_simulated_mlikelihood_diagnostics, "likelihood_algorithm", "parameters",
            "experiment", "simulation"));


    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&>(
            &cmd::calculate_mdlikelihood_predictions, "likelihood_algorithm", "parameters",
            "experiment", "data"));

    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&>(
            &cmd::calculate_simulation_mdlikelihood_predictions, "likelihood_algorithm",
            "parameters", "experiment", "data"));

    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&,
                               const std::vector<Simulated_Recording<var::please_include<> >>&>(
            &cmd::calculate_n_simulation_mdlikelihood_predictions, "likelihood_algorithm",
            "parameters", "experiment", "data_series"));


            cm.push_function(
        "calc_likelihood",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool, bool>(
            &cmd::calculate_likelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation"));

    cm.push_function(
        "calc_likelihood",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool, bool>(
            &cmd::calculate_simulation_likelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation"));

    cm.push_function(
        "calc_dlikelihood",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool, bool>(
            &cmd::calculate_dlikelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation"));

    cm.push_function(
        "calc_dlikelihood",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool, bool>(
            &cmd::calculate_simulation_dlikelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation"));

    cm.push_function(
        "calc_diff_likelihood",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool, bool, double>(
            &cmd::calculate_diff_likelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation",
            "delta_param"));

    cm.push_function(
        "calc_diff_likelihood",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool, bool, double>(
            &cmd::calculate_simulation_diff_likelihood, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation",
            "delta_param"));

    cm.push_function(
        "calc_likelihood_predictions",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool, bool>(
            &cmd::calculate_likelihood_predictions, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation"));

    cm.push_function(
        "calc_likelihood_predictions",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool, bool>(
            &cmd::calculate_simulation_likelihood_predictions, "model", "parameters", "experiment",
            "data", "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation"));
    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool, bool>(
            &cmd::calculate_dlikelihood_predictions, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation"));

    // Numerical (observed) Fisher information per recording. 2·n_params extra
    // dlikelihood evaluations at θ ± h_rel·max(|θ_i|, 1) per coordinate.
    // Variant via likelihood_algorithm_type (matches calc_dlikelihood_predictions
    // script-level usage). Three overloads: Recording, Simulated_Recording,
    // vector<Simulated_Recording>.
    cm.push_function(
        "calc_numerical_fisher_information",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Recording&, double>(
            &cmd::calculate_mnumerical_fisher_information, "likelihood_algorithm", "parameters",
            "experiment", "data", "h_rel"));

    cm.push_function(
        "calc_numerical_fisher_information",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, double>(
            &cmd::calculate_simulation_mnumerical_fisher_information, "likelihood_algorithm",
            "parameters", "experiment", "data", "h_rel"));

    cm.push_function(
        "calc_numerical_fisher_information",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&,
                               const std::vector<Simulated_Recording<var::please_include<>>>&,
                               double, std::size_t>(
            &cmd::calculate_n_simulation_mnumerical_fisher_information, "likelihood_algorithm",
            "parameters", "experiment", "data_series", "h_rel", "decimate"));

    // Per-sample variant: decomposes the global F per replica into per-timestep
    // F_t contributions via per-step FD on the cumulative AD score. Returns a
    // vector<dMacro_State_Ev_per_sample_F> whose Evolution_of<> slot carries
    // one Likelihood_Numerical_Fisher_Information per timestep. By additivity
    // of logL = Σ_t logL_t, summing F_t over the Evolution recovers the global
    // F. Used for diagnostic localization of FD instabilities — find which
    // timestep triggers the kink in score derivatives for specific parameter
    // directions (the bug hunt on Current_Noise / Current_Baseline directions).
    cm.push_function(
        "calc_per_sample_numerical_fisher_information",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&,
                               const std::vector<Simulated_Recording<var::please_include<>>>&,
                               double>(
            &cmd::calculate_per_sample_n_simulation_mnumerical_fisher_information,
            "likelihood_algorithm", "parameters", "experiment", "data_series", "h_rel"));

    // Detailed variant: runs the dlikelihood ONCE at θ and dumps the significant
    // recursion variables (logL, y_mean, y_var, P_mean, P_Cov, trust_coefficient)
    // per sample, each as a Derivative carrying the value X AND the regular AD
    // gradient ∂X/∂θ over ALL parameters (the Jacobian already covers every
    // parameter). The FD-instability perturbation is built OUTSIDE: pass
    // `parameters` as an Indexed over the perturbation axes (see
    // by_parameter_coordinate + apply_relative_perturbation) and the DSL lifts
    // this command over them.
    // sample_min/sample_max window each Evolution INSIDE (memory) — absolute
    // sample indices are preserved into the CSV by the matching write_csv.
    cm.push_function(
        "calc_per_sample_numerical_fisher_information_detailed",
        dsl::to_typed_function<const cmd::likelihood_algorithm_type&,
                               const var::Parameters_transformed&, const Experiment&,
                               const std::vector<Simulated_Recording<var::please_include<>>>&,
                               std::size_t, std::size_t>(
            &cmd::calculate_n_simulation_mdetailed_predictions, "likelihood_algorithm",
            "parameters", "experiment", "data_series", "sample_min", "sample_max"));

    // -------- Per-group MLE refit + figure_2 battery at θ̂_group --------------
    // Builder for the Gauss-Newton inner-loop options (by value, like
    // build_likelihood_function).
    cm.push_function(
        "build_gauss_newton_options",
        dsl::to_typed_function<double, double, double, std::size_t, double, double, bool>(
            &build_gauss_newton_options_dsl, "lambda_kickoff", "lambda_factor", "lambda_max",
            "max_iter", "grad_rtol", "dvalue_tol", "verbose"));

    // Extract plain Recordings from loaded simulations (input adapter for the
    // calc_MLE_per_group_of_replicates std::vector<Recording> argument).
    cm.push_function(
        "get_recordings_from_simulations",
        dsl::to_typed_function<const std::vector<Simulated_Recording<var::please_include<>>>&>(
            &get_recordings_from_simulations_dsl, "simulations"));

    // Per-group Gauss-Newton MLE refit then the figure_2-style diagnostic
    // battery evaluated at θ̂_group. State pinned to Path A
    // (dMacro_State_Hessian_minimal_param); F_h_relative bound to 1e-5.
    cm.push_function(
        "calc_MLE_per_group_of_replicates",
        dsl::to_typed_function<
            const cmd::likelihood_algorithm_type&, const var::Parameters_transformed&,
            const var::Parameters_transformed&, const Experiment&,
            const std::vector<Recording>&, std::size_t, std::size_t,
            std::size_t, std::set<double>, std::set<double>, std::vector<std::string>,
            std::size_t, const macrodr::optimization::gauss_newton_options&>(
            &calc_MLE_per_group_of_replicates_dsl, "likelihood_algorithm", "theta_warmstart",
            "theta_reference", "experiment", "recordings", "group_size", "n_bootstrap_samples",
            "min_groups_for_bootstrap", "probit_cis", "probit_sample_heights",
            "ranking_variables", "seed", "gauss_newton_options"));

    // write_csv for the MLE_Group_Cloud<State> output. The dedicated
    // write_csv(MLE_Group_Cloud<State>) overload (likelihood.h) handles the
    // plain case; the generic write_csv(Indexed<T>) handles the Indexed case.
    // Alias the long type (pattern: the Analisis_derivative_* writers).
    using MLEGroupCloudMinimal =
        cmd::MLE_Group_Cloud<dMacro_State_Hessian_minimal_param>;
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>, const MLEGroupCloudMinimal&,
                                      std::string>(
            &macrodr::cmd::write_csv<dMacro_State_Hessian_minimal_param>, "analysis", "path"));
    // Indexed cloud: dedicated writer that also emits the per-group _runs.csv
    // companion (the generic write_csv(Indexed<T>) would only write the summary).
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>,
                                      const var::Indexed<MLEGroupCloudMinimal>&, std::string>(
            &macrodr::cmd::write_csv_indexed_cloud<dMacro_State_Hessian_minimal_param>,
            "analysis", "path"));

    // θ̄ handoff: rebuild the full-sample mean of θ̂ as a Parameters_transformed
    // from the cloud (for the downstream figure_2 battery at θ̄). State pinned
    // to Path A (dMacro_State_Hessian_minimal_param).
    cm.push_function(
        "get_parameters_mean",
        dsl::to_typed_return_function<Maybe_error<var::Parameters_transformed>,
                                      const MLEGroupCloudMinimal&>(
            &get_parameters_mean_dsl, "cloud"));

    // figure_3 Fase-2 empirical-vs-theoretical capstone. State pinned to Path A;
    // composes the cloud's bare Cov_emp with the figure_2 numerical-Fisher
    // (at θ_sim and θ̄) and score (at θ̄) outputs.
    cm.push_function(
        "calc_empirical_distortion",
        dsl::to_typed_return_function<Maybe_error<cmd::Empirical_Distortion_Bootstrap>,
                                      const MLEGroupCloudMinimal&,
                                      const std::vector<parameter_spd_payload>&,
                                      const std::vector<parameter_spd_payload>&,
                                      const std::vector<dMacro_State_Ev_gradient_all>&,
                                      std::size_t, std::size_t, std::set<double>>(
            &calc_empirical_distortion_dsl, "cloud", "fim_sim", "fim_bar", "dlik_bar",
            "n_bootstrap", "seed", "probit_cis"));

    // write_csv for the empirical-distortion capstone output (point ++ probit CIs;
    // plain + Indexed).
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>,
                                      const cmd::Empirical_Distortion_Bootstrap&, std::string>(
            static_cast<Maybe_error<std::string> (*)(
                cmd::Empirical_Distortion_Bootstrap const&, std::string)>(&macrodr::cmd::write_csv),
            "analysis", "path"));
    cm.push_function(
        "write_csv",
        dsl::to_typed_return_function<Maybe_error<std::string>,
                                      const var::Indexed<cmd::Empirical_Distortion_Bootstrap>&,
                                      std::string>(
            &macrodr::cmd::write_csv<cmd::Empirical_Distortion_Bootstrap>, "analysis", "path"));

    // Raw-model variants (matches calc_likelihood / calc_dlikelihood pattern):
    // takes ModelPtr and individual approximation flags rather than the
    // pre-built likelihood_algorithm_type.
    cm.push_function(
        "calc_numerical_fisher_information",
        dsl::to_typed_function<const ModelPtr&, const var::Parameters_transformed&,
                               const Experiment&, const Recording&, bool, bool, int, bool, bool,
                               bool, double>(
            &cmd::calculate_numerical_fisher_information, "model", "parameters", "experiment",
            "data", "adaptive_approximation", "recursive_approximation",
            "averaging_approximation", "variance_approximation", "taylor_variance_correction",
            "micro_approximation", "h_rel"));

    cm.push_function(
        "calc_numerical_fisher_information",
        dsl::to_typed_function<const ModelPtr&, const var::Parameters_transformed&,
                               const Experiment&,
                               const Simulated_Recording<var::please_include<>>&,
                               bool, bool, int, bool, bool, bool, double>(
            &cmd::calculate_simulation_numerical_fisher_information, "model", "parameters",
            "experiment", "data", "adaptive_approximation", "recursive_approximation",
            "averaging_approximation", "variance_approximation", "taylor_variance_correction",
            "micro_approximation", "h_rel"));

    cm.push_function(
        "calc_numerical_fisher_information",
        dsl::to_typed_function<const ModelPtr&, const var::Parameters_transformed&,
                               const Experiment&,
                               const std::vector<Simulated_Recording<var::please_include<>>>&,
                               bool, bool, int, bool, bool, bool, double, std::size_t>(
            &cmd::calculate_n_simulation_numerical_fisher_information, "model", "parameters",
            "experiment", "data_series", "adaptive_approximation", "recursive_approximation",
            "averaging_approximation", "variance_approximation", "taylor_variance_correction",
            "micro_approximation", "h_rel", "decimate"));

cm.push_function(
        "calc_likelihood_diagnostic",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&, const Recording&,
                               bool, bool, int, bool, bool, bool>(
            &cmd::calculate_likelihood_diagnostics, "model", "parameters", "experiment", "data",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation"));

cm.push_function(
        "calc_likelihood_diagnostic",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const  Simulated_Recording<var::please_include<>>&&,
                               bool, bool, int, bool, bool, bool>(
            &cmd::calculate_simulation_likelihood_diagnostics, "model", "parameters", "experiment", "simulation",
            "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation"));


    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool, bool>(
            &cmd::calculate_simulation_dlikelihood_predictions, "model", "parameters", "experiment",
            "data", "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation"));

    cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const ModelPtr&,
                               const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<Only_Ch_Curent_Sub_Evolution>>&, bool,
                               bool, int, bool, bool, bool>(
            &cmd::calculate_simulation_sub_dlikelihood_predictions, "model", "parameters", "experiment",
            "data", "adaptive_approximation", "recursive_approximation", "averaging_approximation",
            "variance_approximation", "taylor_variance_correction", "micro_approximation"));


            cm.push_function(
        "calc_dlikelihood_predictions",
        dsl::to_typed_function<const std::string&, const var::Parameters_transformed&, const Experiment&,
                               const Simulated_Recording<var::please_include<>>&, bool,
                               bool, int, bool, bool, bool>(
            &cmd::calculate_simulation_dlikelihood_predictions_model, "model", "parameters",
            "experiment", "data", "adaptive_approximation", "recursive_approximation",
            "averaging_approximation", "variance_approximation", "taylor_variance_correction",
            "micro_approximation"));

      
    cm.push_function(
        "likelihood_derivative_basic_diagnostics",
        dsl::to_typed_function<const std::vector<dMacro_State_Ev_gradient_all>&,
                               const std::vector<parameter_spd_payload>&,
                               std::size_t, std::set<double>, std::size_t, std::size_t>(
            &calculate_likelihood_derivative_basic_diagnostics_dsl,
            "dlikelihood_predictions", "numerical_fisher_information",
            "n_boostrap_samples", "probits", "seed", "max_lag"));

    cm.push_function(
        "likelihood_derivative_basic_diagnostics_paired",
        dsl::to_typed_function<const std::vector<dMacro_State_Ev_gradient_all>&,
                               const std::vector<parameter_spd_payload>&,
                               std::size_t, std::set<double>, std::size_t, std::size_t,
                               std::string>(
            &calculate_likelihood_derivative_basic_diagnostics_paired_dsl,
            "dlikelihood_predictions", "numerical_fisher_information",
            "n_boostrap_samples", "probits", "seed", "max_lag", "samples_dir"));

    cm.push_function(
        "likelihood_derivative_series_var_diagnostics",
        dsl::to_typed_function<const std::vector<dMacro_State_Ev_gradient_all>&,
                               const std::vector<parameter_spd_payload>&,
                               std::size_t, std::set<double>, std::size_t, std::size_t>(
            &calculate_likelihood_derivative_series_var_diagnostics_dsl,
            "dlikelihood_predictions", "numerical_fisher_information",
            "n_boostrap_samples", "probits", "seed", "max_lag"));

    cm.push_function(
        "likelihood_derivative_series_cov_diagnostics",
        dsl::to_typed_function<const std::vector<dMacro_State_Ev_gradient_all>&,
                               const std::vector<parameter_spd_payload>&,
                               std::size_t, std::set<double>, std::size_t, std::size_t>(
            &calculate_likelihood_derivative_series_cov_diagnostics_dsl,
            "dlikelihood_predictions", "numerical_fisher_information",
            "n_boostrap_samples", "probits", "seed", "max_lag"));

    cm.push_function(
        "likelihood_derivative_series_kernel_diagnostics",
        dsl::to_typed_function<const std::vector<dMacro_State_Ev_gradient_all>&,
                               const std::vector<parameter_spd_payload>&,
                               std::size_t, std::set<double>, std::size_t, std::size_t>(
            &calculate_likelihood_derivative_series_kernel_diagnostics_dsl,
            "dlikelihood_predictions", "numerical_fisher_information",
            "n_boostrap_samples", "probits", "seed", "max_lag"));

    cm.push_function(
        "likelihood_derivative_series_kernel_full_diagnostics",
        dsl::to_typed_function<const std::vector<dMacro_State_Ev_gradient_all>&,
                               const std::vector<parameter_spd_payload>&,
                               std::size_t, std::set<double>, std::size_t, std::size_t>(
            &calculate_likelihood_derivative_series_kernel_full_diagnostics_dsl,
            "dlikelihood_predictions", "numerical_fisher_information",
            "n_boostrap_samples", "probits", "seed", "max_lag"));


    return cm;
}

}  // namespace macrodr::cli
