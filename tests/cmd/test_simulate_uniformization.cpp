#include <catch_amalgamated.hpp>

#include <macrodr/cmd/load_experiment.h>
#include <macrodr/cmd/load_model.h>
#include <macrodr/cmd/load_parameters.h>
#include <macrodr/cmd/simulate.h>

#include <cmath>
#include <numeric>
#include <string>
#include <vector>

#include <qmodel.h>

namespace {

struct LoadedScheme {
    macrodr::cmd::ModelPtr model;
    var::Parameters_Transformations transformations;
    var::Parameters_values parameter_values;
    var::Parameters_transformed parameter_transformed;
};

LoadedScheme load_scheme_1() {
    auto maybe_model = macrodr::cmd::load_model("scheme_1");
    REQUIRE(maybe_model);

    auto model = std::move(maybe_model.value());
    auto maybe_parameters = macrodr::cmd::load_parameters(model, "../data/scheme_1_PI_par.csv");
    REQUIRE(maybe_parameters);

    LoadedScheme loaded{std::move(model), std::move(maybe_parameters.value()), {}, {}};
    loaded.parameter_values =
        macrodr::cmd::get_standard_parameter_values(loaded.transformations);
    loaded.parameter_transformed =
        macrodr::cmd::get_standard_parameter_transformed_values(loaded.transformations);
    return loaded;
}

double recording_value(
    const macrodr::Simulated_Recording<var::please_include<>>& simulation, std::size_t i) {
    return get<macrodr::Recording>(simulation())()[i]();
}

std::size_t recording_size(const macrodr::Simulated_Recording<var::please_include<>>& simulation) {
    return get<macrodr::Recording>(simulation())().size();
}

std::size_t seed_number(const macrodr::Simulated_Recording<var::please_include<>>& simulation) {
    return get<macrodr::SeedNumber>(simulation())();
}

std::pair<double, double> mean_and_variance(
    const std::vector<macrodr::Simulated_Recording<var::please_include<>>>& simulations) {
    std::vector<double> values;
    values.reserve(simulations.size());
    for (auto const& simulation : simulations) {
        values.push_back(recording_value(simulation, 0));
    }

    const auto mean =
        std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());
    double centered_sum = 0.0;
    for (auto value : values) {
        centered_sum += (value - mean) * (value - mean);
    }
    const auto variance =
        centered_sum / static_cast<double>(values.size() > 1 ? values.size() - 1 : 1);
    return {mean, variance};
}

}  // namespace

TEST_CASE("uniformization single simulation returns experiment-sized recording",
          "[simulate][uniformization]") {
    auto loaded = load_scheme_1();
    auto experiment =
        macrodr::cmd::create_experiment({{1, 10, 10.0}, {1, 8, 0.0}}, 50e3, 0.0, 0.0);
    auto observations = macrodr::cmd::define_recording({0.0, 0.0});

    auto maybe_simulation =
        macrodr::cmd::run_simulations(loaded.model, loaded.parameter_values, experiment,
                                      observations, std::string("uniformization"), 0, 42);
    REQUIRE(maybe_simulation);

    const auto& simulation = maybe_simulation.value();
    CHECK(recording_size(simulation) == observations().size());
    CHECK(seed_number(simulation) == 42);
}

TEST_CASE("uniformization batch simulation returns requested count",
          "[simulate][uniformization]") {
    auto loaded = load_scheme_1();
    auto experiment = macrodr::cmd::create_experiment({{1, 10, 10.0}}, 50e3, 0.0, 0.0);
    auto observations = macrodr::cmd::define_recording({0.0});

    auto maybe_simulations =
        macrodr::cmd::run_n_simulations(loaded.model, loaded.parameter_transformed, 3, experiment,
                                        observations, std::string("uniformization"), 0, 0);
    REQUIRE(maybe_simulations);

    CHECK(maybe_simulations.value().size() == 3);
}

TEST_CASE("string substeps overload dispatches to substep simulation",
          "[simulate][uniformization]") {
    auto loaded = load_scheme_1();
    auto experiment = macrodr::cmd::create_experiment({{1, 10, 10.0}}, 50e3, 0.0, 0.0);
    auto observations = macrodr::cmd::define_recording({0.0});

    auto maybe_simulation =
        macrodr::cmd::run_simulations(loaded.model, loaded.parameter_values, experiment,
                                      observations, std::string("substeps"), 64, 0);
    REQUIRE(maybe_simulation);
    CHECK(recording_size(maybe_simulation.value()) == observations().size());
}

TEST_CASE("simulate_with_sub_intervals rejects uniformization", "[simulate][uniformization]") {
    auto loaded = load_scheme_1();
    auto experiment = macrodr::cmd::create_experiment({{1, 10, 10.0}}, 50e3, 0.0, 0.0);
    auto observations = macrodr::cmd::define_recording({0.0});

    auto maybe_simulation = macrodr::cmd::run_simulations_with_sub_intervals(
        loaded.model, loaded.parameter_values, experiment, observations,
        std::string("uniformization"), 0);
    REQUIRE_FALSE(maybe_simulation);
    CHECK(maybe_simulation.error()().find("does not support") != std::string::npos);
}

TEST_CASE("legacy simulation_algorithm helper accepts uniformization",
          "[simulate][uniformization]") {
    auto maybe_algorithm =
        macrodr::cmd::set_simulation_algorithm(false, std::string("uniformization"));
    REQUIRE(maybe_algorithm);
    CHECK(maybe_algorithm.value().algorithm == "uniformization");
    CHECK(maybe_algorithm.value().number_of_substeps == 0);
}

TEST_CASE("uniformization and dense substeps agree on first-bin moments",
          "[simulate][uniformization]") {
    auto loaded = load_scheme_1();
    auto experiment = macrodr::cmd::create_experiment({{1, 40, 10.0}}, 50e3, 0.0, 0.0);
    auto observations = macrodr::cmd::define_recording({0.0});

    auto maybe_uniform = macrodr::cmd::run_n_simulations(
        loaded.model, loaded.parameter_values, 64, experiment, observations,
        std::string("uniformization"), 0, 12345);
    REQUIRE(maybe_uniform);

    auto maybe_substeps = macrodr::cmd::run_n_simulations(
        loaded.model, loaded.parameter_values, 64, experiment, observations, 512, 12345);
    REQUIRE(maybe_substeps);

    const auto [uniform_mean, uniform_variance] = mean_and_variance(maybe_uniform.value());
    const auto [substeps_mean, substeps_variance] = mean_and_variance(maybe_substeps.value());

    CHECK(substeps_mean == Catch::Approx(uniform_mean).margin(1.0).epsilon(0.5));
    CHECK(substeps_variance == Catch::Approx(uniform_variance).margin(1.0).epsilon(0.75));
}
