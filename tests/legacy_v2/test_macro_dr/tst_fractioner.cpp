
#include <cstddef>
#include <vector>

#include "../catch2/catch.hpp"
#include "../qmodel.h"

TEST_CASE("Another test with Catch2", "[fancy]") {
    using namespace macrodr;
    std::vector<double> yv = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};

    macrodr::Simulated_Recording<macrodr::includes_N_state_evolution(false)> sim;
    std::transform(yv.begin(), yv.end(), std::back_inserter(get<Recording>(sim())()),
                   [](double x) { return Patch_current(x); });

    Experiment x(Recording_conditions{}, Frequency_of_Sampling(50e3),
                 initial_agonist_concentration(Agonist_concentration(0)));
    std::vector<std::tuple<double, double, double>> ev = {
        {0, 50, 0},      {1e-3, 100, 1},   {3e-3, 50, 0},   {0.8e-3, 100, 1},
        {1.4e-3, 50, 0}, {1.8e-3, 100, 1}, {0, 50, 0},      {1e-3, 100, 1},
        {3e-3, 50, 0},   {0.8e-3, 100, 1}, {1.4e-3, 50, 0}, {1.8e-3, 100, 1}};

    std::transform(ev.begin(), ev.end(), std::back_inserter(get<Recording_conditions>(x)()),
                   [](std::tuple<double, double, double> x0) {
                       return Experiment_step(
                           Time(get<0>(x0)),
                           Agonist_evolution(std::vector<Agonist_step>{Agonist_step(
                               number_of_samples(get<1>(x0)), Agonist_concentration(get<2>(x0)))}));
                   });

    auto init_seed = calc_seed(0);

    mt_64i mt(init_seed);

    macrodr::experiment_fractioner frac(std::vector<std::size_t>{}, 0);
    std::vector<std::vector<std::size_t>> indexes{
        {2, 5, 9}, {2, 3, 5, 7, 9, 10}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}};
    auto n_frac = size(indexes);
    deprecated::by_fraction<Simulated_Recording<macrodr::includes_N_state_evolution(false)>> y_out(
        n_frac);
    deprecated::by_fraction<Experiment> x_out(
        n_frac, Experiment(Recording_conditions{}, get<Frequency_of_Sampling>(x),
                           get<initial_agonist_concentration>(x)));

    y_out[n_frac - 1] = sim;
    x_out[n_frac - 1] = x;

    for (std::size_t i = n_frac - 1; i > 0; --i) {
        std::tie(get<Recording_conditions>(x_out[i - 1]), y_out[i - 1]) =
            macrodr::experiment_fractioner::average_Recording(
                get<Recording_conditions>(x_out[i]), y_out[i], indexes[i], indexes[i - 1], false);
    }

    CHECK(x_out == deprecated::by_fraction<Experiment>{});
    CHECK(y_out ==
          std::vector<macrodr::Simulated_Recording<macrodr::includes_N_state_evolution(false)>>{});
}
