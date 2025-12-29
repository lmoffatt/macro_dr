#pragma once
#include "models_linear.h"



namespace macrodr::simple {



static const auto scheme_CO = Model0::Model("scheme_CO", []() {
    auto names_model =
        std::vector<std::string>{"on", "off",  "unitary_current"};
    auto names_other = std::vector<std::string>{"Current_Noise", 
                                                "Current_Baseline", "Num_ch_mean"};

    std::size_t N = 2ul;

    auto v_Q0_formula = Q0_formula(N);
    v_Q0_formula()[1][0] = "off";
   
    auto v_Qa_formula = Qa_formula(N);
    v_Qa_formula()[0][1] = "kon";
    auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
    v_g_formula()[1] = "unitary_current";

    names_model.insert(names_model.end(), names_other.begin(), names_other.end());
    auto p =
        Matrix<double>(6, 1, std::vector<double>{0.1, 100,1, 1e-3, 0, 5000});

    auto tr = std::vector<std::string>(p.size(), "Log10");
    assert(tr.size() == p.size());
    auto tr_param = var::MyTranformations::from_strings(tr).value();

    auto npar = names_model.size();

    return std::tuple(
        [](const auto& p) -> Maybe_error<Transfer_Op_to<std::decay_t<decltype(p)>, Patch_Model>> {
            auto on = p[0];
            auto off = p[1];
            auto v_unitary_current = p[2] * -1.0;
            auto Npar = 3ul;
            auto v_curr_noise = p[Npar];
            auto v_baseline = p[Npar + 1];
            //  auto v_Num_ch_mean=p[Npar+2];
            //  auto v_std_log_Num_ch=p[Npar+3];

            auto v_N0 = p[std::pair(Npar + 2, Npar + 2)];

            auto N = N_St(2);
            return build<Patch_Model>(
                N_St(N()),
                build<Q0>(var::build_<Matrix<double>>(
                    N(), N(), {{1, 0}},
                    {off})),
                build<Qa>(var::build_<Matrix<double>>(N(), N(), {{0, 1}},
                                                      {on})),
                build<P_initial>(
                    var::build_<Matrix<double>>(1, N(), {{0, 0}}, {on * (1.0 / on)})),
                build<g>(var::build_<Matrix<double>>(N(), 1, {{1, 0}}, {v_unitary_current})),
                build<N_Ch_mean>(v_N0),

                build<Current_Noise>(v_curr_noise), build<Pink_Noise>(0.0),
                build<Proportional_Noise>(0.0),
                build<Current_Baseline>(v_baseline), N_Ch_mean_time_segment_duration(121),
                Binomial_magical_number(5.0), min_P(1e-12), Probability_error_tolerance(1e-2),
                Conductance_variance_error_tolerance(1e-2));
        },
        [npar](const auto& patch_model)
            -> Maybe_error<Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>> {
            auto& v_Q0 = get<Q0>(patch_model);
            auto& v_Qa = get<Qa>(patch_model);
            auto& v_g = get<g>(patch_model);
            Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>> out =
                Matrix<double>(1, npar, 0.0);

            assert(get<N_St>(patch_model)() == 2);

            auto on = v_Qa()(0ul, 1ul) ;
            out[0] = on;
            auto off = v_Q0()(1ul, 0ul);
            out[1] = off;
            auto Npar = 3ul;

            auto v_curr_noise = get<Current_Noise>(patch_model);
            out[Npar] = v_curr_noise();

           
            auto v_baseline = get<Current_Baseline>(patch_model);
            out[Npar + 1] = v_baseline();

            auto v_N0 = get<N_Ch_mean>(patch_model);
            out.set(std::pair(Npar + 2, Npar + 2), v_N0());

            return out;
        },
        p, names_model, std::move(v_Q0_formula), std::move(v_Qa_formula), std::move(v_g_formula),
        std::move(tr_param));
});


static const auto scheme_CCO = Model0::Model("scheme_CCO", []() {
    auto names_model =
        std::vector<std::string>{"kon", "koff", "gating_on", "gating_off", "unitary_current"};
    auto names_other = std::vector<std::string>{"Current_Noise",
                                                "Current_Baseline", "Num_ch_mean"};

    std::size_t N = 3ul;

    auto v_Q0_formula = Q0_formula(N);
    v_Q0_formula()[1][0] = "koff";
    v_Q0_formula()[1][2] = "gating_on";
    v_Q0_formula()[2][1] = "gating_off";

    auto v_Qa_formula = Qa_formula(N);
    v_Qa_formula()[0][1] = "kon";
    auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
    v_g_formula()[2] = "unitary_current";

    names_model.insert(names_model.end(), names_other.begin(), names_other.end());
    auto p =
        Matrix<double>(8, 1, std::vector<double>{6.73, 166, 743, 45.3, 1, 1e-3, 0, 5000});

    auto tr = std::vector<std::string>(p.size(), "Log10");
    assert(tr.size() == p.size());
    auto tr_param = var::MyTranformations::from_strings(tr).value();

    auto npar = names_model.size();

    return std::tuple(
        [](const auto& p) -> Maybe_error<Transfer_Op_to<std::decay_t<decltype(p)>, Patch_Model>> {
            auto kon = p[0];
            auto koff = p[1];
            auto gating_on = p[2];
            auto gating_off = p[3];
            auto v_unitary_current = p[4] * -1.0;
            auto Npar = 5ul;
            auto v_curr_noise = p[Npar];
            auto v_baseline = p[Npar + 1];
            //  auto v_Num_ch_mean=p[Npar+2];
            //  auto v_std_log_Num_ch=p[Npar+3];

            auto v_N0 = p[std::pair(Npar + 2, Npar + 2)];

            auto N = N_St(3);
            return build<Patch_Model>(
                N_St(N()),
                build<Q0>(var::build_<Matrix<double>>(
                    N(), N(), {{1, 0}, {1, 2}, {2, 1}},
                    {koff, gating_on, gating_off})),
                build<Qa>(var::build_<Matrix<double>>(N(), N(), {{0, 1}},
                                                      { kon})),
                build<P_initial>(
                    var::build_<Matrix<double>>(1, N(), {{0, 0}}, {kon * (1.0 / kon)})),
                build<g>(var::build_<Matrix<double>>(N(), 1, {{2, 0}}, {v_unitary_current})),
                build<N_Ch_mean>(v_N0),

                 build<Current_Noise>(v_curr_noise), build<Pink_Noise>(0.0),
                build<Proportional_Noise>(0.0),
               build<Current_Baseline>(v_baseline), N_Ch_mean_time_segment_duration(121),
                Binomial_magical_number(5.0), min_P(1e-12), Probability_error_tolerance(1e-2),
                Conductance_variance_error_tolerance(1e-2));
        },
        [npar](const auto& patch_model) 
            -> Maybe_error<Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>> {
            auto& v_Q0 = get<Q0>(patch_model);
            auto& v_Qa = get<Qa>(patch_model);
            auto& v_g = get<g>(patch_model);
            Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>> out =
                Matrix<double>(1, npar, 0.0);

            assert(get<N_St>(patch_model)() == 3);

            auto kon = v_Qa()(0ul, 1ul);
            out[0] = kon;
            auto koff = v_Q0()(1ul, 0ul);
            out[1] = koff;
            auto gating_on = v_Q0()(1ul, 2ul);
            out[2] = gating_on;

            auto gating_off = v_Q0()(2ul, 1ul);
            out[3] = gating_off;
            auto v_unitary_current = v_g()[2ul] * -1.0;
            out[4] = v_unitary_current;
            auto Npar = 5ul;

            auto v_curr_noise = get<Current_Noise>(patch_model);
            out[Npar] = v_curr_noise();


            auto v_baseline = get<Current_Baseline>(patch_model);
            out[Npar + 1] = v_baseline();

            auto v_N0 = get<N_Ch_mean>(patch_model);
            out.set(std::pair(Npar + 2, Npar + 2), v_N0());

            return out;
        },
        p, names_model, std::move(v_Q0_formula), std::move(v_Qa_formula), std::move(v_g_formula),
        std::move(tr_param));
});

static const auto scheme_COC = Model0::Model("scheme_COC", []() {
    auto names_model =
        std::vector<std::string>{"kon", "koff", "gating_on", "gating_off", "unitary_current"};
    auto names_other = std::vector<std::string>{"Current_Noise",
                                                "Current_Baseline", "Num_ch_mean"};

    std::size_t N = 3ul;

    auto v_Q0_formula = Q0_formula(N);
    v_Q0_formula()[1][0] = "off";
    v_Q0_formula()[1][2] = "inactivating_on";
    v_Q0_formula()[2][1] = "inactivating_off";

    auto v_Qa_formula = Qa_formula(N);
    v_Qa_formula()[0][1] = "on";
    auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
    v_g_formula()[2] = "unitary_current";

    names_model.insert(names_model.end(), names_other.begin(), names_other.end());
    auto p =
        Matrix<double>(8, 1, std::vector<double>{6.73, 166, 743, 45.3, 1, 1e-3,  0, 5000});

    auto tr = std::vector<std::string>(p.size(), "Log10");
    assert(tr.size() == p.size());
    auto tr_param = var::MyTranformations::from_strings(tr).value();

    auto npar = names_model.size();

    return std::tuple(
        [](const auto& p) -> Maybe_error<Transfer_Op_to<std::decay_t<decltype(p)>, Patch_Model>> {
            auto on = p[0];
            auto off = p[1];
            auto inactivating_on = p[2];
            auto inactivating_off = p[3];
            auto v_unitary_current = p[4] * -1.0;
            auto Npar = 5ul;
            auto v_curr_noise = p[Npar];
            auto v_baseline = p[Npar + 1];
            //  auto v_Num_ch_mean=p[Npar+2];
            //  auto v_std_log_Num_ch=p[Npar+3];

            auto v_N0 = p[std::pair(Npar + 2, Npar + 2)];

            auto N = N_St(3);
            return build<Patch_Model>(
                N_St(N()),
                build<Q0>(var::build_<Matrix<double>>(
                    N(), N(), {{1, 0}, {1, 2}, {2, 1}},
                    {off, inactivating_on, inactivating_off})),
                build<Qa>(var::build_<Matrix<double>>(N(), N(), {{0, 1}},
                                                      { on})),
                build<P_initial>(
                    var::build_<Matrix<double>>(1, N(), {{0, 0}}, {on * (1.0 / on)})),
                build<g>(var::build_<Matrix<double>>(N(), 1, {{1, 0}}, {v_unitary_current})),
                build<N_Ch_mean>(v_N0),

                build<Current_Noise>(v_curr_noise), build<Pink_Noise>(0.0),
                build<Proportional_Noise>(0.0),
                build<Current_Baseline>(v_baseline), N_Ch_mean_time_segment_duration(121),
                Binomial_magical_number(5.0), min_P(1e-12), Probability_error_tolerance(1e-2),
                Conductance_variance_error_tolerance(1e-2));
        },
        [npar](const auto& patch_model)
            -> Maybe_error<Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>> {
            auto& v_Q0 = get<Q0>(patch_model);
            auto& v_Qa = get<Qa>(patch_model);
            auto& v_g = get<g>(patch_model);
            Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>> out =
                Matrix<double>(1, npar, 0.0);

            assert(get<N_St>(patch_model)() == 3);

            auto kon = v_Qa()(0ul, 1ul);
            out[0] = kon;
            auto koff = v_Q0()(1ul, 0ul);
            out[1] = koff;
            auto gating_on = v_Q0()(1ul, 2ul);
            out[2] = gating_on;

            auto gating_off = v_Q0()(2ul, 1ul);
            out[3] = gating_off;
            auto v_unitary_current = v_g()[1ul] * -1.0;
            out[4] = v_unitary_current;
            auto Npar = 5ul;

            auto v_curr_noise = get<Current_Noise>(patch_model);
            out[Npar] = v_curr_noise();


            auto v_baseline = get<Current_Baseline>(patch_model);
            out[Npar + 1] = v_baseline();

            auto v_N0 = get<N_Ch_mean>(patch_model);
            out.set(std::pair(Npar + 2, Npar + 2), v_N0());

            return out;
        },
        p, names_model, std::move(v_Q0_formula), std::move(v_Qa_formula), std::move(v_g_formula),
        std::move(tr_param));
});



} // namespace macrodr::p2x2


