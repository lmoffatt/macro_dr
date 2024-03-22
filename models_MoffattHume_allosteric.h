#ifndef MODELS_MOFFATTHUME_ALLOSTERIC_H
#define MODELS_MOFFATTHUME_ALLOSTERIC_H

#include "models_MoffattHume_linear.h"

namespace macrodr {

static auto scheme_6 = Allost1::Model("scheme6", []() {
  auto v_binding = Conformational_change_label{"Binding"};
  auto v_rocking = Conformational_change_label{"Rocking"};
  auto v_gating = Conformational_change_label{"Gating"};

  auto mo = make_Conformational_model_standarized(
      Agonist_dependency_map{
          std::map<Conformational_change_label, Agonist_dependency>{
              {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
              {v_rocking, Agonist_dependency{}},
              {v_gating, Agonist_dependency{}}}},
      std::vector<Conformational_change_label>{v_binding, v_binding, v_binding,
                                               v_rocking, v_gating},
      std::vector<Conformational_interaction>{
          {Vector_Space{
               Conformational_interaction_label{"BR"},
               Conformational_interaction_players{{v_binding, v_rocking}},
               Conformational_interaction_positions{{{0, 3}, {1, 3}, {2, 3}}}},
           Vector_Space{
               Conformational_interaction_label{"BG"},
               Conformational_interaction_players{{v_binding, v_gating}},
               Conformational_interaction_positions{{{0, 4}, {1, 4}, {2, 4}}}},
           Vector_Space{
               Conformational_interaction_label{"RG"},
               Conformational_interaction_players{{v_rocking, v_gating}},
               Conformational_interaction_positions{{{3, 4}}}}

          }},
      std::vector<Conductance_interaction>{
          Vector_Space{Conductance_interaction_label{"Gating_Current"},
                       Conductance_interaction_players{{v_gating}},
                       Conductance_interaction_positions{{{{4}}}}}},
      std::map<Conformational_change_label, Conformation_change_standard_state>{
          {v_rocking,
           Conformation_change_standard_state{
               Conformational_interactions_domain_state{std::map<
                   Vector_Space<Conformational_interaction_index,
                                Conformational_interaction_subposition>,
                   int>{
                   {Vector_Space{Conformational_interaction_index{0ul},
                                 Conformational_interaction_subposition{1}},
                    3}}}}},
          {v_gating,
           Conformation_change_standard_state{
               Conformational_interactions_domain_state{std::map<
                   Vector_Space<Conformational_interaction_index,
                                Conformational_interaction_subposition>,
                   int>{
                   {Vector_Space{Conformational_interaction_index{1ul},
                                 Conformational_interaction_subposition{1}},
                    3},
                   {Vector_Space{Conformational_interaction_index{2ul},
                                 Conformational_interaction_subposition{1}},
                    1}}}}}});

  assert(mo);
  auto m = std::move(mo.value());

  auto names = make_ModelNames<Allost1>(m);

  auto names_vec_untransformed = std::vector<std::string>{
      "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
      "Gating_on",  "Gating_off",  "BR",         "BR_0",
      "BR_1",       "BG",          "BG_0",       "BG_1",
      "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
  auto names_vec = std::vector<std::string>{"Binding_Act_on",
                                            "Binding_Act_off",
                                            "Rocking_on_B",
                                            "Rocking_off_B",
                                            "Gating_on_BR",
                                            "Gating_off_BR",
                                            "BR",
                                            "BR_Bon",
                                            "BR_Ron",
                                            "BG",
                                            "BG_Bon",
                                            "BG_Gon",
                                            "RG",
                                            "RG_Ron",
                                            "RG_Gon",
                                            "Gating_Current"};

  auto names_other = std::vector<std::string>{"Current_Noise", "Pink_Noise",
                                              "Proportional_Noise",
                                              "Current_Baseline", "Num_ch"};

  auto p_kinetics =
      std::vector<double>{9.28, 1871, 3875, 1.07, 914, 776, 65.1 * 1.15, 1.15,
                          33.3, 1.77, 0.77, 1.77, 123, 123, 635,         1};
  auto p_other = std::vector<double>{1e-3, 5, 1e-2, 1, 4800};

  p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
  auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

  auto tr = std::vector<std::string>(p.size(), "Log10");
  tr[tr.size() - 2] = "Linear";
  assert(tr.size() == p.size());
  auto tr_param = var::MyTranformations::from_strings(tr).value();

  assert(names() == names_vec_untransformed);

  names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

  auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
  assert(Maybe_modeltyple_formula);
  auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
      std::move(Maybe_modeltyple_formula.value());
  return std::tuple(
      [names, m](const auto &t_p)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(t_p)>, Patch_Model>> {
        using std::pow;
          auto Binding_on = t_p[0];
        auto Binding_off = t_p[1];
          auto Rocking_on_B = t_p[2];
        auto Rocking_off_B = t_p[3];
          
          auto Gating_on_BR = t_p[4];
        auto Gating_off_BR = t_p[5];
          auto BR = t_p[6];
        auto BR_Bon = t_p[7];
          auto BR_Ron = t_p[8];
        auto BG = t_p[9];
          auto BG_Bon = t_p[10];
        auto BG_Gon = t_p[11];
          auto RG = t_p[12];
        auto RG_Ron = t_p[13];
          auto RG_Gon = t_p[14];
        auto Gating_Current = t_p[15] * (-1.0);

        auto Rocking_on = Rocking_on_B;
        auto Rocking_off = Rocking_off_B;

        auto Gating_on = Gating_off_BR;
        auto Gating_off = Gating_off_BR;

        auto BR_0 = log(BR_Bon) / log(BR);
        auto BR_1 = log(BR_Ron) / log(BR);

        auto BG_0 = log(BG_Bon) / log(BG);
        auto BG_1 = log(BG_Gon) / log(BG);

        auto RG_0 = log(RG_Ron) / log(RG);
        auto RG_1 = log(RG_Gon) / log(RG);

        auto p = t_p;
        p[2] = Rocking_on;
        p[3] = Rocking_off;
        
        p[4] = Gating_on;
        p[5] = Gating_off;
        
        p[7] = BR_0;
        p[8] = BR_1;
        
        p[10] = BG_0;
        p[11] = BG_1;
        
        p[13] = RG_0;
        p[14] = RG_1;
        p[15] = Gating_Current;
        auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
        // std::cerr<<"parameters\n"<<p();

        assert(Maybe_Q0Qag);
        auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
        auto Npar = names().size();
        
        // auto v_Inac_rate = t_p()[Npar];
        auto v_curr_noise = t_p[Npar];
        auto v_pink_noise = t_p[Npar + 1];
        auto v_prop_noise = t_p[Npar + 2];
        auto v_baseline = t_p[Npar + 3];
        auto v_N0 = t_p[std::pair{Npar + 4, Npar + 4}];
        auto Nst = get<N_St>(m());
        auto Maybe_v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
            a_Q0, a_Qa, ATP_concentration(0.0), Nst);
        if (!Maybe_v_P_initial)
          return Maybe_v_P_initial.error();
        else
          return build<Patch_Model>(
              N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
              std::move(Maybe_v_P_initial.value()), std::move(a_g),
              build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
              build<Pink_Noise>(v_pink_noise),
              build<Proportional_Noise>(v_prop_noise),
              build<Current_Baseline>(v_baseline),
              N_Ch_mean_time_segment_duration(120000),
              Binomial_magical_number(5.0), min_P(1e-7),
              Probability_error_tolerance(1e-2),
              Conductance_variance_error_tolerance(1e-2));
      },
      p, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula,
      std::move(tr_param));
});

// static auto model6 = Allost1::Model("model6", []() {
//   auto v_binding = Conformational_change_label{"Binding"};
//   auto v_rocking = Conformational_change_label{"Rocking"};
//   auto v_gating = Conformational_change_label{"Gating"};

//   auto mo = make_Conformational_model(
//       Agonist_dependency_map{

//           std::map<Conformational_change_label, Agonist_dependency>{
//               {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
//               {v_rocking, Agonist_dependency{}},
//               {v_gating, Agonist_dependency{}}}},
//       std::vector<Conformational_change_label>{v_binding, v_binding,
//       v_binding,
//                                                v_rocking, v_gating},
//       std::vector<Conformational_interaction>{
//           {Vector_Space{
//                Conformational_interaction_label{"BR"},
//                Conformational_interaction_players{{v_binding, v_rocking}},
//                Conformational_interaction_positions{{{0, 3}, {1, 3}, {2,
//                3}}}},
//            Vector_Space{
//                Conformational_interaction_label{"BG"},
//                Conformational_interaction_players{{v_binding, v_gating}},
//                Conformational_interaction_positions{{{0, 4}, {1, 4}, {2,
//                4}}}},
//            Vector_Space{
//                Conformational_interaction_label{"RG"},
//                Conformational_interaction_players{{v_rocking, v_gating}},
//                Conformational_interaction_positions{{{3, 4}}}}

//           }},
//       std::vector<Conductance_interaction>{
//           Vector_Space{Conductance_interaction_label{"Gating_Current"},
//                        Conductance_interaction_players{{v_gating}},
//                        Conductance_interaction_positions{{{{4}}}}}});

//   assert(mo);
//   auto m = std::move(mo.value());

//   auto names = make_ModelNames<Allost1>(m);

//   auto names_vec = std::vector<std::string>{
//       "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
//       "Gating_on",  "Gating_off",  "BR",         "BR_0",
//       "BR_1",       "BG",          "BG_0",       "BG_1",
//       "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
//   auto names_other = std::vector<std::string>{
//       "Inactivation_rate",
//       "Current_Noise","Pink_Noise","Proportional_Noise","Current_Baseline",
//       "Num_ch"};

//   auto p_kinetics = std::vector<double>{
//       10, 10000, 100, 10000, 1, 10000, 10, 1, 1, 10, 1, 1, 10, 1, 1, 1};
//   auto p_other = std::vector<double>{1e-3, 1e-3, 5, 1e-2, 1, 5000};

//   p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
//   auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

//    auto tr=std::vector<std::string>(p.size(),"Log10");
//   tr[tr.size()-2]= "Linear";
//   assert(tr.size()==p.size());
//   auto tr_param=var::MyTranformations::from_strings(tr).value();

//   assert(names() == names_vec);

//   names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

//   auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
//   assert(Maybe_modeltyple_formula);
//   auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
//       std::move(Maybe_modeltyple_formula.value());
//   return std::tuple(
//       [names, m](const auto &pp)
//           -> Maybe_error<
//               Transfer_Op_to<std::decay_t<decltype(pp)>, Patch_Model>> {
//         using std::pow;
//         auto p = build<Parameters_Transformations<Allost1>>(
//             pp.IdName(), pp.names(),
//             pp());

//         p()[names["BR_0"].value()] =
//             p()[names["BR_0"].value()] / (1.0 + p()[names["BR_0"].value()]);
//         p()[names["BR_1"].value()] =
//             p()[names["BR_1"].value()] / (1.0 + p()[names["BR_1"].value()]);

//         p()[names["BG_0"].value()] =
//             p()[names["BG_0"].value()] / (1.0 + p()[names["BG_0"].value()]);
//         p()[names["BG_1"].value()] =
//             p()[names["BG_1"].value()] / (1.0 + p()[names["BG_1"].value()]);

//         p()[names["RG_0"].value()] =
//             p()[names["RG_0"].value()] / (1.0 + p()[names["RG_0"].value()]);
//         p()[names["RG_1"].value()] =
//             p()[names["RG_1"].value()] / (1.0 + p()[names["RG_1"].value()]);

//         p()[names["Gating_Current"].value()] =
//             p()[names["Gating_Current"].value()] * -1.0;
//         auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
//         assert(Maybe_Q0Qag);
//         auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());

//         auto Npar = names().size();

//         auto v_Inac_rate = p()[Npar];
//         auto v_curr_noise = p()[Npar + 1];
//         auto v_pink_noise = p()[Npar+2];
//         auto v_prop_noise = p()[Npar+3];
//         auto v_baseline = p()[Npar + 4];
//         auto v_N0 = p()[std::pair{Npar + 5, Npar + 5}];
//         auto Nst = get<N_St>(m());
//         auto v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
//             a_Q0, a_Qa, ATP_concentration(0.0), Nst);
//         if (v_P_initial.valid())
//           return add_Patch_inactivation(
//               build<Patch_Model>(
//                   N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
//                   std::move(v_P_initial.value()), std::move(a_g),
//                   build<N_Ch_mean>(v_N0),
//                     build<Current_Noise>(v_curr_noise),
//                     build<Pink_Noise>(v_pink_noise),
//                     build<Proportional_Noise>(v_prop_noise),
//                     build<Current_Baseline>(v_baseline),
//                   N_Ch_mean_time_segment_duration(120000),
//                   Binomial_magical_number(5.0), min_P(1e-7),
//                   Probability_error_tolerance(1e-2),
//                   Conductance_variance_error_tolerance(1e-2)),
//               v_Inac_rate);
//         else
//           return v_P_initial.error();
//       },
//       p, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula,
//       std::move(tr_param));
// });

// static auto model6_no_inactivation =
//     Allost1::Model("model6_no_inactivation", []() {
//       auto v_binding = Conformational_change_label{"Binding"};
//       auto v_rocking = Conformational_change_label{"Rocking"};
//       auto v_gating = Conformational_change_label{"Gating"};

//       auto mo = make_Conformational_model(
//           Agonist_dependency_map{

//               std::map<Conformational_change_label, Agonist_dependency>{
//                   {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
//                   {v_rocking, Agonist_dependency{}},
//                   {v_gating, Agonist_dependency{}}}},
//           std::vector<Conformational_change_label>{
//               v_binding, v_binding, v_binding, v_rocking, v_gating},
//           std::vector<Conformational_interaction>{
//               {Vector_Space{
//                    Conformational_interaction_label{"BR"},
//                    Conformational_interaction_players{{v_binding,
//                    v_rocking}}, Conformational_interaction_positions{
//                        {{0, 3}, {1, 3}, {2, 3}}}},
//                Vector_Space{
//                    Conformational_interaction_label{"BG"},
//                    Conformational_interaction_players{{v_binding, v_gating}},
//                    Conformational_interaction_positions{
//                        {{0, 4}, {1, 4}, {2, 4}}}},
//                Vector_Space{
//                    Conformational_interaction_label{"RG"},
//                    Conformational_interaction_players{{v_rocking, v_gating}},
//                    Conformational_interaction_positions{{{3, 4}}}}

//               }},
//           std::vector<Conductance_interaction>{
//               Vector_Space{Conductance_interaction_label{"Gating_Current"},
//                            Conductance_interaction_players{{v_gating}},
//                            Conductance_interaction_positions{{{{4}}}}}});

//       assert(mo);
//       auto m = std::move(mo.value());

//       auto names = make_ModelNames<Allost1>(m);

//       auto names_vec = std::vector<std::string>{
//           "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
//           "Gating_on",  "Gating_off",  "BR",         "BR_0",
//           "BR_1",       "BG",          "BG_0",       "BG_1",
//           "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //-->
//           8
//       auto names_other =
//       std::vector<std::string>{"Current_Noise","Pink_Noise","Proportional_Noise",
//                                                   "Current_Baseline",
//                                                   "Num_ch"};

//       auto p_kinetics = std::vector<double>{
//           10, 1000, 1000, 100000, 1, 100, 100, 1, 1, 1, 1, 1, 100, 1, 1, 1};
//       auto p_Moffatt_Hume_transformed = std::vector<double>{
//           9.28,   1871,      2547.88,  295207, 0.220378,  150.312,
//           74.865, 0.0323846, 0.187903, 1.77,   -0.457748, 1,
//           123,    1,         1.3411,   1};
//       auto p_other = std::vector<double>{1e-3, 1.0, 5000};

//       p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
//       auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

//        auto tr=std::vector<std::string>(p.size(),"Log10");
//   tr[tr.size()-2]= "Linear";
//   assert(tr.size()==p.size());
//   auto tr_param=var::MyTranformations::from_strings(tr).value();

//       assert(names() == names_vec);

//       names_vec.insert(names_vec.end(), names_other.begin(),
//       names_other.end());

//       auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
//       assert(Maybe_modeltyple_formula);
//       auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
//           std::move(Maybe_modeltyple_formula.value());
//       return std::tuple(
//           [names, m](const auto &logp)
//               -> Maybe_error<
//                   Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>>
//                   {
//             using std::pow;
//               auto p = build<Parameters_Transformations<Allost1>>(
//                 logp.IdName(), logp.names(),

//                 apply([](const auto &x) { return pow(10.0, x); }, logp()));

//             p()[names["BR_0"].value()] =
//                 p()[names["BR_0"].value()] / (1.0 +
//                 p()[names["BR_0"].value()]);
//             p()[names["BR_1"].value()] =
//                 p()[names["BR_1"].value()] / (1.0 +
//                 p()[names["BR_1"].value()]);

//             p()[names["BG_0"].value()] =
//                 p()[names["BG_0"].value()] / (1.0 +
//                 p()[names["BG_0"].value()]);
//             p()[names["BG_1"].value()] =
//                 p()[names["BG_1"].value()] / (1.0 +
//                 p()[names["BG_1"].value()]);

//             p()[names["RG_0"].value()] =
//                 p()[names["RG_0"].value()] / (1.0 +
//                 p()[names["RG_0"].value()]);
//             p()[names["RG_1"].value()] =
//                 p()[names["RG_1"].value()] / (1.0 +
//                 p()[names["RG_1"].value()]);

//             p()[names["Gating_Current"].value()] =
//                 p()[names["Gating_Current"].value()] * -1.0;
//             auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
//             assert(Maybe_Q0Qag);
//             auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
//             auto Npar = names().size();

//             // auto v_Inac_rate = p()[Npar];
//             auto v_curr_noise = p()[Npar];
//             auto v_pink_noise = p()[Npar+1];
//             auto v_prop_noise = p()[Npar+2];
//             auto v_baseline = logp()[Npar + 3];
//             auto v_N0 = p()[std::pair{Npar + 4, Npar + 4}];
//             auto Nst = get<N_St>(m());
//             auto v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
//                 a_Q0, a_Qa, ATP_concentration(0.0), Nst);
//             if (v_P_initial.valid())

//               return build<Patch_Model>(
//                   N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
//                   std::move(v_P_initial.value()), std::move(a_g),
//                   build<N_Ch_mean>(v_N0),
//                     build<Current_Noise>(v_curr_noise),
//                     build<Pink_Noise>(v_pink_noise),
//                     build<Proportional_Noise>(v_prop_noise),
//                     build<Current_Baseline>(v_baseline),
//                   N_Ch_mean_time_segment_duration(120000),
//                   Binomial_magical_number(5.0), min_P(1e-7),
//                   Probability_error_tolerance(1e-2),
//                   Conductance_variance_error_tolerance(1e-2));
//             else
//               return v_P_initial.error();
//           },
//           p, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula,
//           std::move(tr_param));
//     });

// static auto model6_Eff_no_inactivation =
//     Allost1::Model("model6_Eff_no_inactivation", []() {
//       auto v_binding = Conformational_change_label{"Binding"};
//       auto v_rocking = Conformational_change_label{"Rocking"};
//       auto v_gating = Conformational_change_label{"Gating"};

//       auto mo = make_Conformational_model(
//           Agonist_dependency_map{
//               std::map<Conformational_change_label, Agonist_dependency>{
//                   {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
//                   {v_rocking, Agonist_dependency{}},
//                   {v_gating, Agonist_dependency{}}}},
//           std::vector<Conformational_change_label>{
//               v_binding, v_binding, v_binding, v_rocking, v_gating},
//           std::vector<Conformational_interaction>{
//               {Vector_Space{
//                    Conformational_interaction_label{"BR"},
//                    Conformational_interaction_players{{v_binding,
//                    v_rocking}}, Conformational_interaction_positions{
//                        {{0, 3}, {1, 3}, {2, 3}}}},
//                Vector_Space{
//                    Conformational_interaction_label{"BG"},
//                    Conformational_interaction_players{{v_binding, v_gating}},
//                    Conformational_interaction_positions{
//                        {{0, 4}, {1, 4}, {2, 4}}}},
//                Vector_Space{
//                    Conformational_interaction_label{"RG"},
//                    Conformational_interaction_players{{v_rocking, v_gating}},
//                    Conformational_interaction_positions{{{3, 4}}}}

//               }},
//           std::vector<Conductance_interaction>{
//               Vector_Space{Conductance_interaction_label{"Gating_Current"},
//                            Conductance_interaction_players{{v_gating}},
//                            Conductance_interaction_positions{{{{4}}}}}});

//       assert(mo);
//       auto m = std::move(mo.value());

//       auto names = make_ModelNames<Allost1>(m);

//       auto names_vec_untransformed = std::vector<std::string>{
//           "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
//           "Gating_on",  "Gating_off",  "BR",         "BR_0",
//           "BR_1",       "BG",          "BG_0",       "BG_1",
//           "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //-->
//           8
//       auto names_vec = std::vector<std::string>{"Binding_Act_on",
//                                                 "Binding_Act_off",
//                                                 "Rocking_on_B",
//                                                 "Rocking_off_B",
//                                                 "Gating_on_BR",
//                                                 "Gating_off_BR",
//                                                 "BR",
//                                                 "BR_Bon",
//                                                 "BR_Ron",
//                                                 "BG",
//                                                 "BG_Bon",
//                                                 "BG_Gon",
//                                                 "RG",
//                                                 "RG_Ron",
//                                                 "RG_Gon",
//                                                 "Gating_Current"};

//       auto names_other =
//       std::vector<std::string>{"Current_Noise","Pink_Noise","Proportional_Noise",
//                                                   "Current_Baseline",
//                                                   "Num_ch"};

//       auto p_kinetics = std::vector<double>{
//           9.28, 1871, 3875, 1.07, 914, 776, 65.1 * 1.15, 1.15,
//           33.3, 1.77, 0.77, 1.77, 123, 123, 635,         1};
//       auto p_other = std::vector<double>{1e-3, 5, 1e-2, 1, 4800};

//       p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
//       auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

//        auto tr=std::vector<std::string>(p.size(),"Log10");
//   tr[tr.size()-2]= "Linear";
//   assert(tr.size()==p.size());
//   auto tr_param=var::MyTranformations::from_strings(tr).value();

//       assert(names() == names_vec_untransformed);

//       names_vec.insert(names_vec.end(), names_other.begin(),
//       names_other.end());

//       auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
//       assert(Maybe_modeltyple_formula);
//       auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
//           std::move(Maybe_modeltyple_formula.value());
//       return std::tuple(
//           [names, m](const auto &logp)
//               -> Maybe_error<
//                   Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>>
//                   {
//             using std::pow;
//             auto tr_p = build<Parameters<Allost1>>(
//                 logp.IdName(), logp.names(),
//                 apply([](const auto &x) { return pow(10.0, x); }, logp()));
//             auto Binding_on = tr_p()[0];
//             auto Binding_off = tr_p()[1];
//             auto Rocking_on_B = tr_p()[2];
//             auto Rocking_off_B = tr_p()[3];

//             auto Gating_on_BR = tr_p()[4];
//             auto Gating_off_BR = tr_p()[5];
//             auto BR = tr_p()[6];
//             auto BR_Bon = tr_p()[7];
//             auto BR_Ron = tr_p()[8];
//             auto BG = tr_p()[9];
//             auto BG_Bon = tr_p()[10];
//             auto BG_Gon = tr_p()[11];
//             auto RG = tr_p()[12];
//             auto RG_Ron = tr_p()[13];
//             auto RG_Gon = tr_p()[14];
//             auto Gating_Current = tr_p()[15] * (-1.0);

//             auto Rocking_on = Rocking_on_B / pow(BR_Bon, 3);
//             auto Rocking_off = Rocking_off_B * pow(BR / BR_Bon, 3);

//             auto Gating_on = Gating_off_BR / pow(BG_Gon, 3) / RG_Gon;
//             auto Gating_off = Gating_off_BR * pow(BG / BG_Gon, 3) * RG /
//             RG_Gon;

//             auto BR_0 = log(BR_Bon) / log(BR);
//             auto BR_1 = log(BR_Ron) / log(BR);

//             auto BG_0 = log(BG_Bon) / log(BG);
//             auto BG_1 = log(BG_Gon) / log(BG);

//             auto RG_0 = log(RG_Ron) / log(RG);
//             auto RG_1 = log(RG_Gon) / log(RG);

//             auto p = tr_p;
//             p()[2] = Rocking_on;
//             p()[3] = Rocking_off;

//             p()[4] = Gating_on;
//             p()[5] = Gating_off;

//             p()[7] = BR_0;
//             p()[8] = BR_1;

//             p()[10] = BG_0;
//             p()[11] = BG_1;

//             p()[13] = RG_0;
//             p()[14] = RG_1;
//             p()[15] = Gating_Current;
//             auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
//             // std::cerr<<"parameters\n"<<p();

//             assert(Maybe_Q0Qag);
//             auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
//             auto Npar = names().size();

//             // auto v_Inac_rate = p()[Npar];
//             auto v_curr_noise = tr_p()[Npar];
//             auto v_pink_noise = tr_p()[Npar+1];
//             auto v_prop_noise = tr_p()[Npar+2];
//             auto v_baseline = logp()[Npar + 3];
//             auto v_N0 = tr_p()[std::pair{Npar + 4, Npar + 4}];
//             auto Nst = get<N_St>(m());
//             auto Maybe_v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
//                 a_Q0, a_Qa, ATP_concentration(0.0), Nst);
//             if (!Maybe_v_P_initial)
//               return Maybe_v_P_initial.error();
//             else
//               return build<Patch_Model>(
//                   N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
//                   std::move(Maybe_v_P_initial.value()), std::move(a_g),
//                   build<N_Ch_mean>(v_N0),
//                     build<Current_Noise>(v_curr_noise),
//                     build<Pink_Noise>(v_pink_noise),
//                     build<Proportional_Noise>(v_prop_noise),
//                     build<Current_Baseline>(v_baseline),
//                   N_Ch_mean_time_segment_duration(120000),
//                   Binomial_magical_number(5.0), min_P(1e-7),
//                   Probability_error_tolerance(1e-2),
//                   Conductance_variance_error_tolerance(1e-2));
//           },
//           p, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula,
//           std::move(tr_param));
//     });

// static auto model6_Eff_std = Allost1::Model("model6_Eff_std", []() {
//   auto v_binding = Conformational_change_label{"Binding"};
//   auto v_rocking = Conformational_change_label{"Rocking"};
//   auto v_gating = Conformational_change_label{"Gating"};

//   auto mo = make_Conformational_model_standarized(
//       Agonist_dependency_map{
//           std::map<Conformational_change_label, Agonist_dependency>{
//               {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
//               {v_rocking, Agonist_dependency{}},
//               {v_gating, Agonist_dependency{}}}},
//       std::vector<Conformational_change_label>{v_binding, v_binding,
//       v_binding,
//                                                v_rocking, v_gating},
//       std::vector<Conformational_interaction>{
//           {Vector_Space{
//                Conformational_interaction_label{"BR"},
//                Conformational_interaction_players{{v_binding, v_rocking}},
//                Conformational_interaction_positions{{{0, 3}, {1, 3}, {2,
//                3}}}},
//            Vector_Space{
//                Conformational_interaction_label{"BG"},
//                Conformational_interaction_players{{v_binding, v_gating}},
//                Conformational_interaction_positions{{{0, 4}, {1, 4}, {2,
//                4}}}},
//            Vector_Space{
//                Conformational_interaction_label{"RG"},
//                Conformational_interaction_players{{v_rocking, v_gating}},
//                Conformational_interaction_positions{{{3, 4}}}}

//           }},
//       std::vector<Conductance_interaction>{
//           Vector_Space{Conductance_interaction_label{"Gating_Current"},
//                        Conductance_interaction_players{{v_gating}},
//                        Conductance_interaction_positions{{{{4}}}}}},
//       std::map<Conformational_change_label,
//       Conformation_change_standard_state>{
//           {v_rocking,
//            Conformation_change_standard_state{
//                Conformational_interactions_domain_state{std::map<
//                    Vector_Space<Conformational_interaction_index,
//                                 Conformational_interaction_subposition>,
//                    int>{
//                    {Vector_Space{Conformational_interaction_index{0ul},
//                                  Conformational_interaction_subposition{1}},
//                     3}}}}},
//           {v_gating,
//            Conformation_change_standard_state{
//                Conformational_interactions_domain_state{std::map<
//                    Vector_Space<Conformational_interaction_index,
//                                 Conformational_interaction_subposition>,
//                    int>{
//                    {Vector_Space{Conformational_interaction_index{1ul},
//                                  Conformational_interaction_subposition{1}},
//                     3},
//                    {Vector_Space{Conformational_interaction_index{2ul},
//                                  Conformational_interaction_subposition{1}},
//                     1}}}}}});

//   assert(mo);
//   auto m = std::move(mo.value());

//   auto names = make_ModelNames<Allost1>(m);

//   auto names_vec_untransformed = std::vector<std::string>{
//       "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
//       "Gating_on",  "Gating_off",  "BR",         "BR_0",
//       "BR_1",       "BG",          "BG_0",       "BG_1",
//       "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
//   auto names_vec = std::vector<std::string>{"Binding_Act_on",
//                                             "Binding_Act_off",
//                                             "Rocking_on_B",
//                                             "Rocking_off_B",
//                                             "Gating_on_BR",
//                                             "Gating_off_BR",
//                                             "BR",
//                                             "BR_Bon",
//                                             "BR_Ron",
//                                             "BG",
//                                             "BG_Bon",
//                                             "BG_Gon",
//                                             "RG",
//                                             "RG_Ron",
//                                             "RG_Gon",
//                                             "Gating_Current"};

//   auto names_other =
//       std::vector<std::string>{"Current_Noise","Pink_Noise","Proportional_Noise",
//       "Current_Baseline", "Num_ch"};

//   auto p_kinetics =
//       std::vector<double>{9.28, 1871, 3875, 1.07, 914, 776, 65.1
//       * 1.15, 1.15,
//                           33.3, 1.77, 0.77, 1.77, 123, 123, 635,         1};
//   auto p_other = std::vector<double>{1e-3,  5, 1e-2,1, 4800};

//   p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
//   auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

//    auto tr=std::vector<std::string>(p.size(),"Log10");
//   tr[tr.size()-2]= "Linear";
//   assert(tr.size()==p.size());
//   auto tr_param=var::MyTranformations::from_strings(tr).value();

//   assert(names() == names_vec_untransformed);

//   names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

//   auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
//   assert(Maybe_modeltyple_formula);
//   auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
//       std::move(Maybe_modeltyple_formula.value());
//   return std::tuple(
//       [names, m](const auto &logp)
//           -> Maybe_error<
//               Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
//         using std::pow;
//         auto tr_p = build<Parameters<Allost1>>(
//             logp.IdName(), logp.names(),
//             apply([](const auto &x) { return pow(10.0, x); }, logp()));
//         auto Binding_on = tr_p()[0];
//         auto Binding_off = tr_p()[1];
//         auto Rocking_on_B = tr_p()[2];
//         auto Rocking_off_B = tr_p()[3];

//         auto Gating_on_BR = tr_p()[4];
//         auto Gating_off_BR = tr_p()[5];
//         auto BR = tr_p()[6];
//         auto BR_Bon = tr_p()[7];
//         auto BR_Ron = tr_p()[8];
//         auto BG = tr_p()[9];
//         auto BG_Bon = tr_p()[10];
//         auto BG_Gon = tr_p()[11];
//         auto RG = tr_p()[12];
//         auto RG_Ron = tr_p()[13];
//         auto RG_Gon = tr_p()[14];
//         auto Gating_Current = tr_p()[15] * (-1.0);

//         auto Rocking_on = Rocking_on_B;
//         auto Rocking_off = Rocking_off_B;

//         auto Gating_on = Gating_off_BR;
//         auto Gating_off = Gating_off_BR;

//         auto BR_0 = log(BR_Bon) / log(BR);
//         auto BR_1 = log(BR_Ron) / log(BR);

//         auto BG_0 = log(BG_Bon) / log(BG);
//         auto BG_1 = log(BG_Gon) / log(BG);

//         auto RG_0 = log(RG_Ron) / log(RG);
//         auto RG_1 = log(RG_Gon) / log(RG);

//         auto p = tr_p;
//         p()[2] = Rocking_on;
//         p()[3] = Rocking_off;

//         p()[4] = Gating_on;
//         p()[5] = Gating_off;

//         p()[7] = BR_0;
//         p()[8] = BR_1;

//         p()[10] = BG_0;
//         p()[11] = BG_1;

//         p()[13] = RG_0;
//         p()[14] = RG_1;
//         p()[15] = Gating_Current;
//         auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
//         // std::cerr<<"parameters\n"<<p();

//         assert(Maybe_Q0Qag);
//         auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
//         auto Npar = names().size();

//         // auto v_Inac_rate = p()[Npar];
//         auto v_curr_noise = tr_p()[Npar];
//         auto v_pink_noise = tr_p()[Npar+1];
//         auto v_prop_noise = tr_p()[Npar+2];
//         auto v_baseline = logp()[Npar + 3];
//         auto v_N0 = tr_p()[std::pair{Npar + 4, Npar + 4}];
//         auto Nst = get<N_St>(m());
//         auto Maybe_v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
//             a_Q0, a_Qa, ATP_concentration(0.0), Nst);
//         if (!Maybe_v_P_initial)
//           return Maybe_v_P_initial.error();
//         else
//           return build<Patch_Model>(
//               N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
//               std::move(Maybe_v_P_initial.value()), std::move(a_g),
//               build<N_Ch_mean>(v_N0),
//                 build<Current_Noise>(v_curr_noise),
//                 build<Pink_Noise>(v_pink_noise),
//                 build<Proportional_Noise>(v_prop_noise),
//                 build<Current_Baseline>(v_baseline),
//               N_Ch_mean_time_segment_duration(120000),
//               Binomial_magical_number(5.0), min_P(1e-7),
//               Probability_error_tolerance(1e-2),
//               Conductance_variance_error_tolerance(1e-2));
//       },
//       p, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula,
//       std::move(tr_param));
// });

// static auto model7 = Allost1::Model("model7", []() {
//   auto v_rocking = Conformational_change_label{"Rocking"};
//   auto v_binding = Conformational_change_label{"Binding"};
//   auto v_gating = Conformational_change_label{"Gating"};

//   auto mo = make_Conformational_model(
//       Agonist_dependency_map{

//           std::map<Conformational_change_label, Agonist_dependency>{
//               {v_rocking, Agonist_dependency{}},
//               {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
//               {v_gating, Agonist_dependency{}}}},
//       std::vector<Conformational_change_label>{v_binding, v_rocking,
//       v_binding,
//                                                v_rocking, v_binding,
//                                                v_rocking, v_gating},
//       std::vector<Conformational_interaction>{
//           {Vector_Space{
//                Conformational_interaction_label{"BR"},
//                Conformational_interaction_players{{v_binding, v_rocking}},
//                Conformational_interaction_positions{{{0, 1}, {2, 3}, {4,
//                5}}}},
//            Vector_Space{
//                Conformational_interaction_label{"BG"},
//                Conformational_interaction_players{{v_binding, v_gating}},
//                Conformational_interaction_positions{{{0, 6}, {2, 6}, {4,
//                6}}}},
//            Vector_Space{
//                Conformational_interaction_label{"RG"},
//                Conformational_interaction_players{{v_rocking, v_gating}},
//                Conformational_interaction_positions{{{1, 6}, {3, 6}, {5,
//                6}}}}

//           }},
//       std::vector<Conductance_interaction>{
//           Vector_Space{Conductance_interaction_label{"Gating_Current"},
//                        Conductance_interaction_players{{v_gating}},
//                        Conductance_interaction_positions{{{{6}}}}}});

//   assert(mo);
//   auto m = std::move(mo.value());

//   auto names = make_ModelNames<Allost1>(m);

//   auto names_vec = std::vector<std::string>{
//       "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
//       "Gating_on",  "Gating_off",  "BR",         "BR_0",
//       "BR_1",       "BG",          "BG_0",       "BG_1",
//       "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
//   auto names_other = std::vector<std::string>{
//       "Inactivation_rate", "Current_Noise","Pink_Noise","Proportional_Noise",
//       "Current_Baseline", "Num_ch"};

//   auto p_kinetics = std::vector<double>{
//       10, 10000, 100, 10000, 1, 10000, 10, 1, 1, 10, 1, 1, 10, 1, 1, 1};
//   auto p_other = std::vector<double>{1, 1e-3, 5, 1e-2, 1, 5000};

//   p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());

//   auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

//    auto tr=std::vector<std::string>(p.size(),"Log10");
//   tr[tr.size()-2]= "Linear";
//   assert(tr.size()==p.size());
//   auto tr_param=var::MyTranformations::from_strings(tr).value();

//   assert(names() == names_vec);

//   names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

//   auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
//   assert(Maybe_modeltyple_formula);
//   auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
//       std::move(Maybe_modeltyple_formula.value());
//   return std::tuple(
//       [names, m](const auto &logp)
//           -> Maybe_error<
//               Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
//         using std::pow;
//         auto p = build<Parameters<Allost1>>(
//             logp.IdName(), logp.names(),

//             apply([](const auto &x) { return pow(10.0, x); }, logp()));

//         p()[names["BR_0"].value()] =
//             p()[names["BR_0"].value()] / (1.0 + p()[names["BR_0"].value()]);
//         p()[names["BR_1"].value()] =
//             p()[names["BR_1"].value()] / (1.0 + p()[names["BR_1"].value()]);

//         p()[names["BG_0"].value()] =
//             p()[names["BG_0"].value()] / (1.0 + p()[names["BG_0"].value()]);
//         p()[names["BG_1"].value()] =
//             p()[names["BG_1"].value()] / (1.0 + p()[names["BG_1"].value()]);

//         p()[names["RG_0"].value()] =
//             p()[names["RG_0"].value()] / (1.0 + p()[names["RG_0"].value()]);
//         p()[names["RG_1"].value()] =
//             p()[names["RG_1"].value()] / (1.0 + p()[names["RG_1"].value()]);

//         p()[names["Gating_Current"].value()] =
//             p()[names["Gating_Current"].value()] * -1.0;
//         auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
//         assert(Maybe_Q0Qag);
//         auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
//         auto Npar = names().size();

//         auto v_Inac_rate = p()[Npar];
//         auto v_curr_noise = p()[Npar + 1];
//         auto v_pink_noise = p()[Npar+2];
//         auto v_prop_noise = p()[Npar+3];
//         auto v_baseline = logp()[Npar + 4];
//         auto v_N0 = p[std::pair(Npar + 5, Npar + 5)];
//         auto Nst = get<N_St>(m());
//         auto v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
//             a_Q0, a_Qa, ATP_concentration(0.0), Nst);

//         if (v_P_initial.valid())
//           return add_Patch_inactivation(
//               build<Patch_Model>(
//                   N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
//                   std::move(v_P_initial.value()), std::move(a_g),
//                   build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
//                   build<Current_Baseline>(v_baseline),
//                   N_Ch_mean_time_segment_duration(120000),
//                   Binomial_magical_number(5.0), min_P(1e-7),
//                   Probability_error_tolerance(1e-2),
//                   Conductance_variance_error_tolerance(1e-2)),
//               v_Inac_rate);
//         else
//           return v_P_initial.error();
//       },
//       p, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula,
//       std::move(tr_param));
// });

// static auto model8 = Allost1::Model("model8", []() {
//   auto v_binding = Conformational_change_label{"Binding"};
//   auto v_rocking = Conformational_change_label{"Rocking"};
//   auto mo = make_Conformational_model(
//       Agonist_dependency_map{
//           std::map<Conformational_change_label, Agonist_dependency>{
//               {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
//               {v_rocking, Agonist_dependency{}}}},
//       std::vector<Conformational_change_label>{v_binding, v_rocking,
//       v_binding,
//                                                v_rocking, v_binding,
//                                                v_rocking},
//       std::vector<Conformational_interaction>{{Vector_Space{
//           Conformational_interaction_label{"RBR"},
//           Conformational_interaction_players{{v_rocking, v_binding,
//           v_rocking}}, Conformational_interaction_positions{{{5, 0, 1},
//                                                 {1, 0, 5},
//                                                 {1, 2, 3},
//                                                 {3, 2, 1},
//                                                 {3, 4, 5},
//                                                 {5, 4, 3}}}}}},
//       std::vector<Conductance_interaction>{Vector_Space{
//           Conductance_interaction_label{"Rocking_Current_factor"},
//           Conductance_interaction_players{{v_rocking}},
//           Conductance_interaction_positions{{{{1}}, {{3}}, {{5}}}}}});

//   assert(mo);
//   auto m = std::move(mo.value());

//   auto names = make_ModelNames<Allost1>(m);

//   auto names_vec =
//       std::vector<std::string>{
//           "Binding_on",  "Binding_off", "Rocking_on",
//           "Rocking_off", "RBR",         "RBR_0",
//           "RBR_1",       "RBR_2",       "Rocking_Current_factor"}; //--> 8
//   auto names_other =
//       std::vector<std::string>{"Inactivation_rate", "Leaking_current",
//                                "Current_Noise","Pink_Noise","Proportional_Noise",
//                                "Current_Baseline", "Num_ch"};

//   auto p_kinetics =
//       std::vector<double>{10, 10000, 100, 10000, 100, 1.0, 1e-2, 1.0, 100};
//   auto p_other = std::vector<double>{1e-4, 0.01, 1e-3, 5, 1e-2, 1, 5000};

//   p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
//   auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

//    auto tr=std::vector<std::string>(p.size(),"Log10");
//   tr[tr.size()-2]= "Linear";
//   assert(tr.size()==p.size());
//   auto tr_param=var::MyTranformations::from_strings(tr).value();

//   assert(names() == names_vec);

//   names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

//   auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
//   assert(Maybe_modeltyple_formula);
//   auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
//       std::move(Maybe_modeltyple_formula.value());
//   return std::tuple(
//       [names, m](const auto &logp)
//           -> Maybe_error<
//               Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
//         using std::pow;
//         auto p = build<Parameters<Allost1>>(
//             logp.IdName(), logp.names(),

//             apply([](const auto &x) { return pow(10.0, x); }, logp()));
//         p()[names["RBR_0"].value()] =
//             p()[names["RBR_0"].value()] / (1.0 +
//             p()[names["RBR_0"].value()]);
//         p()[names["RBR_1"].value()] =
//             p()[names["RBR_1"].value()] / (1.0 +
//             p()[names["RBR_1"].value()]);
//         p()[names["RBR_2"].value()] =
//             p()[names["RBR_2"].value()] / (1.0 +
//             p()[names["RBR_2"].value()]);

//         auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
//         assert(Maybe_Q0Qag);
//         auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
//         auto Npar = names().size();

//         auto v_Inac_rate = p()[Npar];
//         auto v_leaking_current = p()[Npar + 1];
//         auto v_curr_noise = p()[Npar + 2];
//         auto v_pink_noise = p()[Npar+3];
//         auto v_prop_noise = p()[Npar+4];
//         auto v_baseline = logp()[Npar + 5];
//         auto v_N0 = p[std::pair(Npar + 6, Npar + 6)];

//         auto v_g = build<g>(apply(
//             [&v_leaking_current](const auto &x) {
//               return v_leaking_current * pow(10.0, x) * (-1.0);
//             },
//             a_g()));
//         auto Nst = get<N_St>(m());
//         auto v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
//             a_Q0, a_Qa, ATP_concentration(0.0), Nst);

//         if (v_P_initial.valid())
//           return add_Patch_inactivation(
//               build<Patch_Model>(
//                   N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
//                   std::move(v_P_initial.value()), std::move(v_g),
//                   build<N_Ch_mean>(v_N0),
//                     build<Current_Noise>(v_curr_noise),
//                     build<Pink_Noise>(v_pink_noise),
//                     build<Proportional_Noise>(v_prop_noise),
//                     build<Current_Baseline>(v_baseline),
//                   N_Ch_mean_time_segment_duration(120000),
//                   Binomial_magical_number(5.0), min_P(1e-7),
//                   Probability_error_tolerance(1e-2),
//                   Conductance_variance_error_tolerance(1e-2)),
//               v_Inac_rate);
//         else
//           return v_P_initial.error();
//       },
//       p, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula,
//       std::move(tr_param));
// });

// static auto model9 = Allost1::Model("model9", []() {
//   auto v_binding = Conformational_change_label{"Binding"};
//   auto v_rocking = Conformational_change_label{"Rocking"};
//   auto mo = make_Conformational_model(
//       Agonist_dependency_map{
//           std::map<Conformational_change_label, Agonist_dependency>{
//               {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
//               {v_rocking, Agonist_dependency{}}}},
//       std::vector<Conformational_change_label>{v_binding, v_rocking,
//       v_binding,
//                                                v_rocking, v_binding,
//                                                v_rocking},
//       std::vector<Conformational_interaction>{{Vector_Space{
//           Conformational_interaction_label{"RB"},
//           Conformational_interaction_players{{v_rocking, v_binding}},
//           Conformational_interaction_positions{
//               {{5, 0}, {1, 0}, {1, 2}, {3, 2}, {3, 4}, {5, 4}}}}}},
//       std::vector<Conductance_interaction>{Vector_Space{
//           Conductance_interaction_label{"Rocking_Current_factor"},
//           Conductance_interaction_players{{v_rocking}},
//           Conductance_interaction_positions{{{{1}}, {{3}}, {{5}}}}}});

//   assert(mo);
//   auto m = std::move(mo.value());

//   auto names = make_ModelNames<Allost1>(m);

//   auto names_vec =
//       std::vector<std::string>{"Binding_on", "Binding_off",
//                                "Rocking_on", "Rocking_off",
//                                "RB",         "RB_0",
//                                "RB_1",       "Rocking_Current_factor"}; //-->
//                                8
//   auto names_other =
//       std::vector<std::string>{"Inactivation_rate", "leaking_current",
//                                "Current_Noise","Pink_Noise","Proportional_Noise",
//                                "Current_Baseline", "Num_ch"};

//   auto p_kinetics =
//       std::vector<double>{10, 10000, 100, 10000, 100, 1.0, 1e-2, 100};
//   auto p_other = std::vector<double>{1e-3, 1e-2, 1e-3, 5, 1e-2, 1, 5000};

//   p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
//   auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

//    auto tr=std::vector<std::string>(p.size(),"Log10");
//   tr[tr.size()-2]= "Linear";
//   assert(tr.size()==p.size());
//   auto tr_param=var::MyTranformations::from_strings(tr).value();

//   assert(names() == names_vec);

//   names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

//   auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
//   assert(Maybe_modeltyple_formula);
//   auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
//       std::move(Maybe_modeltyple_formula.value());
//   return std::tuple(
//       [names, m](const auto &logp)
//           -> Maybe_error<
//               Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
//         using std::pow;
//         auto p = build<Parameters<Allost1>>(
//             logp.IdName(), logp.names(),
//             apply([](const auto &x) { return pow(10.0, x); }, logp()));
//         p()[names["RB_0"].value()] =
//             p()[names["RB_0"].value()] / (1.0 + p()[names["RB_0"].value()]);
//         p()[names["RB_1"].value()] =
//             p()[names["RB_1"].value()] / (1.0 + p()[names["RB_1"].value()]);

//         auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
//         assert(Maybe_Q0Qag);
//         auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
//         auto Npar = names().size();

//         auto v_Inac_rate = p()[Npar];
//         auto v_leaking_current = p()[Npar + 1];
//         auto v_curr_noise = p()[Npar + 2];
//         auto v_pink_noise = p()[Npar+3];
//         auto v_prop_noise = p()[Npar+4];
//         auto v_baseline = logp()[Npar + 5];
//         auto v_N0 = p()[std::pair(Npar + 6, Npar + 6)];

//         auto v_g = build<g>(apply(
//             [&v_leaking_current](const auto &x) {
//               return v_leaking_current * pow(10.0, x) * (-1.0);
//             },
//             a_g()));
//         auto Nst = get<N_St>(m());
//         auto v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
//             a_Q0, a_Qa, ATP_concentration(0.0), Nst);
//         if (v_P_initial.valid())

//           return add_Patch_inactivation(
//               build<Patch_Model>(
//                   N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
//                   std::move(v_P_initial.value()), std::move(v_g),
//                   build<N_Ch_mean>(v_N0),
//                     build<Current_Noise>(v_curr_noise),
//                     build<Pink_Noise>(v_pink_noise),
//                     build<Proportional_Noise>(v_prop_noise),
//                   build<Current_Baseline>(v_baseline),
//                   N_Ch_mean_time_segment_duration(120000),
//                   Binomial_magical_number(5.0), min_P(1e-7),
//                   Probability_error_tolerance(1e-2),
//                   Conductance_variance_error_tolerance(1e-2)),
//               v_Inac_rate);
//         else
//           return v_P_initial.error();
//       },
//       p, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula,
//       std::move(tr_param));
// });

inline auto get_model(std::string modelName) {
  auto allmodels =
      // Models_Library(&model00, &model00_7, &model01, &model4, &model4_g_lin,
      //                             &model6, &model6_no_inactivation,
      //                             &model6_Eff_no_inactivation, &model7,
      //                             &model8, &model9);
      Models_Library(&scheme_1, &scheme_2, &scheme_3, &scheme_4, &scheme_1_d,
                       &scheme_2_d, &scheme_3_d, &scheme_4_d, &scheme_6);
  //                      &model00, &model00_7, &model01, &model4,
  //                      &model4_g_lin,
  //                    &model6, &model6_no_inactivation,
  //                    &model6_Eff_no_inactivation, &model6_Eff_std, &model7,
  //                    &model8, &model9);
  return allmodels[modelName];
}

inline Maybe_error<std::size_t> get_num_parameters(std::string model) {
  auto maybe_model = get_model(model);
  if (!maybe_model)
    return maybe_model.error();

  return std::visit(
      [&](auto model0ptr) {
        return model0ptr->parameters_transformations().names().size();
      },
      maybe_model.value());
}

inline auto get_model_scheme(std::string modelName) {
  auto allmodels = // Models_Library(&scheme_1);
      Models_Library(&scheme_1, &scheme_2, &scheme_3, &scheme_4, &scheme_1_d,
                     &scheme_2_d, &scheme_3_d, &scheme_4_d, &scheme_6);
  return allmodels[modelName];
}

inline void print_model_Priors(double covar) {
  auto allmodels = // Models_Library(&scheme_1);
      Models_Library(&scheme_1, &scheme_2, &scheme_3, &scheme_4, &scheme_1_d,
                     &scheme_2_d, &scheme_3_d, &scheme_4_d, &scheme_6);
  //,
  //
  //                     &model6, &model6_no_inactivation,
  //                     &model6_Eff_no_inactivation, &model7, &model8,
  //                     &model9);

  std::apply(
      [&covar](auto... ptr_models) {
        (
            [&covar](auto modelp) {
              auto &tr_par = modelp->parameters_transformations();
              //   var::Parameters_Transformations tr_par(par,
              //   modelp->parameters_transformations());
              auto prior = var::prior_around(tr_par, covar);
              var::write_Parameters(tr_par.IdName() + "_par.csv", ",", tr_par);
              write_Prior(tr_par.IdName() + "_prior.csv", ",", prior);
            }(ptr_models),
            ...);
      },
      allmodels());
}

// inline auto get_model_old(std::string modelName) -> std::variant<
//     /*decltype(&model4),*/ decltype(&model6_Eff_no_inactivation)> {
//   using return_type = std::variant<
//       /*decltype(&model4), */ decltype(&model6_Eff_no_inactivation)>;

//   //  if (modelName=="model4")
//   //     return return_type(&model4);
//   // else
//   return return_type(&model6_Eff_no_inactivation);
//   // }

  using Model_v = decltype(get_model(std::string{}));

} // namespace macrodr

#endif // MODELS_MOFFATTHUME_ALLOSTERIC_H
