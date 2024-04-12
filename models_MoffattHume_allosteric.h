#ifndef MODELS_MOFFATTHUME_ALLOSTERIC_H
#define MODELS_MOFFATTHUME_ALLOSTERIC_H

#include "allosteric_models.h"
#include "models_MoffattHume_linear.h"
#include "variables.h"

namespace macrodr {


static auto scheme_6 = Allost1::Model("scheme_6", []() {
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
                    
                    1}}}}}},
      Conductance_interaction_info(Conductance_interaction_kind::additive,""));

  assert(mo);
  auto m = std::move(mo.value());

  auto names = make_ModelNames<Allost1>(m);

  auto names_vec = std::vector<std::string>{
      "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
      "Gating_on",  "Gating_off",  "BR",         "BR_0",
      "BR_1",       "BG",          "BG_0",       "BG_1",
      "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
  
  auto names_other = std::vector<std::string>{ "Current_Noise",    "Pink_Noise",
      "Proportional_Noise", "Current_Baseline", "Num_ch"};

  auto p_kinetics = std::vector<double>{
                                        9.28, 1871, 3875, 1.07, 914, 776,  65.1 * 1.15, 1.15,
      33.3, 1.77, 1.77, 0.77, 123, 1.00, 635,         1};
  auto p_other = std::vector<double>{1e-3, 5, 1e-2, 1, 4800};

  p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
  auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

  auto tr = std::vector<std::string>(p.size(), "Log10");
  tr[tr.size() - 2] = "Linear";
  assert(tr.size() == p.size());
  auto tr_param = var::MyTranformations::from_strings(tr).value();

  assert(names() == names_vec);

  names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

  auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
  assert(Maybe_modeltyple_formula);
  auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
      std::move(Maybe_modeltyple_formula.value());
  return std::tuple(
      [names, m](const auto &t_p)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(t_p)>, Patch_Model>> {
        auto Maybe_Q0Qag = make_Model<Allost1>(m, names, t_p);
        // std::cerr<<"parameters\n"<<p();

        assert(Maybe_Q0Qag);
        auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
        
        auto Npar = names().size();
        
        // auto v_Inac_rate = t_p()[Npar];
       // auto v_unitary_current = t_p[Npar-1] * -1.0;
        
        auto v_curr_noise = t_p[Npar ];
        auto v_pink_noise = t_p[Npar + 1];
        auto v_prop_noise = t_p[Npar + 2];
        auto v_baseline = t_p[Npar + 3];
        auto v_N0 = t_p[std::pair{Npar + 4, Npar + 4}];
        a_g() = a_g() * -1.0;
        
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
      [names, m](const auto &patch_model)
      -> Maybe_error<
          Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>> {
          
          auto par= get_Parameters_from_Q0_Qa_g(m,names,get<Q0>(patch_model),get<Qa>(patch_model),build<g>(get<g>(patch_model)()*(-1.0)));
          
          if(!par)
              return par.error();
          auto Npar = names().size();
          
          auto out=Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>(Npar+5,1);
          
          out.set(std::pair(0ul,Npar-1),par.value());
          
          out[Npar]=get<Current_Noise>(patch_model)();
          out[Npar+1]=get<Pink_Noise>(patch_model)();
          out[Npar+2]=get<Proportional_Noise>(patch_model)();
          out[Npar+3]=get<Current_Baseline>(patch_model)();
          out.set(std::pair{Npar + 4, Npar + 4},get<N_Ch_mean>(patch_model)());
          return out;
          
          
      },
      p, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula,
      std::move(tr_param));
});


static auto scheme_6_d = add_Patch_inactivation_to_model<Allost1>(
    scheme_6, 1e-5, var::MyTranformations::from_string("Log10").value());


static auto scheme_7 = Allost1::Model("scheme_7", []() {
    auto v_binding = Conformational_change_label{"Binding"};
    auto v_rocking = Conformational_change_label{"Rocking"};
    auto v_gating = Conformational_change_label{"Gating"};
    
    auto mo = make_Conformational_model_standarized(
      Agonist_dependency_map{
                               std::map<Conformational_change_label, Agonist_dependency>{
                                                                      {v_rocking, Agonist_dependency{}},
                {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
                {v_gating, Agonist_dependency{}}}},
        std::vector<Conformational_change_label>{v_binding, v_rocking, v_binding,
                                                 v_rocking, v_binding, v_rocking,
                                                 v_gating},
        std::vector<Conformational_interaction>{
                                                {Vector_Space{
                             Conformational_interaction_label{"BR"},
                    Conformational_interaction_players{{v_binding, v_rocking}},
                    Conformational_interaction_positions{{{0, 1}, {2, 3}, {4, 5}}}},
                Vector_Space{
                             Conformational_interaction_label{"BG"},
                    Conformational_interaction_players{{v_binding, v_gating}},
                    Conformational_interaction_positions{{{0, 6}, {2, 6}, {4, 6}}}},
                Vector_Space{
                             Conformational_interaction_label{"RG"},
                    Conformational_interaction_players{{v_rocking, v_gating}},
                    Conformational_interaction_positions{{{1, 6}, {3, 6}, {5, 6}}}}
                
            }},
        std::vector<Conductance_interaction>{
                                             Vector_Space{Conductance_interaction_label{"Gating_Current"},
                         Conductance_interaction_players{{v_gating}},
                         Conductance_interaction_positions{{{{6}}}}}},
        std::map<Conformational_change_label, Conformation_change_standard_state>{
                                                                                  {v_rocking,
             Conformation_change_standard_state{
                                                Conformational_interactions_domain_state{std::map<
                     Vector_Space<Conformational_interaction_index,
                                  Conformational_interaction_subposition>,
                     int>{
                          {Vector_Space{Conformational_interaction_index{0ul},
                                   Conformational_interaction_subposition{1}},
                      1}}}}},
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
                      
                      3}}}}}},
        
        Conductance_interaction_info(Conductance_interaction_kind::additive,""));
        
    assert(mo);
    auto m = std::move(mo.value());
    
    auto names = make_ModelNames<Allost1>(m);
    
    auto names_vec = std::vector<std::string>{
                                              "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
        "Gating_on",  "Gating_off",  "BR",         "BR_0",
        "BR_1",       "BG",          "BG_0",       "BG_1",
        "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
    
    auto names_other = std::vector<std::string>{ "Current_Noise",    "Pink_Noise",
        "Proportional_Noise", "Current_Baseline", "Num_ch"};
    
    auto p_kinetics = std::vector<double>{
                                          12.0, 6731, 9035, 242,  1380,       110,  4295, 4295,
        4295, 0.39, 1.00, 0.39, 118 * 5.31, 5.31, 62.1, 1};
    auto p_other = std::vector<double>{1e-3, 5, 1e-2, 1, 4800};
    
    p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
    auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);
    
    auto tr = std::vector<std::string>(p.size(), "Log10");
    tr[tr.size() - 2] = "Linear";
    assert(tr.size() == p.size());
    auto tr_param = var::MyTranformations::from_strings(tr).value();
    
    assert(names() == names_vec);
    
    names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());
    
    auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
    assert(Maybe_modeltyple_formula);
    auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
        std::move(Maybe_modeltyple_formula.value());
    return std::tuple(
        [names, m](const auto &t_p)
        -> Maybe_error<
                Transfer_Op_to<std::decay_t<decltype(t_p)>, Patch_Model>> {
            auto Maybe_Q0Qag = make_Model<Allost1>(m, names, t_p);
            // std::cerr<<"parameters\n"<<p();
            
            assert(Maybe_Q0Qag);
            auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
            
            auto Npar = names().size();
            
            // auto v_Inac_rate = t_p()[Npar];
            
            auto v_curr_noise = t_p[Npar ];
            auto v_pink_noise = t_p[Npar + 1];
            auto v_prop_noise = t_p[Npar + 2];
            auto v_baseline = t_p[Npar + 3];
            auto v_N0 = t_p[std::pair{Npar + 4, Npar + 4}];
            a_g() = a_g() * -1.0;
            
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
        [names, m](const auto &patch_model)
        -> Maybe_error<
            Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>> {
            
            auto par= get_Parameters_from_Q0_Qa_g(m,names,get<Q0>(patch_model),get<Qa>(patch_model),build<g>(get<g>(patch_model)()*(-1.0)));
            
            if(!par)
                return par.error();
            auto Npar = names().size();
            
            auto out=Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>(Npar+5,1);
            
            out.set(std::pair(0ul,Npar-1),par.value());
            
            out[Npar]=get<Current_Noise>(patch_model)();
            out[Npar+1]=get<Pink_Noise>(patch_model)();
            out[Npar+2]=get<Proportional_Noise>(patch_model)();
            out[Npar+3]=get<Current_Baseline>(patch_model)();
            out.set(std::pair{Npar + 4, Npar + 4},get<N_Ch_mean>(patch_model)());
            return out;
            
            
        },
        p, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula,
        std::move(tr_param));
});

static auto scheme_7_d = add_Patch_inactivation_to_model<Allost1>(
    scheme_7, 1e-5, var::MyTranformations::from_string("Log10").value());


static auto scheme_8 = Allost1::Model("scheme_8", []() {
    auto v_binding = Conformational_change_label{"Binding"};
    auto v_rocking = Conformational_change_label{"Rocking"};
    
    auto mo = make_Conformational_model_standarized(
        Agonist_dependency_map{
                               std::map<Conformational_change_label, Agonist_dependency>{
                                                                      {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
                {v_rocking, Agonist_dependency{}}}},
        std::vector<Conformational_change_label>{v_binding, v_rocking, v_binding,
                                                 v_rocking, v_binding, v_rocking},
        std::vector<Conformational_interaction>{{Vector_Space{
                                                              Conformational_interaction_label{"RBR"},
            Conformational_interaction_players{{v_rocking, v_binding, v_rocking}},
            Conformational_interaction_positions{{{5, 0, 1},
                                               //   {1, 0, 5},
                                                  {1, 2, 3},
                                               //   {3, 2, 1},
                                                  {3, 4, 5}/*,
                                                  {5, 4, 3}*/}}}}},
        std::vector<Conductance_interaction>{Vector_Space{
                                                          Conductance_interaction_label{"Rocking_Current_factor"},
            Conductance_interaction_players{{v_rocking}},
            Conductance_interaction_positions{{{{1}}, {{3}}, {{5}}}}}},
        std::map<Conformational_change_label, Conformation_change_standard_state>{
                                                                                  {v_rocking,
             Conformation_change_standard_state{
                                                Conformational_interactions_domain_state{std::map<
                     Vector_Space<Conformational_interaction_index,
                                  Conformational_interaction_subposition>,
                     int>{
                          {Vector_Space{Conformational_interaction_index{0ul},
                                   Conformational_interaction_subposition{0}},
                      1},
                     {Vector_Space{Conformational_interaction_index{0ul},
                                   Conformational_interaction_subposition{2}},
                      1}}}}}},
        
        Conductance_interaction_info(Conductance_interaction_kind::equilibrium,"Leakeage_current_ratio"));
        
    assert(mo);
    auto m = std::move(mo.value());
    
    auto names = make_ModelNames<Allost1>(m);
    
    auto names_vec =
        std::vector<std::string>{
                                              "Binding_on",  "Binding_off", "Rocking_on",
        "Rocking_off", "RBR",         "RBR_0",
        "RBR_1",       "RBR_2",       "Rocking_Current_factor","Leakeage_current_ratio"}; //--> 8
    
    auto names_other = std::vector<std::string>{
                                                "Gating_Current",     "Current_Noise",    "Pink_Noise",
        "Proportional_Noise", "Current_Baseline", "Num_ch"};
    
    auto p_kinetics =
        std::vector<double>{10.0, 10000, 10000, 100, 10, 1, 1, 1, 10,1e-6};
    auto p_other = std::vector<double>{1, 1e-3, 5, 1e-2, 1, 4800};
    
    p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
    auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);
    
    auto tr = std::vector<std::string>(p.size(), "Log10");
    tr[tr.size() - 2] = "Linear";
    assert(tr.size() == p.size());
    auto tr_param = var::MyTranformations::from_strings(tr).value();
    
    assert(names() == names_vec);
    
    names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());
    
    auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
    assert(Maybe_modeltyple_formula);
    auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
        std::move(Maybe_modeltyple_formula.value());
    return std::tuple(
        [names, m](const auto &t_p)
        -> Maybe_error<
                Transfer_Op_to<std::decay_t<decltype(t_p)>, Patch_Model>> {
            auto Maybe_Q0Qag = make_Model<Allost1>(m, names, t_p);
            // std::cerr<<"parameters\n"<<p();
            
            assert(Maybe_Q0Qag);
            auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
            
            auto Npar = names().size();
            
            // auto v_Inac_rate = t_p()[Npar];
            auto v_unitary_current = t_p[Npar] * -1.0;
            
            auto v_curr_noise = t_p[Npar + 1];
            auto v_pink_noise = t_p[Npar + 2];
            auto v_prop_noise = t_p[Npar + 3];
            auto v_baseline = t_p[Npar + 4];
            auto v_N0 = t_p[std::pair{Npar + 5, Npar + 5}];
            a_g() = a_g() /var::max(a_g())* v_unitary_current;
            
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
        [names, m](const auto &patch_model)
        -> Maybe_error<
            Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>> {
            
            auto v_g=get<g>(patch_model);
            auto v_unitary_current=var::min(v_g());
            v_g()=v_g()/v_unitary_current;
            
            auto par= get_Parameters_from_Q0_Qa_g(m,names,get<Q0>(patch_model),get<Qa>(patch_model),v_g);
            
            if(!par)
                return par.error();
            auto Npar = names().size();
            
            auto out=Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>(Npar+6 ,1);
            
            out.set(std::pair(0ul,Npar-1),par.value());
            
            out[Npar]=v_unitary_current* (-1.0);
            out[Npar+1]=get<Current_Noise>(patch_model)();
            out[Npar+2]=get<Pink_Noise>(patch_model)();
            out[Npar+3]=get<Proportional_Noise>(patch_model)();
            out[Npar+4]=get<Current_Baseline>(patch_model)();
            out.set(std::pair{Npar + 5, Npar + 5},get<N_Ch_mean>(patch_model)());
            return out;
            
            
        },
        p, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula,
        std::move(tr_param));
});

static auto scheme_8_d = add_Patch_inactivation_to_model<Allost1>(
    scheme_8, 1e-5, var::MyTranformations::from_string("Log10").value());

static auto scheme_9 = Allost1::Model("scheme_9", []() {
    auto v_binding = Conformational_change_label{"Binding"};
    auto v_rocking = Conformational_change_label{"Rocking"};
    
    auto mo = make_Conformational_model_standarized(
        Agonist_dependency_map{
                               std::map<Conformational_change_label, Agonist_dependency>{
                                                                      {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
                {v_rocking, Agonist_dependency{}}}},
        std::vector<Conformational_change_label>{v_binding, v_rocking, v_binding,
                                                 v_rocking, v_binding, v_rocking},
        std::vector<Conformational_interaction>{{Vector_Space{
                                                              Conformational_interaction_label{"RB"},
            Conformational_interaction_players{{v_rocking, v_binding}},
            Conformational_interaction_positions{
                                                 {{5, 0}, {1, 0}, {1, 2}, {3, 2}, {3, 4}, {5, 4}}}}}},
        std::vector<Conductance_interaction>{Vector_Space{
                                                          Conductance_interaction_label{"Rocking_Current_factor"},
            Conductance_interaction_players{{v_rocking}},
            Conductance_interaction_positions{{{{1}}, {{3}}, {{5}}}}}},
        std::map<Conformational_change_label, Conformation_change_standard_state>{
                                                                                  {v_rocking,
             Conformation_change_standard_state{
                                                Conformational_interactions_domain_state{std::map<
                     Vector_Space<Conformational_interaction_index,
                                  Conformational_interaction_subposition>,
                     int>{
                          {Vector_Space{Conformational_interaction_index{0ul},
                                   Conformational_interaction_subposition{0}},
                      2}}}}}},
        
        Conductance_interaction_info(Conductance_interaction_kind::equilibrium,"Leakeage_current_ratio"));
        
    assert(mo);
    auto m = std::move(mo.value());
    
    auto names = make_ModelNames<Allost1>(m);
    
    auto names_vec =
        std::vector<std::string>{"Binding_on", "Binding_off",
                                              "Rocking_on", "Rocking_off",
                                              "RB",         "RB_0",
                                              "RB_1",       "Rocking_Current_factor","Leakeage_current_ratio"}; //--> 8
    
    auto names_other = std::vector<std::string>{
                                                "Gating_Current",     "Current_Noise",    "Pink_Noise",
        "Proportional_Noise", "Current_Baseline", "Num_ch"};
    
    auto p_kinetics =
        std::vector<double>{10.0, 10000, 10000, 100, 10, 1, 1, 10,1e-6};
    auto p_other = std::vector<double>{1, 1e-3, 5, 1e-2, 1, 4800};
    
    p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
    auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);
    
    auto tr = std::vector<std::string>(p.size(), "Log10");
    tr[tr.size() - 2] = "Linear";
    assert(tr.size() == p.size());
    auto tr_param = var::MyTranformations::from_strings(tr).value();
    
    assert(names() == names_vec);
    
    names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());
    
    auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
    assert(Maybe_modeltyple_formula);
    auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
        std::move(Maybe_modeltyple_formula.value());
    return std::tuple(
        [names, m](const auto &t_p)
        -> Maybe_error<
                Transfer_Op_to<std::decay_t<decltype(t_p)>, Patch_Model>> {
            auto Maybe_Q0Qag = make_Model<Allost1>(m, names, t_p);
            // std::cerr<<"parameters\n"<<p();
            
            assert(Maybe_Q0Qag);
            auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
            
            auto Npar = names().size();
            
            // auto v_Inac_rate = t_p()[Npar];
            auto v_unitary_current = t_p[Npar] * -1.0;
            
            auto v_curr_noise = t_p[Npar + 1];
            auto v_pink_noise = t_p[Npar + 2];
            auto v_prop_noise = t_p[Npar + 3];
            auto v_baseline = t_p[Npar + 4];
            auto v_N0 = t_p[std::pair{Npar + 5, Npar + 5}];
            a_g() = a_g() /var::max(a_g())* v_unitary_current;
            
            
            
            
            
            
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
        
        [names, m](const auto &patch_model)
        -> Maybe_error<
            Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>> {
            
            auto v_g=get<g>(patch_model);
            auto v_unitary_current=var::min(v_g());
            v_g()=v_g()/v_unitary_current;
            
            auto par= get_Parameters_from_Q0_Qa_g(m,names,get<Q0>(patch_model),get<Qa>(patch_model),v_g);
            
            if(!par)
                return par.error();
            auto Npar = names().size();
            
            auto out=Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>(Npar+6 ,1);
            
            out.set(std::pair(0ul,Npar-1),par.value());
            
            out[Npar]=v_unitary_current* (-1.0);
            out[Npar+1]=get<Current_Noise>(patch_model)();
            out[Npar+2]=get<Pink_Noise>(patch_model)();
            out[Npar+3]=get<Proportional_Noise>(patch_model)();
            out[Npar+4]=get<Current_Baseline>(patch_model)();
            out.set(std::pair{Npar + 5, Npar + 5},get<N_Ch_mean>(patch_model)());
            return out;
            
            
        },
        p, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula,
        std::move(tr_param));
});

static auto scheme_9_d = add_Patch_inactivation_to_model<Allost1>(
    scheme_9, 1e-5, var::MyTranformations::from_string("Log10").value());




inline auto get_model(std::string modelName) {
  auto allmodels =
      // Models_Library(&model00, &model00_7, &model01, &model4, &model4_g_lin,
      //                             &model6, &model6_no_inactivation,
      //                             &model6_Eff_no_inactivation, &model7,
      //                             &model8, &model9);
      Models_Library(&scheme_1  ,&scheme_2, &scheme_3, &scheme_4 , &scheme_1_d,
                       &scheme_2_d, &scheme_3_d, &scheme_4_d, &scheme_6 ,
                       &scheme_6_d , &scheme_7, &scheme_7_d, &scheme_8,
                       &scheme_8_d, &scheme_9, &scheme_9_d);
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
      Models_Library(&scheme_1 ,&scheme_2 ,&scheme_3, &scheme_4 , &scheme_1_d,
                     &scheme_2_d, &scheme_3_d, &scheme_4_d , &scheme_6 ,
                       &scheme_6_d, &scheme_7, &scheme_7_d, &scheme_8,
                       &scheme_8_d, &scheme_9, &scheme_9_d);
  return allmodels[modelName];
}

inline void print_model_Priors(double covar) {
  auto allmodels = // Models_Library(&scheme_1);
      Models_Library(&scheme_1 , &scheme_2 , &scheme_3, &scheme_4 , &scheme_1_d,
                     &scheme_2_d, &scheme_3_d, &scheme_4_d , &scheme_6 ,
                       &scheme_6_d , &scheme_7, &scheme_7_d, &scheme_8,
                       &scheme_8_d, &scheme_9, &scheme_9_d);
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
