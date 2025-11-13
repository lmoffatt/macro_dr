#pragma once

#include "models_MoffattHume_allosteric.h"

namespace macrodr{

   
    inline auto get_model(std::string modelName) {
    auto allmodels =
        // Models_Library(&model00, &model00_7, &model01, &model4, &model4_g_lin,
        //                             &model6, &model6_no_inactivation,
        //                             &model6_Eff_no_inactivation, &model7,
        //                             &model8, &model9);
        //Models_Library(&p2x2::scheme_1_d);
        Models_Library(&p2x2::scheme_1, /*&p2x2::scheme_2, &p2x2::scheme_3, &p2x2::scheme_4,&p2x2::scheme_5, &p2x2::scheme_6, &p2x2::scheme_7,
                        &p2x2::scheme_8, &p2x2::scheme_9, &p2x2::scheme_10, &p2x2::scheme_11, 
                       /* &p2x2::scheme_12, &p2x2::scheme_13,  &p2x2::scheme_14,  &p2x2::scheme_15,*/ /*&p2x2::scheme_1_d,
                       &p2x2::scheme_2_d, &p2x2::scheme_3_d, &p2x2::scheme_4_d, 
                       &p2x2::scheme_5_d, &p2x2::scheme_6_d,
                       &p2x2::scheme_7_d,  &p2x2::scheme_8_d,
                       &p2x2::scheme_9_d, */ &p2x2::scheme_10_d /*, 
                       &p2x2::scheme_11_d/*,  &p2x2::scheme_12_d,  &p2x2::scheme_13_d, &p2x2::scheme_14_d,  &p2x2::scheme_15_d*/);
    return allmodels[modelName];
}

inline Maybe_error<std::size_t> get_num_parameters(std::string model) {
    auto maybe_model = get_model(model);
    if (!maybe_model)
        return maybe_model.error();

    return std::visit(
        [&](auto model0ptr) { return model0ptr->parameters_transformations().names().size(); },
        maybe_model.value());
}

inline auto get_model_scheme(std::string modelName) {
    auto allmodels =  // Models_Library(&p2x2::scheme_1_d);
        Models_Library(
            &p2x2::scheme_1, /*&p2x2::scheme_2, &p2x2::scheme_3, &p2x2::scheme_4, &p2x2::scheme_5, &p2x2::scheme_6, &p2x2::scheme_7, &p2x2::scheme_8,
            &p2x2::scheme_9, &p2x2::scheme_10, &p2x2::scheme_11,
            /* &p2x2::scheme_12, &p2x2::scheme_13,  &p2x2::scheme_14,  &p2x2::scheme_15,*//* &p2x2::scheme_1_d, &p2x2::scheme_2_d,
            &p2x2::scheme_3_d, &p2x2::scheme_4_d, &p2x2::scheme_5_d, &p2x2::scheme_6_d, &p2x2::scheme_7_d, &p2x2::scheme_8_d,
            &p2x2::scheme_9_d, */&p2x2::scheme_10_d/*,
            &p2x2::scheme_11_d /*,  &p2x2::scheme_12_d,  &p2x2::scheme_13_d, &p2x2::scheme_14_d,  &p2x2::scheme_15_d*/);
    return allmodels[modelName];
}

inline void print_model_Priors(double covar) {
    auto allmodels =  //Models_Library(&p2x2::scheme_1_d);
        Models_Library(
            &p2x2::scheme_1, /*&p2x2::scheme_2, &p2x2::scheme_3, &p2x2::scheme_4, &p2x2::scheme_5, &p2x2::scheme_6, &p2x2::scheme_7, &p2x2::scheme_8,
            &p2x2::scheme_9, &p2x2::scheme_10, &p2x2::scheme_11,
            /* &p2x2::scheme_12, &p2x2::scheme_13,  &p2x2::scheme_14,  &p2x2::scheme_15,*//* &p2x2::scheme_1_d, &p2x2::scheme_2_d,
            &p2x2::scheme_3_d, &p2x2::scheme_4_d, &p2x2::scheme_5_d, &p2x2::scheme_6_d, &p2x2::scheme_7_d, &p2x2::scheme_8_d,
            &p2x2::scheme_9_d, */&p2x2::scheme_10_d/*,
            &p2x2::scheme_11_d /*,  &p2x2::scheme_12_d,  &p2x2::scheme_13_d, &p2x2::scheme_14_d,  &p2x2::scheme_15_d*/);
    //,
    //
    //                     &model6, &model6_no_inactivation,
    //                     &model6_Eff_no_inactivation, &model7, &model8,
    //                     &model9);

    std::apply(
        [&covar](auto... ptr_models) {
            (
                [&covar](auto modelp) {
                    auto& tr_par = modelp->parameters_transformations();
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



}// namespace macrodr
