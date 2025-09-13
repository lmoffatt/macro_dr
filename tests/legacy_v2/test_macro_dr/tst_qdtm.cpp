#include "../CLI_function_table.h"
#include "../catch2/catch.hpp"
#include "../models_MoffattHume_linear.h"
#include "../qmodel.h"
TEST_CASE("Test Qdtm ", "[Qdtm]") {
    auto& model0 = macrodr::scheme_4_d;
    using MyModel = typename std::decay_t<decltype(model0)>::my_Id;

    std::string parfile = "../macro_dr/models/scheme_4_inact_par.csv";

    std::string ModelName = model0.model_name();
    auto Maybe_param1 =
        var::load_Parameters<MyModel>(parfile, ",", model0.model_name(), model0.names());
    REQUIRE(Maybe_param1);
    auto param1 = std::move(Maybe_param1.value());

    auto Maybe_m = model0(param1.standard_parameter());

    REQUIRE(Maybe_m);
    auto m = std::move(Maybe_m.value());

    auto func_table = macrodr::cmd::get_function_Table_maker_St(
        "../macro_dr/test_macro_dr/examples/tst.txt", 0ul)();

    auto Qx = macrodr::Macro_DMR{}.calc_Qx(m, macrodr::ATP_concentration(0.0));

    auto Maybe_Qx_eig = macrodr::Macro_DMR{}.calc_eigen(Qx);
    REQUIRE(Maybe_Qx_eig);
    auto Qx_eig = std::move(Maybe_Qx_eig.value());

    double dt = 1e-3;
    auto Maybe_Qdt = macrodr::Macro_DMR{}.calc_Qdt_eig(func_table, m, Qx_eig,
                                                       macrodr::number_of_samples(10), dt);
    auto Maybe_Qdt_Taylor =
        macrodr::Macro_DMR{}.calc_Qdt_taylor(m, Qx, macrodr::number_of_samples(10), 1e-3, 10);
    REQUIRE(Maybe_Qdt);
    auto r_Qdt = std::move(Maybe_Qdt.value());

    auto Maybe_Qdtm = macrodr::Macro_DMR{}.calc_Qdtm_eig(func_table, m, Qx_eig,
                                                         macrodr::number_of_samples(10), dt);

    REQUIRE(Maybe_Qdtm);
    auto r_Qdtm = std::move(Maybe_Qdtm.value());

    auto r_Qn = macrodr::Macro_DMR{}.get_Qn(r_Qdt);

    r_Qn = macrodr::Macro_DMR{}.sum_Qn(std::move(r_Qn), r_Qn);

    auto r_Qdt2 = macrodr::Macro_DMR{}.Qn_to_Qdt(r_Qn);

    auto Maybe_Qdt2 = macrodr::Macro_DMR{}.calc_Qdt_eig(func_table, m, Qx_eig,
                                                        macrodr::number_of_samples(20), 2 * dt);
    auto Qdt22 = std::move(Maybe_Qdt2.value());
    CHECK(compare_contents(Qdt22, r_Qdt2, std::sqrt(eps), std::sqrt(eps), 10));

    auto Maybe_Qdtm_2 = macrodr::Macro_DMR{}.calc_Qdtm_eig(func_table, m, Qx_eig,
                                                           macrodr::number_of_samples(10), 2e-3);

    CHECK(compare_contents(r_Qdtm, extract<macrodr::Qdtm>(r_Qdt), eps * 100, 10));

    CHECK(compare_contents(r_Qdt, Maybe_Qdt_Taylor.value(), eps * 100, 10));

    REQUIRE(0 == 0);
}
