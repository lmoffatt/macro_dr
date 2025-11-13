#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "allosteric_models.h"
#include "general_algorithm_on_containers.h"
#include "maybe_error.h"
#include "parameters.h"
#include "parameters_distribution.h"
#include "qmodel.h"
#include "variables.h"
namespace macrodr {

template <class... Ts, class S>
std::variant<std::monostate, Ts*..., S*> operator||(std::variant<std::monostate, Ts*...> one,
                                                    std::variant<std::monostate, S*> next) {
    if (!std::holds_alternative<std::monostate>(one))
        return std::visit([](auto x) { return std::variant<std::monostate, Ts*..., S*>(x); }, one);
    else
        return std::visit([](auto x) { return std::variant<std::monostate, Ts*..., S*>(x); }, next);
}

template <class Id>
struct Model_Patch {
    template <class F, class Finv>
    class Model {
        std::string m_name;
        std::tuple<F, Finv, Matrix<double>, std::vector<std::string>, Q0_formula, Qa_formula,
                   g_formula, var::transformations_vector>
            m_f;
        var::Parameters_Transformations m_param;

       public:
        using my_Id = Id;
        static constexpr bool is_Model_Patch = true;
        template <class G>
        Model(std::string t_name, G&& t_g)
            : m_name{t_name},
              m_f{std::forward<G>(t_g)()},
              m_param(t_name, std::get<3>(m_f), std::get<7>(m_f), std::get<Matrix<double>>(m_f)) {}

        auto& model_name() const { return m_name; }

        std::variant<std::monostate, Model const*> operator[](const std::string& name) const {
            if (name == m_name)
                return this;
            else
                return std::monostate();
        }

        auto& names() const { return std::get<std::vector<std::string>>(m_f); }
        auto& parameters_transformations() const { return m_param; }
        auto number_of_states() const { return get_Q0_formula()().size(); }

        auto& get_Q0_formula() const { return std::get<Q0_formula>(m_f); }
        auto& get_Qa_formula() const { return std::get<Qa_formula>(m_f); }
        auto& get_g_formula() const { return std::get<g_formula>(m_f); }

        template <class P>
            requires std::is_same_v<var::untransformed_type_t<P>, Patch_Model>
        auto operator()(const P& t_p) const -> Maybe_error<var::Parameters_values> {
            auto p = std::invoke(std::get<Finv>(m_f), t_p);
            if (p)
                return var::Parameters_values(parameters_transformations(), p.value());
            else {
                std::cerr << p.error()();
                return p.error();
            }
        }

        template <class P>
            requires std::is_same_v<var::untransformed_type_t<P>, var::Parameters_values>
        auto operator()(const P& t_p) const -> std::invoke_result_t<F, P> {
            if (!isfinite(var::fullsum(t_p())))
                return error_message("nan parameter value");
            return std::invoke(std::get<F>(m_f), t_p);
            //     auto result=std::invoke(std::get<F>(m_f), t_p);
            //         assert((
            //           [&t_p, this,&result](){
            //  auto res= var::compare_contents(t_p,(*this)(result.value()).value());
            // if (!res)
            //          std::cerr<<res.error()();
            //  return res;}()));
            //   return std::move(result);
        }

        template <class P>
        friend void report_model(save_Parameter<P>& s, const Model& m) {
            std::ofstream f(s.fname + "_&resumodel.csv");
            f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
            f << "Model Name\t" << m.model_name() << "\n";
            f << "Parameters Names\n";
            f << m.names();
            f << "<<\n---------------------\n";
            f << "Q0_formula\n";
            f << m.get_Q0_formula();
            f << "<<\n---------------------\n";
            f << "Qa_formula\n";
            f << m.get_Qa_formula();
            f << "<<\n---------------------\n";
            f << "g_formula\n";
            f << m.get_g_formula();
            f << "<<\n---------------------\n";
        }
    };

    template <class F>
        requires(!is_of_this_template_type_v<F, Model>)
    Model(std::string, F&& f) -> Model<std::tuple_element_t<0, decltype(std::declval<F&&>()())>,
                                       std::tuple_element_t<1, decltype(std::declval<F&&>()())>>;
};

struct Model0 : public Model_Patch<Model0> {};
struct Model1 : public Model_Patch<Model1> {};

struct Allost1 : public Model_Patch<Allost1> {};


template <class Id, class F, class Finv>
auto add_Patch_inactivation_to_model(typename Model_Patch<Id>::template Model<F, Finv> const& model,
                                     double inactivation_value,
                                     std::unique_ptr<var::base_transformation> transformation) {
    return typename Model_Patch<Id>::Model(
        model.model_name() + "_inact", [model, inactivation_value, &transformation]() {
            auto names_model = model.names();
            auto names_other = std::vector<std::string>{"inactivation_rate"};
            names_model.insert(names_model.end(), names_other.begin(), names_other.end());

            auto v_Q0_formula = model.get_Q0_formula();
            std::size_t N = v_Q0_formula().size();
            auto v_Q0_inactivation_formula =
                insert_new_formula(v_Q0_formula, 0ul, N, "inactivation_rate");
            auto v_Qa_formula = model.get_Qa_formula();
            auto v_Qa_inactivation_formula = change_states_number(v_Qa_formula, N + 1);

            auto v_g_formula = model.get_g_formula();
            auto v_g_inactiavation_formula = change_states_number(v_g_formula, N + 1);

            auto tr_par = model.parameters_transformations();
            tr_par.push_back("inactivation_rate", transformation->clone(), inactivation_value);
            auto p_inactivation = tr_par.standard_values();
            auto tr_param = tr_par.transf();

            return std::tuple(
                [model](const auto& p)
                    -> Maybe_error<Transfer_Op_to<std::decay_t<decltype(p)>, Patch_Model>> {
                    auto n = p.size();
                    auto inactivation = p[n - 1];
                    auto p_no_inactivation =
                        Transfer_Op_to<std::decay_t<decltype(p)>, var::Parameters_values>(
                            model.parameters_transformations(), p[std::pair(0ul, n - 2)]);
                    auto mo = model(p_no_inactivation);
                    using std::pow;
                    if (!mo)
                        return mo.error();
                    else
                        return add_Patch_inactivation(std::move(mo.value()), inactivation);
                },
                [model](const auto& patch_model)
                    -> Maybe_error<
                        Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>> {
                    auto [patch_model_no_inactivation, inactivation_rate] =
                        remove_Patch_inactivation(patch_model);

                    auto p = model(patch_model_no_inactivation);
                    if (!p)
                        return p.error();
                    else {
                        auto n = p.value().size() + 1;
                        auto out =
                            Transfer_Op_to<std::decay_t<decltype(patch_model)>, Matrix<double>>(1,
                                                                                                n);
                        out.set(std::pair(0ul, n - 2), p.value()());
                        out[n - 1] = inactivation_rate;
                        return out;
                    }
                },
                p_inactivation, names_model, std::move(v_Q0_inactivation_formula),
                std::move(v_Qa_inactivation_formula), std::move(v_g_inactiavation_formula),
                std::move(tr_param));
        });
}


template <class... Ms>
class Models_Library {
    std::tuple<Ms*...> m_models;

   public:
    using Model_type = std::variant<Ms*...>;
    Models_Library(Ms*... m) : m_models{m...} {}

    auto operator()() const { return m_models; }
    Maybe_error<Model_type> operator[](const std::string& name) {
        auto out =
            std::apply([&name](auto... models) { return (... || (*models)[name]); }, m_models);
        return std::visit(
            [&name](auto x) -> Maybe_error<Model_type> {
                if constexpr (std::is_same_v<std::monostate, decltype(x)>)
                    return error_message(name + " is not a model");
                else

                    return Maybe_error<Model_type>(x);
            },
            out);
    }
};

template <class... Ms>
Models_Library(Ms*... m) -> Models_Library<Ms...>;

}  // namespace macrodr

