#ifndef MODELS_MOFFATTHUME_LINEAR_H
#define MODELS_MOFFATTHUME_LINEAR_H
//#include "allosteric_models.h"
#include "allosteric_models.h"
#include "maybe_error.h"
#include "parameters.h"
#include "parameters_distribution.h"
#include "qmodel.h"
#include "variables.h"
#include <cstddef>
#include <optional>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>
namespace macrodr {

template <class... Ts, class S>
std::variant<std::monostate, Ts *..., S *>
operator||(std::variant<std::monostate, Ts *...> one,
           std::variant<std::monostate, S *> next) {
  if (!std::holds_alternative<std::monostate>(one))
    return std::visit(
        [](auto x) { return std::variant<std::monostate, Ts *..., S *>(x); },
        one);
  else
    return std::visit(
        [](auto x) { return std::variant<std::monostate, Ts *..., S *>(x); },
        next);
}

template <class Id> struct Model_Patch {
  template <class F> class Model {
    std::string m_name;
    std::tuple<F, Matrix<double>, std::vector<std::string>, Q0_formula,
               Qa_formula, g_formula>
        m_f;
    Parameters<Id> m_par;

  public:
    using my_Id = Id;
    static constexpr bool is_Model_Patch = true;
    template <class G>
    Model(std::string t_name, G &&t_g)
        : m_name{t_name}, m_f{std::forward<G>(t_g)()},
          m_par{Parameters<Id>(m_name, std::get<std::vector<std::string>>(m_f),
                               std::get<Matrix<double>>(m_f))} {}

    Model(const Model &x)
        : m_name{x.m_name}, m_f{x.m_f},
          m_par{Parameters<Id>(m_name, std::get<std::vector<std::string>>(m_f),
                               std::get<Matrix<double>>(m_f))} {}
    Model() {}
    auto &model_name() const { return m_name; }

    std::variant<std::monostate, Model *> operator[](const std::string &name) {
      if (name == m_name)
        return this;
      else
        return std::monostate();
    }

    auto &names() const { return std::get<std::vector<std::string>>(m_f); }
    auto &parameters() const { return m_par; }

    auto &get_Q0_formula() const { return std::get<Q0_formula>(m_f); }
    auto &get_Qa_formula() const { return std::get<Qa_formula>(m_f); }
    auto &get_g_formula() const { return std::get<g_formula>(m_f); }

    template <class P>
      requires std::is_same_v<var::untransformed_type_t<P>, Parameters<Id>>
    auto operator()(const P &t_p) const {
      return std::invoke(std::get<F>(m_f), t_p);
    }

    template <class P>
    friend void report_model(save_Parameter<P> &s, const Model &m) {
      std::ofstream f(s.fname + "_model.csv");
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
  Model(std::string, F &&f)
      ->Model_Patch<Id>::Model<
          std::tuple_element_t<0, decltype(std::declval<F &&>()())>>;
  template <class F>
    requires(std::is_same_v<F, Model<F>>)
  Model(Model<F> &&f)->Model_Patch<Id>::Model<F>;
};

template <class Id, class F>
auto add_Patch_inactivation_to_model(
    typename Model_Patch<Id>::template Model<F> const &model,
    double inactivation_value) {
  return typename Model_Patch<Id>::Model(
      model.model_name() + "_inact", [model, inactivation_value]() {
        auto names_model = model.names();
        auto names_other = std::vector<std::string>{"inactivation_rate"};
        names_model.insert(names_model.end(), names_other.begin(),
                           names_other.end());

        auto v_Q0_formula = model.get_Q0_formula();
        std::size_t N = v_Q0_formula().size();
        auto v_Q0_inactivation_formula =
            insert_new_formula(v_Q0_formula, 0ul, N, "inactivation_rate");
        auto v_Qa_formula = model.get_Qa_formula();
        auto v_Qa_inactivation_formula =
            change_states_number(v_Qa_formula, N + 1);

        auto v_g_formula = model.get_g_formula();
        auto v_g_inactiavation_formula =
            change_states_number(v_g_formula, N + 1);

        auto logp = model.parameters()();

        Matrix<double> logp_inactivation(logp.size()+1, 1, 0.0);
        for (std::size_t i = 0; i < logp.size(); ++i)
          logp_inactivation[i] = logp[i];
        using std::log10;
        logp_inactivation[logp.size()] = log10(inactivation_value);

        return std::tuple(
            [model](const auto &logp)
                -> Maybe_error<
                    Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
              auto n = logp.size();
              auto &inactivation = logp[n - 1];
              auto mo=model(logp);
              if (!mo)
                  return mo.error();
              else
                  return add_Patch_inactivation(std::move(mo.value()), inactivation);
            },
            logp_inactivation, names_model,
            std::move(v_Q0_inactivation_formula),
            std::move(v_Qa_inactivation_formula),
            std::move(v_g_inactiavation_formula));
      });
}

struct Model0 : public Model_Patch<Model0> {};
struct Model1 : public Model_Patch<Model1> {};

struct Allost1 : public Model_Patch<Allost1> {};

class Model00_7;

static auto scheme_1 = Model0::Model("scheme_1", []() {
  auto names_model = std::vector<std::string>{"kon", "koff", "gating_on",
                                              "gating_off", "unitary_current"};
  auto names_other = std::vector<std::string>{
      "Current_Noise", "Current_Baseline", "Num_ch_mean"};

  std::size_t N = 5ul;

  auto v_Q0_formula = Q0_formula(N);
  v_Q0_formula()[1][0] = "koff";
  v_Q0_formula()[2][1] = "2*koff";
  v_Q0_formula()[3][2] = "3*koff";
  v_Q0_formula()[3][4] = "gating_on";
  v_Q0_formula()[4][3] = "gating_off";

  auto v_Qa_formula = Qa_formula(N);
  v_Qa_formula()[0][1] = "3*kon";
  v_Qa_formula()[1][2] = "2*kon";
  v_Qa_formula()[2][3] = "kon";
  auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
  v_g_formula()[4] = "unitary_current";

  names_model.insert(names_model.end(), names_other.begin(), names_other.end());
  auto p = Matrix<double>(
      8, 1, std::vector<double>{6.73, 166, 743, 45.3, 1, 1e-3, 1, 5000});

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  return std::tuple(
      [](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
        auto kon = p[0];
        auto koff = p[1];
        auto gating_on = p[2];
        auto gating_off = p[3];
        auto v_unitary_current = p[4] * -1.0;
        auto Npar = 5ul;
        auto v_curr_noise = p[Npar];
        auto v_baseline = logp()[Npar + 1];
        //  auto v_Num_ch_mean=p[Npar+2];
        //  auto v_std_log_Num_ch=p[Npar+3];

        auto v_N0 = p[std::pair(Npar + 2, Npar + 2)];

        auto N = N_St(5);
        return build<Patch_Model>(
            N_St(N()),
            build<Q0>(var::build_<Matrix<double>>(
                N(), N(), {{1, 0}, {2, 1}, {3, 2}, {3, 4}, {4, 3}},
                {koff, koff * 2.0, koff * 3.0, gating_on, gating_off})),
            build<Qa>(var::build_<Matrix<double>>(N(), N(),
                                                  {{0, 1}, {1, 2}, {2, 3}},
                                                  {kon * 3.0, kon * 2.0, kon})),
            build<P_initial>(
                var::build_<Matrix<double>>(1, N(), {{0, 0}}, {1.0})),
            build<g>(var::build_<Matrix<double>>(N(), 1, {{4, 0}},
                                                 {v_unitary_current})),
            build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
            build<Current_Baseline>(v_baseline),
            N_Ch_mean_time_segment_duration(121), Binomial_magical_number(5.0),
            min_P(1e-7), Probability_error_tolerance(1e-2),
            Conductance_variance_error_tolerance(1e-2));
      },
      logp, names_model, std::move(v_Q0_formula), std::move(v_Qa_formula),
      std::move(v_g_formula));
});


static auto scheme_2 = Model0::Model("scheme_2", []() {
  auto names_model = std::vector<std::string>{
      "kon",       "koff",       "flipping_on",    "flipping_off",
      "gating_on", "gating_off", "unitary_current"};
  auto names_other = std::vector<std::string>{
      "Current_Noise", "Current_Baseline", "Num_ch_mean"};

  std::size_t N = 6ul;

  auto v_Q0_formula = Q0_formula(N);
  v_Q0_formula()[1][0] = "koff";
  v_Q0_formula()[2][1] = "2*koff";
  v_Q0_formula()[3][2] = "3*koff";
  v_Q0_formula()[3][4] = "flipping_on";
  v_Q0_formula()[4][3] = "flipping_off";
  v_Q0_formula()[4][5] = "gating_on";
  v_Q0_formula()[5][4] = "gating_off";

  auto v_Qa_formula = Qa_formula(N);
  v_Qa_formula()[0][1] = "3*kon";
  v_Qa_formula()[1][2] = "2*kon";
  v_Qa_formula()[2][3] = "kon";
  auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
  v_g_formula()[5] = "unitary_current";

  names_model.insert(names_model.end(), names_other.begin(), names_other.end());
  auto p = Matrix<double>(
      10, 1,
      std::vector<double>{7.90, 121, 632, 42.1, 382, 1100, 1, 1e-3, 1, 5000});

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  return std::tuple(
      [N](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
        auto kon = p[0];
        auto koff = p[1];
        auto flipping_on = p[2];
        auto flipping_off = p[3];
        auto gating_on = p[4];
        auto gating_off = p[5];
        auto v_unitary_current = p[6] * -1.0;
        auto Npar = 7ul;
            auto v_curr_noise = p[Npar];
            auto v_baseline = logp()[Npar + 1];
            //  auto v_Num_ch_mean=p[Npar+2];
            //  auto v_std_log_Num_ch=p[Npar+3];
            
            auto v_N0 = p[std::pair(Npar + 2, Npar + 2)];
            
            return build<Patch_Model>(
                N_St(N),
                build<Q0>(var::build_<Matrix<double>>(
                    N, N, {{1, 0}, {2, 1}, {3, 2}, {3, 4}, {4, 3},{4,5},{5,4}},
                    {koff, koff * 2.0, koff * 3.0, flipping_on,flipping_off,gating_on, gating_off})),
                build<Qa>(var::build_<Matrix<double>>(
                    N, N, {{0, 1}, {1, 2}, {2, 3}}, {kon * 3.0, kon * 2.0, kon})),
                build<P_initial>(
                    var::build_<Matrix<double>>(1, N, {{0, 0}}, {1.0})),
                build<g>(var::build_<Matrix<double>>(N, 1, {{5, 0}},
                                                     {v_unitary_current})),
                build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
                build<Current_Baseline>(v_baseline),
                N_Ch_mean_time_segment_duration(121),
                Binomial_magical_number(5.0), min_P(1e-7),
                Probability_error_tolerance(1e-2),
                Conductance_variance_error_tolerance(1e-2));
      },
      logp, names_model, std::move(v_Q0_formula), std::move(v_Qa_formula),
      std::move(v_g_formula));
});

static auto scheme_3 = Model0::Model("scheme_3", []() {
    auto names_model = std::vector<std::string>{"kon_0",
                                                "koff_0",
                                                "kon_1",
                                                "koff_1",
                                                "kon_2",
                                                "koff_2",
                                                "gating_on_0",
                                                "gating_off_0",
                                                "gating_on_1",
                                                "gating_off_1",
                                                "desensitization_on_0",
                                                "desensitization_off_0",
                                                "desensitization_on_1",
                                                "unitary_current"};
    auto names_other = std::vector<std::string>{
                                                "Current_Noise", "Current_Baseline", "Num_ch_mean"};
    
    std::size_t N = 7ul;
    
    auto v_Q0_formula = Q0_formula(N);
    v_Q0_formula()[1][0] = "koff_0";
    v_Q0_formula()[2][1] = "koff_1";
    v_Q0_formula()[3][2] = "koff_2";
    v_Q0_formula()[3][4] = "gating_on_0";
    v_Q0_formula()[4][3] = "gating_off_0";
    v_Q0_formula()[3][5] = "gating_on_1";
    v_Q0_formula()[5][3] = "gating_off_1";
    
    v_Q0_formula()[4][6] = "desensitization_on_0";
    v_Q0_formula()[6][4] = "desensitization_off_0";
    v_Q0_formula()[5][6] = "desensitization_on_1";
    v_Q0_formula()[6][5] = "(desensitization_off_0 * gating_off_0 * gating_on_1 "
                           "* desensitization_on_1)/ (gating_off_1 * gating_on_0 "
                           "* desensitization_on_0)";
    
    auto v_Qa_formula = Qa_formula(N);
    v_Qa_formula()[0][1] = "kon_0";
    v_Qa_formula()[1][2] = "kon_1";
    v_Qa_formula()[2][3] = "kon_2";
    auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
    v_g_formula()[5] = "unitary_current";
    
    names_model.insert(names_model.end(), names_other.begin(), names_other.end());
    
    auto p_k_MH2007 =
        std::vector<double>{18, 0.017, 16.8,  175,   8.42,  4541, 
                            24.08,   29,  3198, 1189, 0.298, 1936, 5.9e7};
    auto p_other = std::vector<double>{1, 1e-3, 1, 5000};
    
    p_k_MH2007.insert(p_k_MH2007.end(), p_other.begin(), p_other.end());
    auto p = Matrix<double>(p_k_MH2007.size(), 1, p_k_MH2007);
    
    auto logp = apply([](auto x) { return std::log10(x); }, p);
    
    return std::tuple(
        [N](const auto &logp)
        -> Maybe_error<
            Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
            using std::pow;
            auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
            auto kon_0 = p[0];
            auto koff_0 = p[1];
            auto kon_1 = p[2];
            auto koff_1 = p[3];
            auto kon_2 = p[4];
            auto koff_2 = p[5];
            auto gating_on_0 = p[6];
            auto gating_off_0 = p[7];
            auto gating_on_1 = p[8];
            auto gating_off_1 = p[9];
            auto desensitization_on_0 = p[10];
            auto desensitization_off_0 = p[11];
            auto desensitization_on_1 = p[12];
            auto desensitization_off_1 =
                (desensitization_off_0 * gating_off_0 * gating_on_1 *
                 desensitization_on_1) /
                (gating_off_1 * gating_on_0 * desensitization_on_0);
            auto v_unitary_current = p[13] * -1.0;
            auto Npar = 14ul;
            auto v_curr_noise = p[Npar];
            auto v_baseline = logp()[Npar + 1];
            //  auto v_Num_ch_mean=p[Npar+2];
            //  auto v_std_log_Num_ch=p[Npar+3];
            
            auto v_N0 = p[std::pair(Npar + 2, Npar + 2)];
            
            return build<Patch_Model>(
                N_St(N),
                build<Q0>(var::build_<Matrix<double>>(
                    N, N, {{1, 0},
                     {2, 1},
                     {3, 2},
                     {3, 4},
                     {4, 3},
                     {3, 5},
                     {5, 3},
                     {4, 6},
                     {6, 4},
                     {5, 6},
                     {6, 5}},
                    {koff_0, koff_1, koff_2,  gating_on_0, gating_off_0,gating_on_1, gating_off_1,desensitization_on_0,desensitization_off_0,desensitization_on_1,desensitization_off_1})),
                build<Qa>(var::build_<Matrix<double>>(
                    N, N, {{0, 1}, {1, 2}, {2, 3}}, {kon_0, kon_1, kon_2})),
                build<P_initial>(
                    var::build_<Matrix<double>>(1, N, {{0, 0}}, {1.0})),
                build<g>(var::build_<Matrix<double>>(N, 1, {{4, 0}, {5, 0}},
                                                     {v_unitary_current, v_unitary_current})),
                build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
                build<Current_Baseline>(v_baseline),
                N_Ch_mean_time_segment_duration(121), Binomial_magical_number(5.0),
                min_P(1e-7), Probability_error_tolerance(1e-2),
                Conductance_variance_error_tolerance(1e-2));
        },
        logp, names_model, std::move(v_Q0_formula), std::move(v_Qa_formula),
        std::move(v_g_formula));
});


static auto scheme_4 = Model0::Model("scheme_4", []() {
  auto names_model = std::vector<std::string>{"kon_0",
                                                "koff_0",
                                                "kon_1",
                                                "koff_1",
                                                "kon_2",
                                                "koff_2",
                                                "flipping_on",
                                              "flipping_off",
                                              "gating_on_0",
                                              "gating_off_0",
                                              "gating_on_1",
                                              "gating_off_1",
                                              "desensitization_on_0",
                                              "desensitization_off_0",
                                              "desensitization_on_1",
                                              "unitary_current"};
  auto names_other = std::vector<std::string>{
      "Current_Noise", "Current_Baseline", "Num_ch_mean"};

  std::size_t N = 8ul;

  auto v_Q0_formula = Q0_formula(N);
  v_Q0_formula()[1][0] = "koff_0";
  v_Q0_formula()[2][1] = "koff_1";
  v_Q0_formula()[3][2] = "koff_2";
  v_Q0_formula()[3][4] = "flipping_on";
  v_Q0_formula()[4][3] = "flipping_off";
  v_Q0_formula()[4][5] = "gating_on_0";
  v_Q0_formula()[5][4] = "gating_off_0";
  v_Q0_formula()[4][6] = "gating_on_1";
  v_Q0_formula()[6][4] = "gating_off_1";

  v_Q0_formula()[5][7] = "desensitization_on_0";
  v_Q0_formula()[7][5] = "desensitization_off_0";
  v_Q0_formula()[6][7] = "desensitization_on_1";
  v_Q0_formula()[7][6] = "(desensitization_off_0 * gating_off_0 * gating_on_1 "
                         "* desensitization_on_1)/ (gating_off_1 * gating_on_0 "
                         "* desensitization_on_0)";

  auto v_Qa_formula = Qa_formula(N);
  v_Qa_formula()[0][1] = "kon_0";
  v_Qa_formula()[1][2] = "kon_1";
  v_Qa_formula()[2][3] = "kon_2";
  auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
  v_g_formula()[5] = "unitary_current";

  names_model.insert(names_model.end(), names_other.begin(), names_other.end());

  auto p_k_MH2007 =
      std::vector<double>{15.98, 0.019, 16.3,  380,   11.6,  6822, 3718, 43.54,
                          540,   1088,  0.033, 0.246, 31.16, 79.0, 4.53};
  auto p_other = std::vector<double>{1, 1e-3, 1, 5000};

  p_k_MH2007.insert(p_k_MH2007.end(), p_other.begin(), p_other.end());
  auto p = Matrix<double>(p_k_MH2007.size(), 1, p_k_MH2007);

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  return std::tuple(
      [N](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
        auto kon_0 = p[0];
        auto koff_0 = p[1];
        auto kon_1 = p[2];
        auto koff_1 = p[3];
        auto kon_2 = p[4];
        auto koff_2 = p[5];
        auto flipping_on = p[6];
        auto flipping_off = p[7];
        auto gating_on_0 = p[8];
        auto gating_off_0 = p[9];
        auto gating_on_1 = p[10];
        auto gating_off_1 = p[11];
        auto desensitization_on_0 = p[12];
        auto desensitization_off_0 = p[13];
        auto desensitization_on_1 = p[14];
        auto desensitization_off_1 =
            (desensitization_off_0 * gating_off_0 * gating_on_1 *
             desensitization_on_1) /
            (gating_off_1 * gating_on_0 * desensitization_on_0);
        auto v_unitary_current = p[15] * -1.0;
        auto Npar = 16ul;
        auto v_curr_noise = p[Npar];
        auto v_baseline = logp()[Npar + 1];
        //  auto v_Num_ch_mean=p[Npar+2];
        //  auto v_std_log_Num_ch=p[Npar+3];

        auto v_N0 = p[std::pair(Npar + 2, Npar + 2)];

        return build<Patch_Model>(
            N_St(N),
            build<Q0>(var::build_<Matrix<double>>(
                N, N, {{1, 0},
                       {2, 1},
                       {3, 2},
                       {3, 4},
                       {4, 3},
                       {4, 5},
                       {5, 4},
                       {4, 6},
                       {6, 4},
                       {5, 7},
                       {7, 5},
                       {6, 7},
                       {7, 6}},
                {koff_0, koff_1, koff_2, flipping_on, flipping_off,
                 gating_on_0, gating_off_0,gating_on_1, gating_off_1,desensitization_on_0,desensitization_off_0,desensitization_on_1,desensitization_off_1})),
            build<Qa>(var::build_<Matrix<double>>(
                N, N, {{0, 1}, {1, 2}, {2, 3}}, {kon_0, kon_1, kon_2})),
            build<P_initial>(
                var::build_<Matrix<double>>(1, N, {{0, 0}}, {1.0})),
            build<g>(var::build_<Matrix<double>>(N, 1, {{5, 0}, {6, 0}},
                                                 {v_unitary_current, v_unitary_current})),
            build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
            build<Current_Baseline>(v_baseline),
            N_Ch_mean_time_segment_duration(121), Binomial_magical_number(5.0),
            min_P(1e-7), Probability_error_tolerance(1e-2),
            Conductance_variance_error_tolerance(1e-2));
      },
      logp, names_model, std::move(v_Q0_formula), std::move(v_Qa_formula),
      std::move(v_g_formula));
});


static auto scheme_1_d = add_Patch_inactivation_to_model<Model0>(scheme_1, 1e-5);
static auto scheme_2_d = add_Patch_inactivation_to_model<Model0>(scheme_2, 1e-5);
static auto scheme_3_d = add_Patch_inactivation_to_model<Model0>(scheme_3, 1e-5);
static auto scheme_4_d = add_Patch_inactivation_to_model<Model0>(scheme_4, 1e-5);


static auto model00_7 = Model0::Model("model00_7", []() {
  auto names_model = std::vector<std::string>{"kon",
                                              "koff",
                                              "gatin_on",
                                              "gating_off",
                                              "inactivation_rate",
                                              "unitary_current"};
  auto names_other = std::vector<std::string>{
      "Current_Noise", "Current_Baseline", "Num_ch_mean"};

  std::size_t N = 6ul;

  auto v_Q0_formula = Q0_formula(N);
  v_Q0_formula()[1][0] = "koff";
  v_Q0_formula()[2][1] = "2*koff";
  v_Q0_formula()[3][2] = "3*koff";
  v_Q0_formula()[3][4] = "gatin_on";
  v_Q0_formula()[4][3] = "gating_off";
  v_Q0_formula()[0][5] = "inactivation_rate";

  auto v_Qa_formula = Qa_formula(N);
  v_Qa_formula()[0][1] = "3*kon";
  v_Qa_formula()[1][2] = "2*kon";
  v_Qa_formula()[2][3] = "kon";
  auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
  v_g_formula()[4] = "unitary_current";

  names_model.insert(names_model.end(), names_other.begin(), names_other.end());
  auto p = Matrix<double>(
      9, 1, std::vector<double>{10, 200, 1500, 50, 1e-5, 1, 1e-3, 1, 5000});

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  return std::tuple(
      [](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
        auto kon = p[0];
        auto koff = p[1];
        auto gating_on = p[2];
        auto gating_off = p[3];
        auto inactivation_rate = p[4];
        auto v_unitary_current = p[5] * -1.0;
        auto Npar = 6ul;
        auto v_curr_noise = p[Npar];
        auto v_baseline = logp()[Npar + 1];
        //  auto v_Num_ch_mean=p[Npar+2];
        //  auto v_std_log_Num_ch=p[Npar+3];

        auto v_N0 = p[std::pair(Npar + 2, Npar + 2)];

        return build<Patch_Model>(
            N_St(6),
            build<Q0>(var::build_<Matrix<double>>(
                6, 6, {{0, 5}, {1, 0}, {2, 1}, {3, 2}, {3, 4}, {4, 3}},
                {inactivation_rate, koff, koff * 2.0, koff * 3.0, gating_on,
                 gating_off})),
            build<Qa>(var::build_<Matrix<double>>(
                6, 6, {{0, 1}, {1, 2}, {2, 3}}, {kon * 3.0, kon * 2.0, kon})),
            build<P_initial>(
                var::build_<Matrix<double>>(1, 6, {{0, 0}}, {1.0})),
            build<g>(var::build_<Matrix<double>>(6, 1, {{4, 0}},
                                                 {v_unitary_current})),
            build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
            build<Current_Baseline>(v_baseline),
            N_Ch_mean_time_segment_duration(12100),
            Binomial_magical_number(5.0), min_P(1e-7),
            Probability_error_tolerance(1e-2),
            Conductance_variance_error_tolerance(1e-2));
      },
      logp, names_model, std::move(v_Q0_formula), std::move(v_Qa_formula),
      std::move(v_g_formula));
});





static auto prior_model00_7 = Custom_Distribution(
    9ul,
    [](mt_64i &mt) {
      double kon = 10.0;
      double koff = 100;
      double gating_on = 100;
      double gating_off = 100;
      double inactivation_rate = 1e-4;
      double unitary_current = 0.5;
      double current_noise = 1e-4;
      double current_baseline = 1; // zero after log10
      double Num_ch_mean = 1000;
      double global_std_log10_par = 2.0;

      auto parv = std::vector<double>{kon,
                                      koff,
                                      gating_on,
                                      gating_off,
                                      inactivation_rate,
                                      unitary_current,
                                      current_noise,
                                      current_baseline,
                                      Num_ch_mean};

      auto n_parv = parv.size();

      auto par = Matrix<double>(parv.size(), 1ul, parv);
      auto logpar = apply([](auto x) { return std::log10(x); }, par);

      auto s_logpar = std::vector<double>(n_parv, global_std_log10_par);

      auto covPar = DiagPosDetMatrix<double>(s_logpar);

      auto distPar = make_multivariate_normal_distribution(logpar, covPar);

      auto sample_par = distPar.value()(mt);

      // now recalculate the N_ch_numbers according to segements

      return Parameters<Model0>(model00_7.model_name(), model00_7.names(),
                                sample_par);
    },
    [](Parameters<Model0> const &logp) {
      double kon = 10.0;
      double koff = 100;
      double gating_on = 100;
      double gating_off = 100;
      double inactivation_rate = 1e-4;
      double unitary_current = 0.5;
      double current_noise = 1e-3;
      double current_baseline = 1; // zero after log10
      double Num_ch_mean = 4000;
      double global_std_log10_par = 2.0;

      auto parv = std::vector<double>{kon,
                                      koff,
                                      gating_on,
                                      gating_off,
                                      inactivation_rate,
                                      unitary_current,
                                      current_noise,
                                      current_baseline,
                                      Num_ch_mean};

      auto n_parv = parv.size();
      auto par = Matrix<double>(parv.size(), 1ul, parv);
      auto logpar = apply([](auto x) { return std::log10(x); }, par);

      auto s_logpar = std::vector<double>(n_parv, global_std_log10_par);

      auto covPar = DiagPosDetMatrix<double>(s_logpar);

      auto distPar = make_multivariate_normal_distribution(logpar, covPar);

      return distPar.value().logP(logp());
    });

static auto model00 = Model0::Model("model00", []() {
  auto names_model = std::vector<std::string>{"kon",
                                              "koff",
                                              "gatin_on",
                                              "gating_off",
                                              "inactivation_rate",
                                              "unitary_current"};
  auto names_other = std::vector<std::string>{
      "Current_Noise", "Current_Baseline", // "Num_ch_mean", "Num_ch_stddev",
      "Num_ch_0",      "Num_ch_1",         "Num_ch_2", "Num_ch_3",
      "Num_ch_4",      "Num_ch_5",         "Num_ch_6"};

  std::size_t N = 6ul;

  auto v_Q0_formula = Q0_formula(std::vector<std::vector<std::string>>(
      N, std::vector<std::string>(N, "")));
  v_Q0_formula()[1][0] = "koff";
  v_Q0_formula()[2][1] = "2*koff";
  v_Q0_formula()[3][2] = "3*koff";
  v_Q0_formula()[3][4] = "gatin_on";
  v_Q0_formula()[4][3] = "gating_off";
  v_Q0_formula()[0][5] = "inactivation_rate";

  auto v_Qa_formula = Qa_formula(N);
  v_Qa_formula()[0][1] = "3*kon";
  v_Qa_formula()[1][2] = "2*kon";
  v_Qa_formula()[2][3] = "kon";
  auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
  v_g_formula()[4] = "unitary_current";

  names_model.insert(names_model.end(), names_other.begin(), names_other.end());
  auto p = Matrix<double>(15, 1,
                          std::vector<double>{10, 200, 1500, 50, 1e-5, 1, 1e-3,
                                              1, 5000, 5000, 5000, 5000, 5000,
                                              5000, 5000});

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  return std::tuple(
      [](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
        auto kon = p[0];
        auto koff = p[1];
        auto gating_on = p[2];
        auto gating_off = p[3];
        auto inactivation_rate = p[4];
        auto v_unitary_current = p[5] * -1.0;
        auto Npar = 6ul;
        auto v_curr_noise = p[Npar];
        auto v_baseline = logp()[Npar + 1];
        //  auto v_Num_ch_mean=p[Npar+2];
        //  auto v_std_log_Num_ch=p[Npar+3];

        auto Npar2 = Npar + 1;
        auto v_N0 = p[std::pair(Npar2 + 1, Npar2 + 7)];
        return build<Patch_Model>(
            N_St(6),
            build<Q0>(var::build_<Matrix<double>>(
                6, 6, {{0, 5}, {1, 0}, {2, 1}, {3, 2}, {3, 4}, {4, 3}},
                {inactivation_rate, koff, koff * 2.0, koff * 3.0, gating_on,
                 gating_off})),
            build<Qa>(var::build_<Matrix<double>>(
                6, 6, {{0, 1}, {1, 2}, {2, 3}}, {kon * 3.0, kon * 2.0, kon})),
            build<P_initial>(
                var::build_<Matrix<double>>(1, 6, {{0, 0}}, {1.0})),
            build<g>(var::build_<Matrix<double>>(6, 1, {{4, 0}},
                                                 {v_unitary_current})),
            build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
            build<Current_Baseline>(v_baseline),
            N_Ch_mean_time_segment_duration(12100),
            Binomial_magical_number(5.0), min_P(1e-7),
            Probability_error_tolerance(1e-2),
            Conductance_variance_error_tolerance(1e-2));
      },
      logp, names_model, std::move(v_Q0_formula), std::move(v_Qa_formula),
      std::move(v_g_formula));
});

static auto prior_model00 = Custom_Distribution(
    17ul,
    [](mt_64i &mt) {
      double kon = 10.0;
      double koff = 100;
      double gating_on = 100;
      double gating_off = 100;
      double inactivation_rate = 1e-4;
      double unitary_current = 0.5;
      double current_noise = 1e-4;
      double current_baseline = 1; // zero after log10
      double Num_ch_mean = 1000;
      double std_log_Num_ch = 5e-2;
      double Num_ch_factor_n = 1.0;
      std::size_t N_ch_segments = 7;
      double global_std_log10_par = 2.0;

      auto parv = std::vector<double>{kon,
                                      koff,
                                      gating_on,
                                      gating_off,
                                      inactivation_rate,
                                      unitary_current,
                                      current_noise,
                                      current_baseline,
                                      Num_ch_mean,
                                      std_log_Num_ch};

      auto n_parv = parv.size();

      parv.insert(parv.end(), N_ch_segments, Num_ch_factor_n);

      auto par = Matrix<double>(parv.size(), 1ul, parv);
      auto logpar = apply([](auto x) { return std::log10(x); }, par);

      auto s_logpar = std::vector<double>(n_parv - 1, global_std_log10_par);
      s_logpar.insert(s_logpar.end(), N_ch_segments + 1, 1);

      auto covPar = DiagPosDetMatrix<double>(s_logpar);

      auto distPar = make_multivariate_normal_distribution(logpar, covPar);

      auto sample_par = distPar.value()(mt);

      // now recalculate the N_ch_numbers according to segements

      auto s_log_Num_ch_mean = sample_par[n_parv - 2];
      auto s_std_log_Num_ch_cv = std::pow(10.0, sample_par[n_parv - 1]);

      for (std::size_t i = n_parv; i < par.size(); ++i)
        sample_par[i] =
            sample_par[i] / global_std_log10_par * s_std_log_Num_ch_cv +
            s_log_Num_ch_mean;

      return Parameters<Model0>(model00.model_name(), model00.names(),
                                sample_par);
    },
    [](Parameters<Model0> const &logp) {
      double kon = 10.0;
      double koff = 100;
      double gating_on = 100;
      double gating_off = 100;
      double inactivation_rate = 1e-4;
      double unitary_current = 0.5;
      double current_noise = 1e-4;
      double current_baseline = 1; // zero after log10
      double Num_ch_mean = 1000;
      double std_logNum_ch = 5e-2;
      double Num_ch_factor_n = 1.0; //
      std::size_t N_ch_segments = 7;
      double global_std_log10_par = 2.0;

      auto parv = std::vector<double>{kon,
                                      koff,
                                      gating_on,
                                      gating_off,
                                      inactivation_rate,
                                      unitary_current,
                                      current_noise,
                                      current_baseline,
                                      Num_ch_mean,
                                      std_logNum_ch};

      auto n_parv = parv.size();
      auto p_Num_ch_mean = std::pow(10.0, logp()[n_parv - 2]);
      auto p_std_logNum_ch = std::pow(10.0, logp()[n_parv - 1]);

      parv.insert(parv.end(), N_ch_segments, p_Num_ch_mean);

      auto par = Matrix<double>(parv.size(), 1ul, parv);
      auto logpar = apply([](auto x) { return std::log10(x); }, par);

      auto s_logpar = std::vector<double>(n_parv, global_std_log10_par);
      s_logpar.insert(s_logpar.end(), N_ch_segments, p_std_logNum_ch);

      auto covPar = DiagPosDetMatrix<double>(s_logpar);

      auto distPar = make_multivariate_normal_distribution(logpar, covPar);

      return distPar.value().logP(logp());
    });

static auto model01 = Model0::Model("model01", []() {
  auto names_model = std::vector<std::string>{"kon",
                                              "koff",
                                              "gatin_on",
                                              "gating_off",
                                              "inactivation_rate",
                                              "unitary_current"};
  auto names_other =
      std::vector<std::string>{"Current_Noise", "Current_Baseline", "Num_ch"};

  std::size_t N = 6ul;

  auto v_Q0_formula = Q0_formula(std::vector<std::vector<std::string>>(
      N, std::vector<std::string>(N, "")));
  v_Q0_formula()[1][0] = "koff";
  v_Q0_formula()[2][1] = "2*koff";
  v_Q0_formula()[3][2] = "3*koff";
  v_Q0_formula()[3][4] = "gatin_on";
  v_Q0_formula()[4][3] = "gating_off";
  v_Q0_formula()[0][5] = "inactivation_rate";

  auto v_Qa_formula = Qa_formula(std::vector<std::vector<std::string>>(
      N, std::vector<std::string>(N, "")));
  v_Qa_formula()[0][1] = "3*kon";
  v_Qa_formula()[1][2] = "2*kon";
  v_Qa_formula()[2][3] = "kon";
  auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
  v_g_formula()[4] = "unitary_current";

  names_model.insert(names_model.end(), names_other.begin(), names_other.end());
  auto p = Matrix<double>(
      9, 1, std::vector<double>{10, 200, 1500, 50, 1e-5, 1, 1e-3, 1, 5000});

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  return std::tuple(
      [](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
        auto kon = p[0];
        auto koff = p[1];
        auto gating_on = p[2];
        auto gating_off = p[3];
        auto inactivation_rate = p[4];
        auto v_unitary_current = p[5] * -1.0;
        auto Npar = 6ul;
        auto v_curr_noise = p[Npar];
        auto v_baseline = logp()[Npar + 1];
        auto v_N0 = p[std::pair(Npar + 2, Npar + 2)];
        return build<Patch_Model>(
            N_St(6),
            build<Q0>(var::build_<Matrix<double>>(
                6, 6, {{0, 5}, {1, 0}, {2, 1}, {3, 2}, {3, 4}, {4, 3}},
                {inactivation_rate, koff, koff * 2.0, koff * 3.0, gating_on,
                 gating_off})),
            build<Qa>(var::build_<Matrix<double>>(
                6, 6, {{0, 1}, {1, 2}, {2, 3}}, {kon * 3.0, kon * 2.0, kon})),
            build<P_initial>(
                var::build_<Matrix<double>>(1, 6, {{0, 0}}, {1.0})),
            build<g>(var::build_<Matrix<double>>(6, 1, {{4, 0}},
                                                 {v_unitary_current})),
            build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
            build<Current_Baseline>(v_baseline),
            N_Ch_mean_time_segment_duration(12100),
            Binomial_magical_number(5.0), min_P(1e-7),
            Probability_error_tolerance(1e-2),
            Conductance_variance_error_tolerance(1e-2));
      },
      logp, names_model, std::move(v_Q0_formula), std::move(v_Qa_formula),
      std::move(v_g_formula));
});

static auto model4 = Model0::Model("model4", []() {
  auto names_model = std::vector<std::string>{"k01",
                                              "k10",
                                              "k12",
                                              "k21",
                                              "k23",
                                              "k32",
                                              "k34",
                                              "k43",
                                              "k45",
                                              "k54",
                                              "k46",
                                              "k64",
                                              "k57",
                                              "k75",
                                              "k67",
                                              "k08",
                                              "unitary_current"};
  auto names_other =
      std::vector<std::string>{"Current_Noise", "Current_Baseline", "Num_ch"};

  std::size_t N = 9ul;

  auto v_Q0_formula = Q0_formula(std::vector<std::vector<std::string>>(
      N, std::vector<std::string>(N, "")));
  v_Q0_formula()[0][8] = "k08";
  v_Q0_formula()[1][0] = "k10";
  v_Q0_formula()[2][1] = "k21";
  v_Q0_formula()[3][2] = "k32";
  v_Q0_formula()[3][4] = "k34";
  v_Q0_formula()[4][3] = "k43";
  v_Q0_formula()[4][5] = "k45";
  v_Q0_formula()[5][4] = "k54";
  v_Q0_formula()[4][6] = "k46";
  v_Q0_formula()[6][4] = "k64";
  v_Q0_formula()[5][7] = "k57";
  v_Q0_formula()[7][5] = "k75";
  v_Q0_formula()[6][7] = "k67";
  v_Q0_formula()[7][6] = "(k75 * k54 * k46 * k67)/ (k64 * k45 * k57)";
  auto v_Qa_formula = Qa_formula(std::vector<std::vector<std::string>>(
      N, std::vector<std::string>(N, "")));
  v_Qa_formula()[0][1] = "k01";
  v_Qa_formula()[1][2] = "k12";
  v_Qa_formula()[2][3] = "k23";
  auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
  v_g_formula()[5] = "unitary_current";
  v_g_formula()[6] = "unitary_current";

  names_model.insert(names_model.end(), names_other.begin(), names_other.end());

  auto p_kinetics = std::vector<double>(16, 100.0);
  p_kinetics[0] = 10;
  p_kinetics[2] = 10;
  p_kinetics[4] = 10;
  p_kinetics[12] = 100;
  p_kinetics[15] = 1e-3;

  auto p_k_MH2007 =
      std::vector<double>{15.98, 0.019, 16.3,  380,   11.6,  6822, 3718, 43.54,
                          540,   1088,  0.033, 0.246, 31.16, 79.0, 4.53, 1e-5};
  auto p_other = std::vector<double>{1, 1e-3, 1, 5000};

  p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
  p_k_MH2007.insert(p_k_MH2007.end(), p_other.begin(), p_other.end());

  //    auto p = Parameters<Model0>(p_kinetics);
  auto p = Matrix<double>(p_k_MH2007.size(), 1, p_k_MH2007);

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  return std::tuple(
      [](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
        auto Nst = 9ul;
        auto k01 = p[0];
        auto k10 = p[1];
        auto k12 = p[2];
        auto k21 = p[3];
        auto k23 = p[4];
        auto k32 = p[5];
        auto k34 = p[6];
        auto k43 = p[7];
        auto k45 = p[8];
        auto k54 = p[9];
        auto k46 = p[10];
        auto k64 = p[11];
        auto k57 = p[12];
        auto k75 = p[13];
        auto k67 = p[14];
        auto k76 = (k75 * k54 * k46 * k67) / (k64 * k45 * k57);
        auto k08 = p[15];
        auto v_g = p[16] * -1.0;
        auto Npar = 17ul;
        auto v_curr_noise = p[Npar];
        auto v_baseline = logp()[Npar + 1];
        auto v_N0 = p[std::pair(Npar + 2, Npar + 2)];

        return build<Patch_Model>(
            N_St(Nst),
            build<Q0>(var::build_<Matrix<double>>(Nst, Nst,
                                                  {{1, 0},
                                                   {2, 1},
                                                   {3, 2},
                                                   {3, 4},
                                                   {4, 3},
                                                   {4, 5},
                                                   {5, 4},
                                                   {4, 6},
                                                   {6, 4},
                                                   {5, 7},
                                                   {7, 5},
                                                   {6, 7},
                                                   {7, 6},
                                                   {0, 8}},
                                                  {k10, k21, k32, k34, k43, k45,
                                                   k54, k46, k64, k57, k75, k67,
                                                   k76, k08})),
            build<Qa>(var::build_<Matrix<double>>(
                Nst, Nst, {{0, 1}, {1, 2}, {2, 3}}, {k01, k12, k23})),
            build<P_initial>(
                var::build_<Matrix<double>>(1, Nst, {{0, 0}}, {1.0})),
            build<g>(var::build_<Matrix<double>>(Nst, 1, {{5, 0}, {6, 0}},
                                                 {v_g, v_g})),
            build<N_Ch_mean>(v_N0),

            build<Current_Noise>(v_curr_noise),
            build<Current_Baseline>(v_baseline),
            N_Ch_mean_time_segment_duration(121000),
            Binomial_magical_number(5.0), min_P(1e-14),
            Probability_error_tolerance(1e-2),
            Conductance_variance_error_tolerance(1e-2));
      },
      std::move(logp), std::move(names_model), v_Q0_formula, v_Qa_formula,
      v_g_formula);
});

static auto model4_g_lin = Model0::Model("model4_g_lin", []() {
  auto names_model = std::vector<std::string>{"k01",
                                              "k10",
                                              "k12",
                                              "k21",
                                              "k23",
                                              "k32",
                                              "k34",
                                              "k43",
                                              "k45",
                                              "k54",
                                              "k46",
                                              "k64",
                                              "k57",
                                              "k75",
                                              "k67",
                                              "k08",
                                              "unitary_current"};
  auto names_other =
      std::vector<std::string>{"Num_ch", "Current_Noise", "Current_Baseline"};

  std::size_t N = 9ul;

  auto v_Q0_formula = Q0_formula(std::vector<std::vector<std::string>>(
      N, std::vector<std::string>(N, "")));
  v_Q0_formula()[0][8] = "k08";
  v_Q0_formula()[1][0] = "k10";
  v_Q0_formula()[2][1] = "k21";
  v_Q0_formula()[3][2] = "k32";
  v_Q0_formula()[3][4] = "k34";
  v_Q0_formula()[4][3] = "k43";
  v_Q0_formula()[4][5] = "k45";
  v_Q0_formula()[5][4] = "k54";
  v_Q0_formula()[4][6] = "k46";
  v_Q0_formula()[6][4] = "k64";
  v_Q0_formula()[5][7] = "k57";
  v_Q0_formula()[7][5] = "k75";
  v_Q0_formula()[6][7] = "k67";
  v_Q0_formula()[7][6] = "(k75 * k54 * k46 * k67)/ (k64 * k45 * k57)";
  auto v_Qa_formula = Qa_formula(std::vector<std::vector<std::string>>(
      N, std::vector<std::string>(N, "")));
  v_Qa_formula()[0][1] = "k01";
  v_Qa_formula()[1][2] = "k12";
  v_Qa_formula()[2][3] = "k23";
  auto v_g_formula = g_formula(std::vector<std::string>(N, ""));
  v_g_formula()[5] = "unitary_current";
  v_g_formula()[6] = "unitary_current";

  names_model.insert(names_model.end(), names_other.begin(), names_other.end());

  auto p_kinetics = std::vector<double>(16, 100.0);
  p_kinetics[0] = 10;
  p_kinetics[2] = 10;
  p_kinetics[4] = 10;
  p_kinetics[12] = 100;
  p_kinetics[15] = 1e-3;

  auto p_other = std::vector<double>{1, 1e-3, 1, 5000};

  p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
  auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  return std::tuple(
      [](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto p = apply([](const auto &x) { return pow(10.0, x); }, logp());
        auto Nst = 9ul;
        auto k01 = p[0];
        auto k10 = p[1];
        auto k12 = p[2];
        auto k21 = p[3];
        auto k23 = p[4];
        auto k32 = p[5];
        auto k34 = p[6];
        auto k43 = p[7];
        auto k45 = p[8];
        auto k54 = p[9];
        auto k46 = p[10];
        auto k64 = p[11];
        auto k57 = p[12];
        auto k75 = p[13];
        auto k67 = p[14];
        auto k76 = (k75 * k54 * k46 * k67) / (k64 * k45 * k57);
        auto k08 = p[15];
        auto v_g = logp()[16];
        auto Npar = 17ul;
        auto v_curr_noise = p[Npar];
        auto v_baseline = logp()[Npar + 1];
        auto v_N0 = p[std::pair(Npar + 2, Npar + 2)];

        return build<Patch_Model>(
            N_St(Nst),
            build<Q0>(var::build_<Matrix<double>>(Nst, Nst,
                                                  {{1, 0},
                                                   {2, 1},
                                                   {3, 2},
                                                   {3, 4},
                                                   {4, 3},
                                                   {4, 5},
                                                   {5, 4},
                                                   {4, 6},
                                                   {6, 4},
                                                   {5, 7},
                                                   {7, 5},
                                                   {6, 7},
                                                   {7, 6},
                                                   {0, 8}},
                                                  {k10, k21, k32, k34, k43, k45,
                                                   k54, k46, k64, k57, k75, k67,
                                                   k76, k08})),
            build<Qa>(var::build_<Matrix<double>>(
                Nst, Nst, {{0, 1}, {1, 2}, {2, 3}}, {k01, k12, k23})),
            build<P_initial>(
                var::build_<Matrix<double>>(1, Nst, {{0, 0}}, {1.0})),
            build<g>(var::build_<Matrix<double>>(Nst, 1, {{5, 0}, {6, 0}},
                                                 {v_g, v_g})),
            build<N_Ch_mean>(v_N0),

            build<Current_Noise>(v_curr_noise),
            build<Current_Baseline>(v_baseline),
            N_Ch_mean_time_segment_duration(121000),
            Binomial_magical_number(5.0), min_P(1e-7),
            Probability_error_tolerance(1e-2),
            Conductance_variance_error_tolerance(1e-2));
      },
      std::move(logp), std::move(names_model), v_Q0_formula, v_Qa_formula,
      v_g_formula);
});

static auto model6 = Allost1::Model("model6", []() {
  auto v_binding = Conformational_change_label{"Binding"};
  auto v_rocking = Conformational_change_label{"Rocking"};
  auto v_gating = Conformational_change_label{"Gating"};

  auto mo = make_Conformational_model(
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
                       Conductance_interaction_positions{{{{4}}}}}});

  assert(mo);
  auto m = std::move(mo.value());

  auto names = make_ModelNames<Allost1>(m);

  auto names_vec = std::vector<std::string>{
      "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
      "Gating_on",  "Gating_off",  "BR",         "BR_0",
      "BR_1",       "BG",          "BG_0",       "BG_1",
      "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
  auto names_other = std::vector<std::string>{
      "Inactivation_rate", "Current_Noise", "Current_Baseline", "Num_ch"};

  auto p_kinetics = std::vector<double>{
      10, 10000, 100, 10000, 1, 10000, 10, 1, 1, 10, 1, 1, 10, 1, 1, 1};
  auto p_other = std::vector<double>{1e-3, 1e-3, 1, 5000};

  p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
  auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  assert(names() == names_vec);

  names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

  auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
  assert(Maybe_modeltyple_formula);
  auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
      std::move(Maybe_modeltyple_formula.value());
  return std::tuple(
      [names, m](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto p = build<Parameters<Allost1>>(
            logp.IdName(), logp.names(),
            apply([](const auto &x) { return pow(10.0, x); }, logp()));

        p()[names["BR_0"].value()] =
            p()[names["BR_0"].value()] / (1.0 + p()[names["BR_0"].value()]);
        p()[names["BR_1"].value()] =
            p()[names["BR_1"].value()] / (1.0 + p()[names["BR_1"].value()]);

        p()[names["BG_0"].value()] =
            p()[names["BG_0"].value()] / (1.0 + p()[names["BG_0"].value()]);
        p()[names["BG_1"].value()] =
            p()[names["BG_1"].value()] / (1.0 + p()[names["BG_1"].value()]);

        p()[names["RG_0"].value()] =
            p()[names["RG_0"].value()] / (1.0 + p()[names["RG_0"].value()]);
        p()[names["RG_1"].value()] =
            p()[names["RG_1"].value()] / (1.0 + p()[names["RG_1"].value()]);

        p()[names["Gating_Current"].value()] =
            p()[names["Gating_Current"].value()] * -1.0;
        auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
        assert(Maybe_Q0Qag);
        auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());

        auto Npar = names().size();

        auto v_Inac_rate = p()[Npar];
        auto v_curr_noise = p()[Npar + 1];
        auto v_baseline = logp()[Npar + 2];
        auto v_N0 = p()[std::pair{Npar + 3, Npar + 3}];
        auto Nst = get<N_St>(m());
        auto v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
            a_Q0, a_Qa, ATP_concentration(0.0), Nst);
        if (v_P_initial.valid())
          return add_Patch_inactivation(
              build<Patch_Model>(
                  N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
                  std::move(v_P_initial.value()), std::move(a_g),
                  build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
                  build<Current_Baseline>(v_baseline),
                  N_Ch_mean_time_segment_duration(120000),
                  Binomial_magical_number(5.0), min_P(1e-7),
                  Probability_error_tolerance(1e-2),
                  Conductance_variance_error_tolerance(1e-2)),
              v_Inac_rate);
        else
          return v_P_initial.error();
      },
      logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
});

static auto model6_no_inactivation =
    Allost1::Model("model6_no_inactivation", []() {
      auto v_binding = Conformational_change_label{"Binding"};
      auto v_rocking = Conformational_change_label{"Rocking"};
      auto v_gating = Conformational_change_label{"Gating"};

      auto mo = make_Conformational_model(
          Agonist_dependency_map{

              std::map<Conformational_change_label, Agonist_dependency>{
                  {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
                  {v_rocking, Agonist_dependency{}},
                  {v_gating, Agonist_dependency{}}}},
          std::vector<Conformational_change_label>{
              v_binding, v_binding, v_binding, v_rocking, v_gating},
          std::vector<Conformational_interaction>{
              {Vector_Space{
                   Conformational_interaction_label{"BR"},
                   Conformational_interaction_players{{v_binding, v_rocking}},
                   Conformational_interaction_positions{
                       {{0, 3}, {1, 3}, {2, 3}}}},
               Vector_Space{
                   Conformational_interaction_label{"BG"},
                   Conformational_interaction_players{{v_binding, v_gating}},
                   Conformational_interaction_positions{
                       {{0, 4}, {1, 4}, {2, 4}}}},
               Vector_Space{
                   Conformational_interaction_label{"RG"},
                   Conformational_interaction_players{{v_rocking, v_gating}},
                   Conformational_interaction_positions{{{3, 4}}}}

              }},
          std::vector<Conductance_interaction>{
              Vector_Space{Conductance_interaction_label{"Gating_Current"},
                           Conductance_interaction_players{{v_gating}},
                           Conductance_interaction_positions{{{{4}}}}}});

      assert(mo);
      auto m = std::move(mo.value());

      auto names = make_ModelNames<Allost1>(m);

      auto names_vec = std::vector<std::string>{
          "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
          "Gating_on",  "Gating_off",  "BR",         "BR_0",
          "BR_1",       "BG",          "BG_0",       "BG_1",
          "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
      auto names_other = std::vector<std::string>{"Current_Noise",
                                                  "Current_Baseline", "Num_ch"};

      auto p_kinetics = std::vector<double>{
          10, 1000, 1000, 100000, 1, 100, 100, 1, 1, 1, 1, 1, 100, 1, 1, 1};
      auto p_Moffatt_Hume_transformed = std::vector<double>{
          9.28,   1871,      2547.88,  295207, 0.220378,  150.312,
          74.865, 0.0323846, 0.187903, 1.77,   -0.457748, 1,
          123,    1,         1.3411,   1};
      auto p_other = std::vector<double>{1e-3, 1.0, 5000};

      p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
      auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

      auto logp = apply([](auto x) { return std::log10(x); }, p);

      assert(names() == names_vec);

      names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

      auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
      assert(Maybe_modeltyple_formula);
      auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
          std::move(Maybe_modeltyple_formula.value());
      return std::tuple(
          [names, m](const auto &logp)
              -> Maybe_error<
                  Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
            using std::pow;
            auto p = build<Parameters<Allost1>>(
                logp.IdName(), logp.names(),

                apply([](const auto &x) { return pow(10.0, x); }, logp()));

            p()[names["BR_0"].value()] =
                p()[names["BR_0"].value()] / (1.0 + p()[names["BR_0"].value()]);
            p()[names["BR_1"].value()] =
                p()[names["BR_1"].value()] / (1.0 + p()[names["BR_1"].value()]);

            p()[names["BG_0"].value()] =
                p()[names["BG_0"].value()] / (1.0 + p()[names["BG_0"].value()]);
            p()[names["BG_1"].value()] =
                p()[names["BG_1"].value()] / (1.0 + p()[names["BG_1"].value()]);

            p()[names["RG_0"].value()] =
                p()[names["RG_0"].value()] / (1.0 + p()[names["RG_0"].value()]);
            p()[names["RG_1"].value()] =
                p()[names["RG_1"].value()] / (1.0 + p()[names["RG_1"].value()]);

            p()[names["Gating_Current"].value()] =
                p()[names["Gating_Current"].value()] * -1.0;
            auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
            assert(Maybe_Q0Qag);
            auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
            auto Npar = names().size();

            // auto v_Inac_rate = p()[Npar];
            auto v_curr_noise = p()[Npar];
            auto v_baseline = logp()[Npar + 1];
            auto v_N0 = p()[std::pair{Npar + 2, Npar + 2}];
            auto Nst = get<N_St>(m());
            auto v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
                a_Q0, a_Qa, ATP_concentration(0.0), Nst);
            if (v_P_initial.valid())

              return build<Patch_Model>(
                  N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
                  std::move(v_P_initial.value()), std::move(a_g),
                  build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
                  build<Current_Baseline>(v_baseline),
                  N_Ch_mean_time_segment_duration(120000),
                  Binomial_magical_number(5.0), min_P(1e-7),
                  Probability_error_tolerance(1e-2),
                  Conductance_variance_error_tolerance(1e-2));
            else
              return v_P_initial.error();
          },
          logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
    });

static auto model6_Eff_no_inactivation =
    Allost1::Model("model6_Eff_no_inactivation", []() {
      auto v_binding = Conformational_change_label{"Binding"};
      auto v_rocking = Conformational_change_label{"Rocking"};
      auto v_gating = Conformational_change_label{"Gating"};

      auto mo = make_Conformational_model(
          Agonist_dependency_map{
              std::map<Conformational_change_label, Agonist_dependency>{
                  {v_binding, Agonist_dependency{Agonist_label{"ATP"}}},
                  {v_rocking, Agonist_dependency{}},
                  {v_gating, Agonist_dependency{}}}},
          std::vector<Conformational_change_label>{
              v_binding, v_binding, v_binding, v_rocking, v_gating},
          std::vector<Conformational_interaction>{
              {Vector_Space{
                   Conformational_interaction_label{"BR"},
                   Conformational_interaction_players{{v_binding, v_rocking}},
                   Conformational_interaction_positions{
                       {{0, 3}, {1, 3}, {2, 3}}}},
               Vector_Space{
                   Conformational_interaction_label{"BG"},
                   Conformational_interaction_players{{v_binding, v_gating}},
                   Conformational_interaction_positions{
                       {{0, 4}, {1, 4}, {2, 4}}}},
               Vector_Space{
                   Conformational_interaction_label{"RG"},
                   Conformational_interaction_players{{v_rocking, v_gating}},
                   Conformational_interaction_positions{{{3, 4}}}}

              }},
          std::vector<Conductance_interaction>{
              Vector_Space{Conductance_interaction_label{"Gating_Current"},
                           Conductance_interaction_players{{v_gating}},
                           Conductance_interaction_positions{{{{4}}}}}});

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

      auto names_other = std::vector<std::string>{"Current_Noise",
                                                  "Current_Baseline", "Num_ch"};

      auto p_kinetics = std::vector<double>{
          9.28, 1871, 3875, 1.07, 914, 776, 65.1 * 1.15, 1.15,
          33.3, 1.77, 0.77, 1.77, 123, 123, 635,         1};
      auto p_other = std::vector<double>{1e-3, 1, 4800};

      p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
      auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

      auto logp = apply([](auto x) { return std::log10(x); }, p);

      assert(names() == names_vec_untransformed);

      names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

      auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
      assert(Maybe_modeltyple_formula);
      auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
          std::move(Maybe_modeltyple_formula.value());
      return std::tuple(
          [names, m](const auto &logp)
              -> Maybe_error<
                  Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
            using std::pow;
            auto tr_p = build<Parameters<Allost1>>(
                logp.IdName(), logp.names(),
                apply([](const auto &x) { return pow(10.0, x); }, logp()));
            auto Binding_on = tr_p()[0];
            auto Binding_off = tr_p()[1];
            auto Rocking_on_B = tr_p()[2];
            auto Rocking_off_B = tr_p()[3];

            auto Gating_on_BR = tr_p()[4];
            auto Gating_off_BR = tr_p()[5];
            auto BR = tr_p()[6];
            auto BR_Bon = tr_p()[7];
            auto BR_Ron = tr_p()[8];
            auto BG = tr_p()[9];
            auto BG_Bon = tr_p()[10];
            auto BG_Gon = tr_p()[11];
            auto RG = tr_p()[12];
            auto RG_Ron = tr_p()[13];
            auto RG_Gon = tr_p()[14];
            auto Gating_Current = tr_p()[15] * (-1.0);

            auto Rocking_on = Rocking_on_B / pow(BR_Bon, 3);
            auto Rocking_off = Rocking_off_B * pow(BR / BR_Bon, 3);

            auto Gating_on = Gating_off_BR / pow(BG_Gon, 3) / RG_Gon;
            auto Gating_off = Gating_off_BR * pow(BG / BG_Gon, 3) * RG / RG_Gon;

            auto BR_0 = log(BR_Bon) / log(BR);
            auto BR_1 = log(BR_Ron) / log(BR);

            auto BG_0 = log(BG_Bon) / log(BG);
            auto BG_1 = log(BG_Gon) / log(BG);

            auto RG_0 = log(RG_Ron) / log(RG);
            auto RG_1 = log(RG_Gon) / log(RG);

            auto p = tr_p;
            p()[2] = Rocking_on;
            p()[3] = Rocking_off;

            p()[4] = Gating_on;
            p()[5] = Gating_off;

            p()[7] = BR_0;
            p()[8] = BR_1;

            p()[10] = BG_0;
            p()[11] = BG_1;

            p()[13] = RG_0;
            p()[14] = RG_1;
            p()[15] = Gating_Current;
            auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
            // std::cerr<<"parameters\n"<<p();

            assert(Maybe_Q0Qag);
            auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
            auto Npar = names().size();

            // auto v_Inac_rate = p()[Npar];
            auto v_curr_noise = tr_p()[Npar];
            auto v_baseline = logp()[Npar + 1];
            auto v_N0 = tr_p()[std::pair{Npar + 2, Npar + 2}];
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
                  build<Current_Baseline>(v_baseline),
                  N_Ch_mean_time_segment_duration(120000),
                  Binomial_magical_number(5.0), min_P(1e-7),
                  Probability_error_tolerance(1e-2),
                  Conductance_variance_error_tolerance(1e-2));
          },
          logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
    });

static auto model6_Eff_std = Allost1::Model("model6_Eff_std", []() {
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

  auto names_other =
      std::vector<std::string>{"Current_Noise", "Current_Baseline", "Num_ch"};

  auto p_kinetics =
      std::vector<double>{9.28, 1871, 3875, 1.07, 914, 776, 65.1 * 1.15, 1.15,
                          33.3, 1.77, 0.77, 1.77, 123, 123, 635,         1};
  auto p_other = std::vector<double>{1e-3, 1, 4800};

  p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
  auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  assert(names() == names_vec_untransformed);

  names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

  auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
  assert(Maybe_modeltyple_formula);
  auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
      std::move(Maybe_modeltyple_formula.value());
  return std::tuple(
      [names, m](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto tr_p = build<Parameters<Allost1>>(
            logp.IdName(), logp.names(),
            apply([](const auto &x) { return pow(10.0, x); }, logp()));
        auto Binding_on = tr_p()[0];
        auto Binding_off = tr_p()[1];
        auto Rocking_on_B = tr_p()[2];
        auto Rocking_off_B = tr_p()[3];

        auto Gating_on_BR = tr_p()[4];
        auto Gating_off_BR = tr_p()[5];
        auto BR = tr_p()[6];
        auto BR_Bon = tr_p()[7];
        auto BR_Ron = tr_p()[8];
        auto BG = tr_p()[9];
        auto BG_Bon = tr_p()[10];
        auto BG_Gon = tr_p()[11];
        auto RG = tr_p()[12];
        auto RG_Ron = tr_p()[13];
        auto RG_Gon = tr_p()[14];
        auto Gating_Current = tr_p()[15] * (-1.0);

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

        auto p = tr_p;
        p()[2] = Rocking_on;
        p()[3] = Rocking_off;

        p()[4] = Gating_on;
        p()[5] = Gating_off;

        p()[7] = BR_0;
        p()[8] = BR_1;

        p()[10] = BG_0;
        p()[11] = BG_1;

        p()[13] = RG_0;
        p()[14] = RG_1;
        p()[15] = Gating_Current;
        auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
        // std::cerr<<"parameters\n"<<p();

        assert(Maybe_Q0Qag);
        auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
        auto Npar = names().size();

        // auto v_Inac_rate = p()[Npar];
        auto v_curr_noise = tr_p()[Npar];
        auto v_baseline = logp()[Npar + 1];
        auto v_N0 = tr_p()[std::pair{Npar + 2, Npar + 2}];
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
              build<Current_Baseline>(v_baseline),
              N_Ch_mean_time_segment_duration(120000),
              Binomial_magical_number(5.0), min_P(1e-7),
              Probability_error_tolerance(1e-2),
              Conductance_variance_error_tolerance(1e-2));
      },
      logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
});

static auto model7 = Allost1::Model("model7", []() {
  auto v_rocking = Conformational_change_label{"Rocking"};
  auto v_binding = Conformational_change_label{"Binding"};
  auto v_gating = Conformational_change_label{"Gating"};

  auto mo = make_Conformational_model(
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
                       Conductance_interaction_positions{{{{6}}}}}});

  assert(mo);
  auto m = std::move(mo.value());

  auto names = make_ModelNames<Allost1>(m);

  auto names_vec = std::vector<std::string>{
      "Binding_on", "Binding_off", "Rocking_on", "Rocking_off",
      "Gating_on",  "Gating_off",  "BR",         "BR_0",
      "BR_1",       "BG",          "BG_0",       "BG_1",
      "RG",         "RG_0",        "RG_1",       "Gating_Current"}; //--> 8
  auto names_other = std::vector<std::string>{
      "Inactivation_rate", "Current_Noise", "Current_Baseline", "Num_ch"};

  auto p_kinetics = std::vector<double>{
      10, 10000, 100, 10000, 1, 10000, 10, 1, 1, 10, 1, 1, 10, 1, 1, 1};
  auto p_other = std::vector<double>{1, 1e-3, 1, 5000};

  p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());

  auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  assert(names() == names_vec);

  names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

  auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
  assert(Maybe_modeltyple_formula);
  auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
      std::move(Maybe_modeltyple_formula.value());
  return std::tuple(
      [names, m](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto p = build<Parameters<Allost1>>(
            logp.IdName(), logp.names(),

            apply([](const auto &x) { return pow(10.0, x); }, logp()));

        p()[names["BR_0"].value()] =
            p()[names["BR_0"].value()] / (1.0 + p()[names["BR_0"].value()]);
        p()[names["BR_1"].value()] =
            p()[names["BR_1"].value()] / (1.0 + p()[names["BR_1"].value()]);

        p()[names["BG_0"].value()] =
            p()[names["BG_0"].value()] / (1.0 + p()[names["BG_0"].value()]);
        p()[names["BG_1"].value()] =
            p()[names["BG_1"].value()] / (1.0 + p()[names["BG_1"].value()]);

        p()[names["RG_0"].value()] =
            p()[names["RG_0"].value()] / (1.0 + p()[names["RG_0"].value()]);
        p()[names["RG_1"].value()] =
            p()[names["RG_1"].value()] / (1.0 + p()[names["RG_1"].value()]);

        p()[names["Gating_Current"].value()] =
            p()[names["Gating_Current"].value()] * -1.0;
        auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
        assert(Maybe_Q0Qag);
        auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
        auto Npar = names().size();

        auto v_Inac_rate = p()[Npar];
        auto v_curr_noise = p()[Npar + 1];
        auto v_baseline = logp()[Npar + 2];
        auto v_N0 = p[std::pair(Npar + 3, Npar + 3)];
        auto Nst = get<N_St>(m());
        auto v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
            a_Q0, a_Qa, ATP_concentration(0.0), Nst);

        if (v_P_initial.valid())
          return add_Patch_inactivation(
              build<Patch_Model>(
                  N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
                  std::move(v_P_initial.value()), std::move(a_g),
                  build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
                  build<Current_Baseline>(v_baseline),
                  N_Ch_mean_time_segment_duration(120000),
                  Binomial_magical_number(5.0), min_P(1e-7),
                  Probability_error_tolerance(1e-2),
                  Conductance_variance_error_tolerance(1e-2)),
              v_Inac_rate);
        else
          return v_P_initial.error();
      },
      logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
});

static auto model8 = Allost1::Model("model8", []() {
  auto v_binding = Conformational_change_label{"Binding"};
  auto v_rocking = Conformational_change_label{"Rocking"};
  auto mo = make_Conformational_model(
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
                                                {1, 0, 5},
                                                {1, 2, 3},
                                                {3, 2, 1},
                                                {3, 4, 5},
                                                {5, 4, 3}}}}}},
      std::vector<Conductance_interaction>{Vector_Space{
          Conductance_interaction_label{"Rocking_Current_factor"},
          Conductance_interaction_players{{v_rocking}},
          Conductance_interaction_positions{{{{1}}, {{3}}, {{5}}}}}});

  assert(mo);
  auto m = std::move(mo.value());

  auto names = make_ModelNames<Allost1>(m);

  auto names_vec =
      std::vector<std::string>{
          "Binding_on",  "Binding_off", "Rocking_on",
          "Rocking_off", "RBR",         "RBR_0",
          "RBR_1",       "RBR_2",       "Rocking_Current_factor"}; //--> 8
  auto names_other =
      std::vector<std::string>{"Inactivation_rate", "Leaking_current",
                               "Current_Noise", "Current_Baseline", "Num_ch"};

  auto p_kinetics =
      std::vector<double>{10, 10000, 100, 10000, 100, 1.0, 1e-2, 1.0, 100};
  auto p_other = std::vector<double>{1e-4, 0.01, 1e-3, 1, 5000};

  p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
  auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  assert(names() == names_vec);

  names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

  auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
  assert(Maybe_modeltyple_formula);
  auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
      std::move(Maybe_modeltyple_formula.value());
  return std::tuple(
      [names, m](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto p = build<Parameters<Allost1>>(
            logp.IdName(), logp.names(),

            apply([](const auto &x) { return pow(10.0, x); }, logp()));
        p()[names["RBR_0"].value()] =
            p()[names["RBR_0"].value()] / (1.0 + p()[names["RBR_0"].value()]);
        p()[names["RBR_1"].value()] =
            p()[names["RBR_1"].value()] / (1.0 + p()[names["RBR_1"].value()]);
        p()[names["RBR_2"].value()] =
            p()[names["RBR_2"].value()] / (1.0 + p()[names["RBR_2"].value()]);

        auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
        assert(Maybe_Q0Qag);
        auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
        auto Npar = names().size();

        auto v_Inac_rate = p()[Npar];
        auto v_leaking_current = p()[Npar + 1];
        auto v_curr_noise = p()[Npar + 2];
        auto v_baseline = logp()[Npar + 3];
        auto v_N0 = p[std::pair(Npar + 4, Npar + 4)];

        auto v_g = build<g>(apply(
            [&v_leaking_current](const auto &x) {
              return v_leaking_current * pow(10.0, x) * (-1.0);
            },
            a_g()));
        auto Nst = get<N_St>(m());
        auto v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
            a_Q0, a_Qa, ATP_concentration(0.0), Nst);

        if (v_P_initial.valid())
          return add_Patch_inactivation(
              build<Patch_Model>(
                  N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
                  std::move(v_P_initial.value()), std::move(v_g),
                  build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
                  build<Current_Baseline>(v_baseline),
                  N_Ch_mean_time_segment_duration(120000),
                  Binomial_magical_number(5.0), min_P(1e-7),
                  Probability_error_tolerance(1e-2),
                  Conductance_variance_error_tolerance(1e-2)),
              v_Inac_rate);
        else
          return v_P_initial.error();
      },
      logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
});

static auto model9 = Allost1::Model("model9", []() {
  auto v_binding = Conformational_change_label{"Binding"};
  auto v_rocking = Conformational_change_label{"Rocking"};
  auto mo = make_Conformational_model(
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
          Conductance_interaction_positions{{{{1}}, {{3}}, {{5}}}}}});

  assert(mo);
  auto m = std::move(mo.value());

  auto names = make_ModelNames<Allost1>(m);

  auto names_vec =
      std::vector<std::string>{"Binding_on", "Binding_off",
                               "Rocking_on", "Rocking_off",
                               "RB",         "RB_0",
                               "RB_1",       "Rocking_Current_factor"}; //--> 8
  auto names_other =
      std::vector<std::string>{"Inactivation_rate", "leaking_current",
                               "Current_Noise", "Current_Baseline", "Num_ch"};

  auto p_kinetics =
      std::vector<double>{10, 10000, 100, 10000, 100, 1.0, 1e-2, 100};
  auto p_other = std::vector<double>{1e-3, 1e-2, 1e-3, 1, 5000};

  p_kinetics.insert(p_kinetics.end(), p_other.begin(), p_other.end());
  auto p = Matrix<double>(p_kinetics.size(), 1, p_kinetics);

  auto logp = apply([](auto x) { return std::log10(x); }, p);

  assert(names() == names_vec);

  names_vec.insert(names_vec.end(), names_other.begin(), names_other.end());

  auto Maybe_modeltyple_formula = make_Model_Formulas<Allost1>(m, names);
  assert(Maybe_modeltyple_formula);
  auto [a_Q0_formula, a_Qa_formula, a_g_formula] =
      std::move(Maybe_modeltyple_formula.value());
  return std::tuple(
      [names, m](const auto &logp)
          -> Maybe_error<
              Transfer_Op_to<std::decay_t<decltype(logp)>, Patch_Model>> {
        using std::pow;
        auto p = build<Parameters<Allost1>>(
            logp.IdName(), logp.names(),
            apply([](const auto &x) { return pow(10.0, x); }, logp()));
        p()[names["RB_0"].value()] =
            p()[names["RB_0"].value()] / (1.0 + p()[names["RB_0"].value()]);
        p()[names["RB_1"].value()] =
            p()[names["RB_1"].value()] / (1.0 + p()[names["RB_1"].value()]);

        auto Maybe_Q0Qag = make_Model<Allost1>(m, names, p);
        assert(Maybe_Q0Qag);
        auto [a_Q0, a_Qa, a_g] = std::move(Maybe_Q0Qag.value());
        auto Npar = names().size();

        auto v_Inac_rate = p()[Npar];
        auto v_leaking_current = p()[Npar + 1];
        auto v_curr_noise = p()[Npar + 2];
        auto v_baseline = logp()[Npar + 3];
        auto v_N0 = p()[std::pair(Npar + 4, Npar + 4)];

        auto v_g = build<g>(apply(
            [&v_leaking_current](const auto &x) {
              return v_leaking_current * pow(10.0, x) * (-1.0);
            },
            a_g()));
        auto Nst = get<N_St>(m());
        auto v_P_initial = macrodr::Macro_DMR{}.calc_Pinitial(
            a_Q0, a_Qa, ATP_concentration(0.0), Nst);
        if (v_P_initial.valid())

          return add_Patch_inactivation(
              build<Patch_Model>(
                  N_St(get<N_St>(m())), std::move(a_Q0), std::move(a_Qa),
                  std::move(v_P_initial.value()), std::move(v_g),
                  build<N_Ch_mean>(v_N0), build<Current_Noise>(v_curr_noise),
                  build<Current_Baseline>(v_baseline),
                  N_Ch_mean_time_segment_duration(120000),
                  Binomial_magical_number(5.0), min_P(1e-7),
                  Probability_error_tolerance(1e-2),
                  Conductance_variance_error_tolerance(1e-2)),
              v_Inac_rate);
        else
          return v_P_initial.error();
      },
      logp, names_vec, a_Q0_formula, a_Qa_formula, a_g_formula);
});

template <class... Ms> class Models_Library {
  std::tuple<Ms *...> m_models;

public:
  using Model_type = std::variant<Ms *...>;
  Models_Library(Ms *...m) : m_models{m...} {}

  auto operator()() const { return m_models; }
  Maybe_error<Model_type> operator[](const std::string &name) {
    auto out = std::apply(
        [&name](auto... models) { return (... || (*models)[name]); }, m_models);
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

template <class... Ms> Models_Library(Ms *...m) -> Models_Library<Ms...>;

inline auto get_model(std::string modelName) {
  auto allmodels =
      // Models_Library(&model00, &model00_7, &model01, &model4, &model4_g_lin,
      //                             &model6, &model6_no_inactivation,
      //                             &model6_Eff_no_inactivation, &model7,
      //                             &model8, &model9);
  Models_Library(&scheme_1,&scheme_2,&scheme_4_d);
//,&scheme_3,&scheme_4,&scheme_1_d,&scheme_2_d,&scheme_3_d,&scheme_4_d,
 //                      &model00, &model00_7, &model01, &model4, &model4_g_lin,
 //                    &model6, &model6_no_inactivation,
 //                    &model6_Eff_no_inactivation, &model6_Eff_std, &model7,
 //                    &model8, &model9);
        return allmodels[modelName];
}


inline Maybe_error<std::size_t> get_num_parameters(std::string model)
{
    auto maybe_model = get_model(model);
    if (!maybe_model)
        return maybe_model.error();
    
    return    std::visit([&](auto model0ptr) { return model0ptr->parameters().size(); },
                   maybe_model.value());
    
}


inline auto get_model_scheme(std::string modelName) {
    auto allmodels =//Models_Library(&scheme_1);
        Models_Library(&scheme_1,&scheme_2,&scheme_3,&scheme_4,&scheme_1_d,&scheme_2_d,&scheme_3_d,&scheme_4_d);
    return allmodels[modelName];
}



inline void print_model_Priors(double covar) {
  auto allmodels =//Models_Library(&scheme_1);
      Models_Library(&scheme_1,&scheme_2,&scheme_3,&scheme_4,&scheme_1_d,&scheme_2_d,&scheme_3_d,&scheme_4_d,
                                    &model00, &model00_7, &model01, &model4, &model4_g_lin,
                     &model6, &model6_no_inactivation,
                     &model6_Eff_no_inactivation, &model7, &model8, &model9);

  std::apply(
      [&covar](auto... ptr_models) {
        (
            [&covar](auto modelp) {
              auto &par = modelp->parameters();
              auto prior = var::prior_around(par, covar);
              var::write_Parameters(par.IdName() + "_par.csv", ",", par);
              write_Prior(par.IdName() + "_prior.csv", ",", prior);
            }(ptr_models),
            ...);
      },
      allmodels());
}



inline auto get_model_old(std::string modelName) -> std::variant<
    /*decltype(&model4),*/ decltype(&model6_Eff_no_inactivation)> {
  using return_type = std::variant<
      /*decltype(&model4), */ decltype(&model6_Eff_no_inactivation)>;

  //  if (modelName=="model4")
  //     return return_type(&model4);
  // else
  return return_type(&model6_Eff_no_inactivation);
}

using Model_v = decltype(get_model(std::string{}));

} // namespace macrodr

#endif // MODELS_MOFFATTHUME_LINEAR_H
