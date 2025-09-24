#ifndef ALLOSTERIC_MODELS_H
#define ALLOSTERIC_MODELS_H
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <map>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>

#include "maybe_error.h"
#include "parameters.h"
#include "qmodel.h"
#include "variables.h"
namespace macrodr {

using var::Constant;
struct Conformational_change_label
    : public var::Constant<Conformational_change_label, std::string> {
    using Constant::Constant;
};

struct Agonist_label : public Constant<Agonist_label, std::string> {
    using Constant::Constant;
};

struct Agonist_dependency : public Constant<Agonist_dependency, std::optional<Agonist_label>> {
    friend std::ostream& print(std::ostream& os, Agonist_dependency const& x) {
        if (x())
            os << "(" << (*x())() << ")";
        return os;
    }
    friend std::ostream& report(std::ostream& os, Agonist_dependency const& x) {
        if (x())
            os << (*x())();
        else
            os << "";
        return os;
    }
    friend std::ostream& report_title(std::ostream& os, Agonist_dependency const&, std::string s) {
        os << "Agonist_dependency" << s;
        return os;
    }

    using Constant::Constant;
};

struct Agonist_dependency_map
    : public Constant<Agonist_dependency_map,
                      std::map<Conformational_change_label, Agonist_dependency>> {
    using Constant::Constant;
};

struct Conformational_change
    : public Constant<Conformational_change,
                      Vector_Space<Conformational_change_label, Agonist_dependency>> {
    friend std::ostream& print(std::ostream& os, Conformational_change const& x) {
        os << get<Conformational_change_label>(x())();
        print(os, get<Agonist_dependency>(x()));
        return os;
    }
    friend std::ostream& report(std::ostream& os, Conformational_change const& x) {
        os << get<Conformational_change_label>(x())() << ",";
        report(os, get<Agonist_dependency>(x()));
        return os;
    }
    friend std::ostream& report_title(std::ostream& os, Conformational_change const&,
                                      std::string s) {
        os << "Conformational_change_label" << s << ",";
        report_title(os, Agonist_dependency{}, s);
        return os;
    }

    using Constant::Constant;
};

struct Conformational_change_domain_state
    : public Constant<Conformational_change_domain_state, bool> {
    using Constant::Constant;
};

struct Conformational_change_state_vector
    : public Constant<Conformational_change_state_vector,
                      std::vector<Conformational_change_domain_state>> {
    using Constant::Constant;
    friend std::ostream& print(std::ostream& os, Conformational_change_state_vector const& x) {
        for (auto e : x()) os << e();
        os << " ";
        return os;
    }
    friend std::ostream& report(std::ostream& os, Conformational_change_state_vector const& x) {
        for (std::size_t i = 0; i + 1 < x().size(); ++i) os << x()[i]() << ",";
        if (!x().empty())
            os << x().back();
        return os;
    }
    friend std::ostream& report_title(std::ostream& os, Conformational_change_state_vector const& x,
                                      std::string s) {
        for (std::size_t i = 0; i + 1 < x().size(); ++i)
            os << "conformational_state_" << i << s << ",";
        if (!x().empty())
            os << "conformational_state_" << x().size() - 1 << s;
        return os;
    }
};

struct Conductance_interaction_label : public Constant<Conductance_interaction_label, std::string> {
    using Constant::Constant;
};

struct Conformational_change_scheme
    : public Constant<Conformational_change_scheme, std::vector<Conformational_change>> {
    using Constant::Constant;
    friend std::ostream& print(std::ostream& os, Conformational_change_scheme const& x) {
        os << "scheme: {";
        for (auto e : x()) {
            print(os, e);
            os << ", ";
        }
        os << "}\n";
        return os;
    }
    friend std::ostream& report(std::ostream& os, Conformational_change_scheme const& x) {
        for (std::size_t i = 0; i + 1 < x().size(); ++i) os << x()[i]() << ",";
        if (!x().empty())
            os << x().back()();
        return os;
    }
    friend std::ostream& report_title(std::ostream& os, Conformational_change_scheme const& x,
                                      std::string) {
        for (std::size_t i = 0; i + 1 < x().size(); ++i)
            report_title(os, x()[i](), "_change_" + std::to_string(i)) << ",";
        if (!x().empty())
            report_title(os, x().back()(), "_change_" + std::to_string(x().size() - 1));
        return os;
    }
};

struct Conformational_position : public Constant<Conformational_position, std::size_t> {
    using Constant::Constant;
};

struct Conformational_interaction_positions
    : public Constant<Conformational_interaction_positions,
                      std::vector<std::vector<Conformational_position>>> {
    using Constant::Constant;
    friend std::ostream& print(std::ostream& os, Conformational_interaction_positions const& x) {
        os << "positions: {";
        for (auto e : x()) {
            os << "{";
            for (auto ee : e) {
                print(os, ee());
                os << ",";
            }
            os << "} ";
        }
        os << "}\n";
        return os;
    }
};

struct Conductance_interaction_positions
    : public Constant<Conductance_interaction_positions,
                      std::vector<std::vector<Conformational_position>>> {
    using Constant::Constant;
    friend std::ostream& print(std::ostream& os, Conductance_interaction_positions const& x) {
        os << "positions: {";
        for (auto e : x()) {
            os << "{";
            for (auto ee : e) {
                print(os, ee());
                os << ",";
            }
            os << "} ";
        }
        os << "}\n";
        return os;
    }
};
struct Conductance_interaction_players
    : public Constant<Conductance_interaction_players, std::vector<Conformational_change_label>> {
    using Constant::Constant;
    friend std::ostream& print(std::ostream& os, Conductance_interaction_players const& x) {
        os << "players: {";
        for (auto e : x()) {
            print(os, e());
            os << ",";
        }
        os << "} ";
        return os;
    }
};

struct Conformational_interaction_label
    : public Constant<Conformational_interaction_label, std::string> {
    using Constant::Constant;
};

struct Conformational_interaction_players
    : public Constant<Conformational_interaction_players,
                      std::vector<Conformational_change_label>> {
    using Constant::Constant;
    friend std::ostream& print(std::ostream& os, Conformational_interaction_players const& x) {
        os << "players: {";
        for (auto e : x()) {
            print(os, e());
            os << ",";
        }
        os << "} ";
        return os;
    }
};

struct Conformational_interaction
    : public Constant<
          Conformational_interaction,
          Vector_Space<Conformational_interaction_label, Conformational_interaction_players,
                       Conformational_interaction_positions>> {
    using Constant::Constant;
    friend std::ostream& print(std::ostream& os, Conformational_interaction const& x) {
        os << "conf_inter{";
        print(os, get<Conformational_interaction_label>(x())());
        os << " => ";
        print(os, get<Conformational_interaction_players>(x()));
        print(os, get<Conformational_interaction_positions>(x()));
        return os << "}\n";
    }
};

struct Conductance_interaction
    : public Constant<Conductance_interaction,
                      Vector_Space<Conductance_interaction_label, Conductance_interaction_players,
                                   Conductance_interaction_positions>> {
    using Constant::Constant;
    friend std::ostream& print(std::ostream& os, Conductance_interaction const& x) {
        os << "cond_inter{";
        print(os, get<Conductance_interaction_label>(x())());
        os << " => ";
        print(os, get<Conductance_interaction_players>(x()));
        print(os, get<Conductance_interaction_positions>(x()));
        return os << "}\n";
    }
};

struct Conformational_interaction_scheme
    : public Constant<Conformational_interaction_scheme, std::vector<Conformational_interaction>> {
    using Constant::Constant;
    friend std::ostream& print(std::ostream& os, Conformational_interaction_scheme const& x) {
        os << "\ninteraction scheme begin:{ \n";
        for (auto e : x()) {
            print(os, e);
        }
        os << "}interaction end\n";
        return os;
    }
};

struct Conductance_interaction_scheme
    : public Constant<Conductance_interaction_scheme, std::vector<Conductance_interaction>> {
    friend std::ostream& print(std::ostream& os, Conductance_interaction_scheme const& x) {
        os << "conductance interaction scheme begin\n";
        for (auto e : x()) {
            print(os, e);
        }
        os << "end\n";
        return os;
    }
};

struct Conformational_interaction_index
    : public Constant<Conformational_interaction_index, std::size_t> {};

struct Conductance_interaction_index : public Constant<Conductance_interaction_index, std::size_t> {
    friend bool operator<(const Conductance_interaction_index& one,
                          const Conductance_interaction_index& two) {
        return one() < two();
    }
};

struct Conformational_interaction_subposition
    : public Constant<Conformational_interaction_subposition, std::size_t> {};

struct Conformational_interactions_domain_state
    : public Constant<Conformational_interactions_domain_state,
                      std::map<Vector_Space<Conformational_interaction_index,
                                            Conformational_interaction_subposition>,
                               int>> {
    using Constant<Conformational_interactions_domain_state,
                   std::map<Vector_Space<Conformational_interaction_index,
                                         Conformational_interaction_subposition>,
                            int>>::Constant;

    int operator[](const Vector_Space<Conformational_interaction_index,
                                      Conformational_interaction_subposition>& x) const {
        if (auto it = (*this)().find(x); it == (*this)().end())
            return 0;
        else
            return it->second;
    }

    friend std::ostream& print(std::ostream& os,
                               Conformational_interactions_domain_state const& x) {
        // os<<"Interactions:[";
        for (auto& e : x()) {
            os << get<Conformational_interaction_index>(e.first)() << "_"
               << get<Conformational_interaction_subposition>(e.first)();
            if (e.second != 1)
                os << "^" << e.second;
            os << " ";
        }
        // os<<"]";
        return os;
    }

    void insert(
        Vector_Space<Conformational_interaction_index, Conformational_interaction_subposition> e) {
        (*this)()[e] += 1;
    }

    friend Conformational_interactions_domain_state operator-(
        const Conformational_interactions_domain_state& one,
        const Conformational_interactions_domain_state& subtract) {
        Conformational_interactions_domain_state out;
        for (auto e : one()) {
            out()[e.first] += e.second;
        }
        for (auto e : subtract()) {
            out()[e.first] -= e.second;
        }
        return out;
    }
};

struct Conformational_transition_mulitplicity
    : public Constant<Conformational_transition_mulitplicity, std::size_t> {};

struct Conformational_interactions_transition_state
    : public Constant<Conformational_interactions_transition_state,
                      std::map<Conformational_interactions_domain_state,
                               Conformational_transition_mulitplicity>> {
    using base_type = Constant<
        Conformational_interactions_transition_state,
        std::map<Conformational_interactions_domain_state, Conformational_transition_mulitplicity>>;
    using base_type::Constant;

    Conformational_interactions_transition_state() = default;

    Conformational_interactions_transition_state(
        const Conformational_interactions_domain_state& v_int)
        : base_type{std::map<Conformational_interactions_domain_state,
                             Conformational_transition_mulitplicity>{
              std::pair{v_int, Conformational_transition_mulitplicity(1)}}} {}

    friend std::ostream& print(std::ostream& os,
                               Conformational_interactions_transition_state const& x) {
        // os<<"Interactions:[";
        for (auto& e : x()) {
            os << e.first;
            if (e.second != 1)
                os << "*" << e.second;
            os << " ";
        }
        // os<<"]";
        return os;
    }

    friend Conformational_interactions_transition_state operator-(
        const Conformational_interactions_transition_state& one,
        const Conformational_interactions_domain_state& subtract) {
        Conformational_interactions_transition_state out;
        for (auto& e : one()) {
            out()[e.first - subtract] = e.second;
        }
        return out;
    }

    friend void sum(Conformational_interactions_transition_state& current,
                    const Conformational_interactions_transition_state& other) {
        for (auto& e : other()) {
            current()[e.first]() += e.second();
        }
    }
};

struct Conformational_domain_state
    : public Constant<Conformational_domain_state,
                      Vector_Space<Conformational_change_label, Conformational_change_domain_state,
                                   Conformational_interactions_domain_state>> {
    friend std::ostream& print(std::ostream& os, Conformational_domain_state const& x) {
        print(os, get<Conformational_change_label>(x())());
        os << " ";
        print(os, get<Conformational_change_domain_state>(x())());
        os << " [";
        print(os, get<Conformational_interactions_domain_state>(x()));
        os << "]";

        return os;
    }
};

struct Conformational_interactions_state_vector
    : public Constant<Conformational_interactions_state_vector,
                      std::vector<Conformational_interactions_domain_state>> {
    friend std::ostream& print(std::ostream& os,
                               Conformational_interactions_state_vector const& x) {
        for (auto& e : x()) {
            os << "[";
            print(os, e);
            os << "]";
        }
        return os;
    }
};

struct Conformational_state_vector
    : public Constant<Conformational_state_vector,
                      Vector_Space<Conformational_change_state_vector,
                                   Conformational_interactions_state_vector>> {
    friend std::ostream& print(std::ostream& os, Conformational_state_vector const& x) {
        os << " state: ";
        print(os, get<Conformational_change_state_vector>(x()));
        os << " interactions: ";
        print(os, get<Conformational_interactions_state_vector>(x()));
        return os;
    }

    friend std::ostream& report(std::ostream& os, Conformational_state_vector const& x) {
        os << " state: ";
        print(os, get<Conformational_change_state_vector>(x()));
        os << " interactions: ";
        print(os, get<Conformational_interactions_state_vector>(x()));
        return os;
    }
};

struct Conformational_state_count
    : public Constant<Conformational_state_count,
                      std::map<Conformational_domain_state, std::size_t>> {
    using base_type =
        Constant<Conformational_state_count, std::map<Conformational_domain_state, std::size_t>>;

    using base_type::Constant;

    friend std::ostream& print(std::ostream& os, Conformational_state_count const& x) {
        os << " Conformational_state_count: {";
        for (auto& e : x()) {
            os << "{";
            print(os, e.first);
            os << "}->" << e.second << "  ";
        }
        return os;
    }

    friend std::ostream& report(std::ostream& os, Conformational_state_count const& x) {
        os << " Conformational_state_count: ";
        for (auto& e : x()) {
            os << " ,";
            print(os, e.first);
            os << "," << e.second << ",";
        }
        return os;
    }
};

struct Conductance_state_count
    : public Constant<Conductance_state_count, std::map<Conductance_interaction_index, int>> {
    using Constant<Conductance_state_count, std::map<Conductance_interaction_index, int>>::Constant;

    friend Conductance_state_count operator-(const Conductance_state_count& one,
                                             const Conductance_state_count& subtract) {
        Conductance_state_count out;
        for (auto e : one()) {
            out()[e.first] += e.second;
        }
        for (auto e : subtract()) {
            out()[e.first] -= e.second;
        }
        return out;
    }
    friend std::ostream& print(std::ostream& os, Conductance_state_count const& x) {
        os << " conductance: ";
        for (auto& e : x()) {
            os << "[";
            print(os, e.first());
            os << "]->" << e.second << "  ";
        }
        return os;
    }
};

struct Conformational_model_scheme
    : public Constant<Conformational_model_scheme,
                      Vector_Space<Conformational_change_scheme, Conformational_interaction_scheme,
                                   Conductance_interaction_scheme>> {
    friend std::ostream& print(std::ostream& os, Conformational_model_scheme const& x) {
        os << "\nConformational model scheme\n";
        print(os, get<Conformational_change_scheme>(x()));
        print(os, get<Conformational_interaction_scheme>(x()));
        print(os, get<Conductance_interaction_scheme>(x()));
        os << "\n end of Conformational model scheme\n-------------------------\n";
        return os;
    }
};

struct Conformational_transition_direction
    : public Constant<Conformational_transition_direction, bool> {};

struct Conformational_state_index : public Constant<Conformational_state_index, std::size_t> {};

struct Conformational_transition_initiating_state_index
    : public Constant<Conformational_transition_initiating_state_index,
                      Conformational_state_index> {};
struct Conformational_transition_landing_state_index
    : public Constant<Conformational_transition_landing_state_index, Conformational_state_index> {};

struct Conformational_transition
    : public Constant<
          Conformational_transition,
          Vector_Space<Conformational_transition_initiating_state_index,
                       Conformational_transition_landing_state_index, Agonist_dependency,
                       Conformational_transition_direction, Conformational_change_label,
                       Conformational_interactions_transition_state>> {
    friend std::ostream& print(std::ostream& os, Conformational_transition const& x) {
        os << get<Conformational_transition_initiating_state_index>(x())() << "->";
        os << get<Conformational_transition_landing_state_index>(x())();

        os << ": ";

        os << get<Conformational_change_label>(x())();
        if (get<Conformational_transition_direction>(x())())
            os << "_on";
        else
            os << "_off";

        print(os, get<Agonist_dependency>(x()));

        os << " [";
        print(os, get<Conformational_interactions_transition_state>(x()));
        os << "]";
        return os;
    }
};

struct Conformational_transition_list
    : public Constant<Conformational_transition_list,
                      std::vector<std::vector<Conformational_transition>>> {
    friend std::ostream& print(std::ostream& os, Conformational_transition_list const& x) {
        for (auto i = 0ul; i < x().size(); ++i) {
            for (auto j = 0ul; j < x()[i].size(); ++j) {
                print(os, x()[i][j]);
                os << "\n";
            }
        }
        return os;
    }
};

struct Conformational_states
    : public Constant<
          Conformational_states,
          std::vector<Vector_Space<Conformational_state_count, Conformational_state_vector,
                                   Conductance_state_count>>> {
    friend Maybe_error<std::size_t> find_state_with_conformation(
        const Conformational_states& x, const Conformational_domain_state& e) {
        auto it =
            std::find_if(x().begin(), x().end(),
                         [&e](Vector_Space<Conformational_state_count, Conformational_state_vector,
                                           Conductance_state_count> const& m) {
                             auto& v_st = get<Conformational_state_count>(m)();
                             auto it = v_st.find(e);
                             return (it != v_st.end() && it->second > 0);
                         });
        if (it == x().end()) {
            std::stringstream ss;
            print(ss, e);
            ss << " not found in ";
            print(ss, x);
            return error_message(ss.str());
        } else {
            return it - x().begin();
        }
    }

    friend std::ostream& print(std::ostream& os, Conformational_states const& x) {
        os << "\nConformational states\n";
        for (auto i = 0ul; i < x().size(); ++i) {
            auto& e = x()[i];
            os << "\n" << i << "-->";
            print(os, get<Conformational_state_vector>(e));
            os << "\t";
            print(os, get<Conformational_state_count>(e));
            os << "\t";
            print(os, get<Conductance_state_count>(e));
        }
        os << "\n end of Conformational states\n";
        return os;
    }

    friend std::ostream& report(std::ostream& os, Conformational_states const& x) {
        for (auto i = 0ul; i < x().size(); ++i) {
            auto& e = x()[i];
            os << i << ",";
            report(os, get<Conformational_state_vector>(e)) << ",";
            report(os, get<Conformational_state_count>(e));
            os << "\t";
            print(os, get<Conductance_state_count>(e));
        }
        os << "\n end of Conformational states\n";
        return os;
    }
};

struct Conformational_state_count_to_representative_index
    : public Constant<Conformational_state_count_to_representative_index,
                      std::map<Conformational_state_count, Conformational_state_index>> {
    Maybe_error<Conformational_state_index> operator[](Conformational_state_count i) const {
        if (auto it = (*this)().find(i); it == (*this)().end()) {
            std::stringstream ss;
            ss << "Conformational state count ";
            ss << i;
            ss << " not found";
            ss << "All states:";
            ss << *this;
            ss << "\npppppppppppppppp\n";
            print(ss, i);
            ss << " not found\n";
            print(ss, *this);
            // std::cerr<<i;
            // std::cerr<<*this;
            std::cerr << "\nerror\n" << ss.str() << "\n";
            return error_message(ss.str());
        } else
            return it->second;
    }

    friend std::ostream& print(std::ostream& os,
                               const Conformational_state_count_to_representative_index& x) {
        os << "\nConformational_state_count_to_representative_index:\n";
        for (auto e : x()) {
            print(os, e.first);
            os << "state_number: " << e.second() << "\n";
        }
        return os;
    }
};

struct Conformational_model
    : public Constant<Conformational_model,
                      Vector_Space<N_St, Conformational_model_scheme, Conformational_states,
                                   Conformational_transition_list>> {
    friend std::ostream& print(std::ostream& os, Conformational_model const& x) {
        os << "Conformational model\n";
        print(os, get<Conformational_model_scheme>(x()));
        print(os, get<Conformational_states>(x()));
        print(os, get<Conformational_transition_list>(x()));
        return os;
    }
};

struct Conformation_change_standard_state
    : public Constant<Conformation_change_standard_state,
                      Conformational_interactions_domain_state> {};

struct Conformation_change_standard_map
    : public Constant<Conformation_change_standard_map,
                      std::map<Conformational_change_label, Conformation_change_standard_state>> {
    Conformational_interactions_domain_state operator[](
        const Conformational_change_label& la) const {
        auto it = (*this)().find(la);
        if (it != (*this)().end()) {
            return it->second();
        } else {
            return Conformational_interactions_domain_state{};
        }
    }

    friend std::ostream& print(std::ostream& os, Conformation_change_standard_map const& x) {
        for (auto& e : x()) {
            print(os, e.first);
            print(os, e.second());
            os << "\n";
        }
        return os;
    }
};

struct Conductance_leakeage_ratio_label
    : public Constant<Conductance_leakeage_ratio_label, std::string> {};

enum class Conductance_interaction_kind { additive, multiplicative, equilibrium };

struct Conductance_interaction_type
    : public Constant<Conductance_interaction_type, Conductance_interaction_kind> {
    friend std::ostream& print(std::ostream& os, Conductance_interaction_type const& x) {
        switch (x()) {
            case Conductance_interaction_kind::additive:
                print(os, "additive");
                break;
            case Conductance_interaction_kind::multiplicative:
                print(os, "multiplicative");

                break;
            case Conductance_interaction_kind::equilibrium:
                print(os, "equilibrium");
                break;
        }
        return os;
    }
};

struct Conductance_interaction_info
    : public Constant<
          Conductance_interaction_info,
          Vector_Space<Conductance_interaction_type, Conductance_leakeage_ratio_label>> {
    using base_type =
        Constant<Conductance_interaction_info,
                 Vector_Space<Conductance_interaction_type, Conductance_leakeage_ratio_label>>;
    using base_type::Constant;
    Conductance_interaction_info(Conductance_interaction_kind k,
                                 const std::string& leakeage_ratio_label)
        : base_type{Vector_Space{Conductance_interaction_type{k},
                                 Conductance_leakeage_ratio_label{leakeage_ratio_label}}} {}
};

struct Conformational_model_standarized
    : public Constant<Conformational_model_standarized,
                      Vector_Space<N_St, Conformational_model_scheme, Conformational_states,
                                   Conformational_state_count_to_representative_index,
                                   Conformational_transition_list, Conformation_change_standard_map,
                                   Conductance_interaction_info>> {
    friend std::ostream& print(std::ostream& os, Conformational_model_standarized const& x) {
        os << "Conformational model\n";
        print(os, get<Conformational_model_scheme>(x()));
        print(os, get<Conformational_states>(x()));
        print(os, get<Conformational_transition_list>(x()));
        print(os, get<Conformation_change_standard_map>(x()));
        print(os, get<Conductance_interaction_info>(x()));
        return os;
    }
};

namespace impl {

inline auto make_Conformational_change_state_vector(const Conformational_model_scheme& model,
                                                    std::size_t n) {
    auto number_units = get<Conformational_change_scheme>(model())().size();

    std::vector<Conformational_change_domain_state> out(number_units);
    for (std::size_t i = 0; i < number_units; ++i)
        out[i] = Conformational_change_domain_state((n & (1ul << i)) == (1ul << i));
    return Conformational_change_state_vector(std::move(out));
}

inline Conformational_interactions_state_vector make_Conformational_interaction_state_vector(
    const Conformational_model_scheme& model, const Conformational_change_state_vector& state) {
    auto number_units = get<Conformational_change_scheme>(model())().size();
    auto v_inter = get<Conformational_interaction_scheme>(model());

    auto number_inter = v_inter().size();
    std::vector<Conformational_interactions_domain_state> out(number_units);
    for (std::size_t i = 0; i < number_units; ++i) {
        for (std::size_t j = 0; j < number_inter; ++j) {
            auto v_ipos = get<Conformational_interaction_positions>(v_inter()[j]());
            for (std::size_t k = 0; k < v_ipos().size(); ++k) {
                bool includes_i = false;
                bool includes_all_other = true;
                std::size_t i_sub_position;
                for (std::size_t kk = 0; kk < v_ipos()[k].size(); ++kk) {
                    if (v_ipos()[k][kk]() == i) {
                        includes_i = true;
                        i_sub_position = kk;
                    } else {
                        includes_all_other = includes_all_other && (state()[v_ipos()[k][kk]()]());
                    }
                }
                if (includes_all_other && includes_i)
                    out[i].insert(
                        Vector_Space(Conformational_interaction_index(j),
                                     Conformational_interaction_subposition(i_sub_position)));
            }
        }
    }
    return Conformational_interactions_state_vector(std::move(out));
}

inline auto make_Conformational_state_vector(const Conformational_model_scheme& model,
                                             std::size_t n) {
    auto v_change = make_Conformational_change_state_vector(model, n);
    auto v_inter = make_Conformational_interaction_state_vector(model, v_change);
    return Conformational_state_vector(Vector_Space(std::move(v_change), std::move(v_inter)));
}

inline Maybe_error<Conformational_state_count> to_state_count(
    const Conformational_model_scheme& model, const Conformational_state_vector& state) {
    if (get<Conformational_change_state_vector>(state())().size() !=
        get<Conformational_change_scheme>(model())().size())
        return error_message("unequal sizes");
    else {
        auto& v_state = get<Conformational_change_state_vector>(state());
        auto& v_inter = get<Conformational_interactions_state_vector>(state());
        auto& v_change = get<Conformational_change_scheme>(model());
        std::map<Conformational_domain_state, std::size_t> out;

        for (std::size_t i = 0; i < v_state().size(); ++i) {
            ++out[Conformational_domain_state(Vector_Space(
                get<Conformational_change_label>(v_change()[i]()), v_state()[i], v_inter()[i]))];
        }
        return Conformational_state_count(std::move(out));
    }
}

inline Conductance_state_count to_state_conductance_count(
    const Conformational_model_scheme& model, const Conformational_state_vector& state_vector)

{
    auto state = get<Conformational_change_state_vector>(state_vector());
    auto v_inter = get<Conductance_interaction_scheme>(model());
    auto number_inter = v_inter().size();
    std::map<Conductance_interaction_index, int> out;
    for (std::size_t i = 0; i < number_inter; ++i) {
        auto v_ipos = get<Conductance_interaction_positions>(v_inter()[i]());
        for (std::size_t k = 0; k < v_ipos().size(); ++k) {
            bool includes_all = true;
            for (std::size_t kk = 0; kk < v_ipos()[k].size(); ++kk) {
                includes_all = includes_all && (state()[v_ipos()[k][kk]()]());
            }
            if (includes_all)
                ++out[Conductance_interaction_index(i)];
        }
    }
    return Conductance_state_count(std::move(out));
}

inline Maybe_error<Conformational_state_index> to_state_index(
    const Conformational_state_count_to_representative_index& model,
    const Conformational_state_count& state_count) {
    if (auto it = model().find(state_count); it == model().end())
        return error_message("state count not found");
    else
        return it->second;
}

inline auto change_conformation(const Conformational_change_state_vector& state,
                                std::size_t ith_domain) {
    Conformational_change_state_vector out(state);
    out()[ith_domain]() = !out()[ith_domain]();
    return out;
}

inline auto change_conformation(const Conformational_model_scheme& model,
                                const Conformational_state_vector& state, std::size_t ith_domain) {
    auto v_change =
        change_conformation(get<Conformational_change_state_vector>(state()), ith_domain);
    auto v_inter = make_Conformational_interaction_state_vector(model, v_change);
    return Conformational_state_vector(Vector_Space(std::move(v_change), std::move(v_inter)));
}

// Conformational_state_count_to_representative_index sends to the last inserted
// index of the same State_count
inline Maybe_error<
    std::tuple<Conformational_states, Conformational_state_count_to_representative_index>>
    make_Conformational_states_and_index(const Conformational_model_scheme& model) {
    auto number_units = get<Conformational_change_scheme>(model())().size();
    Conformational_states out;
    Conformational_state_count_to_representative_index map;
    for (std::size_t n = 0; n < (1ul << number_units); ++n) {
        auto v_state_vector = make_Conformational_state_vector(model, n);
        auto v_state_count = to_state_count(model, v_state_vector);
        auto v_state_conductance = to_state_conductance_count(model, v_state_vector);
        if (!v_state_count) {
            return v_state_count.error();
        } else if (map().find(v_state_count.value()) == map().end()) {
            map()[v_state_count.value()] = Conformational_state_index(out().size());
            out().push_back(Vector_Space(std::move(v_state_count.value()),
                                         std::move(v_state_vector),
                                         std::move(v_state_conductance)));
        }
    }
    return std::tuple(std::move(out), std::move(map));
}

inline Maybe_error<Conformational_transition_list> make_Conformational_transition_list(
    const Conformational_model_scheme& model, const Conformational_states states,
    const Conformational_state_count_to_representative_index& map) {
    std::vector<std::vector<Conformational_transition>> out(states().size());

    auto scheme = get<Conformational_change_scheme>(model());

    for (std::size_t i = 0; i < states().size(); ++i) {
        //  auto v_count=get<Conformational_state_count>(states()[i]);
        auto v_state = get<Conformational_state_vector>(states()[i]);
        auto v_change = get<Conformational_change_state_vector>(v_state());
        auto v_interactions = get<Conformational_interactions_state_vector>(v_state());
        auto i_start =
            Conformational_transition_initiating_state_index(Conformational_state_index(i));

        std::vector<Conformational_transition> tran;

        for (std::size_t j = 0; j < v_change().size(); ++j) {
            auto j_change = scheme()[j];
            auto j_inter = Conformational_interactions_transition_state(v_interactions()[j]);
            auto j_state = change_conformation(model, v_state, j);
            auto Maybe_j_count = to_state_count(model, j_state);
            if (!Maybe_j_count)
                return Maybe_j_count.error();
            else {
                auto& j_count = Maybe_j_count.value();
                auto Maybe_j_index = to_state_index(map, j_count);
                if (!Maybe_j_index)
                    return Maybe_j_index.error();
                else {
                    auto i_end =
                        Conformational_transition_landing_state_index(Maybe_j_index.value());
                    bool already = false;
                    for (auto& elem : tran) {
                        if (get<Conformational_transition_landing_state_index>(elem())()() ==
                            i_end()()) {
                            sum(get<Conformational_interactions_transition_state>(elem()), j_inter);
                            already = true;
                        }
                    }
                    if (!already) {
                        auto a = get<Agonist_dependency>(j_change());
                        auto d = Conformational_transition_direction(!v_change()[j]());
                        auto l = get<Conformational_change_label>(j_change());
                        tran.push_back(Conformational_transition(
                            Vector_Space(i_start, i_end, a, d, l, j_inter)));
                    }
                }
            }
        }
        out[i] = tran;
    }

    return Conformational_transition_list(std::move(out));
}

inline Maybe_error<Conformational_change_scheme> make_Conformational_change_scheme(
    Agonist_dependency_map&& t_agonist_map, std::vector<Conformational_change_label>&& t_scheme) {
    std::vector<Conformational_change> out(t_scheme.size());
    for (std::size_t i = 0; i < t_scheme.size(); ++i) {
        if (auto it = t_agonist_map().find(t_scheme[i]); it != t_agonist_map().end()) {
            out[i] = Conformational_change(Vector_Space(t_scheme[i], it->second));

        } else {
            return error_message("Conformational_change " + t_scheme[i]() +
                                 " agonist dependency condition unknown");
        }
    }
    return Conformational_change_scheme(std::move(out));
}

inline Maybe_error<bool> check_Conformational_interaction(
    const Conformational_change_scheme& t_scheme, Conformational_interaction const& t_interaction) {
    auto& v_players = get<Conformational_interaction_players>(t_interaction());
    auto& v_positions = get<Conformational_interaction_positions>(t_interaction());

    Maybe_error<bool> out(true);

    for (auto i = 0ul; i < v_positions().size(); ++i) {
        if (v_players().size() != v_positions()[i].size())
            out = error_message(out.error()() + "  " + std::to_string(i) +
                                "th interaction position size mismatch: " +
                                std::to_string(v_positions()[i].size()) + " vs " +
                                std::to_string(v_players().size()));
        else {
            Maybe_error<bool> outp(true);
            for (auto j = 0ul; j < v_positions()[i].size(); ++j) {
                auto jpos = v_positions()[i][j]();
                if (v_players()[j]() != get<Conformational_change_label>(t_scheme()[jpos]())())
                    outp =
                        error_message(outp.error()() + " at position " + std::to_string(j) + ": " +
                                      get<Conformational_change_label>(t_scheme()[j]())() +
                                      "is not " + v_players()[j]());
            }
            if (!outp)
                out = error_message(out.error()() + "  " + std::to_string(i) +
                                    "th interaction label mismatch: " + out.error()());
        }
    }
    return out;
}

inline Maybe_error<Conformational_model_scheme> make_Conformational_model_scheme(
    Agonist_dependency_map&& t_agonist_map, std::vector<Conformational_change_label>&& t_scheme,
    std::vector<Conformational_interaction>&& t_interactions,
    std::vector<Conductance_interaction>&& t_conductance) {
    auto Maybe_Conformational_change =
        make_Conformational_change_scheme(std::move(t_agonist_map), std::move(t_scheme));
    if (!Maybe_Conformational_change)
        return Maybe_Conformational_change.error();
    else {
        Maybe_error<bool> out(true);
        for (std::size_t i = 0; i < t_interactions.size(); ++i) {
            if (auto outp = check_Conformational_interaction(Maybe_Conformational_change.value(),
                                                             t_interactions[i]);
                !outp)
                out = error_message(out.error()() + " " + std::to_string(i) +
                                    "th interaction: " + outp.error()() + "\n");
        }
        if (!out)
            return out.error();
        else
            return Conformational_model_scheme(
                Vector_Space(std::move(Maybe_Conformational_change.value()),
                             Conformational_interaction_scheme(std::move(t_interactions)),
                             Conductance_interaction_scheme(std::move(t_conductance))));
    }
}

template <class Id, class P>
    requires std::is_same_v<var::untransformed_type_t<P>, var::Parameters_values>
auto calc_Qij(const Conformational_interaction_scheme& inter, const var::Parameters_Names& names,
              const P& par, const Conformational_transition& tr,
              const Conformational_interactions_transition_state& v_int_tr_map)
    -> Maybe_error<var::Op_t<transformation_type_t<P>, double>> {
    auto ag = get<Agonist_dependency>(tr());
    auto d = get<Conformational_transition_direction>(tr());
    auto chla = get<Conformational_change_label>(tr());

    auto Maybe_i_base = d() ? names[chla() + "_on"] : names[chla() + "_off"];

    if (!Maybe_i_base)
        return Maybe_i_base.error();
    auto i_base = Maybe_i_base.value();
    auto k = par()[i_base];
    std::decay_t<decltype(k)> out = 0.0;
    for (auto& e : v_int_tr_map()) {
        auto ki = k * e.second();

        for (auto ii = e.first().begin(); ii != e.first().end(); ++ii) {
            auto factor_la = get<Conformational_interaction_label>(
                inter()[get<Conformational_interaction_index>(ii->first)()]())();
            auto factor_ipos = get<Conformational_interaction_subposition>(ii->first)();

            auto factor_power = ii->second;

            auto Maybe_i_Factor = names[factor_la];
            auto Maybe_i_Factor_pos = names[factor_la + "_" + std::to_string(factor_ipos)];
            if (!Maybe_i_Factor || !Maybe_i_Factor_pos)
                return error_message(Maybe_i_Factor.error()() + Maybe_i_Factor_pos.error()());
            else {
                auto i_Factor = Maybe_i_Factor.value();
                auto i_Factor_pos = Maybe_i_Factor_pos.value();

                using std::pow;
                if (d())
                    ki = ki * pow(par()[i_Factor_pos], factor_power);
                else
                    ki = ki * pow(par()[i_Factor_pos] / par()[i_Factor], factor_power);
            }
        }

        out = out + ki;
    }
    return out;
}

template <class Id, class P>
    requires std::is_same_v<var::untransformed_type_t<P>, var::Parameters_values>
auto calc_Qij(const Conformational_interaction_scheme& inter, const var::Parameters_Names& names,
              const P& par, const Conformational_transition& tr,
              const Conformation_change_standard_map& st)
    -> Maybe_error<var::Op_t<transformation_type_t<P>, double>> {
    auto chla = get<Conformational_change_label>(tr());
    auto v_int = get<Conformational_interactions_transition_state>(tr()) - st[chla];

    return calc_Qij<Id>(inter, names, par, tr, v_int);
}

template <class Id, class P>
    requires std::is_same_v<var::untransformed_type_t<P>, var::Parameters_values>
auto calc_Qij(const Conformational_interaction_scheme& inter, const var::Parameters_Names& names,
              const P& par, const Conformational_transition& tr)
    -> Maybe_error<var::Op_t<transformation_type_t<P>, double>> {
    auto v_int = get<Conformational_interactions_transition_state>(tr());
    return calc_Qij<Id>(inter, names, par, tr, v_int);
}

template <class Id>
auto calc_Qij_formula(const Conformational_interaction_scheme& inter,
                      const var::Parameters_Names& names, const Conformational_transition& tr,
                      const Conformational_interactions_transition_state& v_int)
    -> Maybe_error<std::string> {
    auto ag = get<Agonist_dependency>(tr());
    auto d = get<Conformational_transition_direction>(tr());
    auto chla = get<Conformational_change_label>(tr());

    auto Maybe_i_base = d() ? names[chla() + "_on"] : names[chla() + "_off"];

    if (!Maybe_i_base) {
        return Maybe_i_base.error();
    } else {
        auto i_base = Maybe_i_base.value();
        auto out = names()[i_base];
        if (v_int().size() > 1) {
            out += " * (";
        } else {
            if (v_int().begin()->second != 1)
                out += " * ";
        }

        for (auto it = v_int().begin(); it != v_int().end(); ++it) {
            auto& e = *it;
            if (it != v_int().begin())
                out += " + ";
            if (e.second != 1)
                out += std::to_string(e.second()) + " * ";
            for (auto ii = e.first().begin(); ii != e.first().end(); ++ii) {
                auto factor_la = get<Conformational_interaction_label>(
                    inter()[get<Conformational_interaction_index>(ii->first)()]())();
                auto factor_ipos = get<Conformational_interaction_subposition>(ii->first)();
                auto factor_power = ii->second;

                auto Maybe_i_Factor = names[factor_la];
                auto Maybe_i_Factor_pos = names[factor_la + "_" + std::to_string(factor_ipos)];
                if (!Maybe_i_Factor || !Maybe_i_Factor_pos)
                    return error_message(Maybe_i_Factor.error()() + Maybe_i_Factor_pos.error()());
                else {
                    auto i_Factor = Maybe_i_Factor.value();
                    auto i_Factor_pos = Maybe_i_Factor_pos.value();

                    using std::pow;
                    if (ii != e.first().begin())
                        out = out + " * ";
                    if (d())
                        out = out + "pow(" + names()[i_Factor_pos] + "," +
                              std::to_string(factor_power) + ")";
                    else
                        out = out + "pow(" + names()[i_Factor_pos] + "/" + names()[i_Factor] + "," +
                              std::to_string(factor_power) + ")";
                    //      out = out * pow(par[i_Factor], par[i_Factor_pos] - 1.0);
                }
            }
        }
        if (v_int().size() > 1) {
            out += ")";
        }
        return out;
    }
}

template <class Id>
auto calc_Qij_formula(const Conformational_interaction_scheme& inter,
                      const var::Parameters_Names& names, const Conformational_transition& tr,
                      const Conformation_change_standard_map& st) -> Maybe_error<std::string> {
    auto chla = get<Conformational_change_label>(tr());
    auto v_int = get<Conformational_interactions_transition_state>(tr()) - st[chla];

    return calc_Qij_formula<Id>(inter, names, tr, v_int);
}

template <class Id>
auto calc_Qij_formula(const Conformational_interaction_scheme& inter,
                      const var::Parameters_Names& names, const Conformational_transition& tr)
    -> Maybe_error<std::string> {
    auto v_int = get<Conformational_interactions_transition_state>(tr());
    return calc_Qij_formula<Id>(inter, names, tr, v_int);
}

template <class Id, class Conformational_model_, class P>
    requires(std::is_same_v<var::untransformed_type_t<P>, var::Parameters_values> &&
             (std::is_same_v<Conformational_model_, Conformational_model> ||
              std::is_same_v<Conformational_model_, Conformational_model_standarized>))
auto make_Q0_Qa(const Conformational_model_& model, const typename var::Parameters_Names& names,
                const P& par)
    -> Maybe_error<std::tuple<Transfer_Op_to<P, Q0>, Transfer_Op_to<P, Qa>>> {
    using Trans = transformation_type_t<P>;

    auto N = get<N_St>(model())();
    auto inter =
        get<Conformational_interaction_scheme>(get<Conformational_model_scheme>(model())());
    auto tr = get<Conformational_transition_list>(model());
    assert(tr().size() == N);
    auto v_Q0 = Op_t<Trans, Q0>(Matrix<double>(N, N, 0.0));
    auto v_Qa = Op_t<Trans, Qa>(Matrix<double>(N, N, 0.0));

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < tr()[i].size(); ++j) {
            auto trr = tr()[i][j];
            auto v_i = get<Conformational_transition_initiating_state_index>(trr());
            assert(v_i()() == i);
            auto v_j = get<Conformational_transition_landing_state_index>(trr());
            auto ag = get<Agonist_dependency>(trr());
            auto d = get<Conformational_transition_direction>(trr());

            Maybe_error<var::Op_t<transformation_type_t<P>, double>> Maybe_qij;
            if constexpr (std::is_same_v<Conformational_model_, Conformational_model_standarized>)
                Maybe_qij = calc_Qij<Id>(inter, names, par, trr,
                                         get<Conformation_change_standard_map>(model()));
            else
                Maybe_qij = calc_Qij<Id>(inter, names, par, trr);
            if (!Maybe_qij)
                return Maybe_qij.error();
            else if (!ag() || !d()) {
                // set(v_Q0(), i, i, v_Q0()(i, i) - Maybe_qij.value());  later change it
                // back
                set(v_Q0(), i, v_j()(), std::move(Maybe_qij.value()));
            } else {
                //  set(v_Qa(), i, i, v_Qa()(i, i) - Maybe_qij.value());  same
                set(v_Qa(), i, v_j()(), std::move(Maybe_qij.value()));
            }
        }
    }
    return std::tuple(v_Q0, v_Qa);
}

// template <class Id, class P>
//     requires std::is_same_v<var::untransformed_type_t<P>, var::Parameters_values>
// auto get_calc_Qij(const Conformational_interaction_scheme &inter,
//                   const typename var::Parameters_Names &names, const P &par,
//                   const Conformational_transition &tr,
//                   const Conformational_interactions_domain_state &v_int)
//     -> Maybe_error<var::Op_t<transformation_type_t<P>, double>> {
//     auto ag = get<Agonist_dependency>(tr());
//     auto d = get<Conformational_transition_direction>(tr());
//     auto chla = get<Conformational_change_label>(tr());
//     auto n = get<Conformational_transition_mulitplicity>(tr());

//     auto Maybe_i_base = d() ? names[chla() + "_on"] : names[chla() + "_off"];

//     if (!Maybe_i_base)
//         return Maybe_i_base.error();
//     auto i_base = Maybe_i_base.value();
//     auto out = n() * par()[i_base];
//     for (auto ii = v_int().begin(); ii != v_int().end(); ++ii) {
//         auto factor_la = get<Conformational_interaction_label>(
//             inter()[get<Conformational_interaction_index>(ii->first)()]())();
//         auto factor_ipos = get<Conformational_interaction_subposition>(ii->first)();

//         auto factor_power = ii->second;

//         auto Maybe_i_Factor = names[factor_la];
//         auto Maybe_i_Factor_pos =
//             names[factor_la + "_" + std::to_string(factor_ipos)];
//         if (!Maybe_i_Factor || !Maybe_i_Factor_pos)
//             return error_message(Maybe_i_Factor.error()() +
//                                  Maybe_i_Factor_pos.error()());
//         else {
//             auto i_Factor = Maybe_i_Factor.value();
//             auto i_Factor_pos = Maybe_i_Factor_pos.value();

//             using std::pow;
//             if (d())
//                 out = out * pow(par()[i_Factor_pos], factor_power);
//             else
//                 out = out * pow(par()[i_Factor_pos] / par()[i_Factor], factor_power);
//         }
//     }
//     return out;
// }

template <class Id, class Conformational_model_>
    requires(std::is_same_v<Conformational_model_, Conformational_model> ||
             std::is_same_v<Conformational_model_, Conformational_model_standarized>)
auto make_Q0_Qa_formula(const Conformational_model_& model,
                        const typename var::Parameters_Names& names)
    -> Maybe_error<std::tuple<Q0_formula, Qa_formula>> {
    auto N = get<N_St>(model())();
    auto inter =
        get<Conformational_interaction_scheme>(get<Conformational_model_scheme>(model())());
    auto tr = get<Conformational_transition_list>(model());
    assert(tr().size() == N);
    auto v_Q0 =
        Q0_formula(std::vector<std::vector<std::string>>(N, std::vector<std::string>(N, "")));
    auto v_Qa =
        Qa_formula(std::vector<std::vector<std::string>>(N, std::vector<std::string>(N, "")));

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < tr()[i].size(); ++j) {
            auto trr = tr()[i][j];
            auto v_i = get<Conformational_transition_initiating_state_index>(trr());
            assert(v_i()() == i);
            auto v_j = get<Conformational_transition_landing_state_index>(trr());
            auto ag = get<Agonist_dependency>(trr());
            Maybe_error<std::string> Maybe_qij;
            if constexpr (std::is_same_v<Conformational_model_, Conformational_model_standarized>)
                Maybe_qij = calc_Qij_formula<Id>(inter, names, trr,
                                                 get<Conformation_change_standard_map>(model()));
            else
                Maybe_qij = calc_Qij_formula<Id>(inter, names, trr);
            if (!Maybe_qij)
                return Maybe_qij.error();
            else if (!ag()) {
                v_Q0()[i][i] = v_Q0()[i][i] + " -" + Maybe_qij.value();
                v_Q0()[i][v_j()()] = Maybe_qij.value();
            } else {
                v_Qa()[i][i] = v_Qa()[i][i] + " -" + Maybe_qij.value();
                v_Qa()[i][v_j()()] = Maybe_qij.value();
            }
        }
    }
    return std::tuple(v_Q0, v_Qa);
}

template <class Id, class P>
    requires std::is_same_v<var::untransformed_type_t<P>, var::Parameters_values>
auto calc_gi(const Conductance_interaction_scheme& scheme, const var::Parameters_Names& names,
             const P& par, const Conductance_state_count& count,
             const Conductance_interaction_info& int_info)
    -> Maybe_error<Transfer_Op_to<P, double>> {
    auto mult = get<Conductance_interaction_type>(int_info());
    Transfer_Op_to<P, double> out = 0.0;
    if (mult() != Conductance_interaction_kind::additive)
        out = 1.0;
    for (auto it = count().begin(); it != count().end(); ++it) {
        auto i_lab = it->first();
        auto lab = get<Conductance_interaction_label>(scheme()[i_lab]())();
        auto n = it->second;
        auto Maybe_i = names[lab];
        if (!Maybe_i)
            return Maybe_i.error();
        auto i = Maybe_i.value();
        switch (mult()) {
            case Conductance_interaction_kind::additive:
                out = out + par()[i] * n;
                break;
            case Conductance_interaction_kind::multiplicative:
            case Conductance_interaction_kind::equilibrium:
                using std::pow;
                out = out * pow(par()[i], n);
                break;
        }
    }
    switch (mult()) {
        case Conductance_interaction_kind::additive:
        case Conductance_interaction_kind::multiplicative:
            break;
        case Conductance_interaction_kind::equilibrium: {
            auto le_lab = get<Conductance_leakeage_ratio_label>(int_info());
            auto Maybe_ile = names[le_lab()];
            if (!Maybe_ile)
                return Maybe_ile.error();
            auto ile = Maybe_ile.value();
            out = (par()[ile] * out) / (1 + par()[ile] * out);
        }
    }
    return out;
}

template <class Id>
auto calc_gi_formula(const Conductance_interaction_scheme& scheme,
                     const var::Parameters_Names& names, const Conductance_state_count& count,
                     const Conductance_interaction_info& inter_info) -> Maybe_error<std::string> {
    std::string out = "";
    auto mult = get<Conductance_interaction_type>(inter_info());
    if (mult() != Conductance_interaction_kind::additive)
        out = "1";
    for (auto it = count().begin(); it != count().end(); ++it) {
        auto i_lab = it->first();
        auto lab = get<Conductance_interaction_label>(scheme()[i_lab]())();
        auto n = it->second;
        auto Maybe_i = names[lab];
        if (!Maybe_i)
            return Maybe_i.error();
        auto i = Maybe_i.value();
        switch (mult()) {
            case Conductance_interaction_kind::additive:
                if (n > 1)
                    out = out + "+" + names()[i] + "*" + std::to_string(n);
                else
                    out = out + "+" + names()[i];
                break;

            case Conductance_interaction_kind::multiplicative:
            case Conductance_interaction_kind::equilibrium:
                if (n > 1) {
                    if (out != "1")
                        out = out + "*pow(" + names()[i] + "," + std::to_string(n) + ")";
                    else
                        out = "pow(" + names()[i] + "," + std::to_string(n) + ")";
                } else {
                    if (out != "1")
                        out = out + "*" + names()[i];
                    else
                        out = names()[i];
                }
                break;
        }
    }
    switch (mult()) {
        case Conductance_interaction_kind::additive:
        case Conductance_interaction_kind::multiplicative:
            break;
        case Conductance_interaction_kind::equilibrium: {
            auto le_lab = get<Conductance_leakeage_ratio_label>(inter_info());
            auto Maybe_ile = names[le_lab()];
            if (!Maybe_ile)
                return Maybe_ile.error();
            auto ile = Maybe_ile.value();
            out = "(" + names()[ile] + "*" + out + ")/(1+" + names()[ile] + "*" + out + ")";
        } break;
    }
    return out;
}

template <class Id, class Conformational_model_, class P>
    requires(std::is_same_v<var::untransformed_type_t<P>, var::Parameters_values> &&
             (std::is_same_v<Conformational_model_, Conformational_model> ||
              std::is_same_v<Conformational_model_, Conformational_model_standarized>))
auto make_g(const Conformational_model_& model, const typename var::Parameters_Names& names,
            const P& par) -> Maybe_error<Transfer_Op_to<P, g>> {
    using Trans = transformation_type_t<P>;

    auto& v_states = get<Conformational_states>(model());

    auto& v_scheme = get<Conformational_model_scheme>(model());

    auto v_conductance = get<Conductance_interaction_scheme>(v_scheme());

    auto N = get<N_St>(model())();
    assert(v_states().size() == N);
    auto v_g = Op_t<Trans, g>(Matrix<double>(N, 1ul, 0.0));

    // auto v_g_max = Op_t<Trans, double>(0.0);

    for (std::size_t i = 0; i < N; ++i) {
        auto& v_cond_map = get<Conductance_state_count>(v_states()[i]);

        auto Maybe_gi = calc_gi<Id>(v_conductance, names, par, v_cond_map,
                                    get<Conductance_interaction_info>(model()));
        if (!Maybe_gi)
            return Maybe_gi.error();
        else {
            set(v_g(), i, 0, Maybe_gi.value());
        }
    }
    return v_g;
}

template <class Id, class Conformational_model_>
    requires(std::is_same_v<Conformational_model_, Conformational_model> ||
             std::is_same_v<Conformational_model_, Conformational_model_standarized>)
auto make_g_formula(const Conformational_model_& model, const var::Parameters_Names& names)
    -> Maybe_error<g_formula> {
    auto& v_states = get<Conformational_states>(model());

    auto& v_scheme = get<Conformational_model_scheme>(model());

    auto v_conductance = get<Conductance_interaction_scheme>(v_scheme());

    auto N = get<N_St>(model())();
    assert(v_states().size() == N);
    auto v_g = g_formula(std::vector<std::string>(N, ""));

    for (std::size_t i = 0; i < N; ++i) {
        auto& v_cond_map = get<Conductance_state_count>(v_states()[i]);

        auto Maybe_gi = calc_gi_formula<Id>(v_conductance, names, v_cond_map,
                                            get<Conductance_interaction_info>(model()));
        if (!Maybe_gi)
            return Maybe_gi.error();

        v_g()[i] = Maybe_gi.value();
    }
    return v_g;
}

inline auto get_conformational_change_names(const Conformational_change_scheme& sch) {
    std::vector<std::string> out;
    std::set<std::string> labels;
    for (auto& e : sch()) {
        auto la = get<Conformational_change_label>(e())();
        if (labels.find(la) == labels.end()) {
            out.push_back(la + "_on");
            out.push_back(la + "_off");
            labels.insert(la);
        }
    }
    return out;
}

inline auto get_conformational_interaction_names(const Conformational_interaction_scheme& sch) {
    std::vector<std::string> out;
    std::set<std::string> labels;
    for (auto& e : sch()) {
        auto la = get<Conformational_interaction_label>(e())();
        if (labels.find(la) == labels.end()) {
            out.push_back(la);
            auto n = get<Conformational_interaction_players>(e())().size();
            for (std::size_t i = 0; i < n; ++i) out.push_back(la + "_" + std::to_string(i));
            labels.insert(la);
        }
    }
    return out;
}

inline auto get_conductance_names(const Conductance_interaction_scheme& sch,
                                  const Conductance_interaction_info& int_info) {
    std::vector<std::string> out;
    std::set<std::string> labels;
    for (auto& e : sch()) {
        auto la = get<Conductance_interaction_label>(e())();
        if (labels.find(la) == labels.end()) {
            out.push_back(la);
            labels.insert(la);
        }
    }
    switch (get<Conductance_interaction_type>(int_info())()) {
        case Conductance_interaction_kind::additive:
        case Conductance_interaction_kind::multiplicative:
            break;
        case Conductance_interaction_kind::equilibrium:
            out.push_back(get<Conductance_leakeage_ratio_label>(int_info())());
            break;
    }

    return out;
}

inline auto get_states_structure(const Conformational_states&) {}

}  // namespace impl

template <class Conformational_model_, class Q0_, class Qa_, class g_>
    requires(std::is_same_v<var::untransformed_type_t<Q0_>, Q0> &&
             std::is_same_v<var::untransformed_type_t<Qa_>, Qa> &&
             std::is_same_v<var::untransformed_type_t<g_>, g> &&
             (std::is_same_v<Conformational_model_, Conformational_model> ||
              std::is_same_v<Conformational_model_, Conformational_model_standarized>))
auto get_Parameters_from_Q0_Qa_g(const Conformational_model_& model,
                                 const var::Parameters_Names& names, const Q0_& t_Q0,
                                 const Qa_& t_Qa, const g_& t_g)
    -> Maybe_error<Transfer_Op_to<Q0_, Matrix<double>>> {
    Transfer_Op_to<Q0_, Matrix<double>> out = Matrix<double>(names().size(), 1, 0.0);

    auto N = get<N_St>(model())();
    auto tr = get<Conformational_transition_list>(model());

    auto states = get<Conformational_states>(model());
    auto& t_m_sch = get<Conformational_model_scheme>(model());
    auto inter = get<Conformational_interaction_scheme>(t_m_sch());
    auto cond_inter = get<Conductance_interaction_scheme>(t_m_sch());

    auto& t_std = get<Conformation_change_standard_map>(model());
    auto& t_m_change = get<Conformational_change_scheme>(t_m_sch());

    auto& t_map = get<Conformational_state_count_to_representative_index>(model());

    std::set<Conformational_change_label> conf;

    for (auto i_conf = 0ul; i_conf < t_m_change().size(); ++i_conf) {
        auto ch_lab = t_m_change()[i_conf];
        auto ch_label = get<Conformational_change_label>(ch_lab());
        if (!conf.contains(ch_label)) {
            conf.insert(ch_label);
            auto v_std = t_std[ch_label];
            Conformational_domain_state s0(
                Vector_Space(ch_label, Conformational_change_domain_state(false), v_std));

            auto Maybe_i0 = find_state_with_conformation(states, s0);
            if (!Maybe_i0)
                return Maybe_i0.error();
            auto i0 = Maybe_i0.value();

            auto it1 = std::find_if(
                tr()[i0].begin(), tr()[i0].end(),
                [&ch_label, &v_std](Conformational_transition const& x) {
                    return (get<Conformational_change_label>(x())() == ch_label()) &&
                           (get<Conformational_interactions_transition_state>(x()) == v_std);
                });

            if (it1 == tr()[i0].end()) {
                std::stringstream ss;
                ss << i0 << " has not transition " << v_std << "\n";
                return error_message(ss.str());
            }
            auto i1 = get<Conformational_transition_landing_state_index>((*it1)())()();

            auto it0 = std::find_if(
                tr()[i1].begin(), tr()[i1].end(),
                [&ch_label, &v_std](Conformational_transition const& x) {
                    return (get<Conformational_change_label>(x())() == ch_label()) &&
                           (get<Conformational_interactions_transition_state>(x()) == v_std);
                });

            if (it0 == tr()[i1].end()) {
                std::stringstream ss;
                ss << i1 << " has not retrotransition " << v_std << "\n";
                return error_message(ss.str());
            }

            auto trr0 = *it1;
            auto trr1 = *it0;

            auto ag = get<Agonist_dependency>(trr0());
            auto d0 = get<Conformational_transition_direction>(trr0());

            assert(ch_label() == get<Conformational_change_label>(trr0())());
            auto n0 = get<Conformational_transition_mulitplicity>(trr0());
            auto n1 = get<Conformational_transition_mulitplicity>(trr1());

            auto Maybe_ion = names[ch_label() + "_on"];
            if (!Maybe_ion)
                return Maybe_ion.error();
            auto Maybe_ioff = names[ch_label() + "_off"];
            if (!Maybe_ioff)
                return Maybe_ioff.error();

            auto i_on = Maybe_ion.value();
            auto i_off = Maybe_ioff.value();
            auto kon = t_Q0()(i0, i1) / n0();
            if (ag())
                out[i_on] = t_Qa()(i0, i1) / n0();
            else
                out[i_on] = t_Q0()(i0, i1) / n0();
            out[i_off] = t_Q0()(i1, i0) / n1();
        }
    }
    for (auto i_inter = 0ul; i_inter < inter().size(); ++i_inter) {
        Conformational_interaction inter_i = inter()[i_inter];
        auto inter_lab = get<Conformational_interaction_label>(inter_i());
        auto inter_players = get<Conformational_interaction_players>(inter_i());
        for (auto i_conf_play = 0ul; i_conf_play < inter_players().size(); ++i_conf_play) {
            auto ch_play_lab = inter_players()[i_conf_play];
            auto ij_inter = Vector_Space(Conformational_interaction_index(i_inter),
                                         Conformational_interaction_subposition(i_conf_play));
            auto v_play_std = t_std[ch_play_lab];
            bool negative_interaction = false;

            if (v_play_std[ij_inter] > 0) {
                --v_play_std()[ij_inter];
                if (v_play_std()[ij_inter] == 0) {
                    v_play_std().erase(ij_inter);
                }
                negative_interaction = true;
            } else {
                ++v_play_std()[ij_inter];
                if (v_play_std()[ij_inter] == 0) {
                    v_play_std().erase(ij_inter);
                }
                negative_interaction = false;
            }

            Conformational_domain_state s0(Conformational_domain_state(
                Vector_Space(ch_play_lab, Conformational_change_domain_state(false), v_play_std)));

            auto Maybe_i0 = find_state_with_conformation(states, s0);
            if (!Maybe_i0)
                return Maybe_i0.error();
            auto i0 = Maybe_i0.value();

            auto it1 = std::find_if(
                tr()[i0].begin(), tr()[i0].end(),
                [&ch_play_lab, &v_play_std](Conformational_transition const& x) {
                    return (get<Conformational_change_label>(x())() == ch_play_lab()) &&
                           (get<Conformational_interactions_transition_state>(x()) == v_play_std);
                });

            if (it1 == tr()[i0].end()) {
                std::stringstream ss;
                ss << i0 << " has not transition " << v_play_std << "\n";
                return error_message(ss.str());
            }
            auto i1 = get<Conformational_transition_landing_state_index>((*it1)())()();

            auto it0 = std::find_if(
                tr()[i1].begin(), tr()[i1].end(),
                [&ch_play_lab, &v_play_std](Conformational_transition const& x) {
                    return (get<Conformational_change_label>(x())() == ch_play_lab()) &&
                           (get<Conformational_interactions_transition_state>(x()) == v_play_std);
                });

            if (it0 == tr()[i1].end()) {
                std::stringstream ss;
                ss << i1 << " has not retrotransition " << v_play_std << "\n";
                return error_message(ss.str());
            }

            auto trr0 = *it1;
            auto trr1 = *it0;

            auto ag = get<Agonist_dependency>(trr0());
            auto d0 = get<Conformational_transition_direction>(trr0());

            assert(ch_play_lab() == get<Conformational_change_label>(trr0())());
            auto n0 = get<Conformational_transition_mulitplicity>(trr0());
            auto n1 = get<Conformational_transition_mulitplicity>(trr1());

            auto Maybe_ion = names[ch_play_lab() + "_on"];
            if (!Maybe_ion)
                return Maybe_ion.error();
            auto Maybe_ioff = names[ch_play_lab() + "_off"];
            if (!Maybe_ioff)
                return Maybe_ioff.error();

            auto Maybe_i_Factor = names[inter_lab()];
            auto Maybe_i_Factor_pos = names[inter_lab() + "_" + std::to_string(i_conf_play)];
            if (!Maybe_i_Factor || !Maybe_i_Factor_pos)
                return error_message(Maybe_i_Factor.error()() + Maybe_i_Factor_pos.error()());
            auto i_Factor = Maybe_i_Factor.value();
            auto i_Factor_pos = Maybe_i_Factor_pos.value();

            auto i_on = Maybe_ion.value();
            auto i_off = Maybe_ioff.value();
            auto kon = t_Q0()(i0, i1) / n0();
            if (ag())
                kon = t_Qa()(i0, i1) / n0();
            auto koff = t_Q0()(i1, i0) / n1();

            auto kon0 = out[i_on];
            auto koff0 = out[i_off];

            out[i_Factor_pos] = negative_interaction ? kon0 / kon : kon / kon0;

            if (i_conf_play == 0)
                out[i_Factor] = negative_interaction ? (kon0 / koff0) / (kon / koff)
                                                     : (kon / koff) / (kon0 / koff0);

            assert(out[i_Factor] == negative_interaction ? (kon0 / koff0) / (kon / koff)
                                                         : (kon / koff) / (kon0 / koff0));
        }
    }

    auto& cond_int_info = get<Conductance_interaction_info>(model());

    Conductance_interaction_type mult = get<Conductance_interaction_type>(cond_int_info());
    Conductance_state_count cond_st(std::map<Conductance_interaction_index, int>{});

    std::size_t ii_state = 0;
    bool cond_state_found = false;

    for (auto i_state = 0ul; i_state < states().size(); ++i_state) {
        if (get<Conductance_state_count>(states()[i_state])() == cond_st()) {
            ii_state = i_state;
            cond_state_found = true;
            break;
        }
    }
    if (!cond_state_found)
        return error_message("resting conductance state not found, weird!!");

    auto g0 = t_g()[ii_state];

    assert((mult() != Conductance_interaction_kind::additive) || (g0 == 0.0));

    for (auto i_cond_inter = 0ul; i_cond_inter < cond_inter().size(); ++i_cond_inter) {
        Conductance_interaction cond_inter_i = cond_inter()[i_cond_inter];
        auto cond_inter_lab = get<Conductance_interaction_label>(cond_inter_i());

        auto Maybe_i_lab_cond_inter = names[cond_inter_lab()];
        if (!Maybe_i_lab_cond_inter)
            return Maybe_i_lab_cond_inter.error();
        auto i_lab_cond_inter = Maybe_i_lab_cond_inter.value();

        Conductance_state_count cond_st(std::map<Conductance_interaction_index, int>{
            {Conductance_interaction_index(i_cond_inter), 1}});

        std::size_t ii_state = 0;
        bool cond_state_found = false;

        for (auto i_state = 0ul; i_state < states().size(); ++i_state) {
            if (get<Conductance_state_count>(states()[i_state])() == cond_st()) {
                ii_state = i_state;
                cond_state_found = true;
                break;
            }
        }
        if (!cond_state_found)
            return error_message("resting conductance state not found, weird!!");

        auto g1 = t_g()[ii_state];

        switch (mult()) {
            case Conductance_interaction_kind::additive:
                out[i_lab_cond_inter] = g1;
                break;
            case Conductance_interaction_kind::multiplicative:
                out[i_lab_cond_inter] = g1 / g0;
                break;
            case Conductance_interaction_kind::equilibrium: {
                auto le_lab = get<Conductance_leakeage_ratio_label>(cond_int_info());
                auto Maybe_ile = names[le_lab()];
                if (!Maybe_ile)
                    return Maybe_ile.error();
                auto ile = Maybe_ile.value();
                Conductance_state_count cond_st_2(std::map<Conductance_interaction_index, int>{
                    {Conductance_interaction_index(i_cond_inter), 2}});

                std::size_t ii_state_2 = 0;
                bool cond_state_found_2 = false;

                for (auto i_state = 0ul; i_state < states().size(); ++i_state) {
                    if (get<Conductance_state_count>(states()[i_state])() == cond_st_2()) {
                        ii_state_2 = i_state;
                        cond_state_found_2 = true;
                        break;
                    }
                }
                if (!cond_state_found_2)
                    return error_message("resting conductance state not found, weird!!");

                auto g2 = t_g()[ii_state_2];

                auto gd1 = g1 / g0;
                auto gd2 = g2 / g1;

                auto a = gd2 - 1.0;
                auto b = -gd1 * (gd2 - 1.0) - gd1 + 1.0;
                auto c = gd2 * (gd1 - 1.0);
                using std::sqrt;
                auto beta = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
                auto alpha = (beta - gd2) / ((gd2 - 1.0) * beta * beta);
                out[ile] = alpha;
                out[i_lab_cond_inter] = beta;

            } break;
        };
    }

    return out;
}

inline Maybe_error<Conformational_model> make_Conformational_model(
    Agonist_dependency_map&& t_agonist_map, std::vector<Conformational_change_label>&& t_scheme,
    std::vector<Conformational_interaction>&& t_interactions,
    std::vector<Conductance_interaction>&& t_conductance) {
    auto Maybe_scheme =
        impl::make_Conformational_model_scheme(std::move(t_agonist_map), std::move(t_scheme),
                                               std::move(t_interactions), std::move(t_conductance));
    if (!Maybe_scheme)
        return Maybe_scheme.error();
    else {
        auto& scheme = Maybe_scheme.value();
        auto Maybe_tuple = impl::make_Conformational_states_and_index(scheme);
        if (!Maybe_tuple)
            return Maybe_tuple.error();
        else {
            auto [states, map] = std::move(Maybe_tuple.value());
            auto Maybe_transition = impl::make_Conformational_transition_list(scheme, states, map);
            if (!Maybe_transition)
                return Maybe_transition.error();
            else {
                return Conformational_model(Vector_Space(N_St(states().size()), std::move(scheme),
                                                         std::move(states),
                                                         std::move(Maybe_transition.value())));
            }
        }
    }
}

inline Maybe_error<Conformational_model_standarized> make_Conformational_model_standarized(
    Agonist_dependency_map&& t_agonist_map, std::vector<Conformational_change_label>&& t_scheme,
    std::vector<Conformational_interaction>&& t_interactions,
    std::vector<Conductance_interaction>&& t_conductance,
    std::map<Conformational_change_label, Conformation_change_standard_state>&& t_standard_states,
    Conductance_interaction_info&& t_mult_cond) {
    auto Maybe_scheme =
        impl::make_Conformational_model_scheme(std::move(t_agonist_map), std::move(t_scheme),
                                               std::move(t_interactions), std::move(t_conductance));
    if (!Maybe_scheme)
        return Maybe_scheme.error();
    auto& scheme = Maybe_scheme.value();
    auto Maybe_tuple = impl::make_Conformational_states_and_index(scheme);
    if (!Maybe_tuple)
        return Maybe_tuple.error();
    auto [states, map] = std::move(Maybe_tuple.value());
    auto Maybe_transition = impl::make_Conformational_transition_list(scheme, states, map);
    if (!Maybe_transition) {
        return Maybe_transition.error();
    }

    return Conformational_model_standarized(Vector_Space(
        N_St(states().size()), std::move(scheme), std::move(states), std::move(map),
        std::move(Maybe_transition.value()),
        Conformation_change_standard_map{std::move(t_standard_states)}, std::move(t_mult_cond)));
}

inline auto get_states_structure(const Conformational_model& model) {
    return impl::get_states_structure(get<Conformational_states>(model()));
}

template <class Id, class Conformational_model_, class P>
    requires(std::is_same_v<var::untransformed_type_t<P>, var::Parameters_values> &&
             (std::is_same_v<Conformational_model_, Conformational_model> ||
              std::is_same_v<Conformational_model_, Conformational_model_standarized>))
auto make_Model(const Conformational_model_& model, const typename var::Parameters_Names& names,
                const P& p)
    -> Maybe_error<std::tuple<Transfer_Op_to<P, Q0>, Transfer_Op_to<P, Qa>, Transfer_Op_to<P, g>>>

{
    auto Maybe_Q0Qa = impl::make_Q0_Qa<Id>(model, names, p);
    auto Maybe_g = impl::make_g<Id>(model, names, p);
    if (!Maybe_Q0Qa || !Maybe_g)
        return error_message(Maybe_Q0Qa.error()() + Maybe_g.error()());
    else {
        auto [v_Q0, v_Qa] = std::move(Maybe_Q0Qa.value());
        auto v_g = std::move(Maybe_g.value());

        return std::tuple(std::move(v_Q0), std::move(v_Qa), std::move(v_g));
    }
}

template <class Id, class Conformational_model_>
    requires(std::is_same_v<Conformational_model_, Conformational_model> ||
             std::is_same_v<Conformational_model_, Conformational_model_standarized>)
auto make_Model_Formulas(const Conformational_model_& model, const var::Parameters_Names& names)
    -> Maybe_error<std::tuple<Q0_formula, Qa_formula, g_formula>>

{
    auto Maybe_Q0Qa = impl::make_Q0_Qa_formula<Id>(model, names);
    auto Maybe_g = impl::make_g_formula<Id>(model, names);
    if (!Maybe_Q0Qa || !Maybe_g)
        return error_message(Maybe_Q0Qa.error()() + Maybe_g.error()());
    else {
        auto [v_Q0, v_Qa] = std::move(Maybe_Q0Qa.value());
        auto v_g = std::move(Maybe_g.value());
        return std::tuple(std::move(v_Q0), std::move(v_Qa), std::move(v_g));
    }
}

template <class Id, class Conformational_model_>
    requires(std::is_same_v<Conformational_model_, Conformational_model> ||
             std::is_same_v<Conformational_model_, Conformational_model_standarized>)
auto make_ModelNames(const Conformational_model_& confmodel) {
    auto& model = get<Conformational_model_scheme>(confmodel());
    auto con_names =
        impl::get_conformational_change_names(get<Conformational_change_scheme>(model()));
    auto inter_names =
        impl::get_conformational_interaction_names(get<Conformational_interaction_scheme>(model()));
    auto cond_names = impl::get_conductance_names(get<Conductance_interaction_scheme>(model()),
                                                  get<Conductance_interaction_info>(confmodel()));

    con_names.insert(con_names.end(), inter_names.begin(), inter_names.end());
    con_names.insert(con_names.end(), cond_names.begin(), cond_names.end());

    return var::Parameters_Names(std::move(con_names));
}

///**
// * @brief The Allosteric_Model class build a kinetic rate model starting with
// a
// * set of conformational changes and their interactions.
// *
// *
// */
// class Allosteric_Model {
// public:
//    typedef Allosteric_Model self_type;
//    typedef Model_Parameter_label myParameter_label;
//    struct transitions {
//        bool on;
//        bool agonist;
//        std::string conformation;
//        std::map<std::vector<std::pair<std::string, std::string>>,
//        std::size_t>
//            coupling;
//    };
//    struct model_definition {
//        std::vector<std::string> conformational_changes;
//        std::map<std::pair<std::string, bool>, std::string>
//            conformational_changes_names;
//        std::set<std::size_t> agonist_changes;
//        std::set<std::size_t> conductance_changes;
//        std::map<std::size_t, std::string> conductance_names;
//        std::multimap<std::size_t, std::pair<std::set<std::size_t>,
//                                             std::pair<std::string,
//                                             std::string>>>
//            conformational_inter_unit_cell;
//    };

//    struct new_model_definition {
//        std::size_t number_of_units;
//        std::map<Conformational_change_label, Conformational_change>
//            conformational_changes;
//        std::vector<Conformational_change_label>
//        unit_of_conformational_changes; std::set<Conformational_interaction>
//        conformational_interactions; std::map<std::size_t,
//        Conductance_Parameter_label> conductance_names;
//    };

//    model_definition new_to_old(const new_model_definition &m) {
//        model_definition out;
//        for (auto &e : m.conductance_names)
//            out.conductance_names[e.first] = e.second;
//        auto k = m.unit_of_conformational_changes.size();
//        std::size_t n = m.number_of_units * k;
//        out.conformational_changes.resize(n);
//        for (std::size_t i = 0; i < m.number_of_units; ++i)
//            for (std::size_t j = 0; j <
//            m.unit_of_conformational_changes.size();
//                 ++j) {
//                auto ii = i * k + j;
//                auto &cf = m.conformational_changes.at(
//                    m.unit_of_conformational_changes[j].name());
//                out.conformational_changes[ii] = cf.label();
//                if (cf.change_in_agonist() != 0)
//                    out.agonist_changes.insert(ii);
//                if (cf.change_in_conductance() != 0)
//                    out.conductance_changes.insert(ii);
//            }
//        for (auto &e : m.conformational_changes) {
//            out.conformational_changes_names[std::pair(e.second.label(),
//            true)] =
//                e.second.par_on();
//            out.conformational_changes_names[std::pair(e.second.label(),
//            false)] =
//                e.second.par_off();
//        }
//        for (auto &e : m.conformational_interactions) {
//            std::vector<std::size_t> cc;
//            std::size_t current_i = 0;
//            for (std::size_t j = 0; j <
//            e.interacting_conformational_changes().size();
//                 ++j) {
//                auto j_n = m.conformational_changes.at(
//                    e.interacting_conformational_changes()[j]);
//                std::size_t i = current_i;
//                while (j_n.label().name() != out.conformational_changes[i])
//                    ++i;
//                cc.push_back(i);
//                ++current_i;
//            }

//            for (std::size_t i = 0; i < cc.size(); ++i) {
//                auto x = cc[i];
//                auto ix = x % k;
//                auto nx = x / k;
//                auto shift = (n - nx) * k;
//                std::set<std::size_t> s;
//                for (std::size_t j = 0; j < cc.size(); ++j)
//                    if (j != i)
//                        s.insert(cc[j]);
//                s = rotate(s, n * k, shift);
//                out.conformational_inter_unit_cell.emplace(
//                    ix, std::pair(
//                        s, std::pair(e.factor_label(),
//                        e.coefficient_labels()[i])));
//            }
//        }
//        return out;
//    }

//    constexpr static auto className = my_static_string("Allosteric_Model");

// protected:
//     new_model_definition new_d_;
//     std::vector<std::vector<std::string>> conformer_;
//     std::vector<std::map<std::size_t, transitions>> transitions_;
//     std::set<std::string> paramNames_;
//     std::vector<std::string> conductances_;
//     model_definition d_;

//    /**
//   * @brief getParameterNamesFrom
//   * @param conformational_changes_names
//   * @param conductance_names
//   * @param conformational_inter
//   */
//    static auto getParameterNamesFrom(
//        const std::map<std::pair<std::string, bool>, std::string>
//            conformational_changes_names,
//        const std::map<std::size_t, std::string> &conductance_names,
//        const std::multimap<
//            std::size_t,
//            std::pair<std::set<std::size_t>, std::pair<std::string,
//            std::string>>> &conformational_inter) {
//        std::set<std::string> out;
//        for (auto &e : conformational_changes_names)
//            out.insert(e.second);
//        for (auto &e :s conductance_names)
//            out.insert(e.second);
//        for (auto &e : conformational_inter) {
//            out.insert(e.second.second.first);
//            out.insert(e.second.second.second);
//        }
//        return out;
//    }

//    static auto
//    getConformers(const std::vector<std::string> &conformational_changes,
//                  std::size_t p) {
//        std::vector<std::vector<std::string>> conformer;
//        std::map<std::vector<std::string>, std::size_t> state_to_conformer;
//        for (std::size_t i = 0; i < (1u << conformational_changes.size());
//        ++i) {
//            auto c = index_to_conformational_change(conformational_changes,
//            i); if (state_to_conformer.find(c) == state_to_conformer.end()) {
//                state_to_conformer[c] = conformer.size();
//                std::size_t n = 1;
//                auto permute = rotate(c, p * n);
//                while (c != permute) {
//                    state_to_conformer[permute] = conformer.size();
//                    ++n;
//                    permute = rotate(c, p * n);
//                }
//                conformer.push_back(c);
//            }
//        }
//        return std::make_tuple(conformer, state_to_conformer);
//    }

//    static auto getTransitions(
//        const std::vector<std::string> &conformational_changes,
//        const std::vector<std::vector<std::string>> &conformer,
//        const std::map<std::vector<std::string>, std::size_t>
//        &state_to_conformer, const std::set<std::size_t> &agonist_changes,
//        const std::map<std::pair<std::string, bool>, std::string>
//            conformational_changes_names,
//        std::multimap<std::size_t, std::pair<std::set<std::size_t>,
//                                             std::pair<std::string,
//                                             std::string>>>
//            conformational_interactions) {
//        std::vector<std::map<std::size_t, transitions>> transitions_out;
//        for (std::size_t i = 0; i < conformer.size(); ++i) {
//            std::map<std::size_t, transitions> myTransition;
//            auto c = conformer[i];
//            for (std::size_t k = 0; k < c.size(); ++k) {
//                auto change = c;
//                if (c[k] == "")
//                    change[k] = conformational_changes[k];
//                else
//                    change[k] = "";
//                auto j = state_to_conformer.at(change);
//                if (myTransition.find(j) == myTransition.end()) {

//                    myTransition[j].agonist =
//                        agonist_changes.find(k) != agonist_changes.end();
//                    bool on = change[k] == conformational_changes[k];
//                    myTransition[j].on = on;
//                    auto name = conformational_changes_names.at(
//                        std::pair{conformational_changes[k], on});
//                    myTransition[j].conformation = name;
//                }
//                auto m = conformational_interactions.equal_range(k);
//                std::vector<std::pair<std::string, std::string>> coupling;
//                for (auto it = m.first; it != m.second; ++it) {
//                    std::set<std::size_t> s = it->second.first;
//                    bool all = true;
//                    for (auto e : s) {
//                        if (c[e].empty()) {
//                            all = false;
//                            break;
//                        }
//                    }
//                    if (all)
//                        coupling.push_back(it->second.second);
//                }
//                ++myTransition[j].coupling[coupling];
//            }
//            transitions_out.push_back(std::move(myTransition));
//        }
//        return transitions_out;
//    }

//    static

//        std::string
//        g_name_of_conformer(
//            const std::vector<std::string> &c,
//            const std::set<std::size_t> &conductance_changes,
//            const std::map<std::size_t, std::string> &conductance_names) {
//        std::size_t i = 0;
//        for (auto e : conductance_changes)
//            if (!c.at(e).empty())
//                ++i;
//        return conductance_names.at(i);
//    }

//    static std::vector<std::string>
//    get_g_names(const std::vector<std::vector<std::string>> &conformers,
//                const std::set<std::size_t> &conductance_changes,
//                const std::map<std::size_t, std::string> &conductance_names) {
//        std::vector<std::string> out(conformers.size());
//        for (std::size_t i = 0; i < conformers.size(); ++i)
//            out[i] = g_name_of_conformer(conformers[i], conductance_changes,
//                                         conductance_names);
//        return out;
//    }

// public:
//     static auto get_constructor_fields_old() {
//         return std::make_tuple(
//             grammar::field(C<self_type>{}, "conformational_changes",
//                            &self_type::get_conformational_changes),
//             grammar::field(C<self_type>{}, "conformational_changes_names",
//                            &self_type::get_conformational_changes_names),
//             grammar::field(C<self_type>{}, "agonist_changes",
//                            &self_type::get_agonist_changes),
//             grammar::field(C<self_type>{}, "conductance_changes",
//                            &self_type::get_conductance_changes),
//             grammar::field(C<self_type>{}, "conductance_names",
//                            &self_type::get_conductance_names),
//             grammar::field(C<self_type>{}, "conformational_inter_unit_cell",
//                            &self_type::get_conformational_inter_unit_cell));
//     }

//    const std::vector<std::string> &get_conformational_changes() const {
//        return d_.conformational_changes;
//    }
//    const std::map<std::pair<std::string, bool>, std::string> &
//    get_conformational_changes_names() const {
//        return d_.conformational_changes_names;
//    }
//    const std::set<std::size_t> &get_agonist_changes() const {
//        return d_.agonist_changes;
//    }
//    const std::set<std::size_t> &get_conductance_changes() const {
//        return d_.conductance_changes;
//    }
//    const std::map<std::size_t, std::string> &get_conductance_names() const {
//        return d_.conductance_names;
//    }
//    const std::multimap<
//        std::size_t,
//        std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
//        &
//    get_conformational_inter_unit_cell() const {
//        return d_.conformational_inter_unit_cell;
//    }

//    std::size_t number_of_units() const { return new_d_.number_of_units; }
//    std::map<Conformational_change_label, Conformational_change> const &
//    conformational_changes() const {
//        return new_d_.conformational_changes;
//    }
//    std::vector<Conformational_change_label> const &
//    unit_of_conformational_changes() const {
//        return new_d_.unit_of_conformational_changes;
//    }
//    std::set<Conformational_interaction> const &
//    conformational_interactions() const {
//        return new_d_.conformational_interactions;
//    }
//    std::map<std::size_t, Conductance_Parameter_label> const &
//    conductance_names() const {
//        return new_d_.conductance_names;
//    }

//    static auto get_constructor_fields() {
//        return std::make_tuple(
//            grammar::field(C<self_type>{}, "number_of_units",
//                           &self_type::number_of_units),
//            grammar::field(C<self_type>{}, "conformational_changes",
//                           &self_type::conformational_changes),
//            grammar::field(C<self_type>{}, "unit_of_conformational_changes",
//                           &self_type::unit_of_conformational_changes),
//            grammar::field(C<self_type>{}, "conformational_interactions",
//                           &self_type::conformational_interactions),
//            grammar::field(C<self_type>{}, "conductance_names",
//                           &self_type::conductance_names));
//    }

//    Allosteric_Model() = default;

//    Allosteric_Model(
//        std::size_t number_of_units,
//        std::map<Conformational_change_label, Conformational_change>
//            conformational_changes,
//        std::vector<Conformational_change_label>
//        unit_of_conformational_changes, std::set<Conformational_interaction>
//        conformational_interactions, std::map<std::size_t,
//        Conductance_Parameter_label> conductance_names) :
//        new_d_{number_of_units, conformational_changes,
//                 unit_of_conformational_changes, conformational_interactions,
//                 conductance_names},
//        d_{new_to_old(new_d_)} {
//        init();
//    }

//    Allosteric_Model(
//        const std::vector<std::string> &conformational_changes,
//        const std::map<std::pair<std::string, bool>, std::string>
//            conformational_changes_names,
//        const std::set<std::size_t> &agonist_changes,
//        const std::set<std::size_t> &conductance_changes,
//        const std::map<std::size_t, std::string> &conductance_names,
//        const std::multimap<
//            std::size_t,
//            std::pair<std::set<std::size_t>, std::pair<std::string,
//            std::string>>> &conformational_inter_unit_cell)
//        : d_{conformational_changes, conformational_changes_names,
//             agonist_changes,        conductance_changes,
//             conductance_names,      conformational_inter_unit_cell} {
//        init();
//    }

//    template <class Parameters> auto Qs(const Parameters &p) const {

//        M_Matrix<double> Q0(conformer_.size(), conformer_.size());
//        M_Matrix<double> Qa(conformer_.size(), conformer_.size());
//        for (std::size_t i = 0; i < Q0.nrows(); ++i) {
//            Q0(i, i) = 0;
//            Qa(i, i) = 0;
//            for (auto it = transitions_[i].begin(); it !=
//            transitions_[i].end();
//                 ++it) {
//                std::size_t j = it->first;
//                if ((it->second.agonist) && (it->second.on)) {
//                    Qa(i, j) = rate(it->second, p);
//                    Qa(i, i) -= Qa(i, j);
//                } else {
//                    Q0(i, j) = rate(it->second, p);
//                    Q0(i, i) -= Q0(i, j);
//                }
//            }
//        }
//        auto r = Q0 * ones<double>(Q0.ncols(), 1);
//        auto q = Qa * ones<double>(Q0.ncols(), 1);

//        return std::pair(Q0, Qa);
//    }

//    template <class Parameters> auto g(const Parameters &p) const {
//        M_Matrix<double> out(conformer_.size(), 1);
//        for (std::size_t i = 0; i < conformer_.size(); ++i)
//            out(i, 0) = p.at(conductances_.at(i));
//        return out;
//    }

//    auto getParameterNames() const { return paramNames_; }

// private:
//     void init() {
//         paramNames_ = getParameterNamesFrom(d_.conformational_changes_names,
//                                             d_.conductance_names,
//                                             d_.conformational_inter_unit_cell);
//         std::cout << paramNames_;
//         auto p = periodicity(d_.conformational_changes);

//        auto conformational_interactions = fill_conformational_interactions(
//            d_.conformational_inter_unit_cell,
//            d_.conformational_changes.size(), p);

//        std::map<std::vector<std::string>, std::size_t> state_to_conformer;

//        std::tie(conformer_, state_to_conformer) =
//            getConformers(d_.conformational_changes, p);

//        transitions_ = getTransitions(d_.conformational_changes, conformer_,
//                                      state_to_conformer, d_.agonist_changes,
//                                      d_.conformational_changes_names,
//                                      conformational_interactions);
//        conductances_ =
//            get_g_names(conformer_, d_.conductance_changes,
//            d_.conductance_names);
//    }

//    static std::vector<std::string>
//    index_to_conformational_change(const std::vector<std::string> &cc,
//                                   std::size_t index) {
//        std::vector<std::string> out(cc.size(), "");
//        for (std::size_t i = 0; i < cc.size(); ++i)
//            if (((index >> i) & 1) == 1)
//                out[i] = cc[i];
//        return out;
//    }

//    static std::vector<std::string> rotate(const std::vector<std::string> &x,
//                                           std::size_t n) {
//        std::vector<std::string> out(x.size());
//        for (std::size_t i = 0; i < out.size(); ++i)
//            out[(i + n) % out.size()] = x[i];
//        return out;
//    }

//    static std::size_t periodicity(const std::vector<std::string> &x) {
//        std::size_t p = 1;
//        auto t = rotate(x, p);
//        while (t != x && p < x.size()) {
//            ++p;
//            t = rotate(x, p);
//        }
//        return p;
//    }

//    static std::set<std::size_t> rotate(const std::set<std::size_t> &c,
//                                        std::size_t n, std::size_t i) {
//        std::set<std::size_t> out;
//        for (auto e : c)
//            out.insert((e + i) % n);
//        return out;
//    }

//    static std::multimap<
//        std::size_t,
//        std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
//    fill_conformational_interactions(
//        const std::multimap<
//            std::size_t,
//            std::pair<std::set<std::size_t>, std::pair<std::string,
//            std::string>>> &conformational_interactions,
//        std::size_t n, std::size_t p) {
//        std::multimap<std::size_t, std::pair<std::set<std::size_t>,
//                                             std::pair<std::string,
//                                             std::string>>>
//            out;
//        for (std::size_t i = 0; i < n / p; ++i) {
//            for (auto &e : conformational_interactions) {
//                out.insert({(e.first + i * p) % n,
//                            {rotate(e.second.first, n, i * p),
//                            e.second.second}});
//            }
//        }
//        return out;
//    }

//    template <class P> static double rate(const transitions &tr, const P &p) {
//        double out = 0;
//        if (tr.on) {
//            for (auto &e : tr.coupling) {
//                double b = e.second;
//                for (auto &e2 : e.first)
//                    b *= std::pow(p.at(e2.first), p.at(e2.second));
//                out += b;
//            }
//        } else {
//            for (auto &e : tr.coupling) {
//                double b = e.second;
//                for (auto &e2 : e.first)
//                    b *= std::pow(p.at(e2.first), p.at(e2.second) - 1.0);
//                out += b;
//            }
//        }
//        return out * p.at(tr.conformation);
//    }
//};

}  // namespace macrodr

#endif  // ALLOSTERIC_MODELS_H
