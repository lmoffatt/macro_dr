#ifndef ALLOSTERIC_MODELS_H
#define ALLOSTERIC_MODELS_H
#include "qmodel.h"
#include <optional>
#include <map>
namespace  macrodr
{

struct Conformational_change_label: public Constant<Conformational_change_label, std::string>{};

struct Agonist_label: public Constant<Agonist_label,std::string>{};

struct Agonist_dependency: public Constant<Agonist_dependency, std::optional<Agonist_label>>{};

struct Agonist_dependency_map: public Constant<Agonist_dependency_map, std::map<Conformational_change_label,Agonist_dependency>>{};

struct Conformational_change: public Constant<Conformational_change,Vector_Space<Conformational_change_label,Agonist_dependency>>{};

struct Conformational_change_domain_state: public Constant<Conformational_change_domain_state, bool>{};

struct Conformational_change_state_vector: public Constant<Conformational_change_state_vector, std::vector<Conformational_change_domain_state>>{};


struct Conformational_change_scheme: public Constant<Conformational_change_scheme,std::vector<Conformational_change>>{};

struct Conformational_position: public Constant<Conformational_position, std::size_t>{};


struct Conformational_interaction_positions: public Constant<Conformational_interaction_positions, std::vector<std::vector<Conformational_position>>>{};

struct Conformational_interaction_label: public Constant<Conformational_interaction_label,
                                                            std::string>{};


struct Conformational_interaction_players: public Constant<Conformational_interaction_players,
                                                            std::vector<Conformational_change_label>>{};

struct Conformational_interaction: public Constant<Conformational_interaction,
                                                    Vector_Space<Conformational_interaction_label,Conformational_interaction_players,Conformational_interaction_positions>>{};

struct Conformational_interaction_scheme: public Constant<Conformational_interaction_scheme,std::vector<Conformational_interaction>>{};

struct Conformational_interaction_index: public Constant<Conformational_interaction_index,std::size_t>{};

struct Conformational_interaction_subposition: public Constant<Conformational_interaction_subposition,std::size_t>{};

struct Conformational_interactions_domain_state: public Constant<Conformational_interactions_domain_state,std::vector<Vector_Space<Conformational_interaction_index,Conformational_interaction_subposition>>>{};

struct Conformational_domain_state: public Constant<Conformational_domain_state, std::pair<Conformational_change, Conformational_interactions_domain_state>>{};

struct Conformational_interactions_state_vector: public Constant<Conformational_interactions_state_vector,std::vector<Conformational_interactions_domain_state>>{};

struct Conformational_state_vector: public Constant<Conformational_state_vector,Vector_Space<Conformational_change_state_vector,Conformational_interactions_state_vector>>{};


struct Conformational_state_count: public Constant<Conformational_state_count,std::map<Conformational_domain_state, std::size_t>>{};





struct Conformational_model_scheme: public Constant<Conformational_model_scheme,Vector_Space<Conformational_change_scheme, Conformational_interaction_scheme>>{};

struct Conformational_transition_direction: public Constant<Conformational_transition_direction,bool>{};


struct Conformational_transition_mulitplicity: public Constant<Conformational_transition_mulitplicity, std::size_t>{};

struct Conformational_state_index: public Constant<Conformational_state_index, std::size_t>{};

struct Conformational_transition_initiating_state_index: public Constant<Conformational_transition_initiating_state_index, Conformational_state_index>{};
struct Conformational_transition_landing_state_index: public Constant<Conformational_transition_landing_state_index, Conformational_state_index>{};




struct Conformational_transition: public Constant<Conformational_transition,
                                                   Vector_Space<Conformational_transition_initiating_state_index,
                                                                Conformational_transition_landing_state_index,
                                                                Agonist_dependency,
                                                                Conformational_transition_direction,
                                                                Conformational_change_label,
                                                                Conformational_interactions_domain_state,
                                                                Conformational_transition_mulitplicity>> {};


struct Conformational_transition_list: public Constant<Conformational_transition_list, std::vector<std::vector<Conformational_transition>>>{};


struct Conformational_states: public Constant<Conformational_states,
    std::vector<Vector_Space<Conformational_state_count,
                                                  Conformational_state_vector>>>{};

struct Conformational_state_count_to_index: public Constant<Conformational_state_count_to_index, std::map<Conformational_state_count, Conformational_state_index>>{};


struct Conformational_model: public Constant<Conformational_model,Vector_Space<N_St,Conformational_model_scheme, Conformational_states, Conformational_transition_list>>{};


namespace impl {



auto make_Conformational_change_state_vector(const Conformational_model_scheme& model, std::size_t n)
{
    auto number_units= get<Conformational_change_scheme>(model())().size();

    std::vector<Conformational_change_domain_state> out(number_units);
    for (std::size_t i=0; i<number_units; ++i)
        out[i]=Conformational_change_domain_state((n&(1ul<<i))==(1ul<<i));
    return Conformational_change_state_vector(std::move(out));   
}

Conformational_interactions_state_vector make_Conformational_interaction_state_vector(const Conformational_model_scheme& model, const Conformational_change_state_vector& state)
{
    auto number_units= get<Conformational_change_scheme>(model())().size();
    auto v_inter=get<Conformational_interaction_scheme>(model());
    
    auto number_inter=v_inter().size();
    std::vector<Conformational_interactions_domain_state> out(number_units);
    for (std::size_t i=0; i<number_units; ++i)
    {
        for (std::size_t j=0; j<number_inter; ++j)
        {
            auto v_ipos=get<Conformational_interaction_positions>(v_inter()[j]());
            for (std::size_t k=0; k<v_ipos().size(); ++k)
            {
                bool includes_i=false;
                bool includes_all_other=true;
                std::size_t i_sub_position;
                for (std::size_t kk=0; kk<v_ipos()[k].size(); ++kk)
                {
                    if (v_ipos()[k][kk]()==i)
                    {
                        includes_i=true;
                        i_sub_position=kk;
                    }
                    else
                    {
                        includes_all_other=includes_all_other&&(state()[v_ipos()[k][kk]()]());
                    }
                }
                if(includes_all_other&& includes_i)
                    out[i]().push_back(Vector_Space(Conformational_interaction_index(j),Conformational_interaction_subposition(i_sub_position)));
            }
        }
         
    }
    return Conformational_interactions_state_vector(std::move(out));   
}


auto make_Conformational_state_vector(const Conformational_model_scheme& model, std::size_t n)
{
    auto v_change=make_Conformational_change_state_vector(model,n);
    auto v_inter=make_Conformational_interaction_state_vector(model,v_change);
    return Conformational_state_vector (Vector_Space(std::move(v_change), std::move(v_inter)));
    
    
}



Maybe_error<Conformational_state_count> to_state_count(const Conformational_model_scheme& model, const Conformational_state_vector& state)
{
    if (get<Conformational_change_state_vector>(state())().size()!=get<Conformational_change_scheme>(model())().size() )
        return error_message("unequal sizes");
    else
    {
        auto& v_state=get<Conformational_change_state_vector>(state());
        auto& v_inter=get<Conformational_interactions_state_vector>(state());
        auto& v_change=get<Conformational_change_scheme>(model());
        std::map<Conformational_domain_state, std::size_t> out;
        
        for (std::size_t i=0; i<v_state().size(); ++i)
        {
            if (v_state()[i]())
                ++out[Conformational_domain_state(std::pair(v_change()[i], v_inter()[i]))];
        }
        return Conformational_state_count(std::move(out)); 
    }
}

Maybe_error<Conformational_state_index> to_state_index(const Conformational_state_count_to_index& model, const Conformational_state_count& state_count)
{
    if (auto it=model().find(state_count); it==model().end())
        return error_message("state count not found");
    else
        return it->second; 
}


auto change_conformation(const Conformational_change_state_vector& state, std::size_t ith_domain)
{
    Conformational_change_state_vector out(state);
    out()[ith_domain]()=!out()[ith_domain]();
    return out; 
}

auto change_conformation(const Conformational_model_scheme& model,const Conformational_state_vector& state, std::size_t ith_domain)
{
    auto v_change=change_conformation(get<Conformational_change_state_vector>(state()), ith_domain);
    auto v_inter= make_Conformational_interaction_state_vector(model,v_change);
    return Conformational_state_vector(Vector_Space(std::move(v_change), std::move(v_inter))); 
}


Maybe_error<std::tuple<Conformational_states,Conformational_state_count_to_index>> make_Conformational_states_and_index(const Conformational_model_scheme& model){
    auto number_units= get<Conformational_change_scheme>(model())().size();
    Conformational_states out;
    Conformational_state_count_to_index map;
    for (std::size_t n = 0; n < (1ul << number_units); ++n) {
        auto v_state_vector=make_Conformational_state_vector(model,n);
        auto v_state_count=to_state_count(model,v_state_vector);
        if (!v_state_count){
            return v_state_count.error();}
        else if (map().find(v_state_count.value())==map().end())
        {
            map()[v_state_count.value()]=Conformational_state_index(out().size());
            out().push_back(Vector_Space(std::move(v_state_count.value()), std::move(v_state_vector)));
        }
    }
    return std::tuple(std::move(out),std::move(map));
}


Maybe_error<Conformational_transition_list> make_Conformational_transition_list(const Conformational_model_scheme& model, const Conformational_states states, const Conformational_state_count_to_index& map){
    std::vector<std::vector<Conformational_transition>> out(states().size());
    
    auto scheme=get<Conformational_change_scheme>(model());
    
    for (std::size_t i=0; i<states().size(); ++i)
    {
      //  auto v_count=get<Conformational_state_count>(states()[i]);
        auto v_state=get<Conformational_state_vector>(states()[i]);
        auto v_change=get<Conformational_change_state_vector>(v_state());
        auto v_interactions=get<Conformational_interactions_state_vector>(v_state());
        auto i_start=Conformational_transition_initiating_state_index(Conformational_state_index(i));
        
        std::vector<Conformational_transition> tran;
        
        for (std::size_t j=0; j<v_change().size(); ++j)
        {
            auto j_change=scheme()[j];
            auto j_inter=v_interactions()[j];
            auto j_state=change_conformation(model,v_state,j);
            auto Maybe_j_count=to_state_count(model,j_state);
            if (!Maybe_j_count) return Maybe_j_count.error();            
            else
            {
                auto& j_count=Maybe_j_count.value();
                auto Maybe_j_index= to_state_index(map,j_count);
                if (!Maybe_j_index) return Maybe_j_index.error();
                else
                {
                    auto i_end=Conformational_transition_landing_state_index(Maybe_j_index.value());
                    bool already=false;
                    for (auto& elem:tran)
                    {
                        if (get<Conformational_transition_landing_state_index>(elem())()()==i_end()())
                        {
                            ++get<Conformational_transition_mulitplicity>(elem())();
                            already=true;
                            
                        }
                    }
                    if (!already)
                    {
                            auto a=get<Agonist_dependency>(j_change());
                            auto d=Conformational_transition_direction(!v_change()[j]());
                            auto l=get<Conformational_change_label>(j_change());
                            auto m=Conformational_transition_mulitplicity(1);
                            tran.push_back(Conformational_transition(Vector_Space(i_start,i_end,a,d,l,j_inter,m)));
                    }    
                    
                }
            }
        }
        out[i]=tran;
    }
    
    return Conformational_transition_list(std::move(out));
    
}




Maybe_error<Conformational_change_scheme> make_Conformational_change_scheme(Agonist_dependency_map&& t_agonist_map, std::vector<Conformational_change_label>&& t_scheme )
{
    std::vector<Conformational_change> out(t_scheme.size());
    for (std::size_t i=0; i<t_scheme.size(); ++i)
    {
        if (auto it=t_agonist_map().find(t_scheme[i]); it!=t_agonist_map().end())
        {
            out[i]=Conformational_change(Vector_Space(t_scheme[i],it->second));
        }
        else
        {
            return error_message("Conformational_change "+t_scheme[i]()+" agonist dependency condition unknown");
        }
        
    }
    return Conformational_change_scheme(std::move(out));
}



Maybe_error<bool> check_Conformational_interaction(const Conformational_change_scheme& t_scheme,
                                                   Conformational_interaction const& t_interaction){
    
    auto& v_players=get<Conformational_interaction_players>(t_interaction());
    auto& v_positions=get<Conformational_interaction_positions>(t_interaction());
    
    Maybe_error<bool> out(true);
       
    for (auto i=0ul; i<v_positions().size(); ++i)
    {
        if (v_players().size()!=v_positions()[i].size())
            out=error_message(out.error()()+"  "+std::to_string(i)+"th interaction position size mismatch: "+
                                std::to_string(v_positions()[i].size())+ " vs "+std::to_string(v_players().size()));
        else
        {
            Maybe_error<bool> outp(true);
            for (auto j=0ul; j<v_positions()[i].size(); ++j)
            {
                if (v_players()[j]()!= get<Conformational_change_label>(t_scheme()[j]())())
                    outp=error_message(outp.error()()+" at position "+std::to_string(j)+": "+
                                         get<Conformational_change_label>(t_scheme()[j]())()+"is not "+
                                        v_players()[j]());
            }
            if (!outp)
                out= error_message(out.error()()+"  "+std::to_string(i)+"th interaction label mismatch: "+ out.error()());
        }
    }
    return out;
}


Maybe_error<Conformational_model_scheme>
make_Conformational_model_scheme(Agonist_dependency_map&& t_agonist_map,
                                 std::vector<Conformational_change_label>&& t_scheme,
                                 std::vector<Conformational_interaction>&& t_interactions )
{
    auto Maybe_Conformational_change=make_Conformational_change_scheme(std::move(t_agonist_map), std::move(t_scheme));
    if (!Maybe_Conformational_change)
        return Maybe_Conformational_change.error();
    else
    {
        Maybe_error<bool> out(true);
        for (std::size_t i=0; i<t_interactions.size(); ++i)
        {
            if (auto outp=check_Conformational_interaction(Maybe_Conformational_change.value(),t_interactions[i]); !outp)
                out=error_message(out.error()()+ " "+std::to_string(i)+"th interaction: "+outp.error()()+"\n");
        }
        if (!out)
            return out.error();
        else
            return Conformational_model_scheme(Vector_Space(std::move(Maybe_Conformational_change.value()),
                                                            Conformational_interaction_scheme(std::move(t_interactions))));
    }
}




template <class P, class Id>
    requires std::is_same_v<var::untransformed_type_t<P>, Parameters<Id>>
auto calc_Qij(const Conformational_interaction_scheme& inter,const  typename Parameters<Id>::Names& names,const P& par,const Conformational_transition& tr )->Maybe_error<var::Op_t<transformation_type_t<P>,double>>
{
    auto i=get<Conformational_transition_initiating_state_index>(tr());
    auto j=get<Conformational_transition_landing_state_index>(tr());
    auto ag= get<Agonist_dependency>(tr());
    auto d=get<Conformational_transition_direction>(tr());
    auto chla=get<Conformational_change_label>(tr());
    auto v_int=get<Conformational_interactions_domain_state>(tr());
    auto n=get<Conformational_transition_mulitplicity>(tr());
    
    auto Maybe_i_base=d()?names[chla()+"_on"]:names[chla()+":off"];
    
    if(!Maybe_i_base)
        return Maybe_i_base.error();
    else{
        auto i_base=Maybe_i_base.value();
        auto out=n()*par[i_base];
        for (std::size_t ii=0; ii<v_int().size(); ++ii)
         {
            auto factor_la=get<Conformational_interaction_label>(inter()[get<Conformational_interaction_index>(v_int()[ii])()]())();
            auto factor_ipos=get<Conformational_interaction_subposition>(v_int()[ii])();
            auto Maybe_i_Factor=names[factor_la];
            auto Maybe_i_Factor_pos=names[factor_la+"_"+std::to_string(factor_ipos)];
            if (!Maybe_i_Factor || !Maybe_i_Factor_pos)
                return error_message(Maybe_i_Factor.error()()+Maybe_i_Factor_pos.error()());
            else
            {
                auto i_Factor=Maybe_i_Factor.value();
                auto i_Factor_pos=Maybe_i_Factor_pos.value();
                
                using std::pow;
                if (d())
                    out=out*pow(par[i_Factor],par[i_Factor_pos]);
                else
                    out=out*pow(par[i_Factor],par[i_Factor_pos]-1.0);
            }
         }
        return out;
    }
}

template <class P, class Id>
    requires std::is_same_v<var::untransformed_type_t<P>, Parameters<Id>>
 auto make_Q0_Qa(const Conformational_model& model,
                                          const  typename Parameters<Id>::Names& names,
                const P& par)->Maybe_error<std::tuple<Transfer_Op_to<P, Q0>,Transfer_Op_to<P, Qa>>>
{
    
    using Trans = transformation_type_t<P>;
    
    auto N=get<N_St>(model())();
    auto tr=get<Conformational_transition_list>(model());
    assert(tr().size()==N);
    auto v_Q0=Op_t<Trans,Q0>(Matrix<double>(N,N,0.0));
    auto v_Qa=Op_t<Trans,Qa>(Matrix<double>(N,N,0.0));
    
    for (std::size_t i=0; i<N; ++i)
    {
        for (std::size_t j=0; j<tr()[i].size(); ++j)
        {
            auto trr=tr()[i][j];
            auto v_i=get<Conformational_transition_initiating_state_index>(trr());
            assert(v_i()()==i);
            auto v_j=get<Conformational_transition_landing_state_index>(trr());
            auto ag= get<Agonist_dependency>(trr());
            auto Maybe_qij=calc_Qij(model,names, par,trr);
            if (!Maybe_qij)
                return Maybe_qij.error();
            else if (!ag())
            {
                set(v_Q0(),i,i,v_Q0()(i,i)-Maybe_qij.value());
                set(v_Q0(),i,v_j(),Maybe_qij.value());
            }
            else
            {
                set(v_Qa(),i,i,v_Qa()(i,i)-Maybe_qij.value());
                set(v_Qa(),i,v_j(),Maybe_qij.value());
            }    
            
            
        }
    }
    return std::tuple(v_Q0,v_Qa);
    
}


}

Maybe_error<Conformational_model> make_Conformational_model(Agonist_dependency_map&& t_agonist_map,
                               std::vector<Conformational_change_label>&& t_scheme,
                               std::vector<Conformational_interaction>&& t_interactions )
{
    auto Maybe_scheme=impl::make_Conformational_model_scheme(std::move(t_agonist_map),std::move(t_scheme),std::move(t_interactions));
    if (!Maybe_scheme)
        return Maybe_scheme.error();
    else
    {
        auto& scheme=Maybe_scheme.value();
        auto Maybe_tuple=impl::make_Conformational_states_and_index(scheme);
        if (!Maybe_tuple) return Maybe_tuple.error();
        else
        {
            auto [states,map]=std::move(Maybe_tuple.value());
            auto Maybe_transition=impl::make_Conformational_transition_list(scheme,states,map);
            if (!Maybe_transition)
                return Maybe_transition.error();
            else
            {
                return Conformational_model(Vector_Space(N_St(states().size()),std::move(scheme),std::move(states), std::move(Maybe_transition.value())));
            }
        }
    }
}








template <class P, class Id>
requires std::is_same_v<var::untransformed_type_t<P>, Parameters<Id>>
auto make_Patch_Model(const Conformational_model& model, const P& par)
{
    
}








///**
// * @brief The Allosteric_Model class build a kinetic rate model starting with a
// * set of conformational changes and their interactions.
// *
// *
// */
//class Allosteric_Model {
//public:
//    typedef Allosteric_Model self_type;
//    typedef Model_Parameter_label myParameter_label;
//    struct transitions {
//        bool on;
//        bool agonist;
//        std::string conformation;
//        std::map<std::vector<std::pair<std::string, std::string>>, std::size_t>
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
//                                             std::pair<std::string, std::string>>>
//            conformational_inter_unit_cell;
//    };
    
//    struct new_model_definition {
//        std::size_t number_of_units;
//        std::map<Conformational_change_label, Conformational_change>
//            conformational_changes;
//        std::vector<Conformational_change_label> unit_of_conformational_changes;
//        std::set<Conformational_interaction> conformational_interactions;
//        std::map<std::size_t, Conductance_Parameter_label> conductance_names;
//    };
    
//    model_definition new_to_old(const new_model_definition &m) {
//        model_definition out;
//        for (auto &e : m.conductance_names)
//            out.conductance_names[e.first] = e.second;
//        auto k = m.unit_of_conformational_changes.size();
//        std::size_t n = m.number_of_units * k;
//        out.conformational_changes.resize(n);
//        for (std::size_t i = 0; i < m.number_of_units; ++i)
//            for (std::size_t j = 0; j < m.unit_of_conformational_changes.size();
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
//            out.conformational_changes_names[std::pair(e.second.label(), true)] =
//                e.second.par_on();
//            out.conformational_changes_names[std::pair(e.second.label(), false)] =
//                e.second.par_off();
//        }
//        for (auto &e : m.conformational_interactions) {
//            std::vector<std::size_t> cc;
//            std::size_t current_i = 0;
//            for (std::size_t j = 0; j < e.interacting_conformational_changes().size();
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
//                        s, std::pair(e.factor_label(), e.coefficient_labels()[i])));
//            }
//        }
//        return out;
//    }
    
//    constexpr static auto className = my_static_string("Allosteric_Model");
    
//protected:
//    new_model_definition new_d_;
//    std::vector<std::vector<std::string>> conformer_;
//    std::vector<std::map<std::size_t, transitions>> transitions_;
//    std::set<std::string> paramNames_;
//    std::vector<std::string> conductances_;
//    model_definition d_;
    
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
//            std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
//            &conformational_inter) {
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
//        for (std::size_t i = 0; i < (1u << conformational_changes.size()); ++i) {
//            auto c = index_to_conformational_change(conformational_changes, i);
//            if (state_to_conformer.find(c) == state_to_conformer.end()) {
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
//        const std::map<std::vector<std::string>, std::size_t> &state_to_conformer,
//        const std::set<std::size_t> &agonist_changes,
//        const std::map<std::pair<std::string, bool>, std::string>
//            conformational_changes_names,
//        std::multimap<std::size_t, std::pair<std::set<std::size_t>,
//                                             std::pair<std::string, std::string>>>
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
    
//public:
//    static auto get_constructor_fields_old() {
//        return std::make_tuple(
//            grammar::field(C<self_type>{}, "conformational_changes",
//                           &self_type::get_conformational_changes),
//            grammar::field(C<self_type>{}, "conformational_changes_names",
//                           &self_type::get_conformational_changes_names),
//            grammar::field(C<self_type>{}, "agonist_changes",
//                           &self_type::get_agonist_changes),
//            grammar::field(C<self_type>{}, "conductance_changes",
//                           &self_type::get_conductance_changes),
//            grammar::field(C<self_type>{}, "conductance_names",
//                           &self_type::get_conductance_names),
//            grammar::field(C<self_type>{}, "conformational_inter_unit_cell",
//                           &self_type::get_conformational_inter_unit_cell));
//    }
    
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
//        std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>> &
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
//        std::vector<Conformational_change_label> unit_of_conformational_changes,
//        std::set<Conformational_interaction> conformational_interactions,
//        std::map<std::size_t, Conductance_Parameter_label> conductance_names)
//        : new_d_{number_of_units, conformational_changes,
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
//            std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
//            &conformational_inter_unit_cell)
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
//            for (auto it = transitions_[i].begin(); it != transitions_[i].end();
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
    
//private:
//    void init() {
//        paramNames_ = getParameterNamesFrom(d_.conformational_changes_names,
//                                            d_.conductance_names,
//                                            d_.conformational_inter_unit_cell);
//        std::cout << paramNames_;
//        auto p = periodicity(d_.conformational_changes);
        
//        auto conformational_interactions = fill_conformational_interactions(
//            d_.conformational_inter_unit_cell, d_.conformational_changes.size(), p);
        
//        std::map<std::vector<std::string>, std::size_t> state_to_conformer;
        
//        std::tie(conformer_, state_to_conformer) =
//            getConformers(d_.conformational_changes, p);
        
//        transitions_ = getTransitions(d_.conformational_changes, conformer_,
//                                      state_to_conformer, d_.agonist_changes,
//                                      d_.conformational_changes_names,
//                                      conformational_interactions);
//        conductances_ =
//            get_g_names(conformer_, d_.conductance_changes, d_.conductance_names);
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
//            std::pair<std::set<std::size_t>, std::pair<std::string, std::string>>>
//            &conformational_interactions,
//        std::size_t n, std::size_t p) {
//        std::multimap<std::size_t, std::pair<std::set<std::size_t>,
//                                             std::pair<std::string, std::string>>>
//            out;
//        for (std::size_t i = 0; i < n / p; ++i) {
//            for (auto &e : conformational_interactions) {
//                out.insert({(e.first + i * p) % n,
//                            {rotate(e.second.first, n, i * p), e.second.second}});
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


}


#endif // ALLOSTERIC_MODELS_H
