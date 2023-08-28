#pragma once
#include "matrix.h"
#include <functional>
#include <numeric>
#ifndef QMODEL_H
#define QMODEL_H
#include <map>
#include <string>

#include "maybe_error.h"

namespace macrodr {

class State_Model;

Maybe_error<State_Model>
to_State_Model(std::size_t t_number_of_states,
               std::map<std::pair<std::size_t, std::size_t>, std::string>
                   &t_transition_rates,
               std::map<std::pair<std::size_t, std::size_t>, std::string>
                   &t_agonist_transition_rates,
               std::map<std::size_t, std::string> t_conductances);

class State_Model {
  std::size_t m_number_of_states;

  std::map<std::pair<std::size_t, std::size_t>, std::string> m_transition_rates;
  std::map<std::pair<std::size_t, std::size_t>, std::string>
      m_agonist_transition_rates;
  std::map<std::size_t, std::string> m_conductances;

  State_Model(std::size_t t_number_of_states,
              std::map<std::pair<std::size_t, std::size_t>, std::string>
                  &t_transition_rates,
              std::map<std::pair<std::size_t, std::size_t>, std::string>
                  &t_agonist_transition_rates,
              std::map<std::size_t, std::string> t_conductances)
      : m_number_of_states{t_number_of_states},
        m_transition_rates{t_transition_rates},
        m_agonist_transition_rates{t_agonist_transition_rates},
        m_conductances{t_conductances} {}

public:
  friend Maybe_error<State_Model>
  to_State_Model(std::size_t t_number_of_states,
                 std::map<std::pair<std::size_t, std::size_t>, std::string>
                     &t_transition_rates,
                 std::map<std::pair<std::size_t, std::size_t>, std::string>
                     &t_agonist_transition_rates,
                 std::map<std::size_t, std::string> t_conductances) {
    for (auto &elem : t_transition_rates) {
      if (elem.first.first >= t_number_of_states)
        return error_message(
            "transition start state greater than number of "
            "states; number_of_states = " +
            std::to_string(t_number_of_states) +
            " start state= " + std::to_string(elem.first.first));
      else if (elem.first.second >= t_number_of_states)
        return error_message(
            "transition end state greater than number of "
            "states; number_of_states = " +
            std::to_string(t_number_of_states) +
            " end state= " + std::to_string(elem.first.second));
      else if (elem.first.second == t_number_of_states)
        return error_message(
            "transition start state same as end state;  start state: " +
            std::to_string(elem.first.first) +
            " end state= " + std::to_string(elem.first.second));
    }
    //
    for (auto &elem : t_agonist_transition_rates) {
      if (elem.first.first >= t_number_of_states)
        return error_message(
            "agonist transition start state greater than number of "
            "states; number_of_states = " +
            std::to_string(t_number_of_states) +
            " start state= " + std::to_string(elem.first.first));
      else if (elem.first.second >= t_number_of_states)
        return error_message(
            "agonist transition end state greater than number of "
            "states; number_of_states = " +
            std::to_string(t_number_of_states) +
            " end state= " + std::to_string(elem.first.second));
      else if (elem.first.second == t_number_of_states)
        return error_message(
            "agonist transition start state same as end state;  start state: " +
            std::to_string(elem.first.first) +
            " end state= " + std::to_string(elem.first.second));
    }

    std::size_t i_cond = 0;
    for (auto &elem : t_conductances) {
      if (elem.first >= t_number_of_states)
        return error_message(
            "state conductance number greater than number of states:"
            " proposed= " +
            std::to_string(elem.first) +
            " number_of_states = " + std::to_string(t_number_of_states));

      if (elem.first != i_cond)
        return error_message("state conductance skipped: current is" +
                             std::to_string(i_cond) +
                             " proposed= " + std::to_string(elem.first));
      ++i_cond;
    }
    if (i_cond != t_number_of_states)
      return error_message(
          "state conductance missing: number of proposed states=" +
          std::to_string(i_cond) +
          " number_of_states = " + std::to_string(t_number_of_states));

    return State_Model(t_number_of_states, t_transition_rates,
                       t_agonist_transition_rates, t_conductances);
  }
};

template <typename Id, class T>
class Named: public T {
  
  public:
  using value_type=T;
  
  
  Named(T &&t_x) : value_type(std::move(t_x)) {}
  Named(T const &t_x) : value_type(t_x) {}
  Named() {}
  
  auto &operator()() const& { return static_cast<T const&>(*this) ; }
  auto &&operator()() && { return static_cast<T &&>(*this) ; }
  
  
};


template <typename Id>
class Quantity {
  public:
  constexpr static bool is_Quantity = true;
};

template<class T, int N>
class Power;

template<class... T>
class Product;



namespace quantities {




class Time : public Quantity<Time> {};

class Longitude : public Quantity<Longitude> {};

class Current : public Quantity<Current> {};

template<class Id>
class Count: public Quantity<Count<Id>>{};

template<class Id>
class Probability_Mass: public Quantity<Probability_Mass<Id>>{};

template<class Id, class MyQuantity>
    requires MyQuantity::is_quantity
class Probability_Density: public Quantity<Probability_Density<Id, MyQuantity>>{};

template<class Id>
class Molarity : public Quantity<Molarity<Id>> {};

class Volume : public Quantity<Volume> {};



}
template <typename Id, typename Q>
  requires std::is_base_of_v<Quantity<Q>, Q>
class StandardUnit {
public:
  constexpr static bool is_StandardUnit = true;
  constexpr static bool is_Unit = true;
  using unit_type=Id;  
};

template <typename Id, int dec> class Prefix {
public:
  constexpr static bool is_Prefix = true;
};

class Mili : public Prefix<Mili, -3> {};
class Micro : public Prefix<Micro, -6> {};
class Nano : public Prefix<Nano, -9> {};
class Pico : public Prefix<Pico, -9> {};

class Second : public StandardUnit<Second, quantities::Time> {};

template<class Id>
class Molar : public StandardUnit<Molar<Id>, quantities::Molarity<Id>> {};

template<class Id>
class Probability_Mass : public StandardUnit<Probability_Mass<Id>, quantities::Probability_Mass<Id>> {};

template<class Id, class U>
    requires U::is_unit
class Probability_Density : public StandardUnit<Probability_Density<Id,U>, quantities::Probability_Mass<Id>> {};

class Ampere : public StandardUnit<Ampere, quantities::Current> {};


template<class Id>
class Number_of: public StandardUnit<Number_of<Id>, quantities::Count<Id>> {};


template <typename Id, typename P, typename U>
  requires P::is_Prefix && U::is_StandardUnit
class Unit {
public:
  constexpr static bool is_Unit = true;
    using unit_type=U;  
};

class MiliSecond : public Unit<MiliSecond, Mili, Second> {};

template<class Id>
class MiliMolar : public Unit<MiliMolar<Id>, Mili, Molar<Id>> {};

class PicoAmpere : public Unit<PicoAmpere, Pico, Ampere> {};






template <class T, class U>
    requires U::is_Unit
class Value {
    T m_x;
    
public:
    constexpr static bool is_Value = true;
    Value(T &&t_x) : m_x(std::move(t_x)) {}
    Value(T const &t_x) : m_x(t_x) {}
    Value() {}
    
    auto &value() const { return m_x; }
    auto unit()const {return U{};}
    
    friend auto& value(const Value& x){return x.value();}
};




template <class Id, class V>
    requires V::is_Value
class Variable {
    V m_x;
    
public:
    constexpr static bool is_Variable = true;
    Variable(V &&t_x) : m_x(std::move(t_x)) {}
    Variable(V const &t_x) : m_x(t_x) {}
    Variable() {}
    
    auto &operator()() const { return m_x; }
    auto &value() const { return m_x.value(); }
    a
};
template <class Id, class T>
class Probability; 

template <class Id, class T>
    requires Id::is_Variable&&
             std::is_floating_point_v<T>&&
std::is_integral_v<std::decay_t<decltype(Id{}.value())>>
class Probability<Id,T> {
    Value<T,Probability_Mass<Id>> m_p;
    
public:
    constexpr static bool is_Probability = true;
    Probability(T t_p) : m_p(t_p) {}
    Probability() {}
    
    auto &operator()() const { return m_p; }
    auto &value() const { return m_p.value(); }
};


template <class Id, class T>
    requires Id::is_Variable&&
             std::is_floating_point_v<T>&&
             std::is_floating_point_v<std::decay_t<decltype(Id{}.value())>>
class Probability<Id,T> {
    using U=
    Value<T,Probability_Mass<Id>> m_p;
    
public:
    constexpr static bool is_Probability = true;
    Probability(T t_p) : m_p(t_p) {}
    Probability() {}
    
    auto &operator()() const { return m_p; }
    auto &value() const { return m_p.value(); }
};


template <typename Id, class T, class U>
  requires U::is_Unit
class Position {
  T m_x;

public:
  constexpr static bool is_Position = true;
  Position(T &&t_x) : m_x(std::move(t_x)) {}
  Position(T const &t_x) : m_x(t_x) {}
  Position() {}

  auto &operator()() const { return m_x; }
};




template <typename Id, class T, class U>
    requires U::is_Unit
class Magnitude {
  T m_x;
  
  public:
  constexpr static bool is_Magnitude = true;
  Magnitude(T &&t_x) : m_x(std::move(t_x)) {}
  Magnitude(T const &t_x) : m_x(t_x) {}
  Magnitude() {}
  
  auto &operator()() const { return m_x; }
};


template<typename Id>
struct Chemical_Species{
  constexpr static bool is_Chemical_Species = true;
  
};

struct ATP: public Chemical_Species<ATP>{};
struct Agonist: public Chemical_Species<Agonist>{};

struct Channel{};

template <typename Id, class T>
  requires std::is_unsigned_v<T>
class Size {
  T m_n;

public:
  constexpr static bool is_Size = true;

  Size(T &&t_x) : m_n(std::move(t_x)) {}
  Size(T const &t_x) : m_n(t_x) {}
  Size() {}

  auto &operator()() const { return m_n; }
};

template <typename Id, class S>
  requires S::is_Size
class Index {
  using T = std::decay_t<decltype(S{}())>;
  T m_n;

public:
  constexpr static bool is_Index = true;

  Index(T &&t_x) : m_n(std::move(t_x)) {}
  Index(T const &t_x) : m_n(t_x) {}
  Index() {}

  auto &operator()() const { return m_n; }
};



template<class F,class ResultType, class... ArgumentsTypes>
    requires requires (F f, ArgumentsTypes const&...args){ ResultType(std::invoke(f, args()...));}
class Dependence{
  F m_f;
 public:
  ResultType operator()(ArgumentsTypes const&...args)const 
     {
      return std::invoke(m_f,args...);
     }
};





template <class Numerator, class Denominator>
class Ratio {
  using T = decltype(type(Numerator{}() / Denominator{}()));
  T m_x;

public:
  Ratio(T &&t_x) : m_x(std::move(t_x)) {}
  Ratio(T const &t_x) : m_x(t_x) {}
  Ratio() {}
  auto &operator()() const { return m_x; }

  friend Numerator operator*(const Ratio &x, const Denominator &d) {
    return x() * d();
  }
};


template <class Numerator, class Denominator>
class Derivative_Value {
  using T = decltype(Numerator{}() / Denominator{}());
  T m_x;
  
  public:
  Derivative_Value(T &&t_x) : m_x(std::move(t_x)) {}
  Derivative_Value(T const &t_x) : m_x(t_x) {}
  Derivative_Value() {}
  
  auto &operator()() const { return m_x; }
  
  friend Numerator operator*(const Derivative_Value &x, const Denominator &d) {
    return x() * d();
  }
};


template <class EventType, class ProbabilityType>
    requires ProbabilityType::is_Probability
class Probability_Value {
  ProbabilityType m_logp;
  
  public:
  Probability_Value(ProbabilityType t_logp)
      : m_logp{t_logp} {}
  Probability_Value() {}
  auto &operator()() const { return m_logp; }
};



template <class EventType, class ConditionalType, class ProbabilityType>
class Probability_Conditional_Value {
  ProbabilityType m_p;
 public:
  Probability_Conditional_Value(ProbabilityType t_p)
      : m_p{t_p} {}
  Probability_Conditional_Value() {}
  auto operator()() const { return m_p; }
};

template <class EventType, class ProbabilityType>
class Probability_Event {
  EventType m_x;
  ProbabilityType m_logp;
  
 public:
  Probability_Event(EventType &&t_x, ProbabilityType t_logp)
      : m_x(std::move(t_x)), m_logp{t_logp} {}
  Probability_Event(EventType const &t_x, ProbabilityType t_logp)
      : m_x(t_x), m_logp{t_logp} {}
  Probability_Event() {}
  auto &event() const { return m_x; }
  auto &logP() const { return m_logp; }
};






template <class EventType, class PMF>
    requires std::is_floating_point_v<std::invoke_result_t<PMF,decltype(EventType{}())>>
class Probability_Mass_Function {
  using T=std::decay_t<std::invoke_result_t<PMF,decltype(EventType{}())>>;
  
  using ProbabilityType=Probability_Mass<T>;
  
  
  PMF m_pmf;
  
  public:
  Probability_Mass_Function(PMF &&t_pmf)
          : m_pmf{std::move(t_pmf)} {}
  Probability_Mass_Function(PMF const&t_pmf)
      : m_pmf{std::move(t_pmf)} {}
  Probability_Mass_Function() {}
  
  
  Probability_Value<EventType, ProbabilityType> operator()(EventType const & x)const
  {
    return std::invoke(m_pmf,x);
  }
};



template <class EventType, class Conditional_type,class PMF>
    requires std::is_floating_point_v<std::decay_t<std::invoke_result_t<PMF,decltype(Conditional_type{}()),decltype(EventType{}())>>>
class Conditional_Probability_Mass_Function {
  using T=std::decay_t<std::invoke_result_t<PMF,decltype(Conditional_type{}()),decltype(EventType{}())>>;
  
  using ProbabilityType=Probability_Mass<T>;
  
  
  PMF m_pmf;
  
  public:
  Conditional_Probability_Mass_Function(PMF &&t_pmf)
      : m_pmf{std::move(t_pmf)} {}
  Conditional_Probability_Mass_Function(PMF const&t_pmf)
      : m_pmf{std::move(t_pmf)} {}
  Conditional_Probability_Mass_Function() {}
  
  auto& operator()()const{return m_pmf;}
  
  Probability_Value<EventType, ProbabilityType> operator()(Conditional_type const& X,EventType const & x)const
  {
    return std::invoke(m_pmf,X,x);
  }
};



class Number_of_channel_states
    : public Size<Number_of_channel_states, std::size_t> {};

class Number_of_channels
    : public Variable<Number_of_channels, Value<std::size_t,Number_of<Channel>>> {};


class i_State : public Index<i_State, Number_of_channel_states> {};


class Time: public Variable<Time,Value<double,Second>>{};

class Agonist_concentration : public Variable<Agonist_concentration,Value<double,MiliMolar<Agonist>>> {};



class Markov_Transition_rate
    : public Variable<Markov_Transition_rate, Derivative_Value<Conditional_Probability_Mass_Function<i_State,i_State,Matrix<double>>,Time>> {
};

class Time_of_experiment : public Named<Time_of_experiment,Time> {};


class Markov_Transition_rate_by_Agonist
    : public Named<Markov_Transition_rate_by_Agonist, Derivative_Value<
                                                          Derivative_Value<Conditional_Probability_Mass_Function<i_State,i_State,Matrix<double>>,
                                                                                        Time<Second>>, Agonist_concentration>> {
};

class Markov_Transition_rate_resting
    : public Named<Markov_Transition_rate, Derivative_Value<Conditional_Probability_Mass_Function<i_State,i_State,Matrix<double>>,Time<Second>>> {
};



class Channels_Current:public Named<Channels_Current,Current<PicoAmpere>>{};

auto r=Ratio<Channels_Current,Number_of_channels>{};

class Unitary_Current:public Named<Unitary_Current,Ratio<Channels_Current,Number_of_channels>>{
  using Named<Unitary_Current,Ratio<Channels_Current,Number_of_channels>>::Named;
  
};



class States_Current
    : public Named<States_Current, Dependence<Matrix<double>,Unitary_Current,i_State>> {
};

class Markov_Transition_step
    : public Named<Markov_Transition_step, Matrix<double>> {
  using Named<Markov_Transition_step, Matrix<double>>::Named;
};



class Markov_Model
    : public Named<Markov_Model,
                   std::tuple<Markov_Transition_rate_resting,
                              Markov_Transition_rate_by_Agonist,
                              Instantaneous_conductance, Number_of_channels>> {
};

struct eig_method {
  Maybe_error<Matrix<double>> expm(const Markov_Transition_rate &x, double dt) {
    auto Maybe_exp = eigs(x());
    if (Maybe_exp) {
      auto [V, la, W] = std::move(Maybe_exp.value());
      return V * apply([dt](auto lai) { return std::exp(dt * lai); }, la) * W;
    } else
      return Maybe_exp.error();
  }
};

class Macro_DR {

public:
  auto init(const Markov_Model &model, Agonist_concentration initial) {}

  auto logLikelihood(
      const Markov_Model &model, Agonist_concentration initial,
      std::vector<std::tuple<Time_stamp, Time_interval, Agonist_concentration,
                             Signal_current>> const &data) {
    auto markov_state =
        return std::reduce(std::execution::seq, data.begin(), data.end(), )
  }
};

} // namespace macrodr

#endif // QMODEL_H
