#ifndef VARIABLES_DERIVATIVE_H
#define VARIABLES_DERIVATIVE_H
#include "variables.h"
#include "derivative_operator.h"

namespace var {


template<class Id, class T, class X>
class Derivative<Var<Id,T>,X>{
    Derivative<T,X> m_x;
public:
    
    static constexpr bool is_variable=true;
    template<class S>
        requires (std::constructible_from<Derivative<T,X>,S>)
    constexpr Derivative(S&& t_x):m_x{std::forward<S>(t_x)}{}
    constexpr auto& operator()(){return m_x;}
    constexpr auto& operator()()const{return m_x;}
    constexpr auto& operator[](Var<Id>)const{return *this;}
    constexpr Derivative(){}
    constexpr auto& value()const {return m_x;;}
   // friend auto& print(std::ostream& os, const Var& x){ os<<typeid(Id).name()<<": \t"<<x.value()<<"\t"; return os;}
    
   // friend auto& operator<<(std::ostream& os, const Var& x){ os<<x.value(); return os;}
   // friend auto& put(std::ostream& os, const Var& x){ os<<x.value()<<"\t"; return os;}
};


template<class...Ids, class X>
    requires (Ids::is_variable&&...)
class Derivative<Vector_Space<Ids...>,X>: public Vector_Space<Derivative_t<Ids,X>...>
{
    using base_type=Vector_Space<Derivative_t<Ids,X>...>;
public:
    template<class Id>
    friend auto& get(Derivative const& x){return static_cast<Derivative_t<Id,X> const&>(x);}
    template<class Id>
    friend auto& get(Derivative & x){return static_cast<Derivative_t<Id,X> &>(x);}
    static constexpr bool is_vector_space=true;
    
    Derivative(){}
    Derivative(Derivative_t<Ids,X>&&...t_vars): base_type{std::move(t_vars)...}{}
    Derivative(Derivative_t<Ids,X> const&...t_vars): base_type{t_vars...}{}
    // Vector_Space(std::decay_t <decltype(std::declval<Vars const&>().value())> ... t_vars): Vars{std::move(t_vars)}...{}
    
};



template<class Id, class X>
    requires(Id::is_variable)
class Derivative<Id,X>: public Derivative<typename Id::variable_type,X>{};






}


#endif // VARIABLES_DERIVATIVE_H
