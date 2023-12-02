#ifndef VARIABLES_H
#define VARIABLES_H

#include "maybe_error.h"
#include "general_output_operator.h"
#include <map>
#include <ostream>
#include <functional>
namespace var {

 


template<class Id>
struct Unit{ public: static constexpr bool is_unit=true;};


struct s: public Unit<s>{};
struct kHz: public Unit<kHz>{};

struct ms: public Unit<ms>{};
struct number: public Unit<number>{};
struct prob: public Unit<prob>{};
struct dimensionless: public Unit<dimensionless>{};

struct microMolar: public Unit<microMolar>{};
struct pA: public Unit<pA>{};

template<class Id, int n>
struct Power;

template<class... Ids>
struct Product;


template<class Id, int n>
    requires Id::is_unit
struct Power<Id,n>: public Unit<Power<Id,n>>{};



template<class... Ids>
    requires (Ids::is_unit&&...)
struct Product<Ids...>: public Unit<Product<Ids...>>{};




template<class T, class Unit>
    requires Unit::is_unit
class Value{
    T m_x;
public:
    Value(T t_x):m_x{t_x}{}
    
    auto& value()const {return m_x;}
};


template<class...>
class Var;

template<class Id>
class Var<Id>{};




template<class Id, class T>
class Var<Id,T>{
    T m_x;
public:
    static constexpr bool is_variable=true;
    static constexpr bool is_constant=false;
    
    using variable_type=Var<Id,T>;
    
    template<class S>
        requires std::is_convertible_v<S,T>
    constexpr Var(S&& t_x):m_x{std::forward<S>(t_x)}{}
    
    constexpr Var(T const& t_x):m_x{t_x}{}
    
    constexpr Var(T&& t_x):m_x{std::move(t_x)}{}
    constexpr auto& operator()(){return m_x;}
    constexpr auto const& operator()()const{return m_x;}
    constexpr auto& operator[](Var<Id>){return *this;}
    constexpr auto& operator[](Var<Id>) const{return *this;}
    constexpr Var(){}
    constexpr auto& value()const {return m_x;}
    
    
    
    template<class... Ts>
    constexpr auto operator()(const Ts&...){return Id(*this);}
    
    template<class... Ts>
    constexpr auto operator()(const Ts&...)const {return Id(*this);}
    
    
    friend auto& print(std::ostream& os, const Var& x){
        os<<typeid(Id).name()<<": \n";
        print(os,x.value());
        os<<"\n";
        return os;
    }
    
    friend Id  operator-(const Var& one, const Var& two){return Id(one()-two());}
    
    friend bool operator<(const Var& one,const Var& two){ return one.value()<two.value();}
    friend auto& operator<<(std::ostream& os, const Var& x){ os<<x.value(); return os;}
    friend auto& put(std::ostream& os, const Var& x){ os<<x.value()<<"\t"; return os;}
};

template<class...>
class Constant;


template<class Id>
class Constant<Id>{};


template<class Id, class T>
class Constant<Id,T>{
    T m_x;
public:
    static constexpr bool is_variable=true;
    static constexpr bool is_constant=true;
    
    //operator Id()const {return Id(*this);}
    
    constexpr Constant(T&& t_x):m_x{std::move(t_x)}{}
    constexpr Constant(T const& t_x):m_x{t_x}{}
    constexpr auto& operator()(){return m_x;}
    constexpr auto& operator()()const{return m_x;}
    constexpr auto& operator[](Constant<Id>){return *this;}
    constexpr auto& operator[](Constant<Id>) const{return *this;}
    
    template<class... Ts>
    constexpr auto operator()(const Ts&...){return Id(*this);}
    
    template<class... Ts>
    constexpr auto operator()(const Ts&...)const {return Id(*this);}
    
    constexpr Constant(){}
    constexpr auto& value()const {return m_x;}
    friend auto& print(std::ostream& os, const Constant& x){
        os<<typeid(Id).name()<<": \n";
        print(os,x.value());
        os<<"\n";
        return os;
        
    }
    friend bool operator<(const Constant& one,const Constant& two){ return one.value()<two.value();}
    
    friend auto& operator<<(std::ostream& os, const Constant& x){ os<<x.value(); return os;}
    friend auto& put(std::ostream& os, const Constant& x){ os<<x.value()<<"\t"; return os;}
};




template<class Id, class F, class ...T>
class Fun{
    std::tuple<T...> m_x; 
    F m_f;
public:
    static constexpr bool is_variable=true;
    static constexpr bool is_constant=false;
    
    using variable_type=Fun<Id,F,T...>;
    
    
    constexpr Fun(Var<Id>,F const& t_f,T const&... t_x):m_x{t_x...},m_f{t_f}{}
      
    constexpr Fun(Var<Id>,F && t_f,T && ...t_x):m_x{std::move(t_x)...},m_f{std::move(t_f)}{}
    
    template<class... Ts>
    constexpr auto operator()(const Ts&...ts){return Id(std::apply([this,&ts...](auto&...xs) {return std::invoke(m_f,xs...,ts...);},m_x));}
    
    template<class... Ts>
    constexpr auto operator()(const Ts&...ts)const {return Id(std::apply([this,&ts...](auto&...xs) {return std::invoke(m_f,xs...,ts...);},m_x));}
    
    constexpr Fun(){}
    
    constexpr auto& operator[](Var<Id>){return *this;}
    constexpr auto& operator[](Var<Id>) const{return *this;}
    
   };
   
   
   
   
   template<class Id, class T, class F>
   Fun(Var<Id>,T,F)->Fun<Id,std::decay_t<T>,F>;


template<class var, class... T>
    requires (std::constructible_from<var,T...>)
var build(T...x){return var(std::forward<T>(x)...);}

// template<template<class...> class var, class... T>
//     requires ((std::constructible_from<var<T...>,T...>)&&(!std::is_same_v<var<T...>,Fun<T...>>))
// var<T...> build(T...x){return var(std::forward<T>(x)...);}






template<class...>
class struct_Var;

template<class Id>
class struct_Var<Id>{};

template<class Id, class T>
class struct_Var<Id,T>{
 public:
    T value;
    static constexpr bool is_variable=true;
    constexpr struct_Var(T t_x):value{t_x}{}
    constexpr auto& operator[](struct_Var<Id>)const{return *this;}
    constexpr struct_Var(){}
    
    friend auto& operator<<(std::ostream& os, struct_Var x){ os<<x.value; return os;}
    friend auto& put(std::ostream& os, const struct_Var& x){ os<<x.value<<"\t"; return os;}
};

/*

template<class Id, class T, class Unit>
    requires Unit::is_unit
class Var<Id,T,Unit>{
    Value<T,Unit> m_x;
public:
    static constexpr bool is_variable=true;
    Var(T t_x):m_x{t_x}{}
    auto& operator()()const{return m_x;}
    auto& operator()(){return m_x;}
    template<class Id2>
    auto& operator[](Var<Id>)const{return *this;}
    auto& value()const {return m_x.value();}
    //    friend auto& operator<<(std::ostream& os, Var x){ os<<x.scalar(); return os;}
    friend auto& operator<<(std::ostream& os, const Var& x){ os<<x.value(); return os;}
    friend auto& put(std::ostream& os, const Var& x){ os<<x.value()<<"\t"; return os;}    
};

*/

template<class...Vars>
requires (Vars::is_variable&&...)
class Vector_Space: public Vars...
{
    static bool islessthan(const Vector_Space& ,const Vector_Space& )
    {
        return false;
    }
    template<class V0,class V2,class... Vs>
    static bool islessthan(const Vector_Space& a,const Vector_Space& b)
    {
        if (static_cast<V0 const&>(a).value()<static_cast<V0 const&>(b).value())
        {
            return true;
        } else if (static_cast<V0 const&>(b).value() < static_cast<V0 const&>(a).value())
        {
            return false;
        }
        else return islessthan<V2,Vs...>(a,b);
    }
    template<class V0>
    static bool islessthan(const Vector_Space& a,const Vector_Space& b)
    {
        if (static_cast<V0 const&>(a).value()<static_cast<V0 const&>(b).value())
        {
            return true;
        } else if (static_cast<V0 const&>(b).value() < static_cast<V0 const&>(a).value())
        {
                return false;
        }
        else return islessthan(a,b);
    }
    
public:
    using Vars::operator[]...;
    template<class Id>
    friend auto const& get(Vector_Space const& x){return static_cast<Id const&>(x);}
    template<class Id>
    friend auto& get(Vector_Space & x){return static_cast<Id &>(x);}
    
    
    
    static constexpr bool is_vector_space=true;
    
    Vector_Space(){}
    Vector_Space(Vars&&...t_vars): Vars{std::move(t_vars)}...{}
    Vector_Space(Vars const&...t_vars): Vars{t_vars}...{}
    
    template<class... Vars2>
    friend auto concatenate(Vector_Space&& one, Vector_Space<Vars2...>&& two)
    {
        return Vector_Space<Vars...,Vars2...>(std::move(get<Vars>(one))...,std::move(get<Vars2>(two))...);
    }        
     
   // Vector_Space(std::decay_t <decltype(std::declval<Vars const&>().value())> ... t_vars): Vars{std::move(t_vars)}...{}
    friend std::ostream& operator<<(std::ostream& os, const Vector_Space& tu)
    {
        return ((os<<static_cast<Vars const&>(tu).value()<<"\t"),...);
    }
    
    friend std::ostream& print(std::ostream& os, const Vector_Space& tu)
    {
        (print(os,static_cast<Vars const&>(tu)),...);
        os<<"\n";
        return os;
       //  ((os<<typeid(Vars).name()<<" :=\t"<<static_cast<Vars const&>(tu).value()<<"\t"),...);
       // return os
    }
    
    
    friend bool operator<(const Vector_Space& a,const Vector_Space& b)
    {
        return islessthan<Vars...>(a,b);
    }
    
    friend Vector_Space operator- (const Vector_Space& one, const Vector_Space& two)
    {
        return Vector_Space(Vars(get<Vars>(one)()-get<Vars>(two)())...);
    }
    
};

template<class Id,class...Vars>
auto const& fun(Vector_Space<Vars...> const& x){return x[Var<Id>{}];}

template<class Id,class...Vars>
auto& fun(Vector_Space<Vars...> & x){return x[Var<Id>{}];}


template<class...>
class Vector_Map;

template<class Id>
class Vector_Map<Id>{};


template<class Id, class VS>
    requires VS::is_vector_space
class Vector_Map<Id,VS>{
    std::map<VS,Id> m_x;
public:
    static constexpr bool is_vector_map=true;
    Vector_Map()=default;
    auto& operator[](Vector_Map<Id>)const{return *this;}
    
    Maybe_error<Id const &> operator[](const VS& v)const
    {
        auto it=m_x.find(v);
        if (it!=m_x.end())
                return *it;
        else
                return error_message("");   
    }
    auto& emplace(VS&& v, Id&& x)
    {
        m_x.emplace(std::move(v), std::move(x));
    }
    
};


template<class...VarMaps>
    requires (VarMaps::is_vector_map&&...)
class Vector_Map_Space: public VarMaps...
{
    
public:
    static constexpr bool is_vector_map_space=true;
    
    Vector_Map_Space(){}
};

inline Maybe_error<bool> test_equality(double x, double y, double eps)
{
    if (std::abs(x-y)/(std::abs(x)+std::abs(y)+1.0)>eps)
        return error_message(ToString(x)+" is not equal to "+ToString(y)+
                             "|x-y|/(|x|+|y|) = "+
                             ToString(std::abs(x-y)/(std::abs(x)+std::abs(y)+1.0))+
                             " is greater than "+ToString (eps));
    else return true;
}



template<class Id, class T>
Maybe_error<bool> test_equality(const Var<Id,T>& one, const Var<Id,T>& two, double eps)
{
    auto test=test_equality(one(), two(), eps);
    if (!test)
    {
        return error_message(std::string("\n\n")+typeid(Id).name()+std::string(" error: ")+test.error()() );
    }
    else return true;
}

template<class Id, class T>
Maybe_error<bool> test_equality(const Constant<Id,T>& one, const Constant<Id,T>& two, double eps)
{
    auto test=test_equality(one(), two(), eps);
    if (!test)
    {
        return error_message(std::string("\n\n")+typeid(Id).name()+std::string(" error: ")+test.error()() );
    }
    else return true;
}


template<class...Vars>
    requires (Vars::is_variable&&...)
Maybe_error<bool> test_equality(const Vector_Space<Vars...>& one, const  Vector_Space<Vars...>& two, double eps)
{
    return ((test_equality(get<Vars>(one), get<Vars>(two),eps)&&...));
}






} // namespace var

#endif // VARIABLES_H
