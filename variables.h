#ifndef VARIABLES_H
#define VARIABLES_H

#include "maybe_error.h"
#include <map>
#include <ostream>
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
    constexpr Var(T t_x):m_x{t_x}{}
    constexpr auto& operator()(){return m_x;}
    constexpr auto& operator()()const{return m_x;}
    constexpr auto& operator[](Var<Id>)const{return *this;}
    constexpr Var(){}
    constexpr auto& value()const {return m_x;}
    friend auto& print(std::ostream& os, const Var& x){ os<<typeid(Id).name()<<": \t"<<x.value()<<"\t"; return os;}
    
    friend auto& operator<<(std::ostream& os, const Var& x){ os<<x.value(); return os;}
    friend auto& put(std::ostream& os, const Var& x){ os<<x.value()<<"\t"; return os;}
};


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


template<class...Vars>
requires (Vars::is_variable&&...)
class Vector_Space: public Vars...
{
    static bool islessthan(const Vector_Space& ,const Vector_Space& )
    {
        return false;
    }
    
    template<class V0,class... Vs>
    static bool islessthan(const Vector_Space& a,const Vector_Space& b)
    {
        if (static_cast<V0 const&>(a).value()<static_cast<V0 const&>(b).value())
        {
            return true;
        } else if (static_cast<V0 const&>(b).value() < static_cast<V0 const&>(a).value())
        {
                return false;}
        else return islessthan<Vs...>(a,b);
    }
    
public:
    using Vars::operator[]...;
    template<class Id>
    friend auto& get(Vector_Space const& x){return static_cast<Id const&>(x);}
    template<class Id>
    friend auto& get(Vector_Space & x){return static_cast<Id&>(x);}
    static constexpr bool is_vector_space=true;
    
    Vector_Space(){}
    Vector_Space(Vars&&...t_vars): Vars{std::move(t_vars)}...{}
    Vector_Space(Vars const&...t_vars): Vars{t_vars}...{}
   // Vector_Space(std::decay_t <decltype(std::declval<Vars const&>().value())> ... t_vars): Vars{std::move(t_vars)}...{}
    friend std::ostream& operator<<(std::ostream& os, const Vector_Space& tu)
    {
        return ((os<<static_cast<Vars const&>(tu).value()<<"\t"),...);
    }
    
    friend std::ostream& print(std::ostream& os, const Vector_Space& tu)
    {
        return ((os<<typeid(Vars).name()<<" :=\t"<<static_cast<Vars const&>(tu).value()<<"\t"),...);
    }
    
    
    friend bool operator<(const Vector_Space& a,const Vector_Space& b)
    {
        return islessthan<Vars...>(a,b);
    }
    
};
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


} // namespace var

#endif // VARIABLES_H
