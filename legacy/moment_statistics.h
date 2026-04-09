#ifndef MOMENT_STATISTICS_H
#define MOMENT_STATISTICS_H

#include <maybe_error.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <functional>
#include <iterator>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>
#include <set>
#include "matrix.h"
#include "parameter_indexed.h"
#include "variables.h"

template<class T>
struct value_type_impl {
    using type = std::remove_cvref_t<T>;
};

template<class T>
requires requires (T const& v) { v(); }
struct value_type_impl<T> {
    using type = std::remove_cvref_t<decltype(std::declval<T const&>()())>;
};

template<class T>
using value_type_t = typename value_type_impl<T>::type;






template<typename T>
struct mean_value_type_impl {
    using type = std::conditional_t<std::is_arithmetic_v<std::remove_cvref_t<T>>, double,std::remove_cvref_t<T>>;
};
template<class T>
using mean_value_type_t = typename mean_value_type_impl<T>::type;



template<class T>
requires requires (T const& v) { v(); }
struct mean_value_type_impl<T> {
    using type = mean_value_type_t<value_type_t<T>>;
};













template <class Va>
struct count : public var::Constant<count<Va>, std::size_t> {
    count() : var::Constant<count, std::size_t>{0} {
    }
    count(std::size_t n) : var::Constant<count<Va>, std::size_t>{n} {
    }
};

template <class Va>
struct mean : public var::Var<mean<Va>, mean_value_type_t<Va>>{
   using mean_value_type= mean_value_type_t<Va>;
    using value_type= value_type_t<Va>;

    mean() : var::Var<mean<Va>, mean_value_type>{value_type_t<Va>{}} {
    }
    mean(Va&& x) : var::Var<mean<Va>, mean_value_type>{std::move(x)()} {
    }
    mean(Va const& x) : var::Var<mean<Va>, mean_value_type>{x()} {
    }  
    
    template<typename T>
    requires std::convertible_to<T, mean_value_type>
    mean(T&& x) : var::Var<mean<Va>, mean_value_type>{std::forward<T>(x)} {
    }
 };

template <class Va>
requires std::convertible_to<value_type_t<Va>, double>
struct median : public var::Var<median<Va>, double>{
   using mean_type= double;

    median() : var::Var<median<Va>, mean_type>{value_type_t<Va>{}} {
    }
    median(Va&& x) : var::Var<median<Va>, mean_type>{std::move(x)()} {
    }
    median(Va const& x) : var::Var<median<Va>, mean_type>{x()} {
    }  
    median(value_type_t<Va>&& x) : var::Var<median<Va>, mean_type>{std::move(x)} {
    }
    median(value_type_t<Va> const& x) : var::Var<median<Va>, mean_type>{x} {
    }   
};

template <class Va>
requires std::convertible_to<value_type_t<Va>, double>
struct stddev : public var::Var<stddev<Va>, double>{
   using mean_type= double;

    stddev() : var::Var<stddev<Va>, mean_type>{value_type_t<Va>{}} {
    }
    stddev(Va&& x) : var::Var<stddev<Va>, mean_type>{std::move(x)()} {
    }
    stddev(Va const& x) : var::Var<stddev<Va>, mean_type>{x()} {
    }  
    stddev(value_type_t<Va>&& x) : var::Var<stddev<Va>, mean_type>{std::move(x)} {
    }
    stddev(value_type_t<Va> const& x) : var::Var<stddev<Va>, mean_type>{x} {
    }   
};

template <class Va>
struct Probits: public var::Var<Probits<Va>, std::map<double,mean_value_type_t<Va>>>{
   using mean_value_type= mean_value_type_t<Va>;

    Probits() : var::Var<Probits<Va>, std::map<double,mean_value_type>>{std::map<double,mean_value_type>{}} {
    }
    Probits(std::map<double,mean_value_type>&& x) : var::Var<Probits<Va>, std::map<double,mean_value_type>>{std::move(x)} {
    }
    Probits(std::map<double,mean_value_type>const & x) : var::Var<Probits<Va>, std::map<double,mean_value_type>>{x} {
    }
};


template <class Va>
struct Probit_025 : public var::Var<Probit_025<Va>, value_type_t<Va>>{
   using value_type= value_type_t<Va>;

    Probit_025() : var::Var<Probit_025<Va>, value_type>{value_type_t<Va>{}} {
    }
    Probit_025(Va&& x) : var::Var<Probit_025<Va>, value_type>{std::move(x)()} {
    }
    Probit_025(Va const& x) : var::Var<Probit_025<Va>, value_type>{x()} {
    }  
    Probit_025(value_type&& x) : var::Var<Probit_025<Va>, value_type>{std::move(x)} {
    }
    Probit_025(value_type const& x) : var::Var<Probit_025<Va>, value_type>{x} {
    }   
};

template <class Va>
struct Probit_975 : public var::Var<Probit_975<Va>, value_type_t<Va>>{
   using value_type= value_type_t<Va>;

    Probit_975() : var::Var<Probit_975<Va>, value_type>{value_type_t<Va>{}} {
    }
    Probit_975(Va&& x) : var::Var<Probit_975<Va>, value_type>{std::move(x)()} {
    }
    Probit_975(Va const& x) : var::Var<Probit_975<Va>, value_type>{x()} {
    }  
    Probit_975(value_type&& x) : var::Var<Probit_975<Va>, value_type>{std::move(x)} {
    }
    Probit_975(value_type const& x) : var::Var<Probit_975<Va>, value_type>{x} {
    }   
};


template <class Va>
requires std::convertible_to<value_type_t<Va>, double>
struct Probit_minus_2_sigma : public var::Var<Probit_minus_2_sigma<Va>, double>{
   using mean_type= double;

    Probit_minus_2_sigma() : var::Var<Probit_minus_2_sigma<Va>, mean_type>{value_type_t<Va>{}} {
    }
    Probit_minus_2_sigma(Va&& x) : var::Var<Probit_minus_2_sigma<Va>, mean_type>{std::move(x)()} {
    }
    Probit_minus_2_sigma(Va const& x) : var::Var<Probit_minus_2_sigma<Va>, mean_type>{x()} {
    }  
    Probit_minus_2_sigma(mean_type&& x) : var::Var<Probit_minus_2_sigma<Va>, mean_type>{std::move(x)} {
    }
    Probit_minus_2_sigma(mean_type const& x) : var::Var<Probit_minus_2_sigma<Va>, mean_type>{x} {
    }   
};

template <class Va>
requires std::convertible_to<value_type_t<Va>, double>
struct Probit_minus_1_sigma : public var::Var<Probit_minus_1_sigma<Va>, double>{
   using mean_type= double;

    Probit_minus_1_sigma() : var::Var<Probit_minus_1_sigma<Va>, mean_type>{value_type_t<Va>{}} {
    }
    Probit_minus_1_sigma(Va&& x) : var::Var<Probit_minus_1_sigma<Va>, mean_type>{std::move(x)()} {
    }
    Probit_minus_1_sigma(Va const& x) : var::Var<Probit_minus_1_sigma<Va>, mean_type>{x()} {
    }  
    Probit_minus_1_sigma(mean_type&& x) : var::Var<Probit_minus_1_sigma<Va>, mean_type>{std::move(x)} {
    }
    Probit_minus_1_sigma(mean_type const& x) : var::Var<Probit_minus_1_sigma<Va>, mean_type>{x} {
    }   
};

template <class Va>
requires std::convertible_to<value_type_t<Va>, double>
struct Probit_plus_1_sigma : public var::Var<Probit_plus_1_sigma<Va>, double>{
   using mean_type= double;

    Probit_plus_1_sigma() : var::Var<Probit_plus_1_sigma<Va>, mean_type>{value_type_t<Va>{}} {
    }
    Probit_plus_1_sigma(Va&& x) : var::Var<Probit_plus_1_sigma<Va>, mean_type>{std::move(x)()} {
    }
    Probit_plus_1_sigma(Va const& x) : var::Var<Probit_plus_1_sigma<Va>, mean_type>{x()} {
    }  
    Probit_plus_1_sigma(mean_type&& x) : var::Var<Probit_plus_1_sigma<Va>, mean_type>{std::move(x)} {
    }
    Probit_plus_1_sigma(mean_type const& x) : var::Var<Probit_plus_1_sigma<Va>, mean_type>{x} {
    }   
};

template <class Va>
requires std::convertible_to<value_type_t<Va>, double>
struct Probit_plus_2_sigma : public var::Var<Probit_plus_2_sigma<Va>, double>{
   using mean_type= double;

    Probit_plus_2_sigma() : var::Var<Probit_plus_2_sigma<Va>, mean_type>{value_type_t<Va>{}} {
    }
    Probit_plus_2_sigma(Va&& x) : var::Var<Probit_plus_2_sigma<Va>, mean_type>{std::move(x)()} {
    }
    Probit_plus_2_sigma(Va const& x) : var::Var<Probit_plus_2_sigma<Va>, mean_type>{x()} {
    }  
    Probit_plus_2_sigma(mean_type&& x) : var::Var<Probit_plus_2_sigma<Va>, mean_type>{std::move(x)} {
    }
    Probit_plus_2_sigma(mean_type const& x) : var::Var<Probit_plus_2_sigma<Va>, mean_type>{x} {
    }   
};




template <class Va>
struct Sum : public var::Var<Sum<Va>, value_type_t<Va>>{
   using sum_type= value_type_t<Va>;

    Sum() : var::Var<Sum<Va>, sum_type>{value_type_t<Va>{}} {
    }
    Sum(Va&& x) : var::Var<Sum<Va>, sum_type>{std::move(x)()} {
    }
    Sum(Va const& x) : var::Var<Sum<Va>, sum_type>{x()} {
    }  
    template<class T>
    requires std::convertible_to<T, sum_type>
     Sum(T&& x) : var::Var<Sum<Va>, sum_type>{std::forward<T>(x)} {
    }
       
 template<class VectorSpace, class F>
    requires requires (VectorSpace v, F f) { {f(v)} -> std::convertible_to<value_type_t<Va>>; }
    Sum(std::vector<VectorSpace>const& x,  F&& f)
        : var::Var<Sum<Va>, sum_type>{value_type_t<Va>{}} {
        for (const auto& v : x) {
         (*this)() = (*this)() + sum_type(std::invoke(f, v));
    }
       }
   

};



template <class Va>
struct variance : public var::Var<variance<Va>, std::decay_t<decltype(sqr_X<false>(std::declval<value_type_t<Va>>()))>> {
    using variance_type= std::decay_t<decltype(sqr_X<false>(std::declval<value_type_t<Va>>()))>;
    variance() : var::Var<variance<Va>, variance_type>{variance_type{}} {
    }
    variance(variance_type&& x) : var::Var<variance<Va>, variance_type>{std::move(x)} {
    }
    variance(variance_type const& x) : var::Var<variance<Va>, variance_type>{x} {
    }
    variance(Va const& x) : var::Var<variance<Va>, variance_type>{sqr_X<false>(x())} {
    }
};




template <class Va>
struct covariance : public var::Var<covariance<Va>, decltype(sqr_X<true>(std::declval<value_type_t<Va>>()))> {
    using covariance_type= decltype(sqr_X<true>(std::declval<value_type_t<Va>>()));
    covariance() : var::Var<covariance<Va>, covariance_type>{covariance_type{}} {
    }
    covariance(covariance_type&& x) : var::Var<covariance<Va>, covariance_type>{std::move(x)} {
    }
    covariance(covariance_type const& x) : var::Var<covariance<Va>, covariance_type>{x} {
    }

};

template <bool include_covariance>
struct variance_kind; 

template <>
struct variance_kind<true> {
    template <class Va>
    using type = covariance<Va>;
};  

template <>
struct variance_kind<false  > {
    template <class Va>
    using type = variance<Va>;
};  

template<class Va,bool include_covariance>
using variance_t= typename variance_kind<include_covariance>::template type<Va>;


template <class Id, bool include_covariance=false>
class Moment_statistics
    : public var::Var<Moment_statistics<Id, include_covariance>, var::Vector_Space<count<Id>, mean<Id>,
    variance_t<Id,include_covariance>>>
     {
   public:
    using variance_type = variance_t<Id, include_covariance>;
    using base_type =var::Var<Moment_statistics<Id, include_covariance>, var::Vector_Space<count<Id>, mean<Id>,
    variance_type>> ;
 
    using mean_type=typename mean<Id>::mean_value_type;
    
    Moment_statistics()
        : base_type{
              var::Vector_Space{count<Id>{0}, mean<Id>(Id{}()),
                                variance_type{}}} {
    }
    Moment_statistics(Id x)
        : base_type{var::Vector_Space{count<Id>{1}, mean<Id>(x()), variance_type(sqr_X<include_covariance>(x()-x()))}} {
    }
    Moment_statistics(value_type_t<Id> x)
        : base_type{var::Vector_Space{count<Id>{1}, mean<Id>(x), variance_type(sqr_X<include_covariance>(x-x))}} {
    }
    


    Moment_statistics(std::vector<Id>const& x)
        : base_type{
              var::Vector_Space{count<Id>{x.size()}, mean<Id>(Id{}()),
                                variance_type{}}} {
        if (x.empty()) {
            return;
        }
        auto sum=Id{}();
        auto sum2=variance_type{}();
        for (const auto& xi : x) {
            sum = sum + xi();
            sum2 = sum2 + sqr_X<include_covariance>(xi());
        }
        get<mean<Id>>((*this)())()= sum/x.size();
        auto var= x.size()>1 ? (sum2 - x.size() * sqr_X<include_covariance>(get<mean<Id>>((*this)())()))/(x.size()-1)
                             : variance_type{}();
        get<variance_type>((*this)())()= var;
    }
    
    template<class VectorSpace, class F>
    requires requires (VectorSpace v, F&& f) { {f(v)} -> std::convertible_to<value_type_t<Id>>; }
    Moment_statistics(std::vector<VectorSpace>const& x, const std::vector<std::size_t>& indices, F&& f)
        : base_type{
              var::Vector_Space{count<Id>{indices.size()}, mean<Id>(Id{}()),
                                variance_type{}}} {
        if (indices.empty()) {
            return;
        }
        auto sum=Id{}();
        auto sum2=variance_type{}();
        for (const auto& i : indices) {
            auto xi = value_type_t<Id>(std::forward<F>(f)(x[i]));
            sum = sum + xi;
            sum2 = sum2 + sqr_X<include_covariance>(xi);
        }
        get<mean<Id>>((*this)())()= sum/indices.size();
        auto var= indices.size()>1 ? (sum2 - indices.size() * sqr_X<include_covariance>(get<mean<Id>>((*this)())()))/(indices.size()-1)
                                   : variance_type{}();
        get<variance_type>((*this)())()= var;
    }
    

    Moment_statistics(var::Vector_Space<count<Id>, mean<Id>, variance_type>&& vs) : base_type{std::move(vs)} {
    }

    Moment_statistics(var::Vector_Space<count<Id>, mean<Id>, variance_type>const& vs) : base_type{vs} {
    }


    Moment_statistics(count<Id> n,mean<Id> x, variance_type v) : base_type{var::Vector_Space{n, x, v}} {
    }

    friend Moment_statistics operator&(const Moment_statistics& one,
                                       const Moment_statistics& other) {
        double n0 = get<count<Id>>((one)())();
        double n1 = get<count<Id>>((other)())();

        auto m0 = get<mean<Id>>((one)())();
        auto m1 = get<mean<Id>>((other)())();

        auto v0 = get<variance_type>((one)())();
        auto v1 = get<variance_type>((other)())();
        
        
        double n = n0 + n1;

        auto m = m0 + n1 / n * (m1 - m0);

        auto malt = (m0 * n0 + m1 * n1) / n;

        //  auto v = n>1?v0 + (n1-1)/(n-1)* (v1-v0) + n0/(n-1) * sqr_X(m0-m) + n1/(n-1) * sqr_X(m1-m):0;
        auto v =
            n > 1
                ? ((n0 - 1) * v0 + n0 * sqr_X<include_covariance>(m0) + (n1 - 1) * v1 + 
                n1 * sqr_X<include_covariance>(m1) - n * sqr_X<include_covariance>(m)) /
                      (n - 1)
                : sqr_X<include_covariance>(malt-malt);
        

        return Moment_statistics{count<Id>(n), mean<Id>(Id(m)), variance_type(v)};
    }


 

    Moment_statistics& operator&=(const Moment_statistics& other) {
        *this = *this & other;
        return *this;
    }

    Moment_statistics& operator&=(const Id& x) {
        auto n0 = get<count<Id>>((*this)())();

        auto m0 = get<mean<Id>>((*this)())();

        auto v0 = get<variance_type>((*this)())();
        auto n = n0 + 1;

        auto m = m0 + (x() - m0) / n;
        auto v = n > 1 ? v0 + (sqr_X<include_covariance>(x() - m) - v0) / n0 + sqr_X<include_covariance>(m0 - m) : Id{}();

        get<count<Id>>((*this)())() = n;
        get<variance_type>((*this)())()= std::move(v);
        get<mean<Id>>((*this)())() = std::move(m);
        return *this;
    }
    void reset() {
        get<count<Id>>((*this)())() = 0;
        get<mean<Id>>((*this)())() = Id{}();
        get<variance_type>((*this)())() = sqr_X<include_covariance>(Id{}());
    }

    friend Moment_statistics operator+(const Moment_statistics& one,
                                       const Moment_statistics& other) {
        auto n0 = get<count<Id> >((one)())();
        auto m0 = get<mean<Id>>((one)())();
        auto v0 = get<variance_type>((one)())();
        auto n1 = get<count<Id>>((other)())();
        auto m1 = get<mean<Id>>((other)())();
        auto v1 = get<variance_type>((other)())();
        return Moment_statistics(count<Id>(std::min(n0, n1)), mean<Id>(Id(m0 + m1)), variance_type(Id(v0 + v1)));
                                
    }

    friend Moment_statistics operator*(double a, Moment_statistics const& one) {
        auto n0 = get<count<Id>>((one)())();
        auto m0 = get<mean<Id>>((one)())();
        auto v0 = get<variance_type>((one)())();
        return Moment_statistics(count<Id>(n0),mean<Id>(Id(m0 * a)), variance_type(Id(v0 * a * a)));
    }

    constexpr auto& operator[](var::Var<Id>) {
        return *this;
    }
    constexpr auto& operator[](var::Var<Id>) const {
        return *this;
    }
};





template <class... Ids, bool include_covariance>
struct Moment_statistics<var::Vector_Space<Ids...>, include_covariance>
    : var::Vector_Space<Moment_statistics<Ids, include_covariance>...> {
    using base_type = var::Vector_Space<Moment_statistics<Ids, include_covariance>...>;

    Moment_statistics() : base_type{Moment_statistics<Ids, include_covariance>()...} {}
    explicit Moment_statistics(const var::Vector_Space<Ids...>& sample)
        : base_type{Moment_statistics<Ids, include_covariance>(get<Ids>(sample))...} {}

    template <class... Ms>
        requires(sizeof...(Ms) == sizeof...(Ids))
    explicit Moment_statistics(Ms&&... ms) : base_type{std::forward<Ms>(ms)...} {}

    void reset() { (get<Moment_statistics<Ids, include_covariance>>(*this).reset(), ...); }

    friend auto operator&(const Moment_statistics& a, const Moment_statistics& b) {
        return Moment_statistics((get<Moment_statistics<Ids, include_covariance>>(a) & get<Moment_statistics<Ids, include_covariance>>(b))...);
    }

    Moment_statistics& operator&=(const Moment_statistics& other) {
        ((get<Moment_statistics<Ids, include_covariance>>(*this) &= get<Moment_statistics<Ids, include_covariance>>(other)), ...);
        return *this;
    }

    Moment_statistics& operator&=(const var::Vector_Space<Ids...>& sample) {
        ((get<Moment_statistics<Ids, include_covariance>>(*this) &= get<Ids>(sample)), ...);
        return *this;
    }

    friend auto operator+(const Moment_statistics& a, const Moment_statistics& b) {
        return Moment_statistics((get<Moment_statistics<Ids, include_covariance>>(a) + get<Moment_statistics<Ids, include_covariance>>(b))...);
    }

    friend auto operator*(double k, const Moment_statistics& a) {
        return Moment_statistics((k * get<Moment_statistics<Ids, include_covariance>>(a))...);
    }

    template <class Id>
    friend auto& get(Moment_statistics& x) {
        return get<Moment_statistics<Id, include_covariance>>(static_cast<base_type&>(x));
    }
    template <class Id>
    friend auto const& get(Moment_statistics const& x) {
        return get<Moment_statistics<Id, include_covariance>>(static_cast<base_type const&>(x));
    }
};









#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <limits>

template <typename T, class F = std::identity>
requires (std::is_arithmetic_v<std::decay_t<std::invoke_result_t<F, const T&>>>)
auto get_mean_Probits(std::vector<T> const& bootstrap_estimates, const std::set<double>& cis, F f = {}) {
    assert(!bootstrap_estimates.empty());

    std::vector<double> valid_estimates;
    valid_estimates.reserve(bootstrap_estimates.size()); // Prevent reallocations
    
    double sum = 0.0;
    
    // 1. Evaluate, filter NaNs, and Accumulate (O(N))
    for (const auto& x : bootstrap_estimates) {
        double value = static_cast<double>(std::invoke(f, x));
        if (!std::isnan(value)) {
            valid_estimates.push_back(value);
            sum += value;
        }
    }
    
    std::size_t n = valid_estimates.size();
    
    // Guard against an array entirely filled with NaNs
    if (n == 0) {
        return std::make_tuple(std::numeric_limits<double>::quiet_NaN(), std::map<double, double>{}, 0UL);
    }
    
    double mean = sum / static_cast<double>(n);
    
    // 2. Direct Sort (Cache-friendly, no NaN-logic required)
    std::sort(valid_estimates.begin(), valid_estimates.end());
    
    // 3. Map Percentiles
    std::map<double, double> out;

    for (double ci : cis) {
        double real_idx = std::clamp(ci * (n + 1) - 1.0, 0.0, static_cast<double>(n - 1));
        
        std::size_t low = static_cast<std::size_t>(std::floor(real_idx));
        std::size_t high = static_cast<std::size_t>(std::ceil(real_idx));
        
        double weight = real_idx - static_cast<double>(low);
        
        double v_low  = valid_estimates[low];
        double v_high = valid_estimates[high];

        out[ci] = (1.0 - weight) * v_low + weight * v_high;
    }

    return std::make_tuple(mean, std::move(out), n); 
}

template <typename T, class F = std::identity>
requires requires(T x, F f) { { std::invoke(f, x)() }; }
auto get_mean_Probits(std::vector<T> const& bootstrap_estimates, const std::set<double>& cis, F f = {}) {
    assert(!bootstrap_estimates.empty());

    // Recursive call: Dispatch to scalar, Matrix, or Vector overloads 
    // by unwrapping the Var/Functor layer.
    auto result = get_mean_Probits(bootstrap_estimates, cis, 
        [&f](const auto& x) -> decltype(auto) { 
            return std::invoke(f, x)(); 
        });

    // We can return the pair directly because the inner call already
    // produced the unwrapped value types (double, Matrix, etc.)
    return result; 
}

template<class M, class F = std::identity>
 requires (is_Matrix_v<std::decay_t<std::invoke_result_t<F, const M&>>>)
auto get_mean_Probits(
    std::vector<M> const& bootstrap_estimates, 
    const std::set<double>& cis, F f = {}) 
{
    assert(!bootstrap_estimates.empty());
    
    // Determine exactly what f returns
    using ProjResult = std::invoke_result_t<F, const M&>;
    using RawMatrix  = std::decay_t<ProjResult>;
    
    // 1. The Safe Storage Type: 
    // If f returns an lvalue reference, store a reference_wrapper. 
    // If f returns a value, store the raw matrix.
    using StorageType = std::conditional_t<
        std::is_lvalue_reference_v<ProjResult>,
        std::reference_wrapper<std::remove_reference_t<ProjResult>>,
        RawMatrix
    >;

    std::vector<StorageType> valid_estimates;
    valid_estimates.reserve(bootstrap_estimates.size()); // Prevent reallocations   
    
    // 2. Populate the valid estimates
    for (const auto& x : bootstrap_estimates) {
        decltype(auto) value = std::invoke(f, x);
        if (value.size() > 0) {
            if constexpr (std::is_lvalue_reference_v<ProjResult>) {
                valid_estimates.emplace_back(std::ref(value));
            } else {
                valid_estimates.emplace_back(std::move(value));
            }
        }
    }
    
    if (valid_estimates.empty()) {
        return std::make_tuple(RawMatrix{}, std::map<double, RawMatrix>{}, 0UL);
    }
    
    // 3. Extract the first matrix safely
    const RawMatrix& first_matrix = [] (const StorageType& item) -> const RawMatrix& {
        if constexpr (std::is_lvalue_reference_v<ProjResult>) return item.get();
        else return item;
    }(valid_estimates.front());

    RawMatrix mean_matrix(first_matrix.nrows(), first_matrix.ncols());
    std::map<double, RawMatrix> probits_map;
    for (double ci : cis) {
        probits_map.emplace(ci, RawMatrix(first_matrix.nrows(), first_matrix.ncols()));
    }

    std::size_t n = std::numeric_limits<std::size_t>::max();
    
    // 4. Delegate to the scalar 1D version
    for (std::size_t i = 0; i < mean_matrix.size(); ++i) {
        
        // Notice: `valid_estimates` already holds the matrices. 
        // We just need a lambda that extracts the i-th element.
        auto prob = get_mean_Probits(valid_estimates, cis, [i](const StorageType& mat_ref) { 
            if constexpr (std::is_lvalue_reference_v<ProjResult>) {
                return mat_ref.get()[i]; 
            } else {
                return mat_ref[i];
            }
        });
        
        mean_matrix[i] = std::get<0>(prob);
        for (const auto& [ci, value] : std::get<1>(prob)) {
            probits_map[ci][i] = value;  
        }
        
        if (std::get<2>(prob) < n) {
            n = std::get<2>(prob);
        }
    }
    
    if (n == std::numeric_limits<std::size_t>::max()) n = 0;

    return std::make_tuple(std::move(mean_matrix), std::move(probits_map), n);
}

template<class V, class F = std::identity> requires (is_of_this_template_type_v<std::decay_t<std::invoke_result_t<F, const V&>>, std::vector>) 
auto get_mean_Probits(
    std::vector<V> const& bootstrap_estimates, 
    const std::set<double>& cis, F f = {}) 
{
    assert(!bootstrap_estimates.empty());
    
    // Determine exactly what f returns
    using ProjResult = std::invoke_result_t<F, const V&>;
    using RawVector  = std::decay_t<ProjResult>;
    
    // Fallback: reference_wrapper for lvalues, raw vector for prvalues
    using StorageType = std::conditional_t<
        std::is_lvalue_reference_v<ProjResult>,
        std::reference_wrapper<std::remove_reference_t<ProjResult>>,
        RawVector
    >;

    std::vector<StorageType> valid_estimates;
    valid_estimates.reserve(bootstrap_estimates.size()); // Prevent reallocations   

    // 1. Evaluate, filter empty vectors, and store safely
    for (const auto& x : bootstrap_estimates) {
        decltype(auto) value = std::invoke(f, x);
        if (!value.empty()) {
            if constexpr (std::is_lvalue_reference_v<ProjResult>) {
                valid_estimates.emplace_back(std::ref(value));
            } else {
                valid_estimates.emplace_back(std::move(value));
            }
        }
    }

    // Guard against an array entirely filled with empty vectors
    if (valid_estimates.empty()) {
        return std::make_tuple(RawVector{}, std::map<double, RawVector>{}, 0UL);
    }

    // 2. Extract the first vector safely to get sizes
    const RawVector& first_vector = [] (const StorageType& item) -> const RawVector& {
        if constexpr (std::is_lvalue_reference_v<ProjResult>) return item.get();
        else return item;
    }(valid_estimates.front());

    // Check size consistency only on the valid vectors
    assert(std::all_of(valid_estimates.begin(), valid_estimates.end(), [&](const auto& item) {
        if constexpr (std::is_lvalue_reference_v<ProjResult>) return item.get().size() == first_vector.size();
        else return item.size() == first_vector.size();
    }));

    RawVector mean_vector(first_vector.size());
    std::map<double, RawVector> probits_map;
    for (double ci : cis) {
        probits_map.emplace(ci, RawVector(first_vector.size()));
    }

    std::size_t n = std::numeric_limits<std::size_t>::max();

    // 3. Delegate to the scalar 1D version
    for (std::size_t i = 0; i < mean_vector.size(); ++i) {
        
        auto prob = get_mean_Probits(valid_estimates, cis, [i](const StorageType& vec_ref) { 
            if constexpr (std::is_lvalue_reference_v<ProjResult>) {
                return vec_ref.get()[i]; 
            } else {
                return vec_ref[i];
            }
        });
        
        mean_vector[i] = std::get<0>(prob);
        for (const auto& [ci, value] : std::get<1>(prob)) {
            probits_map[ci][i] = value;
        }

        if (std::get<2>(prob) < n) {
            n = std::get<2>(prob);
        }
    }
    
    if (n == std::numeric_limits<std::size_t>::max()) n = 0;

    return std::make_tuple(std::move(mean_vector), std::move(probits_map), n);
}

template <class ParamIndexed, class F = std::identity>
 requires (is_of_this_template_type_v<std::decay_t<std::invoke_result_t<F, const ParamIndexed&>>, var::ParameterIndexed>)
auto get_mean_Probits(std::vector<ParamIndexed> const& bootstrap_estimates,
                      const std::set<double>& cis, F f = {}) 
{
    assert(!bootstrap_estimates.empty());

    using ProjResult = std::invoke_result_t<F, const ParamIndexed&>;
    using RawParamIdx = std::decay_t<ProjResult>;
    using ValueT = typename RawParamIdx::value_type;
    using Params = typename RawParamIdx::params_type;

    // Fallback: reference_wrapper for lvalues, raw object for prvalues
    using StorageType = std::conditional_t<
        std::is_lvalue_reference_v<ProjResult>,
        std::reference_wrapper<std::remove_reference_t<ProjResult>>,
        RawParamIdx
    >;

    std::vector<StorageType> valid_estimates;
    valid_estimates.reserve(bootstrap_estimates.size());

    const Params* metadata = nullptr;

    // 1. Evaluate and store safely (Zero-Copy)
    for (const auto& estimate : bootstrap_estimates) {
        decltype(auto) projected = std::invoke(f, estimate);
        
        if constexpr (std::is_lvalue_reference_v<ProjResult>) {
            valid_estimates.emplace_back(std::ref(projected));
        } else {
            valid_estimates.emplace_back(std::move(projected));
        }

        // Grab metadata from the first item that has it
        if (metadata == nullptr) {
            const RawParamIdx& raw_item = [&]() -> const RawParamIdx& {
                if constexpr (std::is_lvalue_reference_v<ProjResult>) return valid_estimates.back().get();
                else return valid_estimates.back();
            }();
            
            if (raw_item.has_parameters()) {
                metadata = raw_item.parameters_ptr();
            }
        }
    }

    // 2. Delegate to the underlying ValueT implementation.
    // Instead of copying the values, we project them on the fly!
    auto prob = get_mean_Probits(valid_estimates, cis, [](const StorageType& item) -> const ValueT& {
        if constexpr (std::is_lvalue_reference_v<ProjResult>) {
            return item.get().value();
        } else {
            return item.value();
        }
    });

    // 3. Re-wrap the results with the metadata
    var::ParameterIndexed<ValueT, Params> mean_out(std::move(std::get<0>(prob)), metadata);
    
    std::map<double, var::ParameterIndexed<ValueT, Params>> probits_out;
    for (auto& [level, val] : std::get<1>(prob)) {
        probits_out.emplace(level, var::ParameterIndexed<ValueT, Params>(std::move(val), metadata));
    }

    return std::make_tuple(std::move(mean_out), std::move(probits_out), std::get<2>(prob));
}

template<class F=std::identity, class VecSpace>
requires (is_of_this_template_type_v<std::decay_t<std::invoke_result_t<F, const VecSpace&>>, var::Vector_Space>)
auto get_mean_Probits(std::vector<VecSpace> const& bootstrap_estimates, const std::set<double>& cis, F f = {}) ;


namespace old{
template <class Id>
class Probit_statistics
    : public var::Var<Probit_statistics<Id>, 
    var::Vector_Space<mean<Id>,Probit_025<Id>, Probit_975<Id>>>
     {
   public:
    using base_type =var::Var<Probit_statistics<Id>, 
    var::Vector_Space<mean<Id>,Probit_025<Id>, Probit_975<Id>>>
    ;

    using value_type= value_type_t<mean<Id>>;
    
    static constexpr double LOWER_PERCENTILE = 0.025;
    static constexpr double UPPER_PERCENTILE = 0.975;

    Probit_statistics()
        : base_type{var::Vector_Space{mean<Id>{}, Probit_025<Id>{}, Probit_975<Id>{}}} {}
    
    Probit_statistics(value_type&& m, value_type&& p025, value_type&& p975):
        base_type{var::Vector_Space{mean<Id>(std::move( m)), Probit_025<Id>(std::move(p025)), Probit_975<Id>(std::move(p975))}} {}    

    Probit_statistics(std::vector<value_type>& values) {
        auto [mean_probit, probits] = get_mean_Probits(values, {LOWER_PERCENTILE, UPPER_PERCENTILE});
        *this = Probit_statistics(std::move(mean_probit), std::move(probits[0]), std::move(probits[1]));


    } 
    template<class VectorSpace, class F>
    requires requires (VectorSpace v, F f) { {f(v)} -> std::convertible_to<value_type_t<Id>>; }
    Probit_statistics(std::vector<VectorSpace>const& xs,  F&& f)
    {
        std::vector<value_type > values;
        values.reserve(xs.size());
        for (const auto& x : xs) {
            values.push_back(std::forward<F>(f)(x));
        }
        auto [mean_probit, probbits] = get_mean_Probits(values, {LOWER_PERCENTILE, UPPER_PERCENTILE});
        *this = Probit_statistics(std::move(mean_probit), std::move(probbits[0]), std::move(probbits[1]));
    }


    
    constexpr auto& operator[](var::Var<Id>) {
        return *this;
    }
    constexpr auto& operator[](var::Var<Id>) const {
        return *this;
    }
};
}

template <class Id>
class Probit_statistics
    : public var::Var<Probit_statistics<Id>, 
                      var::Vector_Space<mean<Id>, Probits<Id>, count<Id>>> 
{
   public:
    using base_type = var::Var<Probit_statistics<Id>, 
                               var::Vector_Space<mean<Id>, Probits<Id>, count<Id>>>;

    using value_type = value_type_t<mean<Id>>;
    using raw_value_type = value_type_t<Id>;
    
private:
    // Helper constructor to unpack the tuple directly into the main constructor
    template <typename TTuple>
    Probit_statistics(TTuple&& data)
        : Probit_statistics(std::move(std::get<0>(data)), 
                            std::move(std::get<1>(data)), 
                            std::get<2>(data)) {}

public:
    Probit_statistics()
        : base_type{var::Vector_Space{mean<Id>{}, Probits<Id>{}, count<Id>{}}} {}
    
    Probit_statistics(value_type&& m, std::map<double,value_type>&& pbits, std::size_t n)
        : base_type{var::Vector_Space{mean<Id>(std::move(m)), Probits<Id>(std::move(pbits)), count<Id>(n)}} {}

    template <class Dummy = void>
    requires (!std::is_same_v<raw_value_type, value_type> &&
              std::constructible_from<Id, raw_value_type>)
    Probit_statistics(raw_value_type&& m, std::map<double, raw_value_type>&& pbits, std::size_t n)
        : base_type{var::Vector_Space{
              mean<Id>(Id(std::move(m))),
              Probits<Id>([&pbits]() {
                  std::map<double, value_type> wrapped;
                  for (auto& [ci, value] : pbits) {
                      wrapped.emplace(ci, value_type(Id(std::move(value))));
                  }
                  return wrapped;
              }()),
              count<Id>(n)}} {}
      
    // Zero-copy delegation directly to the tuple constructor
    Probit_statistics(std::vector<value_type> const& values, const std::set<double>& cis) 
        : Probit_statistics(get_mean_Probits(values, cis)) {} 

    // Zero-copy functional pipeline! No intermediate vector allocations.
    template<class VectorSpace, class F>
    requires requires (VectorSpace v, F f) { { std::invoke(f, v) } -> std::convertible_to<value_type_t<Id>>; }
    Probit_statistics(std::vector<VectorSpace> const& xs, F f, const std::set<double>& cis)
        : Probit_statistics(get_mean_Probits(xs, cis, std::move(f))) {}

    constexpr auto& operator[](var::Var<Id>) {
        return *this;
    }
    constexpr auto& operator[](var::Var<Id>) const {
        return *this;
    }
};

template<class F, class VecSpace, class... Vs>
requires (std::is_same_v<std::decay_t<std::invoke_result_t<F, const VecSpace&>>, var::Vector_Space<Vs...>>)
auto get_mean_Probits_impl(std::vector<VecSpace> const& bootstrap_estimates, 
                           const std::set<double>& cis, F f, 
                           std::type_identity<var::Vector_Space<Vs...>>) 
{
    using ProjResult = std::invoke_result_t<F, const VecSpace&>;
    using RawVecSpace = std::decay_t<ProjResult>;
    
    using StorageType = std::conditional_t<
        std::is_lvalue_reference_v<ProjResult>,
        std::reference_wrapper<std::remove_reference_t<ProjResult>>,
        RawVecSpace
    >;

    std::vector<StorageType> valid_estimates;
    valid_estimates.reserve(bootstrap_estimates.size());
    
    // Evaluate 'f' EXACTLY ONCE per estimate
    for (const auto& v : bootstrap_estimates) {
        decltype(auto) val = std::invoke(f, v);
        if constexpr (std::is_lvalue_reference_v<ProjResult>) {
            valid_estimates.emplace_back(std::ref(val));
        } else {
            valid_estimates.emplace_back(std::move(val));
        }
    }

    auto probit = var::Vector_Space<Probit_statistics<Vs>...>(
        Probit_statistics<Vs>(
            valid_estimates,
            [](const StorageType& st) -> decltype(auto) { 
                if constexpr (std::is_lvalue_reference_v<ProjResult>) {
                    return get<Vs>(st.get())();
                } else {
                    return get<Vs>(st)();
                }
            },
            cis)...
    );

    auto meanout = var::Vector_Space<Vs...>(
        Vs(std::move(get<mean<Vs>>(get<Probit_statistics<Vs>>(probit)())()))...
    );
    
    // FIX: Extract the minimum count as a std::size_t instead of a Vector_Space
    std::size_t min_n = valid_estimates.size();
    (([&]{
        std::size_t current_n = get<count<Vs>>(get<Probit_statistics<Vs>>(probit)())();
        if (current_n < min_n) {
            min_n = current_n;
        }
    }()), ...);
    
    std::map<double, var::Vector_Space<Vs...>> probitsout;
    for (double ci : cis) {
        probitsout.emplace(
            ci,
            var::Vector_Space<Vs...>(
                Vs([&]() {
                    const auto& component = get<Probit_statistics<Vs>>(probit)();
                    const auto& probits_map = get<Probits<Vs>>(component)();
                    if (auto it = probits_map.find(ci); it != probits_map.end()) {
                        return it->second;
                    }
                    return get<mean<Vs>>(component)();
                }())...
            )
        );
    }
     
    // Return min_n so it perfectly matches the Probit_statistics tuple constructor
    return std::make_tuple(std::move(meanout), std::move(probitsout), min_n);
}



template<class F, class VecSpace>
requires (is_of_this_template_type_v<std::decay_t<std::invoke_result_t<F, const VecSpace&>>, var::Vector_Space>)
auto get_mean_Probits(std::vector<VecSpace> const& bootstrap_estimates, const std::set<double>& cis, F f) {
   using RVecSp = std::decay_t<std::invoke_result_t<F, const VecSpace&>>;
   return get_mean_Probits_impl(bootstrap_estimates, cis, std::move(f), std::type_identity<RVecSp>{});
}

#endif  // MOMENT_STATISTICS_H
