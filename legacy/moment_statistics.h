#ifndef MOMENT_STATISTICS_H
#define MOMENT_STATISTICS_H

#include <maybe_error.h>
#include <algorithm>
#include <cassert>
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









template <typename T, class F = std::identity>
requires (std::is_arithmetic_v<std::decay_t<std::invoke_result_t<F, T>>>)
auto get_mean_Probits(std::vector<T> const& bootstrap_estimates, const std::set<double>& cis, F f = {}) {
    assert(!bootstrap_estimates.empty());

    // We use a local copy of 'f' to ensure it's available for all steps.
    // Projections are usually cheap to copy.
    
    // 1. Calculate Mean
    double sum = std::accumulate(bootstrap_estimates.begin(), bootstrap_estimates.end(), 0.0, 
        [&f](double acc, const T& x) {
            return acc + static_cast<double>(std::invoke(f, x));
        });
    double mean = sum / static_cast<double>(bootstrap_estimates.size());

    // 2. Indirect Sort
    std::vector<std::size_t> i_index(bootstrap_estimates.size());
    std::iota(i_index.begin(), i_index.end(), 0);

    std::sort(i_index.begin(), i_index.end(), [&](std::size_t i1, std::size_t i2) {
        return std::invoke(f, bootstrap_estimates[i1]) < std::invoke(f, bootstrap_estimates[i2]);
    });

    // 3. Map Percentiles
    std::size_t n = bootstrap_estimates.size();
    std::map<double, double> out;

    for (double ci : cis) {
        double real_idx = std::clamp(ci * (n + 1) - 1.0, 0.0, static_cast<double>(n - 1));
        
        std::size_t low = static_cast<std::size_t>(std::floor(real_idx));
        std::size_t high = static_cast<std::size_t>(std::ceil(real_idx));
        
        double weight = real_idx - static_cast<double>(low);
        
        double v_low  = static_cast<double>(std::invoke(f, bootstrap_estimates[i_index[low]]));
        double v_high = static_cast<double>(std::invoke(f, bootstrap_estimates[i_index[high]]));

        out[ci] = (1.0 - weight) * v_low + weight * v_high;
    }

    return std::make_pair(mean, std::move(out)); 
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

template<class M, class F>
requires (is_Matrix_v<std::decay_t<std::invoke_result_t<F, M>>>)
 auto get_mean_Probits(
    std::vector<M> const& bootstrap_estimates, 
    const std::set<double>& cis, F f=std::identity{}) 
{
    assert(!bootstrap_estimates.empty());
    assert( std::all_of(bootstrap_estimates.begin(), bootstrap_estimates.end(), [&bootstrap_estimates, &f](const auto& matrix) {
        return std::invoke(f, matrix).size() == std::invoke(f, bootstrap_estimates.front()).size();
    }));

    decltype (auto) first_matrix = std::invoke(f, bootstrap_estimates.front());
    using MatrixType= std::decay_t<decltype(std::invoke(f, bootstrap_estimates.front()))>;
    MatrixType mean_matrix(first_matrix.nrows(), first_matrix.ncols());
    std::map<double, MatrixType> probits_map;
    for (double ci : cis) {
        probits_map.emplace(ci, MatrixType(first_matrix.nrows(), first_matrix.ncols()));
    }

    for (std::size_t i = 0; i < mean_matrix.size(); ++i) {
            auto prob=get_mean_Probits( bootstrap_estimates, cis,[i, &f](const M& matrix) { return std::invoke(f, matrix)[i]; });
            mean_matrix[i] = prob.first;
            for (const auto& [ci, value] : prob.second) {
                probits_map[ci][i] = value;  
            }

        }
    return std::make_pair(std::move(mean_matrix), std::move(probits_map));
}

template<class V, class F>
requires (is_of_this_template_type_v<std::decay_t<std::invoke_result_t<F, V>>, std::vector>) 
auto get_mean_Probits(
    std::vector<V>const & bootstrap_estimates, 
    const std::set<double>& cis, F f=std::identity{}) 
{
    assert(!bootstrap_estimates.empty());
    assert(std::all_of(bootstrap_estimates.begin(), bootstrap_estimates.end(), [&](const auto& v) {
    return std::invoke(f, v).size() == std::invoke(f, bootstrap_estimates.front()).size();}));
    using VectorType= std::decay_t<std::invoke_result_t<F, V>>;
    decltype (auto) first_vector = std::invoke(f, bootstrap_estimates.front());
    VectorType mean_vector(first_vector.size());
    std::map<double, VectorType> probits_map;
    for (double ci : cis) {
        probits_map.emplace(ci, VectorType(first_vector.size()));
    }

    for (std::size_t i = 0; i < mean_vector.size(); ++i) {
        auto prob=get_mean_Probits( bootstrap_estimates, cis,[i, &f](const V& vec) { return std::invoke(f, vec)[i]; });
        mean_vector[i] = prob.first;
        for (const auto& [ci, value] : prob.second) {
            probits_map[ci][i] = value;
        }

    }
    

    return std::make_pair(std::move(mean_vector), std::move(probits_map));
}

template <class ParamIndexed, class F>
requires (is_of_this_template_type_v<std::decay_t<std::invoke_result_t<F, ParamIndexed>>, var::ParameterIndexed>)
auto get_mean_Probits(std::vector<ParamIndexed> const& bootstrap_estimates,
                             const std::set<double>& cis, F f=std::identity{}) {
    assert(!bootstrap_estimates.empty());

    using Projected = std::decay_t<std::invoke_result_t<F, ParamIndexed>>;
    using ValueT = typename Projected::value_type;
    using Params = typename Projected::params_type;

    std::vector<ValueT> values;
    values.reserve(bootstrap_estimates.size());

    const Params* metadata = nullptr;
    for (const auto& estimate : bootstrap_estimates) {
        const auto& projected = std::invoke(f, estimate);
        values.push_back(projected.value());
        if (metadata == nullptr && projected.has_parameters()) {
            metadata = projected.parameters_ptr();
        }
    }

    auto [mean_value, probits] = get_mean_Probits(values, cis, std::identity{});

    var::ParameterIndexed<ValueT, Params> mean_out(std::move(mean_value), metadata);
    std::map<double, var::ParameterIndexed<ValueT, Params>> probits_out;
    for (auto& [level, value] : probits) {
        probits_out.emplace(level,
                            var::ParameterIndexed<ValueT, Params>(std::move(value), metadata));
    }

    return std::make_pair(std::move(mean_out), std::move(probits_out));
}


template<class F, class VS>
requires (is_of_this_template_type_v<std::decay_t<std::invoke_result_t<F, VS>>, var::Vector_Space>)
auto get_mean_Probits(std::vector<VS> const& bootstrap_estimates, const std::set<double>& cis, F=std::identity{});


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
    var::Vector_Space<mean<Id>,Probits<Id>>>
     {
   public:
    using base_type =var::Var<Probit_statistics<Id>, 
    var::Vector_Space<mean<Id>,Probits<Id>>>
    ;

    using value_type= value_type_t<mean<Id>>;
    using raw_value_type = value_type_t<Id>;
    
    Probit_statistics()
        : base_type{var::Vector_Space{mean<Id>{}, Probits<Id>{}}} {}
    
    Probit_statistics(value_type&& m, std::map<double,value_type>&& pbits):
         base_type{var::Vector_Space{mean<Id>(std::move( m)), Probits<Id>(std::move(pbits))}} {}

    template <class Dummy = void>
    requires (!std::is_same_v<raw_value_type, value_type> &&
              std::constructible_from<Id, raw_value_type>)
    Probit_statistics(raw_value_type&& m, std::map<double, raw_value_type>&& pbits)
        : base_type{var::Vector_Space{
              mean<Id>(Id(std::move(m))),
              Probits<Id>([&pbits]() {
                  std::map<double, value_type> wrapped;
                  for (auto& [ci, value] : pbits) {
                      wrapped.emplace(ci, value_type(Id(std::move(value))));
                  }
                  return wrapped;
              }())}} {}
      
    Probit_statistics(std::vector<value_type> const& values, const std::set<double>& cis) {
        auto [mean_probit, probits] = get_mean_Probits(values, cis);
        *this = Probit_statistics(std::move(mean_probit), std::move(probits));
   } 
    template<class VectorSpace, class F>
    requires requires (VectorSpace v, F f) { {f(v)} -> std::convertible_to<value_type_t<Id>>; }
    Probit_statistics(std::vector<VectorSpace>const& xs,  F&& f, const std::set<double>& cis)
    {
        std::vector<value_type> values;
        values.reserve(xs.size());
        for (const auto& x : xs) {
            values.push_back(std::invoke(f, x));
        }
        auto [mean_probit, probbits] = get_mean_Probits(values, cis, std::identity{});
        *this = Probit_statistics(std::move(mean_probit), std::move(probbits));
    }


    
    constexpr auto& operator[](var::Var<Id> ) {
        return *this;
    }
    constexpr auto& operator[](var::Var<Id>) const {
        return *this;
    }
};


template<class F, class VecSpace, class...Vs>
requires (std::is_same_v<std::decay_t<std::invoke_result_t<F, VecSpace>>, var::Vector_Space<Vs...>>)
auto get_mean_Probits_impl(std::vector<VecSpace> const& bootstrap_estimates, const std::set<double>& cis, F f, 
    std::type_identity<var::Vector_Space<Vs...>>/**/){
     auto probit = var::Vector_Space<Probit_statistics<Vs>...>(
         Probit_statistics<Vs>(
             bootstrap_estimates,
             [&f](const auto& v) -> decltype(auto) { return get<Vs>(std::invoke(f, v))(); },
             cis)...);

     auto meanout=var::Vector_Space<Vs...>(Vs(std::move(get<mean<Vs>>(get<Probit_statistics<Vs>>(probit)())()))...);
     std::map<double,var::Vector_Space<Vs...>> probitsout;
     for (double ci : cis) {
         probitsout.emplace(
             ci,
             var::Vector_Space<Vs...>(Vs(
                 std::move(get<Probits<Vs>>(get<Probit_statistics<Vs>>(probit)())().at(ci)))...));
     }
     
     return std::make_pair(std::move(meanout), std::move(probitsout));
}

template<class F, class VecSpace>
requires (is_of_this_template_type_v<std::decay_t<std::invoke_result_t<F, VecSpace>>, var::Vector_Space>)
auto get_mean_Probits(std::vector<VecSpace> const& bootstrap_estimates, const std::set<double>& cis, F f){
   using RVecSp=std::decay_t<std::invoke_result_t<F, VecSpace>>;
   return get_mean_Probits_impl(bootstrap_estimates, cis, f, std::type_identity<RVecSp>{});
}



#endif  // MOMENT_STATISTICS_H
