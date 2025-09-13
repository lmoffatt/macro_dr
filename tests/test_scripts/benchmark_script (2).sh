#!/bin/bash

# benchmark_compilation.sh
# Script para medir tiempos de compilación de templates

echo "=== BENCHMARK DE COMPILACIÓN - TEMPLATE ALGEBRA ==="
echo "Fecha: $(date)"
echo "Compilador: $(g++ --version | head -1)"
echo "Sistema: $(uname -a)"
echo ""

# Crear archivos de prueba separados
cat > original_test.cpp << 'EOF'
#include <type_traits>

template <typename C, template <typename...> class V>
struct is_of_this_template_type : std::false_type {};
template <template <typename...> class V, typename... Ts>
struct is_of_this_template_type<V<Ts...>, V> : std::true_type {};
template <typename C, template <typename...> class V>
inline constexpr bool is_of_this_template_type_v =
    is_of_this_template_type<std::decay_t<C>, V>::value;

template <template <class T0, class T1, class... Ts> class constexpr_P,
          template <class T0, class T1, class... Ts> class constexpr_S>
struct algebra_2 {
    template<class T>
    inline static constexpr bool is_regular = !is_of_this_template_type_v<T,constexpr_S>&&
                                               !is_of_this_template_type_v<T,constexpr_P>;
    template <class X, class Y, class... Us>
    struct P_constexpr_impl;
    template <class X, class Y, class... Us>
    using P_constexpr=typename P_constexpr_impl<X,Y,Us...>::type;
    
    template <class X, class Y, class... Us>
    struct S_constexpr_impl;
    template <class X, class Y, class... Us>
    using S_constexpr=typename S_constexpr_impl<X,Y,Us...>::type;
    
    template <class X, class Y, class Z,class... Us>
    struct P_constexpr_impl<X,Y,Z,Us...>{
        using type=  P_constexpr<P_constexpr<X, Y>, Z, Us...>;
    };   
    
    template<class X, class Y, class Z, class... Us>
    struct  S_constexpr_impl<X, Y, Z, Us...> {
        using type= S_constexpr<S_constexpr<X, Y>, Z, Us...>;
    };
    
    template <class X, class Y>
    requires (is_regular<X> && is_regular<Y>)
    struct P_constexpr_impl<X, Y> {
        using type= constexpr_P<X, Y>;
    };
    
    template <class X, class Y0, class Y1, class... Ys>
    struct  P_constexpr_impl<X, constexpr_P<Y0, Y1, Ys...>> {
        using type= constexpr_P<X, Y0, Y1, Ys...>;
    };
    
    template <class X0, class X1, class... Xs, class Y>
    struct P_constexpr_impl<constexpr_P<X0, X1, Xs...>, Y> {
        using type= constexpr_P<X0, X1, Xs..., Y>;
    };
    
    template <class X0, class X1, class... Xs, class Y0, class Y1, class... Ys>
    struct P_constexpr_impl<constexpr_P<X0, X1, Xs...>, constexpr_P<Y0, Y1, Ys...>> {
        using type =constexpr_P<X0, X1, Xs..., Y0, Y1, Ys...>;
    };
    
    template <class X, class Y>
    struct S_constexpr_impl<X, Y> {
        using type =constexpr_S<X, Y>;
    };
    
    template <class X, class Y0, class Y1, class... Ys>
    struct S_constexpr_impl<X, constexpr_S<Y0, Y1, Ys...>> {
        using type =constexpr_S<X, Y0, Y1, Ys...>;
    };
    template <class X0, class X1, class... Xs, class Y>
    struct S_constexpr_impl<constexpr_S<X0, X1, Xs...>, Y> {
        using type =constexpr_S<X0, X1, Xs..., Y>;
    };
    
    template <class X0, class X1, class... Xs, class Y0, class Y1, class... Ys>
    struct S_constexpr_impl<constexpr_S<X0, X1, Xs...>, constexpr_S<Y0, Y1, Ys...>> {
        using type =constexpr_S<X0, X1, Xs..., Y0, Y1, Ys...>;
    };
    
    template <class X0, class X1, class... Xs, class Y>
    struct P_constexpr_impl<constexpr_S<X0, X1, Xs...>, Y> {
        using type =S_constexpr<P_constexpr<X0, Y>, P_constexpr<X1, Y>,
                           P_constexpr<Xs, Y>...>;
    };
    
    template <class X, class Y0, class Y1, class... Ys>
        requires (!is_of_this_template_type_v<X, constexpr_S>)
    struct P_constexpr_impl<X, constexpr_S<Y0, Y1, Ys...>> {
        using type =S_constexpr<P_constexpr<X, Y0>, P_constexpr<X, Y1>,
                           P_constexpr<X, Ys>...>;
    };
};

template<class T0, class T1, class... Ts> struct TestP {};
template<class T0, class T1, class... Ts> struct TestS {};
struct A {}; struct B {}; struct C {}; struct D {}; struct E {};

using Algebra = algebra_2<TestP, TestS>;
using Test1 = Algebra::P_constexpr<A, B>;
using Test2 = Algebra::S_constexpr<A, B, C>;
using Test3 = Algebra::P_constexpr<TestP<A, B>, C>;
using Test4 = Algebra::P_constexpr<TestS<A, B>, C>;
using Test5 = Algebra::P_constexpr<A, TestS<B, C, D>>;
using Complex = Algebra::P_constexpr<Algebra::S_constexpr<A, B, C>, Algebra::P_constexpr<D, E>>;

int main() { return 0; }
EOF

cat > optimized_test.cpp << 'EOF'
#include <type_traits>

template <typename C, template <typename...> class V>
struct is_template_type : std::false_type {};
template <template <typename...> class V, typename... Ts>
struct is_template_type<V<Ts...>, V> : std::true_type {};
template <typename C, template <typename...> class V>
inline constexpr bool is_template_type_v = is_template_type<std::decay_t<C>, V>::value;

template <template <class T0, class T1, class... Ts> class P_tmpl,
          template <class T0, class T1, class... Ts> class S_tmpl>
struct algebra_2_opt {
    template<class T>
    static constexpr bool is_regular = !is_template_type_v<T, S_tmpl> && !is_template_type_v<T, P_tmpl>;

    template <class X, class Y, class... Us> struct P_impl;
    template <class X, class Y, class... Us> struct S_impl;
    
    template <class X, class Y, class... Us>
    using P = typename P_impl<X, Y, Us...>::type;
    template <class X, class Y, class... Us>
    using S = typename S_impl<X, Y, Us...>::type;
    
    template <class X, class Y>
    requires (is_regular<X> && is_regular<Y>)
    struct P_impl<X, Y> { using type = P_tmpl<X, Y>; };
    
    template <class X, class Y>
    struct S_impl<X, Y> { using type = S_tmpl<X, Y>; };
    
    template <class X, class Y, class Z, class... Us>
    struct P_impl<X, Y, Z, Us...> { using type = P<P<X, Y>, Z, Us...>; };
    
    template<class X, class Y, class Z, class... Us>
    struct S_impl<X, Y, Z, Us...> { using type = S<S<X, Y>, Z, Us...>; };
    
    template <class X, class Y0, class Y1, class... Ys>
    struct P_impl<X, P_tmpl<Y0, Y1, Ys...>> { using type = P_tmpl<X, Y0, Y1, Ys...>; };
    
    template <class X0, class X1, class... Xs, class Y>
    struct P_impl<P_tmpl<X0, X1, Xs...>, Y> { using type = P_tmpl<X0, X1, Xs..., Y>; };
    
    template <class X0, class X1, class... Xs, class Y0, class Y1, class... Ys>
    struct P_impl<P_tmpl<X0, X1, Xs...>, P_tmpl<Y0, Y1, Ys...>> { 
        using type = P_tmpl<X0, X1, Xs..., Y0, Y1, Ys...>; 
    };
    
    template <class X, class Y0, class Y1, class... Ys>
    struct S_impl<X, S_tmpl<Y0, Y1, Ys...>> { using type = S_tmpl<X, Y0, Y1, Ys...>; };
    
    template <class X0, class X1, class... Xs, class Y>
    struct S_impl<S_tmpl<X0, X1, Xs...>, Y> { using type = S_tmpl<X0, X1, Xs..., Y>; };
    
    template <class X0, class X1, class... Xs, class Y0, class Y1, class... Ys>
    struct S_impl<S_tmpl<X0, X1, Xs...>, S_tmpl<Y0, Y1, Ys...>> { 
        using type = S_tmpl<X0, X1, Xs..., Y0, Y1, Ys...>; 
    };
    
    template <class X0, class X1, class... Xs, class Y>
    struct P_impl<S_tmpl<X0, X1, Xs...>, Y> {
        using type = S<P<X0, Y>, P<X1, Y>, P<Xs, Y>...>;
    };
    
    template <class X, class Y0, class Y1, class... Ys>
    requires (!is_template_type_v<X, S_tmpl>)
    struct P_impl<X, S_tmpl<Y0, Y1, Ys...>> {
        using type = S<P<X, Y0>, P<X, Y1>, P<X, Ys>...>;
    };
};

template<class T0, class T1, class... Ts> struct TestP {};
template<class T0, class T1, class... Ts> struct TestS {};
struct A {}; struct B {}; struct C {}; struct D {}; struct E {};

using Algebra = algebra_2_opt<TestP, TestS>;
using Test1 = Algebra::P<A, B>;
using Test2 = Algebra::S<A, B, C>;
using Test3 = Algebra::P<TestP<A, B>, C>;
using Test4 = Algebra::P<TestS<A, B>, C>;
using Test5 = Algebra::P<A, TestS<B, C, D>>;
using Complex = Algebra::P<Algebra::S<A, B, C>, Algebra::P<D, E>>;

int main() { return 0; }
EOF

# Función para medir tiempo de compilación
measure_compilation() {
    local file=$1
    local description=$2
    local flags=$3
    
    echo "Compilando $description con flags: $flags"
    
    # Limpiar archivos anteriores
    rm -f a.out
    
    # Medir tiempo de compilación
    local start_time=$(date +%s.%N)
    g++ $flags -std=c++20 "$file" 2>/dev/null
    local end_time=$(date +%s.%N)
    
    # Calcular tiempo transcurrido
    local elapsed=$(echo "$end_time - $start_time" | bc -l 2>/dev/null || python3 -c "print($end_time - $start_time)")
    
    if [ $? -eq 0 ]; then
        printf "  ✓ Tiempo: %.3f segundos\n" "$elapsed"
        return 0
    else
        echo "  ✗ Error de compilación"
        return 1
    fi
}

# Probar diferentes niveles de optimización
FLAGS_SETS=(
    "-O0"
    "-O1"
    "-O2"
    "-O3"
    "-O0 -ftime-report"
    "-O2 -ftime-report"
)

echo "--- COMPILACIÓN VERSIÓN ORIGINAL ---"
for flags in "${FLAGS_SETS[@]}"; do
    measure_compilation "original_test.cpp" "Original" "$flags"
done

echo ""
echo "--- COMPILACIÓN VERSIÓN OPTIMIZADA ---"
for flags in "${FLAGS_SETS[@]}"; do
    measure_compilation "optimized_test.cpp" "Optimizada" "$flags"
done

# Test de stress con templates más complejos
echo ""
echo "--- TEST DE STRESS ---"

cat > stress_test.cpp << 'EOF'
#include <type_traits>

template <typename C, template <typename...> class V>
struct is_template_type : std::false_type {};
template <template <typename...> class V, typename... Ts>
struct is_template_type<V<Ts...>, V> : std::true_type {};
template <typename C, template <typename...> class V>
inline constexpr bool is_template_type_v = is_template_type<std::decay_t<C>, V>::value;

template <template <class T0, class T1, class... Ts> class P_tmpl,
          template <class T0, class T1, class... Ts> class S_tmpl>
struct algebra_2_opt {
    template<class T>
    static constexpr bool is_regular = !is_template_type_v<T, S_tmpl> && !is_template_type_v<T, P_tmpl>;

    template <class X, class Y, class... Us> struct P_impl;
    template <class X, class Y, class... Us> struct S_impl;
    
    template <class X, class Y, class... Us>
    using P = typename P_impl<X, Y, Us...>::type;
    template <class X, class Y, class... Us>
    using S = typename S_impl<X, Y, Us...>::type;
    
    template <class X, class Y>
    requires (is_regular<X> && is_regular<Y>)
    struct P_impl<X, Y> { using type = P_tmpl<X, Y>; };
    
    template <class X, class Y>
    struct S_impl<X, Y> { using type = S_tmpl<X, Y>; };
    
    template <class X, class Y, class Z, class... Us>
    struct P_impl<X, Y, Z, Us...> { using type = P<P<X, Y>, Z, Us...>; };
    
    template<class X, class Y, class Z, class... Us>
    struct S_impl<X, Y, Z, Us...> { using type = S<S<X, Y>, Z, Us...>; };
    
    template <class X, class Y0, class Y1, class... Ys>
    struct P_impl<X, P_tmpl<Y0, Y1, Ys...>> { using type = P_tmpl<X, Y0, Y1, Ys...>; };
    
    template <class X0, class X1, class... Xs, class Y>
    struct P_impl<P_tmpl<X0, X1, Xs...>, Y> { using type = P_tmpl<X0,