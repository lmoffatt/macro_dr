#include <catch_amalgamated.hpp>

#include <macrodr/dsl/function_builder.h>
#include <macrodr/dsl/grammar_typed.h>
#include <macrodr/dsl/grammar_untyped.h>
#include <macrodr/dsl/lexer_untyped.h>

#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace dsl = macrodr::dsl;

namespace {

auto make_identifier(std::string name) {
    auto maybe = dsl::to_Identifier<dsl::Lexer>(std::move(name));
    REQUIRE(maybe);
    return maybe.value();
}

template <class T>
auto* as_typed(
    const dsl::base_typed_expression<dsl::Lexer, dsl::Compiler>* expr) {
    return dynamic_cast<const dsl::typed_expression<dsl::Lexer, dsl::Compiler, T>*>(expr);
}

template <class T>
T read_value(const dsl::Environment<dsl::Lexer, dsl::Compiler>& env,
             const std::string& name) {
    auto maybe_expr = env.get_RunValue(make_identifier(name));
    REQUIRE(maybe_expr);
    auto typed = as_typed<T>(maybe_expr.value());
    REQUIRE(typed != nullptr);
    auto result = typed->run(env);
    REQUIRE(result);
    return result.value();
}

template <class T>
void bind_literal(dsl::Environment<dsl::Lexer, dsl::Compiler>& env, const std::string& name,
                  T value) {
    auto id = make_identifier(name);
    env.insert(id, new dsl::typed_literal<dsl::Lexer, dsl::Compiler, T>(value));
    env.push_back(id, new dsl::Identifier_compiler<dsl::Lexer, dsl::Compiler, T>(
                          new dsl::typed_literal<dsl::Lexer, dsl::Compiler, T>(std::move(value))));
}

int echo_int(int x) { return x; }
int read_const(const int& x) { return x; }
int bump(int& x) { return ++x; }
std::vector<int> echo_vector(std::vector<int> xs) { return xs; }
std::set<int> echo_set(std::set<int> xs) { return xs; }
std::tuple<int, int> echo_tuple(std::tuple<int, int> xs) { return xs; }

}  // namespace

TEST_CASE("DSL shared argument adaptation and homogeneous containers remain behavior-preserving") {
    dsl::Compiler compiler;
    REQUIRE(compiler.push_function("echo_int", dsl::to_typed_function<int>(&echo_int, "x")));
    REQUIRE(compiler.push_function("read_const",
                                   dsl::to_typed_function<const int&>(&read_const, "x")));
    REQUIRE(compiler.push_function("bump", dsl::to_typed_function<int&>(&bump, "x")));
    REQUIRE(compiler.push_function("echo_vector",
                                   dsl::to_typed_function<std::vector<int>>(&echo_vector, "xs")));
    REQUIRE(compiler.push_function("echo_set",
                                   dsl::to_typed_function<std::set<int>>(&echo_set, "xs")));
    REQUIRE(compiler.push_function("echo_tuple",
                                   dsl::to_typed_function<std::tuple<int, int>>(&echo_tuple,
                                                                                "xs")));

    dsl::Environment<dsl::Lexer, dsl::Compiler> env(compiler);
    bind_literal(env, "source", 9);
    bind_literal(env, "counter", 4);

    const std::string program_text =
        "literal_result = echo_int(x=7)\n"
        "identifier_result = echo_int(x=source)\n"
        "const_ref_result = read_const(x=counter)\n"
        "mut_ref_result = bump(x=counter)\n"
        "vector_result = echo_vector(xs=[3,1,3])\n"
        "set_result = echo_set(xs=[3,1,3])\n"
        "tuple_result = echo_tuple(xs={2,5})\n";

    auto parsed = dsl::extract_program(program_text);
    REQUIRE(parsed);

    auto compiled = dsl::compile_program(env, parsed.value());
    REQUIRE(compiled);

    auto executed = compiled.value().run(env);
    REQUIRE(executed);

    CHECK(read_value<int>(env, "literal_result") == 7);
    CHECK(read_value<int>(env, "identifier_result") == 9);
    CHECK(read_value<int>(env, "const_ref_result") == 4);
    CHECK(read_value<int>(env, "mut_ref_result") == 5);
    CHECK(read_value<int>(env, "counter") == 5);
    CHECK(read_value<std::vector<int>>(env, "vector_result") == std::vector<int>{3, 1, 3});
    CHECK(read_value<std::set<int>>(env, "set_result") == std::set<int>{1, 3});
    CHECK(read_value<std::tuple<int, int>>(env, "tuple_result") == std::tuple<int, int>{2, 5});
}
