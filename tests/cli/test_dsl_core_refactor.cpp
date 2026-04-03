#include <catch_amalgamated.hpp>

#include <macrodr/cmd/indexed_construction.h>
#include <macrodr/dsl/function_builder.h>
#include <macrodr/dsl/grammar_typed.h>
#include <macrodr/dsl/grammar_untyped.h>
#include <macrodr/dsl/lexer_untyped.h>

#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <indexed.h>

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
bool echo_bool(bool x) { return x; }
int read_const(const int& x) { return x; }
int bump(int& x) { return ++x; }
std::vector<int> echo_vector(std::vector<int> xs) { return xs; }
std::set<int> echo_set(std::set<int> xs) { return xs; }
std::tuple<int, int> echo_tuple(std::tuple<int, int> xs) { return xs; }
std::size_t echo_size(std::size_t x) { return x; }
std::size_t indexed_size(const var::Indexed<std::size_t>& xs) { return xs.size(); }

}  // namespace

TEST_CASE("DSL shared argument adaptation and homogeneous containers remain behavior-preserving") {
    dsl::Compiler compiler;
    REQUIRE(compiler.push_function("echo_int", dsl::to_typed_function<int>(&echo_int, "x")));
    REQUIRE(compiler.push_function("echo_bool", dsl::to_typed_function<bool>(&echo_bool, "x")));
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
        "bool_true_result = echo_bool(x=true)\n"
        "bool_false_result = echo_bool(x=false)\n"
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
    CHECK(read_value<bool>(env, "bool_true_result"));
    CHECK_FALSE(read_value<bool>(env, "bool_false_result"));
    CHECK(read_value<int>(env, "identifier_result") == 9);
    CHECK(read_value<int>(env, "const_ref_result") == 4);
    CHECK(read_value<int>(env, "mut_ref_result") == 5);
    CHECK(read_value<int>(env, "counter") == 5);
    CHECK(read_value<std::vector<int>>(env, "vector_result") == std::vector<int>{3, 1, 3});
    CHECK(read_value<std::set<int>>(env, "set_result") == std::set<int>{1, 3});
    CHECK(read_value<std::tuple<int, int>>(env, "tuple_result") == std::tuple<int, int>{2, 5});
}

TEST_CASE("DSL keeps identifiers beginning with true or false distinct from boolean literals") {
    dsl::Compiler compiler;
    REQUIRE(compiler.push_function("echo_bool", dsl::to_typed_function<bool>(&echo_bool, "x")));

    dsl::Environment<dsl::Lexer, dsl::Compiler> env(compiler);
    bind_literal(env, "true_value", true);
    bind_literal(env, "false_value", false);

    const std::string program_text =
        "from_identifier = echo_bool(x=true_value)\n"
        "from_literal = echo_bool(x=true)\n"
        "from_false_identifier = echo_bool(x=false_value)\n"
        "from_false_literal = echo_bool(x=false)\n";

    auto parsed = dsl::extract_program(program_text);
    REQUIRE(parsed);

    auto compiled = dsl::compile_program(env, parsed.value());
    REQUIRE(compiled);

    auto executed = compiled.value().run(env);
    REQUIRE(executed);

    CHECK(read_value<bool>(env, "from_identifier"));
    CHECK(read_value<bool>(env, "from_literal"));
    CHECK_FALSE(read_value<bool>(env, "from_false_identifier"));
    CHECK_FALSE(read_value<bool>(env, "from_false_literal"));
}

TEST_CASE("DSL lifts scalar functions over indexed arguments and prefers exact indexed overloads") {
    dsl::Compiler compiler;
    REQUIRE(compiler.push_function(
        "axis", dsl::to_typed_function<std::string, std::vector<std::string>>(
                    &macrodr::cmd::axis, "name", "labels")));
    REQUIRE(compiler.push_function(
        "indexed_by", dsl::to_typed_function<var::Axis, std::vector<std::size_t>>(
                          &macrodr::cmd::indexed_by<std::size_t>, "axis", "values")));
    REQUIRE(compiler.push_function("identity_size",
                                   dsl::to_typed_function<std::size_t>(&echo_size, "x")));
    REQUIRE(compiler.push_function("shape", dsl::to_typed_function<std::size_t>(&echo_size, "x")));
    REQUIRE(compiler.push_function(
        "shape", dsl::to_typed_function<const var::Indexed<std::size_t>&>(&indexed_size, "x")));

    dsl::Environment<dsl::Lexer, dsl::Compiler> env(compiler);

    const std::string program_text =
        "model_axis = axis(name=\"models\", labels=[\"scheme_CO\",\"scheme_CCO\"])\n"
        "n_steps_axis = axis(name=\"n_substeps\", labels=[\"10\",\"30\",\"100\"])\n"
        "weights = indexed_by(axis=model_axis, values=[1,2])\n"
        "n_steps = indexed_by(axis=n_steps_axis, values=[10,30,100])\n"
        "inline_steps = indexed_by(axis=axis(name=\"inline\", labels=[\"A\",\"B\"]), values=[4,8])\n"
        "lifted = identity_size(x=weights)\n"
        "shape_result = shape(x=weights)\n";

    auto parsed = dsl::extract_program(program_text);
    REQUIRE(parsed);

    auto compiled = dsl::compile_program(env, parsed.value());
    REQUIRE(compiled);

    auto executed = compiled.value().run(env);
    REQUIRE(executed);

    auto model_axis = read_value<var::Axis>(env, "model_axis");
    CHECK(model_axis.m_id.idName == "models");
    CHECK(model_axis.m_labels == std::vector<std::string>{"scheme_CO", "scheme_CCO"});

    auto lifted = read_value<var::Indexed<std::size_t>>(env, "lifted");
    REQUIRE(lifted.values().size() == 2);
    CHECK(lifted.index_space().m_axes.size() == 1);
    CHECK(lifted.index_space().m_axes[0].m_id.idName == "models");
    CHECK(lifted.index_space().m_axes[0].m_labels == std::vector<std::string>{"scheme_CO",
                                                                               "scheme_CCO"});
    CHECK(lifted.values() == std::vector<std::size_t>{1, 2});

    auto n_steps = read_value<var::Indexed<std::size_t>>(env, "n_steps");
    REQUIRE(n_steps.values().size() == 3);
    CHECK(n_steps.index_space().m_axes.size() == 1);
    CHECK(n_steps.index_space().m_axes[0].m_id.idName == "n_substeps");
    CHECK(n_steps.index_space().m_axes[0].m_labels ==
          std::vector<std::string>{"10", "30", "100"});
    CHECK(n_steps.values() == std::vector<std::size_t>{10, 30, 100});

    auto inline_steps = read_value<var::Indexed<std::size_t>>(env, "inline_steps");
    REQUIRE(inline_steps.values().size() == 2);
    CHECK(inline_steps.index_space().m_axes[0].m_id.idName == "inline");
    CHECK(inline_steps.index_space().m_axes[0].m_labels == std::vector<std::string>{"A", "B"});
    CHECK(inline_steps.values() == std::vector<std::size_t>{4, 8});

    CHECK(read_value<std::size_t>(env, "shape_result") == 2);
}

TEST_CASE("DSL typed indexed constructors force bool and int payloads") {
    dsl::Compiler compiler;
    REQUIRE(compiler.push_function(
        "axis", dsl::to_typed_function<std::string, std::vector<std::string>>(
                    &macrodr::cmd::axis, "name", "labels")));
    REQUIRE(compiler.push_function(
        "indexed_bool_by", dsl::to_typed_function<var::Axis, std::vector<bool>>(
                               &macrodr::cmd::indexed_by<bool>, "axis", "values")));
    REQUIRE(compiler.push_function(
        "indexed_int_by", dsl::to_typed_function<var::Axis, std::vector<int>>(
                              &macrodr::cmd::indexed_by<int>, "axis", "values")));
    REQUIRE(compiler.push_function("echo_bool", dsl::to_typed_function<bool>(&echo_bool, "x")));
    REQUIRE(compiler.push_function("echo_int", dsl::to_typed_function<int>(&echo_int, "x")));

    dsl::Environment<dsl::Lexer, dsl::Compiler> env(compiler);

    const std::string program_text =
        "algorithm_axis = axis(name=\"algorithm\", labels=[\"macro_MRV\",\"macro_IRV\",\"macro_RV\"])\n"
        "recursive = indexed_bool_by(axis=algorithm_axis, values=[false,true,true])\n"
        "averaging = indexed_int_by(axis=algorithm_axis, values=[1,1,2])\n"
        "recursive_echo = echo_bool(x=recursive)\n"
        "averaging_echo = echo_int(x=averaging)\n";

    auto parsed = dsl::extract_program(program_text);
    REQUIRE(parsed);

    auto compiled = dsl::compile_program(env, parsed.value());
    REQUIRE(compiled);

    auto executed = compiled.value().run(env);
    REQUIRE(executed);

    auto recursive = read_value<var::Indexed<bool>>(env, "recursive");
    CHECK(recursive.index_space().m_axes[0].m_id.idName == "algorithm");
    CHECK(recursive.index_space().m_axes[0].m_labels ==
          std::vector<std::string>{"macro_MRV", "macro_IRV", "macro_RV"});
    CHECK(recursive.values() == std::vector<bool>{false, true, true});

    auto averaging = read_value<var::Indexed<int>>(env, "averaging");
    CHECK(averaging.index_space().m_axes[0].m_id.idName == "algorithm");
    CHECK(averaging.values() == std::vector<int>{1, 1, 2});

    auto recursive_echo = read_value<var::Indexed<bool>>(env, "recursive_echo");
    CHECK(recursive_echo.values() == std::vector<bool>{false, true, true});

    auto averaging_echo = read_value<var::Indexed<int>>(env, "averaging_echo");
    CHECK(averaging_echo.values() == std::vector<int>{1, 1, 2});
}

TEST_CASE("DSL indexed_by rejects string axis lookup") {
    dsl::Compiler compiler;
    REQUIRE(compiler.push_function(
        "axis", dsl::to_typed_function<std::string, std::vector<std::string>>(
                    &macrodr::cmd::axis, "name", "labels")));
    REQUIRE(compiler.push_function(
        "indexed_by", dsl::to_typed_function<var::Axis, std::vector<std::size_t>>(
                          &macrodr::cmd::indexed_by<std::size_t>, "axis", "values")));

    dsl::Environment<dsl::Lexer, dsl::Compiler> env(compiler);

    const std::string program_text =
        "weights = indexed_by(axis=\"models\", values=[1,2])\n";

    auto parsed = dsl::extract_program(program_text);
    REQUIRE(parsed);

    auto compiled = dsl::compile_program(env, parsed.value());
    REQUIRE_FALSE(compiled);
    CHECK(compiled.error()().find("indexed_by") != std::string::npos);
}

TEST_CASE("DSL no longer exposes indexed(name, labels, values)") {
    dsl::Compiler compiler;
    REQUIRE(compiler.push_function(
        "axis", dsl::to_typed_function<std::string, std::vector<std::string>>(
                    &macrodr::cmd::axis, "name", "labels")));
    REQUIRE(compiler.push_function(
        "indexed_by", dsl::to_typed_function<var::Axis, std::vector<std::size_t>>(
                          &macrodr::cmd::indexed_by<std::size_t>, "axis", "values")));

    dsl::Environment<dsl::Lexer, dsl::Compiler> env(compiler);

    const std::string program_text =
        "weights = indexed(name=\"models\", labels=[\"scheme_CO\",\"scheme_CCO\"], values=[1,2])\n";

    auto parsed = dsl::extract_program(program_text);
    REQUIRE(parsed);

    auto compiled = dsl::compile_program(env, parsed.value());
    REQUIRE_FALSE(compiled);
    CHECK(compiled.error()().find("indexed") != std::string::npos);
}
