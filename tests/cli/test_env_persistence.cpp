#include <catch_amalgamated.hpp>

#include <macrodr/cli/app/execute_program.h>
#include <macrodr/cli/app/workspace_persistence.h>
#include <macrodr/dsl/grammar_typed.h>
#include <macrodr/io/json/environment_io.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <variables.h>
#include <matrix.h>

namespace cli = macrodr::cli;

namespace {

std::filesystem::path make_unique_temp_dir() {
    auto base = std::filesystem::temp_directory_path();
    auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    for (int attempt = 0; attempt < 1000; ++attempt) {
        auto candidate = base /
                         ("macrodr_env_test-" + std::to_string(stamp + attempt));
        if (!std::filesystem::exists(candidate)) {
            return candidate;
        }
    }
    return base / "macrodr_env_test-fallback";
}

struct TempDirGuard {
    std::filesystem::path original{std::filesystem::current_path()};
    std::filesystem::path dir{make_unique_temp_dir()};

    TempDirGuard() {
        std::filesystem::create_directories(dir);
        std::filesystem::current_path(dir);
    }

    ~TempDirGuard() {
        std::error_code ec;
        std::filesystem::current_path(original, ec);
        std::filesystem::remove_all(dir, ec);
    }
};

auto make_identifier(std::string name) {
    auto maybe = macrodr::dsl::to_Identifier<macrodr::dsl::Lexer>(std::move(name));
    REQUIRE(maybe);
    return maybe.value();
}

template <class T>
auto* as_typed(const macrodr::dsl::base_typed_expression<macrodr::dsl::Lexer, macrodr::dsl::Compiler>* expr) {
    return dynamic_cast<const macrodr::dsl::typed_expression<macrodr::dsl::Lexer, macrodr::dsl::Compiler, T>*>(expr);
}

struct FooVar : public var::Var<FooVar, double> {
    using var::Var<FooVar, double>::Var;
    friend std::string className(FooVar) { return "FooVar"; }
};

struct BarVar : public var::Var<BarVar, Matrix<double>> {
    using var::Var<BarVar, Matrix<double>>::Var;
    friend std::string className(BarVar) { return "BarVar"; }
};

using SpaceVar = var::Vector_Space<FooVar, BarVar>;

static Matrix<double> make_matrix(std::size_t rows, std::size_t cols) {
    Matrix<double> m(rows, cols, false);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            m(i, j) = static_cast<double>(i * cols + j + 1);
        }
    }
    return m;
}

}  // namespace

TEST_CASE("environment IO round-trips supported primitives") {
    macrodr::dsl::Compiler compiler;
    macrodr::dsl::Environment<macrodr::dsl::Lexer, macrodr::dsl::Compiler> env(compiler);

    env.insert(make_identifier("d"),
               new macrodr::dsl::typed_literal<macrodr::dsl::Lexer, macrodr::dsl::Compiler, double>(1.5));
    env.insert(make_identifier("i"),
               new macrodr::dsl::typed_literal<macrodr::dsl::Lexer, macrodr::dsl::Compiler, int64_t>(-2));
    env.insert(make_identifier("b"),
               new macrodr::dsl::typed_literal<macrodr::dsl::Lexer, macrodr::dsl::Compiler, bool>(true));
    env.insert(make_identifier("s"),
               new macrodr::dsl::typed_literal<macrodr::dsl::Lexer, macrodr::dsl::Compiler, std::string>("hi"));

    TempDirGuard guard;
    auto output = std::filesystem::current_path() / "env.json";
    auto saved = macrodr::io::json::envio::save_environment_json(output, env, compiler);
    REQUIRE(saved);

    macrodr::dsl::Environment<macrodr::dsl::Lexer, macrodr::dsl::Compiler> loaded_env(compiler);
    auto loaded = macrodr::io::json::envio::load_environment_json(output, loaded_env);
    REQUIRE(loaded);

    auto check_value = [&](const char* name, auto expected) {
        auto maybe_expr = loaded_env.get(make_identifier(name));
        REQUIRE(maybe_expr);
        using ValueType = decltype(expected);
        auto typed = as_typed<ValueType>(maybe_expr.value());
        REQUIRE(typed != nullptr);
        auto result = typed->run(loaded_env);
        REQUIRE(result);
        CHECK(result.value() == expected);
    };

    check_value("d", 1.5);
    check_value("i", static_cast<int64_t>(-2));
    check_value("b", true);
    check_value("s", std::string{"hi"});
}

TEST_CASE("load replace clears prior environment") {
    TempDirGuard guard;
    const auto json_path = std::filesystem::current_path() / "load.json";
    std::ofstream(json_path)
        << R"({"variables":{"z":{"type":"double","value":3.0}}})";

    macrodr::dsl::Compiler compiler;
    macrodr::dsl::Environment<macrodr::dsl::Lexer, macrodr::dsl::Compiler> env(compiler);
    env.insert(make_identifier("x"),
               new macrodr::dsl::typed_literal<macrodr::dsl::Lexer, macrodr::dsl::Compiler, double>(1.0));

    auto result = macrodr::io::json::envio::load_environment_json(json_path, env, "replace");
    REQUIRE(result);

    auto missing = env.get(make_identifier("x"));
    CHECK_FALSE(missing);

    auto maybe_z = env.get(make_identifier("z"));
    REQUIRE(maybe_z);
    auto typed_z = as_typed<double>(maybe_z.value());
    REQUIRE(typed_z != nullptr);
    auto z_value = typed_z->run(env);
    REQUIRE(z_value);
    CHECK(z_value.value() == Catch::Approx(3.0));
}

TEST_CASE("load reports skipped unsupported types") {
    TempDirGuard guard;
    const auto json_path = std::filesystem::current_path() / "unsupported.json";
    std::ofstream(json_path)
        << R"({"variables":{"u":{"type":"custom","value":42}}})";

    macrodr::dsl::Compiler compiler;
    macrodr::dsl::Environment<macrodr::dsl::Lexer, macrodr::dsl::Compiler> env(compiler);
    auto loaded = macrodr::io::json::envio::load_environment_json(json_path, env);
    REQUIRE(loaded);
    CHECK(loaded.value().find("skipped 1 unsupported variable") != std::string::npos);
    auto maybe = env.get(make_identifier("u"));
    CHECK_FALSE(maybe);
}

TEST_CASE("environment IO preserves matrix and vector space literals") {
    TempDirGuard guard;

    macrodr::dsl::Compiler compiler;
    compiler.ensure_type_registered<Matrix<double>>();
    compiler.ensure_type_registered<FooVar>();
    compiler.ensure_type_registered<BarVar>();
    compiler.ensure_type_registered<SpaceVar>();

    macrodr::dsl::Environment<macrodr::dsl::Lexer, macrodr::dsl::Compiler> env(compiler);

    Matrix<double> mat = make_matrix(2, 2);
    env.insert(make_identifier("matrix"),
               new macrodr::dsl::typed_literal<macrodr::dsl::Lexer, macrodr::dsl::Compiler, Matrix<double>>(mat));

    SpaceVar space;
    static_cast<FooVar&>(space)() = 3.5;
    static_cast<BarVar&>(space)() = mat;
    env.insert(make_identifier("space"),
               new macrodr::dsl::typed_literal<macrodr::dsl::Lexer, macrodr::dsl::Compiler, SpaceVar>(space));

    const auto json_path = std::filesystem::current_path() / "complex_env.json";
    auto saved = macrodr::io::json::envio::save_environment_json(json_path, env, compiler);
    REQUIRE(saved);

    std::ifstream json_stream(json_path);
    REQUIRE(json_stream.is_open());
    std::string json_contents((std::istreambuf_iterator<char>(json_stream)), std::istreambuf_iterator<char>());
    auto parsed = macrodr::io::json::parse(json_contents);
    REQUIRE(parsed);
    const auto& root = parsed.value();
    const auto* vars = root.find("variables");
    REQUIRE(vars != nullptr);
    REQUIRE(vars->type == macrodr::io::json::Json::Type::Object);
    auto matrix_it = vars->obj.find("matrix");
    REQUIRE(matrix_it != vars->obj.end());
    CHECK(matrix_it->second.obj.find("value") != matrix_it->second.obj.end());
    CHECK(matrix_it->second.obj.at("value").type == macrodr::io::json::Json::Type::Array);
    auto space_it = vars->obj.find("space");
    REQUIRE(space_it != vars->obj.end());
    CHECK(space_it->second.obj.find("value") != space_it->second.obj.end());
    CHECK(space_it->second.obj.at("value").type == macrodr::io::json::Json::Type::Object);

    macrodr::dsl::Environment<macrodr::dsl::Lexer, macrodr::dsl::Compiler> loaded_env(compiler);
    auto loaded = macrodr::io::json::envio::load_environment_json(json_path, loaded_env);
    REQUIRE(loaded);

    auto loaded_matrix_expr = loaded_env.get(make_identifier("matrix"));
    REQUIRE(loaded_matrix_expr);
    auto matrix_typed = as_typed<Matrix<double>>(loaded_matrix_expr.value());
    REQUIRE(matrix_typed != nullptr);
    auto matrix_value = matrix_typed->run(loaded_env);
    REQUIRE(matrix_value);
    const auto& matrix_decoded = matrix_value.value();
    CHECK(matrix_decoded.nrows() == mat.nrows());
    CHECK(matrix_decoded.ncols() == mat.ncols());
    for (std::size_t i = 0; i < mat.nrows(); ++i)
        for (std::size_t j = 0; j < mat.ncols(); ++j)
            CHECK(matrix_decoded(i, j) == Catch::Approx(mat(i, j)));

    auto loaded_space_expr = loaded_env.get(make_identifier("space"));
    REQUIRE(loaded_space_expr);
    auto space_typed = as_typed<SpaceVar>(loaded_space_expr.value());
    REQUIRE(space_typed != nullptr);
    auto space_value = space_typed->run(loaded_env);
    REQUIRE(space_value);
    const auto& decoded_space = space_value.value();
    CHECK(static_cast<const FooVar&>(decoded_space)() == Catch::Approx(static_cast<FooVar&>(space)()));
    const auto& loaded_bar_mat = static_cast<const BarVar&>(decoded_space)();
    REQUIRE(loaded_bar_mat.nrows() == mat.nrows());
    REQUIRE(loaded_bar_mat.ncols() == mat.ncols());
    for (std::size_t i = 0; i < mat.nrows(); ++i)
        for (std::size_t j = 0; j < mat.ncols(); ++j)
            CHECK(loaded_bar_mat(i, j) == Catch::Approx(mat(i, j)));
}

TEST_CASE("CLI env-save end writes environment.json") {
    TempDirGuard guard;

    cli::CliOptions opts;
    opts.env_save = cli::CliOptions::EnvSaveMode::End;
    opts.has_env_save_path = true;
    opts.env_save_path = std::filesystem::current_path().string();

    const std::string script = "alpha = 1.0\n";
    std::vector<std::string> argv{"macro_dr"};

    auto run_dir = cli::app::persist_run_workspace(script, opts, argv);
    macrodr::dsl::Compiler compiler;
    auto exit_code = cli::app::execute_program(script, opts, compiler, run_dir);
    REQUIRE(exit_code == 0);

    const auto env_json = std::filesystem::current_path() / "environment.json";
    REQUIRE(std::filesystem::exists(env_json));

    std::ifstream env_stream(env_json);
    REQUIRE(env_stream.is_open());
    std::string env_contents((std::istreambuf_iterator<char>(env_stream)), std::istreambuf_iterator<char>());
    auto parsed = macrodr::io::json::parse(env_contents);
    REQUIRE(parsed);
    const auto& root = parsed.value();
    const auto* vars = root.find("variables");
    REQUIRE(vars != nullptr);
    REQUIRE(vars->type == macrodr::io::json::Json::Type::Object);
    auto it = vars->obj.find("alpha");
    REQUIRE(it != vars->obj.end());
    REQUIRE(it->second.type == macrodr::io::json::Json::Type::Object);
    auto type_it = it->second.obj.find("type");
    REQUIRE(type_it != it->second.obj.end());
    CHECK(type_it->second.str == "double");
    auto value_it = it->second.obj.find("value");
    REQUIRE(value_it != it->second.obj.end());
    REQUIRE(value_it->second.type == macrodr::io::json::Json::Type::Number);
    CHECK(value_it->second.num == Catch::Approx(1.0));
}

TEST_CASE("CLI env-save step writes snapshots") {
    TempDirGuard guard;

    cli::CliOptions opts;
    opts.env_save = cli::CliOptions::EnvSaveMode::Step;

    const std::string script = "x = 1.0\ny = 2.0\n";
    std::vector<std::string> argv{"macro_dr"};

    auto run_dir = cli::app::persist_run_workspace(script, opts, argv);
    macrodr::dsl::Compiler compiler;
    auto exit_code = cli::app::execute_program(script, opts, compiler, run_dir);
    REQUIRE(exit_code == 0);

    const auto snapshots = run_dir / "snapshots";
    REQUIRE(std::filesystem::exists(snapshots));
    REQUIRE(std::filesystem::exists(snapshots / "step_0.json"));
    REQUIRE(std::filesystem::exists(snapshots / "step_1.json"));
    REQUIRE(std::filesystem::exists(run_dir / "environment.json"));
}
