#include <catch_amalgamated.hpp>

#include <macrodr/cli/app/arguments.h>
#include <macrodr/cli/app/immediate_flags.h>
#include <macrodr/cli/app/workspace_persistence.h>
#include <macrodr/cli/app/working_directory.h>
#include <macrodr/cli/cli_parser.h>
#include <macrodr/cli/command_manager.h>
#include <macrodr/dsl/lexer_typed.h>

#include <filesystem>
#include <string>
#include <vector>

namespace app = macrodr::cli::app;

struct CwdGuard {
    std::filesystem::path original{std::filesystem::current_path()};
    ~CwdGuard() {
        std::error_code ec;
        std::filesystem::current_path(original, ec);
    }
};

TEST_CASE("collect_arguments copies argv content") {
    const char* argv[] = {"macro_dr", "--foo", "bar"};
    auto args = app::collect_arguments(3, const_cast<char**>(argv));
    REQUIRE(args.size() == 3);
    CHECK(args[0] == "macro_dr");
    CHECK(args[1] == "--foo");
    CHECK(args[2] == "bar");
}

TEST_CASE("handle_immediate_flags recognises help/version") {
    macrodr::cli::CliOptions opts;
    opts.help = true;
    auto exit_code = app::handle_immediate_flags(opts);
    CHECK(exit_code.has_value());
    CHECK(*exit_code == 0);

    opts.help = false;
    opts.version = true;
    exit_code = app::handle_immediate_flags(opts);
    CHECK(exit_code.has_value());
    CHECK(*exit_code == 0);

    opts.version = false;
    exit_code = app::handle_immediate_flags(opts);
    CHECK_FALSE(exit_code.has_value());
}

TEST_CASE("apply_working_directory succeeds and fails gracefully") {
    macrodr::cli::CliOptions opts;
    opts.has_chdir = false;
    CHECK(app::apply_working_directory(opts));

    opts.has_chdir = true;
    opts.chdir = std::filesystem::current_path().string();
    CHECK(app::apply_working_directory(opts));

    opts.chdir = "/path/that/does/not/exist";
    auto result = app::apply_working_directory(opts);
    CHECK_FALSE(result);
    CHECK(result.error()() == "Error: cannot chdir to '/path/that/does/not/exist'");
}

TEST_CASE("persist_run_workspace writes artefacts") {
    CwdGuard guard;

    auto temp_root = std::filesystem::temp_directory_path() / "macrodr_cli_helper_test";
    std::filesystem::remove_all(temp_root);
    std::filesystem::create_directories(temp_root);
    std::filesystem::current_path(temp_root);

    macrodr::cli::CliOptions opts;
    opts.verbosity = 1;

    const std::string script = "// test script\n";
    std::vector<std::string> args{"macro_dr", "--eval", "version()"};

    auto run_dir = app::persist_run_workspace(script, opts, args);
    REQUIRE(!run_dir.empty());
    CHECK(std::filesystem::exists(run_dir / "script.macroir"));
    CHECK(std::filesystem::exists(run_dir / "meta.json"));

    std::filesystem::remove_all(temp_root);
}

