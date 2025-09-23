#include <macrodr/cli/app/workspace_persistence.h>

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace macrodr::cli::app {

std::filesystem::path persist_run_workspace(const std::string& script,
                                            const macrodr::cli::CliOptions& opts,
                                            const std::vector<std::string>& arguments) {
    const auto now = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &now);
#else
    tm = *std::localtime(&now);
#endif
    std::ostringstream oss;
    oss << std::put_time(&tm, "run-%Y%m%d-%H%M%S");
    const auto run_dir = std::filesystem::path("runs") / oss.str();

    std::error_code ec;
    std::filesystem::create_directories(run_dir, ec);
    if (ec) {
        std::cerr << "Warning: could not create run directory '" << run_dir.string()
                  << "' (" << ec.message() << ")\n";
        return {};
    }

    if (opts.verbosity > 0) {
        std::cout << "[macro_dr] run dir: " << run_dir.string() << '\n';
    }

    {
        std::ofstream script_out(run_dir / "script.macroir");
        script_out << script;
    }

    {
        std::ofstream meta(run_dir / "meta.json");
        meta << "{\n";
        meta << "  \"cwd\": \"" << std::filesystem::current_path().string() << "\",\n";
        meta << "  \"timestamp\": " << static_cast<long long>(now) << ",\n";
        meta << "  \"args\": [";
        for (std::size_t i = 0; i < arguments.size(); ++i) {
            meta << (i == 0 ? "" : ", ") << "\"" << arguments[i] << "\"";
        }
        meta << "]\n";
        meta << "}\n";
    }

    return run_dir;
}

}  // namespace macrodr::cli::app

