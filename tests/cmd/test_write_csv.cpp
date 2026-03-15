#include <catch_amalgamated.hpp>

#include <macrodr/cmd/likelihood.h>
#include <macrodr/cmd/load_experiment.h>

#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct TempDirGuard {
    std::filesystem::path path;

    explicit TempDirGuard(std::string_view stem)
        : path(std::filesystem::temp_directory_path() /
               std::filesystem::path(std::string(stem) + "_macrodr")) {
        std::error_code ec;
        std::filesystem::remove_all(path, ec);
        std::filesystem::create_directories(path);
    }

    ~TempDirGuard() {
        std::error_code ec;
        std::filesystem::remove_all(path, ec);
    }
};

std::vector<std::string> read_lines(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("cannot open test csv: " + path.string());
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }
    return lines;
}

std::vector<std::string> split_csv_row(const std::string& row) {
    std::vector<std::string> columns;
    std::stringstream ss(row);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
        columns.push_back(cell);
    }
    if (!row.empty() && row.back() == ',') {
        columns.emplace_back();
    }
    return columns;
}

constexpr std::string_view kUnifiedHeader =
    "scope,simulation_index,sample_index,segment_index,sub_index,n_step,"
    "time_start,time_end,time_middle,agonist,patch_current,"
    "component_path,value_row,value_col,"
    "probit,calculus,statistic,quantile_level,"
    "param_index,param_col,param_name,"
    "value";

}  // namespace

TEST_CASE("generic write_csv serializes moment and probit statistics", "[write_csv]") {
    using namespace macrodr;
    using namespace macrodr::cmd;

    TempDirGuard dir("write_csv_summary");
    const auto base = dir.path / "summary";

    auto summary = var::Vector_Space{
        Moment_statistics<logL>(std::vector<logL>{logL(-2.0), logL(-1.0)}),
        Probit_statistics<logL>(-1.5, std::map<double, double>{{0.025, -2.5}, {0.975, -0.5}})};

    auto out = write_csv(summary, base.string());
    REQUIRE(out);

    const auto lines = read_lines(base.string() + ".csv");
    REQUIRE(lines.size() == 7);
    CHECK(lines.front() == kUnifiedHeader);

    CHECK(lines[1].find("summary") == 0);
    CHECK(lines[1].find(",point,primitive,count,") != std::string::npos);
    CHECK(lines[2].find(",point,primitive,mean,") != std::string::npos);
    CHECK(lines[3].find(",point,primitive,variance,") != std::string::npos);
    CHECK(lines[4].find(",mean,primitive,value,") != std::string::npos);

    bool saw_low_quantile = false;
    bool saw_high_quantile = false;
    for (const auto& line : lines) {
        saw_low_quantile =
            saw_low_quantile || line.find(",quantile,primitive,value,0.025,") != std::string::npos;
        saw_high_quantile =
            saw_high_quantile || line.find(",quantile,primitive,value,0.975,") != std::string::npos;
    }
    CHECK(saw_low_quantile);
    CHECK(saw_high_quantile);
}

TEST_CASE("recording write_csv uses the unified schema", "[write_csv]") {
    using namespace macrodr::cmd;

    TempDirGuard dir("write_csv_recording");
    const auto base = dir.path / "recording";

    auto experiment = create_experiment({{1, 2, 5.0}, {1, 3, 7.0}}, 10.0, 0.0, 1.0);
    auto recording = define_recording({0.1, 0.2});

    auto out = write_csv(experiment, recording, base.string());
    REQUIRE(out);

    const auto lines = read_lines(base.string() + ".csv");
    REQUIRE(lines.size() == 3);
    CHECK(lines.front() == kUnifiedHeader);

    const auto first = split_csv_row(lines[1]);
    const auto second = split_csv_row(lines[2]);
    REQUIRE(first.size() == 22);
    REQUIRE(second.size() == 22);

    CHECK(first[0] == "recording");
    CHECK(first[1].empty());
    CHECK(first[2] == "0");
    CHECK(first[3] == "0");
    CHECK(first[4].empty());
    CHECK(first[5] == "0");
    CHECK(first[9] == "5");
    CHECK(first[10] == "0.1");
    CHECK(first[11] == "patch_current");
    CHECK(first[14] == "point");
    CHECK(first[15] == "primitive");
    CHECK(first[16] == "value");
    CHECK(first[21] == "0.1");

    CHECK(second[0] == "recording");
    CHECK(second[2] == "1");
    CHECK(second[3] == "0");
    CHECK(second[5] == "1");
    CHECK(second[9] == "7");
    CHECK(second[10] == "0.2");
    CHECK(second[11] == "patch_current");
    CHECK(second[14] == "point");
    CHECK(second[15] == "primitive");
    CHECK(second[16] == "value");
    CHECK(second[21] == "0.2");
}
