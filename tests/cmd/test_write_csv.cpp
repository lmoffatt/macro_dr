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
    "step_start,step_end,step_middle,agonist,patch_current,"
    "component_path,value_row,value_col,"
    "probit,calculus,statistic,quantile_level,"
    "param_index,param_col,param_name,"
    "value";

SymPosDefMatrix<double> make_diag_spd(std::initializer_list<double> diag_vals) {
    const std::size_t n = diag_vals.size();
    SymPosDefMatrix<double> out(n, n, false);
    std::size_t i = 0;
    for (auto v : diag_vals) {
        out.set(i, i, v);
        ++i;
    }
    return out;
}

const var::Parameters_Transformations& test_parameter_transformations() {
    static const auto transformations = [] {
        auto maybe_tr = var::MyTranformations::from_strings({"Linear", "Linear"});
        if (!maybe_tr) {
            throw std::runtime_error(maybe_tr.error()());
        }

        return var::Parameters_Transformations("test_model", {"kon", "koff"}, maybe_tr.value(),
                                               std::vector<double>{1.0, 2.0});
    }();
    return transformations;
}

var::Parameters_transformed make_test_parameters() {
    return test_parameter_transformations().standard_parameter_transformed();
}

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

TEST_CASE("generic write_csv emits rows for filtered non-empty matrix probits and skips empty matrices",
          "[write_csv]") {
    using namespace macrodr;
    using namespace macrodr::cmd;

    TempDirGuard dir("write_csv_matrix_shapes");

    {
        const auto base = dir.path / "filtered_matrix_probit";
        std::vector<Information_Distortion_Matrix> samples{
            Information_Distortion_Matrix(SymPosDefMatrix<double>{}),
            Information_Distortion_Matrix(make_diag_spd({2.0, 6.0})),
            Information_Distortion_Matrix(make_diag_spd({4.0, 10.0}))};

        auto summary = var::Vector_Space{
            Probit_statistics<Information_Distortion_Matrix>(
                samples, [](const auto& x) { return x(); }, std::set<double>{0.5})};

        auto out = write_csv(summary, base.string());
        REQUIRE(out);

        const auto lines = read_lines(base.string() + ".csv");
        CHECK(lines.front() == kUnifiedHeader);
        REQUIRE(lines.size() == 9);
    }

    {
        const auto base = dir.path / "empty_matrix";
        auto summary = var::Vector_Space{Information_Distortion_Matrix(SymPosDefMatrix<double>{})};

        auto out = write_csv(summary, base.string());
        REQUIRE(out);

        const auto lines = read_lines(base.string() + ".csv");
        REQUIRE(lines.size() == 1);
        CHECK(lines.front() == kUnifiedHeader);
    }
}

TEST_CASE("generic write_csv emits parameter names for wrapped vectors and matrices", "[write_csv]") {
    using namespace macrodr;
    using namespace macrodr::cmd;

    TempDirGuard dir("write_csv_parameter_indexed");
    const auto base = dir.path / "parameter_indexed";

    auto params = make_test_parameters();
    Matrix<double> score(2, 1);
    score[0] = 1.25;
    score[1] = -0.75;

    auto summary = var::Vector_Space{dlogL(std::move(score), params), FIM(make_diag_spd({2.0, 5.0}), params)};

    auto out = write_csv(summary, base.string());
    REQUIRE(out);

    const auto lines = read_lines(base.string() + ".csv");
    REQUIRE(lines.size() == 7);
    CHECK(lines.front() == kUnifiedHeader);

    const auto dlogl0 = split_csv_row(lines[1]);
    const auto dlogl1 = split_csv_row(lines[2]);
    const auto fim00 = split_csv_row(lines[3]);
    const auto fim01 = split_csv_row(lines[4]);
    const auto fim10 = split_csv_row(lines[5]);
    const auto fim11 = split_csv_row(lines[6]);

    REQUIRE(dlogl0.size() == 22);
    REQUIRE(fim00.size() == 22);

    CHECK(dlogl0[11] == "dlogL");
    CHECK(dlogl0[18] == "0");
    CHECK(dlogl0[19].empty());
    CHECK(dlogl0[20] == "kon");

    CHECK(dlogl1[11] == "dlogL");
    CHECK(dlogl1[18] == "1");
    CHECK(dlogl1[19].empty());
    CHECK(dlogl1[20] == "koff");

    CHECK(fim00[11] == "FIM");
    CHECK(fim00[18] == "0");
    CHECK(fim00[19] == "0");
    CHECK(fim00[20] == "kon_kon");

    CHECK(fim01[20] == "kon_koff");
    CHECK(fim10[20] == "koff_kon");
    CHECK(fim11[20] == "koff_koff");
}
