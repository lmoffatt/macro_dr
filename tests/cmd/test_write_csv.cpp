#include <catch_amalgamated.hpp>

#include <macrodr/cmd/detail/write_csv_common.h>
#include <macrodr/cmd/likelihood.h>
#include <macrodr/cmd/load_experiment.h>
#include <macrodr/cmd/simulate.h>

#include <algorithm>
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
    std::string cell;
    bool in_quotes = false;
    for (std::size_t i = 0; i < row.size(); ++i) {
        const char ch = row[i];
        if (in_quotes) {
            if (ch == '"') {
                if (i + 1 < row.size() && row[i + 1] == '"') {
                    cell.push_back('"');
                    ++i;
                } else {
                    in_quotes = false;
                }
            } else {
                cell.push_back(ch);
            }
        } else if (ch == ',') {
            columns.push_back(cell);
            cell.clear();
        } else if (ch == '"') {
            in_quotes = true;
        } else {
            cell.push_back(ch);
        }
    }
    columns.push_back(cell);
    return columns;
}

constexpr std::string_view kBaseHeader =
    "scope,simulation_index,sample_index,segment_index,sub_index,n_step,"
    "step_start,step_end,step_middle,agonist,patch_current,"
    "component_path,value_row,value_col,"
    "probit,calculus,statistic,quantile_level,"
    "param_index,param_col,param_name";

std::string expected_header(std::initializer_list<std::string_view> axis_names = {}) {
    std::string out(kBaseHeader);
    for (const auto axis_name : axis_names) {
        out += ",";
        out += axis_name;
    }
    out += ",value";
    return out;
}

std::size_t column_index(const std::vector<std::string>& header, std::string_view name) {
    const auto it = std::find(header.begin(), header.end(), name);
    if (it == header.end()) {
        throw std::runtime_error("missing CSV column: " + std::string(name));
    }
    return static_cast<std::size_t>(std::distance(header.begin(), it));
}

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

Matrix<double> make_col_vector(std::initializer_list<double> vals) {
    Matrix<double> out(vals.size(), 1, false);
    std::size_t i = 0;
    for (auto v : vals) {
        out(i, 0) = v;
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
        Probit_statistics<logL>(-1.5, std::map<double, double>{{0.025, -2.5}, {0.975, -0.5}}, 5)};

    auto out = write_csv(summary, base.string());
    REQUIRE(out);

    const auto lines = read_lines(base.string() + ".csv");
    REQUIRE(lines.size() == 8);
    CHECK(lines.front() == expected_header());

    CHECK(lines[1].find("summary") == 0);
    CHECK(lines[1].find(",point,primitive,count,") != std::string::npos);
    CHECK(lines[2].find(",point,primitive,mean,") != std::string::npos);
    CHECK(lines[3].find(",point,primitive,variance,") != std::string::npos);
    CHECK(lines[4].find(",mean,primitive,bootstrap_count,") != std::string::npos);
    CHECK(lines[5].find(",mean,primitive,value,") != std::string::npos);

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
    CHECK(lines.front() == expected_header());

    const auto header = split_csv_row(lines.front());

    const auto first = split_csv_row(lines[1]);
    const auto second = split_csv_row(lines[2]);
    REQUIRE(first.size() == header.size());
    REQUIRE(second.size() == header.size());

    CHECK(first[column_index(header, "scope")] == "recording");
    CHECK(first[column_index(header, "simulation_index")].empty());
    CHECK(first[column_index(header, "sample_index")] == "0");
    CHECK(first[column_index(header, "segment_index")] == "0");
    CHECK(first[column_index(header, "sub_index")].empty());
    CHECK(first[column_index(header, "n_step")] == "0");
    CHECK(first[column_index(header, "agonist")] == "5");
    CHECK(first[column_index(header, "patch_current")] == "0.1");
    CHECK(first[column_index(header, "component_path")] == "patch_current");
    CHECK(first[column_index(header, "probit")] == "point");
    CHECK(first[column_index(header, "calculus")] == "primitive");
    CHECK(first[column_index(header, "statistic")] == "value");
    CHECK(first[column_index(header, "value")] == "0.1");

    CHECK(second[column_index(header, "scope")] == "recording");
    CHECK(second[column_index(header, "sample_index")] == "1");
    CHECK(second[column_index(header, "segment_index")] == "0");
    CHECK(second[column_index(header, "n_step")] == "1");
    CHECK(second[column_index(header, "agonist")] == "7");
    CHECK(second[column_index(header, "patch_current")] == "0.2");
    CHECK(second[column_index(header, "component_path")] == "patch_current");
    CHECK(second[column_index(header, "probit")] == "point");
    CHECK(second[column_index(header, "calculus")] == "primitive");
    CHECK(second[column_index(header, "statistic")] == "value");
    CHECK(second[column_index(header, "value")] == "0.2");
}

TEST_CASE("generic write_csv serializes one-dimensional indexed values with axis metadata",
          "[write_csv]") {
    using namespace macrodr::cmd;

    TempDirGuard dir("write_csv_indexed");
    const auto base = dir.path / "indexed_summary";

    var::Axis axis{var::AxisId{"models"}, std::vector<std::string>{"scheme_CO", "scheme_CCO"}};
    var::Indexed<std::size_t> indexed(var::IndexSpace{{axis}}, {1, 2});

    auto out = write_csv(indexed, base.string());
    REQUIRE(out);

    const auto lines = read_lines(base.string() + ".csv");
    REQUIRE(lines.size() == 3);
    CHECK(lines.front() == expected_header({"models"}));

    const auto header = split_csv_row(lines.front());

    const auto first = split_csv_row(lines[1]);
    const auto second = split_csv_row(lines[2]);
    REQUIRE(first.size() == header.size());
    REQUIRE(second.size() == header.size());

    CHECK(first[column_index(header, "scope")] == "summary");
    CHECK(first[column_index(header, "component_path")].empty());
    CHECK(first[column_index(header, "models")] == "scheme_CO");
    CHECK(first[column_index(header, "value")] == "1");

    CHECK(second[column_index(header, "models")] == "scheme_CCO");
    CHECK(second[column_index(header, "value")] == "2");
}

TEST_CASE("simulation write_csv merges indexed simulation batches and emits n_substeps",
          "[write_csv]") {
    using namespace macrodr;
    using namespace macrodr::cmd;

    TempDirGuard dir("write_csv_indexed_simulations");
    const auto base = dir.path / "indexed_simulations";

    auto experiment = create_experiment({{1, 2, 5.0}, {1, 3, 7.0}}, 10.0, 0.0, 1.0);
    auto recording_a = define_recording({0.1, 0.2});
    auto recording_b = define_recording({0.3, 0.4});
    auto recording_c = define_recording({0.5, 0.6});

    var::Axis axis{var::AxisId{"n_substeps"}, std::vector<std::string>{"10", "30"}};
    var::Indexed<std::vector<Simulated_recording>> simulations(
        var::IndexSpace{{axis}},
        {{Simulated_recording{{SeedNumber{1}, recording_a}}},
         {Simulated_recording{{SeedNumber{2}, recording_b}},
          Simulated_recording{{SeedNumber{3}, recording_c}}}});

    auto out = write_csv(experiment, simulations, base.string());
    REQUIRE(out);

    const auto lines = read_lines(base.string() + ".csv");
    REQUIRE(lines.size() == 7);
    CHECK(lines.front() == expected_header({"n_substeps"}));

    const auto header = split_csv_row(lines.front());

    const auto first = split_csv_row(lines[1]);
    const auto third = split_csv_row(lines[3]);
    const auto last = split_csv_row(lines[6]);
    REQUIRE(first.size() == header.size());
    REQUIRE(third.size() == header.size());
    REQUIRE(last.size() == header.size());

    CHECK(first[column_index(header, "scope")] == "simulation");
    CHECK(first[column_index(header, "simulation_index")] == "0");
    CHECK(first[column_index(header, "n_substeps")] == "10");
    CHECK(first[column_index(header, "value")] == "0.1");

    CHECK(third[column_index(header, "simulation_index")] == "0");
    CHECK(third[column_index(header, "n_substeps")] == "30");
    CHECK(third[column_index(header, "value")] == "0.3");

    CHECK(last[column_index(header, "simulation_index")] == "1");
    CHECK(last[column_index(header, "n_substeps")] == "30");
    CHECK(last[column_index(header, "value")] == "0.6");
}

TEST_CASE("simulation write_csv supports indexed single simulations", "[write_csv]") {
    using namespace macrodr;
    using namespace macrodr::cmd;

    TempDirGuard dir("write_csv_indexed_single_simulation");
    const auto base = dir.path / "indexed_single_simulation";

    auto experiment = create_experiment({{1, 2, 5.0}, {1, 3, 7.0}}, 10.0, 0.0, 1.0);
    auto recording_a = define_recording({0.1, 0.2});
    auto recording_b = define_recording({0.3, 0.4});

    var::Axis axis{var::AxisId{"n_substeps"}, std::vector<std::string>{"10", "30"}};
    var::Indexed<Simulated_recording> simulations(
        var::IndexSpace{{axis}},
        {Simulated_recording{{SeedNumber{1}, recording_a}},
         Simulated_recording{{SeedNumber{2}, recording_b}}});

    auto out = write_csv(experiment, simulations, base.string());
    REQUIRE(out);

    const auto lines = read_lines(base.string() + ".csv");
    REQUIRE(lines.size() == 5);
    CHECK(lines.front() == expected_header({"n_substeps"}));

    const auto header = split_csv_row(lines.front());

    const auto first = split_csv_row(lines[1]);
    const auto last = split_csv_row(lines[4]);
    REQUIRE(first.size() == header.size());
    REQUIRE(last.size() == header.size());

    CHECK(first[column_index(header, "scope")] == "simulation");
    CHECK(first[column_index(header, "simulation_index")].empty());
    CHECK(first[column_index(header, "n_substeps")] == "10");
    CHECK(first[column_index(header, "value")] == "0.1");

    CHECK(last[column_index(header, "simulation_index")].empty());
    CHECK(last[column_index(header, "n_substeps")] == "30");
    CHECK(last[column_index(header, "value")] == "0.4");
}

TEST_CASE("indexed write_csv helpers merge spaces for broadcast and cartesian product",
          "[write_csv]") {
    using namespace macrodr;
    using namespace macrodr::cmd;

    var::Axis n_substeps_a{var::AxisId{"n_substeps"}, std::vector<std::string>{"10", "30"}};
    var::Axis n_substeps_b{var::AxisId{"n_substeps"}, std::vector<std::string>{"10", "30"}};
    var::Axis models{var::AxisId{"models"}, std::vector<std::string>{"scheme_CO", "scheme_CCO"}};
    var::Axis conflicting_n_substeps{var::AxisId{"n_substeps"},
                                     std::vector<std::string>{"10", "100"}};

    auto same = macrodr::cmd::detail::merge_indexed_write_csv_spaces(
        var::IndexSpace{{n_substeps_a}}, "simulation", var::IndexSpace{{n_substeps_b}},
        "likelihood");
    REQUIRE(same);
    CHECK(same.value().m_axes.size() == 1);
    CHECK(same.value().m_axes[0].m_id.idName == "n_substeps");

    auto cartesian = macrodr::cmd::detail::merge_indexed_write_csv_spaces(
        var::IndexSpace{{n_substeps_a}}, "simulation", var::IndexSpace{{models}}, "likelihood");
    REQUIRE(cartesian);
    CHECK(cartesian.value().m_axes.size() == 2);
    CHECK(cartesian.value().m_axes[0].m_id.idName == "n_substeps");
    CHECK(cartesian.value().m_axes[1].m_id.idName == "models");

    var::Axis algorithms{var::AxisId{"algorithm"}, std::vector<std::string>{"macro_MRV", "macro_IRV"}};
    auto broadcast = macrodr::cmd::detail::merge_indexed_write_csv_spaces(
        var::IndexSpace{{n_substeps_a, algorithms}}, "simulation",
        var::IndexSpace{{n_substeps_b}}, "likelihood");
    REQUIRE(broadcast);
    CHECK(broadcast.value().m_axes.size() == 2);
    CHECK(broadcast.value().m_axes[0].m_id.idName == "n_substeps");
    CHECK(broadcast.value().m_axes[1].m_id.idName == "algorithm");

    auto reordered = macrodr::cmd::detail::merge_indexed_write_csv_spaces(
        var::IndexSpace{{n_substeps_a, algorithms}}, "simulation",
        var::IndexSpace{{algorithms, n_substeps_b}}, "likelihood");
    REQUIRE(reordered);
    CHECK(reordered.value().m_axes.size() == 2);
    CHECK(reordered.value().m_axes[0].m_id.idName == "n_substeps");
    CHECK(reordered.value().m_axes[1].m_id.idName == "algorithm");

    auto conflict = macrodr::cmd::detail::merge_indexed_write_csv_spaces(
        var::IndexSpace{{n_substeps_a}}, "simulation",
        var::IndexSpace{{conflicting_n_substeps}}, "likelihood");
    REQUIRE_FALSE(conflict);
    CHECK(conflict.error()().find("could not be merged") != std::string::npos);

    auto multi_axis = macrodr::cmd::detail::validate_indexed_write_csv_space(
        var::IndexSpace{{n_substeps_a, models}}, "simulation");
    REQUIRE(multi_axis);

    var::Axis datasets{var::AxisId{"dataset"}, std::vector<std::string>{"A", "B"}};
    auto three_axis = macrodr::cmd::detail::validate_indexed_write_csv_space(
        var::IndexSpace{{n_substeps_a, models, datasets}}, "simulation");
    REQUIRE(three_axis);
}

TEST_CASE("generic write_csv serializes two-dimensional indexed values with one column per axis",
          "[write_csv]") {
    using namespace macrodr::cmd;

    TempDirGuard dir("write_csv_indexed_2d");
    const auto base = dir.path / "indexed_summary_2d";

    var::Axis algorithm{var::AxisId{"algorithm"},
                        std::vector<std::string>{"macro_MRV", "macro_IRV"}};
    var::Axis n_substeps{var::AxisId{"n_substeps"}, std::vector<std::string>{"10", "30"}};
    var::Indexed<std::size_t> indexed(var::IndexSpace{{algorithm, n_substeps}}, {1, 2, 3, 4});

    auto out = write_csv(indexed, base.string());
    REQUIRE(out);

    const auto lines = read_lines(base.string() + ".csv");
    REQUIRE(lines.size() == 5);
    CHECK(lines.front() == expected_header({"algorithm", "n_substeps"}));

    const auto header = split_csv_row(lines.front());

    const auto first = split_csv_row(lines[1]);
    const auto second = split_csv_row(lines[2]);
    const auto third = split_csv_row(lines[3]);
    const auto fourth = split_csv_row(lines[4]);
    REQUIRE(first.size() == header.size());
    REQUIRE(second.size() == header.size());
    REQUIRE(third.size() == header.size());
    REQUIRE(fourth.size() == header.size());

    CHECK(first[column_index(header, "algorithm")] == "macro_MRV");
    CHECK(first[column_index(header, "n_substeps")] == "10");
    CHECK(first[column_index(header, "value")] == "1");

    CHECK(second[column_index(header, "algorithm")] == "macro_IRV");
    CHECK(second[column_index(header, "n_substeps")] == "10");
    CHECK(second[column_index(header, "value")] == "2");

    CHECK(third[column_index(header, "algorithm")] == "macro_MRV");
    CHECK(third[column_index(header, "n_substeps")] == "30");
    CHECK(third[column_index(header, "value")] == "3");

    CHECK(fourth[column_index(header, "algorithm")] == "macro_IRV");
    CHECK(fourth[column_index(header, "n_substeps")] == "30");
    CHECK(fourth[column_index(header, "value")] == "4");
}

TEST_CASE("generic write_csv serializes three-dimensional indexed values with one column per axis",
          "[write_csv]") {
    using namespace macrodr::cmd;

    TempDirGuard dir("write_csv_indexed_3d");
    const auto base = dir.path / "indexed_summary_3d";

    var::Axis algorithm{var::AxisId{"algorithm"},
                        std::vector<std::string>{"macro_MRV", "macro_IRV"}};
    var::Axis n_substeps{var::AxisId{"n_substeps"}, std::vector<std::string>{"10", "30"}};
    var::Axis dataset{var::AxisId{"dataset"}, std::vector<std::string>{"A", "B"}};
    var::Indexed<std::size_t> indexed(
        var::IndexSpace{{algorithm, n_substeps, dataset}}, {1, 2, 3, 4, 5, 6, 7, 8});

    auto out = write_csv(indexed, base.string());
    REQUIRE(out);

    const auto lines = read_lines(base.string() + ".csv");
    REQUIRE(lines.size() == 9);
    CHECK(lines.front() == expected_header({"algorithm", "n_substeps", "dataset"}));

    const auto header = split_csv_row(lines.front());

    const auto first = split_csv_row(lines[1]);
    const auto second = split_csv_row(lines[2]);
    const auto last = split_csv_row(lines[8]);
    REQUIRE(first.size() == header.size());
    REQUIRE(second.size() == header.size());
    REQUIRE(last.size() == header.size());

    CHECK(first[column_index(header, "algorithm")] == "macro_MRV");
    CHECK(first[column_index(header, "n_substeps")] == "10");
    CHECK(first[column_index(header, "dataset")] == "A");
    CHECK(first[column_index(header, "value")] == "1");

    CHECK(second[column_index(header, "algorithm")] == "macro_IRV");
    CHECK(second[column_index(header, "n_substeps")] == "10");
    CHECK(second[column_index(header, "dataset")] == "A");
    CHECK(second[column_index(header, "value")] == "2");

    CHECK(last[column_index(header, "algorithm")] == "macro_IRV");
    CHECK(last[column_index(header, "n_substeps")] == "30");
    CHECK(last[column_index(header, "dataset")] == "B");
    CHECK(last[column_index(header, "value")] == "8");
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
        CHECK(lines.front() == expected_header());
        REQUIRE(lines.size() == 9);
    }

    {
        const auto base = dir.path / "empty_matrix";
        auto summary = var::Vector_Space{Information_Distortion_Matrix(SymPosDefMatrix<double>{})};

        auto out = write_csv(summary, base.string());
        REQUIRE(out);

        const auto lines = read_lines(base.string() + ".csv");
        REQUIRE(lines.size() == 1);
        CHECK(lines.front() == expected_header());
    }
}

TEST_CASE("generic write_csv emits parameter-indexed rows for DIB and skips empty vectors",
          "[write_csv]") {
    using namespace macrodr;
    using namespace macrodr::cmd;

    TempDirGuard dir("write_csv_dib_shapes");
    auto params = make_test_parameters();

    {
        const auto base = dir.path / "dib_probit";
        std::vector<Distortion_Induced_Bias> samples{
            Distortion_Induced_Bias(Matrix<double>{}, params),
            Distortion_Induced_Bias(make_col_vector({1.0, -2.0}), params),
            Distortion_Induced_Bias(make_col_vector({3.0, -4.0}), params)};

        auto summary = var::Vector_Space{
            Probit_statistics<Distortion_Induced_Bias>(
                samples, [](const auto& x) { return x(); }, std::set<double>{0.5})};

        auto out = write_csv(summary, base.string());
        REQUIRE(out);

        const auto lines = read_lines(base.string() + ".csv");
        CHECK(lines.front() == expected_header());
        REQUIRE(lines.size() == 5);

        const auto header = split_csv_row(lines.front());
        const auto first = split_csv_row(lines[1]);
        const auto second = split_csv_row(lines[2]);
        REQUIRE(first.size() == header.size());
        REQUIRE(second.size() == header.size());
        CHECK(first[column_index(header, "component_path")] ==
              "Probit_statistics_Distortion_Induced_Bias");
        CHECK(first[column_index(header, "param_name")] == "kon");
        CHECK(second[column_index(header, "param_name")] == "koff");
    }

    {
        const auto base = dir.path / "empty_dib";
        auto summary = var::Vector_Space{Distortion_Induced_Bias(Matrix<double>{}, params)};

        auto out = write_csv(summary, base.string());
        REQUIRE(out);

        const auto lines = read_lines(base.string() + ".csv");
        REQUIRE(lines.size() == 1);
        CHECK(lines.front() == expected_header());
    }
}

TEST_CASE("generic write_csv keeps Gaussian Fisher variance matrix-shaped", "[write_csv]") {
    using namespace macrodr;
    using namespace macrodr::cmd;

    TempDirGuard dir("write_csv_gfi_variance");
    const auto base = dir.path / "gfi_variance";
    auto params = make_test_parameters();

    std::vector<Gaussian_Fisher_Information> samples{
        Gaussian_Fisher_Information(make_diag_spd({2.0, 5.0}), params),
        Gaussian_Fisher_Information(make_diag_spd({4.0, 9.0}), params)};

    auto summary =
        var::Vector_Space{Moment_statistics<Gaussian_Fisher_Information, false>(samples)};

    auto out = write_csv(summary, base.string());
    REQUIRE(out);

    const auto lines = read_lines(base.string() + ".csv");
    CHECK(lines.front() == expected_header());
    REQUIRE(lines.size() == 10);

    const auto header = split_csv_row(lines.front());
    std::vector<std::vector<std::string>> variance_rows;
    for (std::size_t i = 1; i < lines.size(); ++i) {
        auto row = split_csv_row(lines[i]);
        REQUIRE(row.size() == header.size());
        if (row[column_index(header, "component_path")] ==
                "Moment_statistics_Gaussian_Fisher_Information_false" &&
            row[column_index(header, "statistic")] == "variance") {
            variance_rows.push_back(std::move(row));
        }
    }

    REQUIRE(variance_rows.size() == 4);
    CHECK(variance_rows[0][column_index(header, "param_index")] == "0");
    CHECK(variance_rows[0][column_index(header, "param_col")] == "0");
    CHECK(variance_rows[0][column_index(header, "param_name")] == "kon_kon");
    CHECK(variance_rows[1][column_index(header, "param_index")] == "0");
    CHECK(variance_rows[1][column_index(header, "param_col")] == "1");
    CHECK(variance_rows[1][column_index(header, "param_name")] == "kon_koff");
    CHECK(variance_rows[2][column_index(header, "param_index")] == "1");
    CHECK(variance_rows[2][column_index(header, "param_col")] == "0");
    CHECK(variance_rows[2][column_index(header, "param_name")] == "koff_kon");
    CHECK(variance_rows[3][column_index(header, "param_index")] == "1");
    CHECK(variance_rows[3][column_index(header, "param_col")] == "1");
    CHECK(variance_rows[3][column_index(header, "param_name")] == "koff_koff");
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
    CHECK(lines.front() == expected_header());

    const auto header = split_csv_row(lines.front());

    const auto dlogl0 = split_csv_row(lines[1]);
    const auto dlogl1 = split_csv_row(lines[2]);
    const auto fim00 = split_csv_row(lines[3]);
    const auto fim01 = split_csv_row(lines[4]);
    const auto fim10 = split_csv_row(lines[5]);
    const auto fim11 = split_csv_row(lines[6]);

    REQUIRE(dlogl0.size() == header.size());
    REQUIRE(fim00.size() == header.size());

    CHECK(dlogl0[column_index(header, "component_path")] == "dlogL");
    CHECK(dlogl0[column_index(header, "param_index")] == "0");
    CHECK(dlogl0[column_index(header, "param_col")].empty());
    CHECK(dlogl0[column_index(header, "param_name")] == "kon");

    CHECK(dlogl1[column_index(header, "component_path")] == "dlogL");
    CHECK(dlogl1[column_index(header, "param_index")] == "1");
    CHECK(dlogl1[column_index(header, "param_col")].empty());
    CHECK(dlogl1[column_index(header, "param_name")] == "koff");

    CHECK(fim00[column_index(header, "component_path")] == "FIM");
    CHECK(fim00[column_index(header, "param_index")] == "0");
    CHECK(fim00[column_index(header, "param_col")] == "0");
    CHECK(fim00[column_index(header, "param_name")] == "kon_kon");

    CHECK(fim01[column_index(header, "param_name")] == "kon_koff");
    CHECK(fim10[column_index(header, "param_name")] == "koff_kon");
    CHECK(fim11[column_index(header, "param_name")] == "koff_koff");
}
