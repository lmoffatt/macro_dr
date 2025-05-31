#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <cstddef>
#include <fstream>
#include <map>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "derivative_operator.h"
#include "general_output_operator.h"
#include "matrix.h"
#include "maybe_error.h"
#include "variables.h"

namespace var {
class base_transformation {
   public:
    virtual std::string to_string() const = 0;

    virtual Maybe_unique<base_transformation> from_string(const std::string&) const = 0;
    virtual std::unique_ptr<base_transformation> clone() const = 0;

    virtual double tr(double x) const = 0;

    virtual Derivative<double, double> tr(Derivative<double, double> x) const = 0;

    virtual std::string tr(std::string name) const = 0;

    virtual double inv(double tr_x) const = 0;

    virtual Derivative<double, double> inv(Derivative<double, double> tr_x) const = 0;

    virtual Derivative<double, Matrix<double>> inv(
        Derivative<double, Matrix<double>> tr_x) const = 0;

    virtual bool is_fixed() const = 0;
};

template <class Transformation>
class template_transformation : public base_transformation {
   public:
    template_transformation() = default;
    virtual std::string to_string() const override {
        return Transformation{}.name();
    }
    virtual Maybe_unique<base_transformation> from_string(const std::string& m) const override {
        if (m == to_string())
            return new template_transformation();
        else
            return error_message("expected: " + to_string() + " found: " + m);
    }
    virtual std::unique_ptr<base_transformation> clone() const override {
        return std::unique_ptr<base_transformation>(new template_transformation());
    }
    virtual double tr(double x) const override {
        return Transformation{}.tr(x);
    }
    virtual Derivative<double, double> tr(Derivative<double, double> x) const override {
        return Transformation{}.tr(x);
    }
    virtual std::string tr(std::string x) const override {
        return Transformation{}.tr(x);
    }
    virtual double inv(double tr_x) const override {
        return Transformation{}.inv(tr_x);
    }
    virtual Derivative<double, double> inv(Derivative<double, double> tr_x) const override {
        return Transformation{}.inv(tr_x);
    }

    virtual Derivative<double, Matrix<double>> inv(
        Derivative<double, Matrix<double>> tr_x) const override {
        return Transformation{}.inv(tr_x);
    }

    virtual bool is_fixed() const override {
        return Transformation{}.isfixed();
    }
};

inline auto clone(const std::vector<std::unique_ptr<base_transformation>>& p) {
    std::vector<std::unique_ptr<base_transformation>> out;
    out.reserve(p.size());
    for (auto& e : p) {
        out.push_back(e->clone());
    }
    return out;
}

using transformations_vector_t = std::vector<std::unique_ptr<base_transformation>>;

class transformations_vector : public transformations_vector_t {
    std::vector<std::size_t> m_fixed;
    static auto get_fixed(transformations_vector_t const& tr) {
        std::vector<std::size_t> out;
        for (std::size_t i = 0; i < tr.size(); ++i) {
            if (tr[i]->is_fixed())
                out.push_back(i);
        }
        return out;
    }

   public:
    using transformations_vector_t::operator[];
    using transformations_vector_t::size;
    using base_type = transformations_vector_t;
    transformations_vector(base_type&& other)
        : base_type(std::move(other)), m_fixed{get_fixed(*this)} {
    }
    transformations_vector(base_type const& other)
        : base_type(clone(other)), m_fixed{get_fixed(*this)} {
    }

    transformations_vector(const transformations_vector& other)
        : base_type(clone(static_cast<base_type const&>(other))), m_fixed{get_fixed(*this)} {
    }
    transformations_vector() = default;
    transformations_vector(transformations_vector&&) = default;
    transformations_vector& operator=(transformations_vector&&) = default;
    transformations_vector& operator=(transformations_vector const&) = default;

    auto& fixed_set() const {
        return m_fixed;
    }

    auto get_transformed_names(const std::vector<std::string>& names) const {
        assert(names.size() == size());
        std::vector<std::string> out;
        for (std::size_t i = 0; i < size(); ++i) {
            auto name = (*this)[i]->tr(names[i]);
            if (!name.empty())
                out.push_back(name);
        }
        return out;
    }

    auto remove_fixed(std::vector<double> const& x) const {
        std::vector<double> out(size() - fixed_set().size());
        std::size_t out_i = 0;
        std::size_t fixed_i = 0;
        for (std::size_t i = 0; i < size(); ++i) {
            if ((fixed_i < fixed_set().size()) && (i == fixed_set()[fixed_i])) {
                ++fixed_i;
            } else {
                out[out_i] = x[i];
                ++out_i;
            }
        }
        return out;
    }
};

struct Identity_Tr {
    static constexpr auto name() {
        return "Linear";
    }
    double tr(double x) const {
        return x;
    }
    Derivative<double, double> tr(Derivative<double, double> x) const {
        return x;
    }
    std::string tr(std::string x) const {
        return x;
    }
    double inv(double x) const {
        return x;
    }
    Derivative<double, double> inv(Derivative<double, double> x) const {
        return x;
    }
    Derivative<double, Matrix<double>> inv(Derivative<double, Matrix<double>> x) const {
        return x;
    }

    bool isfixed() {
        return false;
    }
};

struct Log10_Tr {
    static constexpr auto name() {
        return "Log10";
    }
    double tr(double x) const {
        return std::log10(x);
    }
    Derivative<double, double> tr(Derivative<double, double> x) const {
        return log(x);
    }

    std::string tr(std::string x) const {
        return "log10_" + x;
    }
    double inv(double x) const {
        return std::pow(10.0, x);
    }
    Derivative<double, double> inv(Derivative<double, double> x) const {
        return pow(10.0, x);
    }
    Derivative<double, Matrix<double>> inv(Derivative<double, Matrix<double>> x) const {
        return pow(10.0, x);
    }
    bool isfixed() const {
        return false;
    }
};

struct Fixed_Tr {
    static constexpr auto name() {
        return "Fixed";
    }
    double tr(double x) const {
        return x;
    }
    double inv(double x) const {
        return x;
    }
    Derivative<double, double> tr(Derivative<double, double> x) const {
        return x;
    }
    std::string tr(std::string) const {
        return "";
    }
    Derivative<double, double> inv(Derivative<double, double> x) const {
        return x;
    }
    Derivative<double, Matrix<double>> inv(Derivative<double, Matrix<double>> x) const {
        return x;
    }
    bool isfixed() const {
        return true;
    }
};

template <class... Transformations>
struct Transformation_library_impl {
    static auto from_string(const std::string& s) {
        return (template_transformation<Transformations>{}.from_string(s) || ...);
    }

    static Maybe_error<transformations_vector> from_strings(const std::vector<std::string>& s) {
        std::vector<std::unique_ptr<base_transformation>> out(s.size());
        bool result = true;
        std::string error;
        for (std::size_t i = 0; i < s.size(); ++i) {
            auto Maybe_tr = from_string(s[i]);
            if (Maybe_tr)
                out[i] = std::move(Maybe_tr.value());
            else {
                result = false;
                error += std::to_string(i) + "th paramter: " + Maybe_tr.error()() + "\n";
            }
        }
        if (result)
            return out;
        else {
            return error_message(error);
        }
    }
};

using MyTranformations_list = Transformation_library_impl<Identity_Tr, Log10_Tr, Fixed_Tr>;

struct MyTranformations : public MyTranformations_list {
    using MyTranformations_list::from_string;
    using MyTranformations_list::from_strings;
};

template <template <class...> class V>
bool verify(const transformations_vector& p, const V<double>& value, const V<double>& tr_value) {
    if ((p.size() != value.size()) || (value.size() != tr_value.size()))
        return false;
    for (std::size_t i = 0; i < p.size(); ++i) {
        if (p[i]->inv(tr_value[i]) != value[i])
            return false;
    }
    return true;
}

class Parameters_Transformations;

class Parameters_values {
   public:
    using value_type = Matrix<double>;

   private:
    value_type m_values;
    Parameters_Transformations const* m_par;

   public:
    auto& parameters() const {
        return *m_par;
    }
    const std::vector<std::string>& names() const;

    template <class VectorType>
        requires std::is_same_v<std::vector<double>, std::decay_t<VectorType>>
    Parameters_values(const Parameters_Transformations& par, VectorType&& x)
        : m_values{x.size(), 1, std::forward<VectorType>(x)}, m_par{&par} {
    }

    template <class MatrixType>
        requires std::is_same_v<Matrix<double>, std::decay_t<MatrixType>>
    Parameters_values(const Parameters_Transformations& par, MatrixType&& x)
        : m_values{std::forward<MatrixType>(x)}, m_par{&par} {
    }

    template <class MatrixType>
        requires std::is_same_v<Matrix<double>, std::decay_t<MatrixType>>
    Parameters_values create(MatrixType&& x) const {
        return Parameters_values(parameters(), std::forward<MatrixType>(x));
    }

    template <class VectorType>
        requires std::is_same_v<std::vector<double>, std::decay_t<VectorType>>
    Parameters_values create(VectorType&& x) {
        return Parameters_values(parameters(), std::forward<VectorType>(x));
    }

    // Parameters(std::string const &IdName,
    //            std::vector<std::string> const &ParNames)
    //     : ptr_IdName{&IdName}, ptr_ParNames{&ParNames},
    //     m_values{ParNames.size(), 1} {}
    Parameters_values() {
    }

    auto& operator()() const {
        return m_values;
    }
    auto& operator()() {
        return m_values;
    }

    auto size() const {
        return m_values.size();
    }

    auto& operator[](std::size_t i) const {
        return (*this)()[i];
    }
    auto& operator[](std::size_t i) {
        return (*this)()[i];
    }

    Matrix<double> operator[](std::pair<std::size_t, std::size_t> ij) const {
        if (m_values.size() == m_values.ncols()) {
            auto out = Matrix<double>(1ul, ij.second - ij.first + 1);
            for (std::size_t i = 0; i < out.size(); ++i) out[i] = m_values[ij.first + i];
            return out;
        } else {
            auto out = Matrix<double>(ij.second - ij.first + 1, 1ul);
            for (std::size_t i = 0; i < out.size(); ++i) out[i] = m_values[ij.first + i];
            return out;
        }
    }
};

Maybe_error<bool> compare_contents(
    const Parameters_values& s0, const Parameters_values& s1,
    double RelError = std::sqrt(std::numeric_limits<double>::epsilon()) * 100,
    double AbsError = std::sqrt(std::numeric_limits<double>::epsilon()) * 100,
    std::size_t max_errors = 10) {
    if (s0.names() != s1.names())
        return error_message("\n not the same parameters names!!\n");
    if (s0() == s1())
        return true;
    std::size_t n_errors = 0;
    std::string message;
    if (s0.size() != s1.size()) {
        message +=
            "differ in size: " + std::to_string(s0.size()) + " vs " + std::to_string(s1.size());
        ++n_errors;
    }
    std::size_t i = 0;
    auto n = std::min(s0.size(), s1.size());
    while ((n_errors < max_errors) && (i < n)) {
        auto Maybe_equal = compare_contents(s0[i], s1[i], RelError, AbsError);
        if (!Maybe_equal) {
            message += "\n " + std::to_string(i) + "th parameter " + s0.names()[i] + ": " +
                       Maybe_equal.error()() + "\n";
            ++n_errors;
        }
        ++i;
    }
    if (n_errors > 0) {
        return error_message(message);

    } else
        return true;
}

class Parameters_Names {
    std::map<std::string, std::size_t> m_index_map;
    std::vector<std::string> m_ParNames;
    static auto get_index_map(const std::vector<std::string>& names) {
        std::map<std::string, std::size_t> out;
        for (std::size_t i = 0; i < names.size(); ++i) out[names[i]] = i;
        return out;
    }

   public:
    Parameters_Names(std::vector<std::string>&& x)
        : m_index_map{get_index_map(x)}, m_ParNames{std::move(x)} {
    }
    Parameters_Names(std::vector<std::string> const& x)
        : m_index_map{get_index_map(x)}, m_ParNames{x} {
    }
    auto& operator()() const {
        return m_ParNames;
    }
    Maybe_error<std::size_t> operator[](const std::string& name) const {
        auto it = m_index_map.find(name);
        if (it != m_index_map.end())
            return it->second;
        else
            return error_message("name not found");
    }
};

void write_Parameters(std::ostream& f, std::string sep, Parameters_Transformations const& m);

class Parameters_Transformations {
   public:
    using value_type = Matrix<double>;
    static auto get_index_map(const std::vector<std::string>& names) {
        std::map<std::string, std::size_t> out;
        for (std::size_t i = 0; i < names.size(); ++i) out[names[i]] = i;
        return out;
    }
    auto& transf() const {
        return m_tr;
    }
    auto& standard_values() const {
        return m_standard_values;
    }

   private:
    std::string m_IdName;
    std::vector<std::string> m_ParNames;

    std::map<std::string, std::size_t> m_index_map;

    transformations_vector m_tr;
    value_type m_standard_values;
    std::vector<std::string> m_TrParNames;

    auto free_to_all(const Matrix<double>& x) const {
        Matrix<double> out;
        auto& m_fixed = transf().fixed_set();
        if (x.nrows() > x.ncols())
            out = Matrix<double>(x.nrows() + m_fixed.size(), 1);
        else
            out = Matrix<double>(1, x.ncols() + m_fixed.size());
        assert(out.size() == names().size());
        std::size_t i_in = 0;
        std::size_t i_fi = 0;

        for (std::size_t i_out = 0; i_out < out.size(); ++i_out) {
            if ((m_fixed.size() > i_fi) && (i_out == m_fixed[i_fi])) {
                out[i_out] = (*this).standard_values()[i_out];
                ++i_fi;
            } else {
                out[i_out] = (*this).transf()[i_out]->inv(x[i_in]);
                ++i_in;
            }
        }
        return out;
    }

    auto all_to_free(const Matrix<double>& x) const {
        auto& m_fixed = transf().fixed_set();
        Matrix<double> out;
        if (x.nrows() > x.ncols())
            out = Matrix<double>(x.nrows() - m_fixed.size(), 1);
        else
            out = Matrix<double>(1, x.ncols() - m_fixed.size());
        std::size_t i_out = 0;
        std::size_t i_fi = 0;

        for (std::size_t i_in = 0; i_in < x.size(); ++i_in) {
            if ((m_fixed.size() > i_fi) && (i_in == m_fixed[i_fi])) {
                ++i_fi;
            } else {
                out[i_out] = (*this).transf()[i_in]->tr(x[i_in]);
                ++i_out;
            }
        }
        return out;
    }

   public:
    auto standard_parameter() const {
        return Parameters_values(*this, m_standard_values);
    }

    template <class C_Matrix>
        requires U<C_Matrix, Matrix<double>>
    auto inv(const C_Matrix& x) const {
        return free_to_all(x);
    }

    auto tr(const Matrix<double>& x) const {
        return all_to_free(x);
    }

    void push_back(const std::string& name, std::unique_ptr<base_transformation>&& transf,
                   double value) {
        m_ParNames.push_back(name);
        m_tr.push_back(std::move(transf));
        value_type val;
        if (m_standard_values.ncols() > m_standard_values.nrows())
            val = value_type(m_standard_values.nrows(), m_standard_values.ncols() + 1);
        else
            val = value_type(m_standard_values.nrows() + 1, m_standard_values.ncols());
        for (std::size_t i = 0; i < m_standard_values.size(); ++i) val[i] = m_standard_values[i];
        val[m_standard_values.size()] = value;
        m_standard_values = std::move(val);
    }

    class Names {
        std::vector<std::string> m_names;
        std::map<std::string, std::size_t> m_index_map;

       public:
        static constexpr bool is_Parameters = true;
        Names() {
        }
        Names(std::vector<std::string> const& names)
            : m_names{names}, m_index_map{get_index_map(names)} {
        }
        auto& operator()() const {
            return m_names;
        }
        auto& names() const {
            return m_names;
        }
        auto& operator[](std::size_t i) const {
            return m_names[i];
        }

        Maybe_error<std::size_t> operator[](const std::string& name) const {
            auto it = m_index_map.find(name);
            if (it != m_index_map.end())
                return it->second;
            else
                return error_message("name not found");
        }

        friend std::ostream& print(std::ostream& os, Names const& x) {
            ::print(os, x.m_index_map);
            return os;
        }
    };

    template <class MatrixType>
        requires std::is_same_v<Matrix<double>, std::decay_t<MatrixType>>
    Parameters_Transformations(std::string const& IdName, std::vector<std::string> const& ParNames,
                               transformations_vector const& tr, MatrixType&& std_values)
        : m_IdName{IdName},
          m_ParNames{ParNames},
          m_index_map{get_index_map(ParNames)},
          m_tr{tr},
          m_standard_values{std::forward<MatrixType>(std_values)},
          m_TrParNames{tr.get_transformed_names(ParNames)} {
    }

    template <class MatrixType>
        requires std::is_same_v<std::vector<double>, std::decay_t<MatrixType>>
    Parameters_Transformations(std::string const& IdName, std::vector<std::string> const& ParNames,
                               transformations_vector const& tr, MatrixType&& std_values)
        : m_IdName{IdName},
          m_ParNames{ParNames},
          m_tr{tr},
          m_index_map{get_index_map(ParNames)},
          m_standard_values{
              Matrix<double>(std_values.size(), 1, std::forward<MatrixType>(std_values))},
          m_TrParNames{tr.get_transformed_names(ParNames)} {
    }

    // Parameters(std::string const &IdName,
    //            std::vector<std::string> const &ParNames)
    //     : ptr_IdName{&IdName}, ptr_ParNames{&ParNames},
    //     m_values{ParNames.size(), 1} {}
    Parameters_Transformations() {
    }

    std::vector<std::string> const& names() const {
        return m_ParNames;
    }
    auto& transformed_names() const {
        return m_TrParNames;
    }

    auto& IdName() const {
        return m_IdName;
    }
    Maybe_error<std::size_t> operator[](const std::string& name) const {
        auto it = m_index_map.find(name);
        if (it != m_index_map.end())
            return it->second;
        else
            return error_message("name not found");
    }

    auto size() const {
        return m_standard_values.size();
    }

    friend std::ostream& operator<<(std::ostream& os, Parameters_Transformations const& p) {
        write_Parameters(os, "\t", p);
        return os;
    }
};

class Parameters_transformed {
   public:
    using value_type = Matrix<double>;

   private:
    value_type m_values;
    Parameters_Transformations const* m_par;

   public:
    Maybe_error<bool> is_good() const {
        if (m_values.ncols() * m_values.nrows() != m_values.size())
            return error_message("wrong dimensions: ncols=" + std::to_string(m_values.ncols()) +
                                 " nrows=" + std::to_string(m_values.nrows()) +
                                 " size=" + std::to_string(m_values.size()));
        if (m_par == nullptr)
            return error_message("lacks Parameters_Tranformations");
        return true;
    }

    auto& parameters() const {
        return *m_par;
    }
    template <class VectorType>
        requires std::is_same_v<std::vector<double>, std::decay_t<VectorType>>
    Parameters_transformed(const Parameters_Transformations& par, VectorType&& x)
        : m_par{&par}, m_values{x.size(), 1, std::forward<VectorType>(x)} {
    }

    template <class MatrixType>
        requires std::is_same_v<Matrix<double>, std::decay_t<MatrixType>>
    Parameters_transformed(const Parameters_Transformations& par, MatrixType&& x)
        : m_values{std::forward<MatrixType>(x)}, m_par{&par} {
    }

    template <class MatrixType>
        requires std::is_same_v<Matrix<double>, std::decay_t<MatrixType>>
    Parameters_transformed create(MatrixType&& x) const {
        return Parameters_transformed(parameters(), std::forward<MatrixType>(x));
    }

    template <class VectorType>
        requires std::is_same_v<std::vector<double>, std::decay_t<VectorType>>
    Parameters_transformed create(VectorType&& x) {
        return Parameters_transformed(parameters(), std::forward<VectorType>(x));
    }

    // Parameters(std::string const &IdName,
    //            std::vector<std::string> const &ParNames)
    //     : ptr_IdName{&IdName}, ptr_ParNames{&ParNames},
    //     m_values{ParNames.size(), 1} {}
    Parameters_transformed() {
    }

    auto& operator()() const {
        return m_values;
    }
    auto& operator()() {
        return m_values;
    }

    auto size() const {
        return m_values.size();
    }

    auto& operator[](std::size_t i) const {
        return (*this)()[i];
    }
    auto& operator[](std::size_t i) {
        return (*this)()[i];
    }

    auto to_value() const {
        return Parameters_values(parameters(), parameters().inv((*this)()));
    }
};

// template <class Id>
// void report_model(std::string filename, std::string sep,
//                   Parameters<Id> const &m) {
//   std::ofstream f(filename);
//   f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
//   auto n = m.size();
//   f << "model_name" << sep << "i_par" << sep << "parameter_name" << sep
//     << "moment" << sep << "value"
//     << "\n";
//   for (auto i_par = 0ul; i_par < n; ++i_par)
//     f << m.IdName() << sep << i_par << sep << m.names()[i_par] << sep <<
//     "mean"
//       << sep << m[i_par] << "\n";
// }
void write_Parameters(std::ostream& f, std::string sep, Parameters_Transformations const& m) {
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    auto n = m.size();

    f << "model_name" << sep << "i_par" << sep << "parameter_name" << sep
      << "parameter_transformation" << sep << "parameter_value" << sep << "transformed_mean"
      << "\n";

    for (auto i_par = 0ul; i_par < n; ++i_par)
        f << m.IdName() << sep << i_par << sep << m.names()[i_par] << sep
          << m.transf()[i_par]->to_string() << sep << m.standard_values()[i_par] << sep
          << m.transf()[i_par]->tr(m.standard_values()[i_par]) << "\n";
}

void write_Parameters(std::string filename, std::string sep, Parameters_Transformations const& m) {
    std::ofstream f(filename);
    write_Parameters(f, sep, m);
}

Maybe_error<Parameters_Transformations> load_Parameters(
    const std::string filename, const std::string separator, const std::string& ModelName,
    const std::vector<std::string>& ParamNames) {
    std::ifstream f(filename);
    if (!f)
        return error_message("file " + filename + " does not exists or cannot be opened");
    else {
        std::string line;
        std::getline(f, line);
        std::stringstream ss(line);

        ss >> ::septr("model_name") >> ::septr(separator) >> ::septr("i_par") >>
            ::septr(separator) >> ::septr("parameter_name") >> ::septr(separator) >>
            ::septr("parameter_transformation") >> ::septr(separator) >>
            ::septr("parameter_value") >> ::septr(separator) >> ::septr("transformed_mean");

        if (!ss)
            return error_message("file " + filename + " column titles" + line +
                                 " do not correspond ");
        else {
            std::size_t i_par;
            string_and_separator transformation_and_sep(separator);
            double value;
            double mean;
            std::vector<double> v;
            std::vector<double> means;
            transformations_vector trs;
            std::getline(f, line);
            ss = std::stringstream(line);
            while (ss >> ::septr(ModelName) >> ::septr(separator) >> i_par >> ::septr(separator) >>
                   ::septr(ParamNames[i_par]) >> ::septr(separator) >> transformation_and_sep >>
                   value >> ::septr(separator) >> mean) {
                if (i_par != v.size())
                    return error_message("i_par out of order: i_par=" + std::to_string(i_par) +
                                         " size=" + std::to_string(v.size()));

                auto Maybe_tr = MyTranformations::from_string(transformation_and_sep());
                if (!Maybe_tr)
                    return Maybe_tr.error();
                if (Maybe_tr.value()->tr(value) != mean)
                    mean = Maybe_tr.value()->tr(value);
                trs.push_back(std::move(Maybe_tr.value()));
                v.push_back(value);

                means.push_back(mean);
                std::getline(f, line);
                ss = std::stringstream(line);
            }

            return Parameters_Transformations(ModelName, ParamNames, std::move(trs), v);
        }
    }
}

inline std::vector<std::string> const& Parameters_values::names() const {
    return parameters().names();
}

}  // namespace var

#endif  // PARAMETERS_H
