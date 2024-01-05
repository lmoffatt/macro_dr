#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "general_output_operator.h"
#include "matrix.h"
#include "maybe_error.h"
#include <cstddef>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace var {

template <class Id> class Parameters {
public:
    using value_type = Matrix<double>;
    
private:
    std::string const *ptr_IdName;
    std::vector<std::string> const *ptr_ParNames;
    value_type m_values;
    
public:
    class Names {
    std::vector<std::string> m_names;
    std::map<std::string, std::size_t> m_index_map;

    static auto get_index_map(const std::vector<std::string> &names) {
      std::map<std::string, std::size_t> out;
      for (std::size_t i = 0; i < names.size(); ++i)
        out[names[i]] = i;
      return out;
    }

  public:
    static constexpr bool is_Parameters = true;
      Names() {}
      Names(std::vector<std::string> const &names)
        : m_names{names}, m_index_map{get_index_map(names)} {}
    auto &operator()() const { return m_names; }
    auto &names() const { return m_names; }
    auto &operator[](std::size_t i) const { return m_names[i]; }
    
    Maybe_error<std::size_t> operator[](const std::string &name) const {
      auto it = m_index_map.find(name);
      if (it != m_index_map.end())
        return it->second;
      else
        return error_message("name not found");
    }
    
    friend std::ostream &print(std::ostream &os, Names const &x) {
        ::print(os, x.m_index_map);
        return os;
    }
    };
    
    template <class MatrixType>
        requires std::is_same_v<Matrix<double>, std::decay_t<MatrixType>>
    Parameters(std::string const &IdName,
               std::vector<std::string> const &ParNames, MatrixType &&x)
        : ptr_IdName{&IdName}, ptr_ParNames{&ParNames},
        m_values{std::forward<MatrixType>(x)} {}
    
    template <class VectorType>
        requires std::is_same_v<std::vector<double>, std::decay_t<VectorType>>
    Parameters(std::string const &IdName,
               std::vector<std::string> const &ParNames, VectorType &&x)
        : ptr_IdName{&IdName}, ptr_ParNames{&ParNames},
        m_values{x.size(), 1, std::forward<VectorType>(x)} {}
    
    template <class MatrixType>
        requires std::is_same_v<Matrix<double>, std::decay_t<MatrixType>>
    Parameters create(MatrixType &&x) const {
        return Parameters(IdName(), names(), std::forward<MatrixType>(x));
    }
    
    template <class VectorType>
        requires std::is_same_v<std::vector<double>, std::decay_t<VectorType>>
    Parameters create(VectorType &&x) {
        return Parameters(IdName(), names(), std::forward<VectorType>(x));
    }
    
    Parameters(std::string const &IdName,
               std::vector<std::string> const &ParNames)
        : ptr_IdName{&IdName}, ptr_ParNames{&ParNames},
        m_values{ParNames.size(), 1} {}
    Parameters() {}
    
    auto &operator()() const { return m_values; }
    auto &operator()() { return m_values; }
    
    auto &names() const { return *ptr_ParNames; }
    auto &IdName() const { return *ptr_IdName; }
    
    auto size() const { return m_values.size(); }
    
    auto &operator[](std::size_t i) const { return (*this)()[i]; }
    auto &operator[](std::size_t i) { return (*this)()[i]; }
    
    Matrix<double> operator[](std::pair<std::size_t, std::size_t> ij) const {
        if (m_values.size() == m_values.ncols()) {
            auto out = Matrix<double>(1ul, ij.second - ij.first + 1);
            for (std::size_t i = 0; i < out.size(); ++i)
                out[i] = m_values[ij.first + i];
            return out;
        } else {
            auto out = Matrix<double>(ij.second - ij.first + 1, 1ul);
            for (std::size_t i = 0; i < out.size(); ++i)
                out[i] = m_values[ij.first + i];
            return out;
        }
    }
};

template <class Id>
void report_model(std::string filename, std::string sep,
                  Parameters<Id> const &m) {
    std::ofstream f(filename);
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    auto n = m.size();
    f << "model_name" << sep << "i_par" << sep << "parameter_name" << sep
      << "moment" << sep << "value"
      << "\n";
    for (auto i_par = 0ul; i_par < n; ++i_par)
        f << m.IdName() << sep << i_par << sep << m.names()[i_par] << sep << "mean"
          << sep << m[i_par] << "\n";
}

template <class Id>
void write_Parameters(std::string filename, std::string sep,
                      Parameters<Id> const &m) {
    std::ofstream f(filename);
    f << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    auto n = m.size();
    f << "model_name" << sep << "i_par" << sep << "parameter_name" << sep
      << "moment" << sep << "value"
      << "\n";
    for (auto i_par = 0ul; i_par < n; ++i_par)
        f << m.IdName() << sep << i_par << sep << m.names()[i_par] << sep << "mean"
          << sep << m[i_par] << "\n";
}

template <class Id>
Maybe_error<Parameters<Id>>
load_Parameters(const std::string filename, std::string separator,
                std::string const &ModelName,
                const std::vector<std::string> &ParamNames) {
    std::ifstream f(filename);
    if (!f)
        return error_message("file " + filename +
                             " does not exists or cannot be opened");
    else {
        std::string line;
        std::getline(f, line);
        std::stringstream ss(line);
        
        //    f<<"model_name"<<sep<<"i_par"<<sep<<"parameter_name"<<sep<<"moment"<<sep<<"value"<<"\n";
        
        ss >> ::septr("model_name") >> ::septr(separator) >> ::septr("i_par") >>
            ::septr(separator) >> ::septr("parameter_name") >> ::septr(separator) >>
            ::septr("moment") >> ::septr(separator) >> ::septr("value");
        if (!ss)
            return error_message("file " + filename + " column titles" + line +
                                 " do not correspond ");
        else {
            Parameters<Id> out(ModelName, ParamNames);
            std::size_t i_par=0;
            std::getline(f, line);
            ss = std::stringstream(line);
            std::vector<Maybe_error<bool>> has_par(ParamNames.size(),
                                                   error_message("not provided "));
            while (bool((ss >> ::septr(ModelName) >> ::septr(separator))) &&
                   bool((ss >> i_par >> ::septr(separator))) &&
                   bool((ss >> ::septr(ParamNames[i_par]) >> ::septr(separator)) && bool((ss >> ::septr("mean") >> ::septr(separator) >> out[i_par])))) {
        has_par[i_par] = true;
        std::getline(f, line);
        ss = std::stringstream(line);
            }
            auto allPars = promote_Maybe_error(has_par);
            if (!allPars)
                return allPars.error();
            else
                return std::move(out);
    }
  }
}

} // namespace var

#endif // PARAMETERS_H
