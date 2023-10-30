#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "matrix.h"
#include "maybe_error.h"
#include "general_output_operator.h"
#include <map>
#include <string>
#include <vector>

namespace var {

template <class Id> class Parameters {
public:
    using value_type=Matrix<double>;
private:
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
    template <class stringvector>
      requires std::is_same_v<std::decay_t<stringvector>,
                              std::vector<std::string>>
    Names(stringvector &&names)
        : m_names{std::forward<stringvector>(names)},
          m_index_map{get_index_map(names)} {}
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
    
    friend std::ostream& print(std::ostream& os, Names const & x)
    {
        ::print(os,x.m_index_map);
        return os;
    }
  };
  
  template<class MatrixType>
      requires std::is_same_v<Matrix<double>,std::decay_t<MatrixType>>
  Parameters(MatrixType&& x): m_values{std::forward<MatrixType>(x)}{}
  
  Parameters(){}
  
  
  auto& operator()()const {return m_values;}
  auto& operator()() {return m_values;}
  
  auto size()const {return m_values.size();}
  
  auto& operator[](std::size_t i)const {return (*this)()[i];}
  auto& operator[](std::size_t i) {return (*this)()[i];}
  
    
};



} // namespace var

#endif // PARAMETERS_H
