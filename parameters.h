#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "matrix.h"
#include "maybe_error.h"
#include "general_output_operator.h"
#include <cstddef>
#include <fstream>
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
    Names(){}  
    Names(std::vector<std::string> const &names)
        : m_names{names},
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
  
  template<class VectorType>
      requires std::is_same_v<std::vector<double>,std::decay_t<VectorType>>
  Parameters(VectorType&& x):m_values{x.size(),1,std::forward<VectorType>(x)}{}
  
  Parameters(){}
  
  
  auto& operator()()const {return m_values;}
  auto& operator()() {return m_values;}
  
  auto size()const {return m_values.size();}
  
  auto& operator[](std::size_t i)const {return (*this)()[i];}
  auto& operator[](std::size_t i) {return (*this)()[i];}
  
  Matrix<double> operator[](std::pair<std::size_t, std::size_t> ij)const
  {
      if (m_values.size()==m_values.ncols())
      {
          auto out=Matrix<double>(1ul, ij.second-ij.first+1);
          for (std::size_t i=0; i<out.size(); ++i) out[i]=m_values[ij.first+i];
          return out;
      }
      else{
          auto out=Matrix<double>(ij.second-ij.first+1,1ul);
          for (std::size_t i=0; i<out.size(); ++i) out[i]=m_values[ij.first+i];
          return out;
      }
  }
    
};

template<class Id>
void report_model(std::string filename, std::string sep, Parameters<Id> const & m)
{
    std::ofstream f(filename);
    f<<std::setprecision(std::numeric_limits<double>::digits10 + 1);
    auto n=m.size();
    f<<"i_par"<<sep<<"moment"<<sep<<"value"<<"\n";
    for (auto i_par=0ul; i_par<n; ++i_par)
        f<<i_par<<sep<<"mean"<<sep<<m[i_par]<<"\n";
    
}

template<class Id>
Maybe_error<Parameters<Id>> load_Parameters(const std::string filename, std::string separator)
{
    std::ifstream f(filename);
    if (!f)
        return error_message("file "+filename+ " does not exists or cannot be opened");
    else
    {
        f>>::septr("i_par")>>::septr(separator)>>::septr("moment")>>::septr(separator)>>::septr("value");
        if (!f)
            return error_message("file "+filename+ " column titles do not correspond");
        else
        {
        
        std::size_t i_par;
        double value;
        std::vector<double> v;
        while (f>>i_par>>::septr(separator)>>::septr("mean")>>::septr(separator)>>value)
        {
            if (i_par!=v.size())
                return error_message("i_par out of order: i_par="+std::to_string(i_par)+" size="+std::to_string(v.size()));
            else
            {
                v.push_back(value);
                std::string line;
                std::getline(f,line);
            }
        }
        return Parameters<Id>(std::move(v));
    }
        
    }
}


} // namespace var

#endif // PARAMETERS_H
