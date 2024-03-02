#ifndef PARAMETERS_DISTRIBUTION_H
#define PARAMETERS_DISTRIBUTION_H

#include "general_output_operator.h"
#include "matrix.h"
#include "maybe_error.h"
#include "multivariate_normal_distribution.h"
#include "parameters.h"
#include "random_samplers.h"
#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace var {

template <class Id>
class Parameters_Normal_Distribution;

template <class Id>
void write_Prior(const std::string filename, const std::string separator,
                 Parameters_Normal_Distribution<Id> const &d); 

template <class Id>
class Parameters_Normal_Distribution
    : public multivariate_normal_distribution<double,
                                              DiagPosDetMatrix<double>> {
public:
  using base_type =
      multivariate_normal_distribution<double, DiagPosDetMatrix<double>>;

private:
  Parameters_Transformations<Id> m_tr;

public:
  Parameters_Normal_Distribution(Parameters_Transformations<Id> &&t_tr,
                                 base_type &&dist)
      : base_type{std::move(dist)}, m_tr{std::move(t_tr)} {}

  Parameters_Normal_Distribution(Parameters_Transformations<Id> const &t_tr,
                                 base_type &&dist)
      : base_type{std::move(dist)}, m_tr{t_tr} {}

  static Maybe_error<Parameters_Normal_Distribution>
  create(const Parameters_Transformations<Id> &t_par,
         base_type &&dist) {
   
      return Parameters_Normal_Distribution(t_par, 
                                            std::move(dist));
  }

  static Maybe_error<Parameters_Normal_Distribution>
  create(const Parameters_Transformations<Id> &t_par,
         base_type const &dist) {
      return Parameters_Normal_Distribution(t_par, 
                                            dist);
  }

  auto &IdName() const { return m_tr.IdName(); }

  auto &names() const { return m_tr.names(); }

  auto &parameters() const { return m_tr; }
  Parameters_transformed<Id> operator()(mt_64i &mt)const {
    return Parameters_transformed<Id>(parameters(),base_type::operator()(mt));
  }

  Maybe_error<double> logP(const Matrix<double> &x) const {
    return base_type::logP(x);
  }
  Maybe_error<double> logP(const Parameters_transformed<Id> &x) const {
    return base_type::logP(x());
  }

  template <class Parameter>
  friend void report_model(save_Parameter<Parameter> &s,
                           Parameters_Normal_Distribution const &d) {
      
      std::string filename=s.fname + "_prior.csv";
      auto separator=s.sep;
      write_Prior( filename, separator,d);  
      
      }
};

template <class Id>
void write_Prior(const std::string filename, const std::string separator,
                 Parameters_Normal_Distribution<Id> const &d) {
  std::ofstream f(filename);
  using base_type = typename Parameters_Normal_Distribution<Id>::base_type;
  f << std::setprecision(std::numeric_limits<double>::digits10 +3);
  auto tr_m = static_cast<base_type const &>(d).mean();
  auto cov = static_cast<base_type const &>(d).cov();
  
  auto m=d.parameters().inv(tr_m);
  auto n = m.size();
  f << "model_name" << separator << "i_par" << separator << "parameter_name"
    << separator << "parameter_transformation" << separator << "parameter_value"
    << separator << "transformed_mean" << separator << "transformed_variance"
    << "\n";
  
  std::size_t ii_par=0;
  
  for (auto i_par = 0ul; i_par < n; ++i_par) {
    f << d.IdName() << separator << i_par << separator << d.names()[i_par]
      << separator
        
        << d.parameters().transf()[i_par]->to_string() << separator
        << m[i_par];
      if (!d.parameters().transf()[i_par]->is_fixed())
      {
          f << separator << tr_m[ii_par]<< separator << cov(ii_par, ii_par) << "\n";
          ++ii_par;
      } else
      {
          f<< separator<< separator<<"\n";
      }  
  }
}

template <class Id>
Maybe_error<Parameters_Normal_Distribution<Id>>
load_Prior(const std::string filename, const std::string separator,
           const std::string &ModelName,
           const std::vector<std::string> &ParamNames) {
  std::ifstream f(filename);
  if (!f)
    return error_message("file " + filename +
                         " does not exists or cannot be opened");
  else {
    std::string line;
    std::getline(f, line);
    std::stringstream ss(line);

    ss >> ::septr("model_name") >> ::septr(separator) >> ::septr("i_par") >>
        ::septr(separator) >> ::septr("parameter_name") >> ::septr(separator) >>
        ::septr("parameter_transformation") >> ::septr(separator) >>
        ::septr("parameter_value") >> ::septr(separator) >>
        ::septr("transformed_mean") >> ::septr(separator) >>
        ::septr("transformed_variance");

    if (!ss)
      return error_message("file " + filename + " column titles" + line +
                           " do not correspond ");
    else {

      std::size_t i_par;
      string_and_separator transformation_and_sep(separator);
      double value;
      double mean;
      double var;
      std::vector<double> v;
      std::vector<double> means;
      std::vector<double> vars;
      transformations_vector trs;
      std::getline(f, line);
      ss = std::stringstream(line);
      while (ss >> ::septr(ModelName) >> ::septr(separator) >> i_par >>
             ::septr(separator) >> ::septr(ParamNames[i_par]) >>
             ::septr(separator) >> transformation_and_sep >>
             value >> ::septr(separator) >> mean >> ::septr(separator) >> var) {
        if (i_par != v.size())
          return error_message(
              "i_par out of order: i_par=" + std::to_string(i_par) +
              " size=" + std::to_string(v.size()));
        
        auto Maybe_tr = MyTranformations::from_string(transformation_and_sep());
        if (!Maybe_tr)
          return Maybe_tr.error();
        if (Maybe_tr.value()->tr(value)!=mean)
            mean=Maybe_tr.value()->tr(value);
        trs.push_back(std::move(Maybe_tr.value()));
        v.push_back(value);
        
        means.push_back(mean);
        vars.push_back(var);
        std::getline(f, line);
        ss = std::stringstream(line);
      }

      auto Maybe_dist = make_multivariate_normal_distribution(
          Matrix<double>(means.size(), 1, means),
          DiagPosDetMatrix<double>(vars));
      if (!Maybe_dist)
        return Maybe_dist.error();
      else
          return Parameters_Normal_Distribution<Id>(
            Parameters_Transformations<Id>(ModelName, ParamNames, std::move(trs),v),
              std::move(Maybe_dist.value()));
    }
  }
}

template <class Id>
auto prior_around(const Parameters_Transformations<Id> &tr, double error) {
  assert(error > 0);

  auto x = tr.tr(tr.standard_values());
  auto Maybe_dist = make_multivariate_normal_distribution(
      x, DiagPosDetMatrix<double>(x.size(), x.size(), error));

  return Parameters_Normal_Distribution<Id>(tr, std::move(Maybe_dist.value()));
}

template <class Id>
auto prior_around(const Parameters_Transformations<Id> &tr,
                  std::vector<double> values) {
  assert(*std::min(values.begin(), values.end()) > 0);

  auto x = tr.tr(tr.parameter_value());

  auto Maybe_dist = make_multivariate_normal_distribution(
      x, DiagPosDetMatrix<double>(values));

  return Parameters_Normal_Distribution<Id>(tr, std::move(Maybe_dist.value()));
}
} // namespace var

#endif // PARAMETERS_DISTRIBUTION_H
