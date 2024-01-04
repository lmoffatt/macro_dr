#ifndef PARAMETERS_DISTRIBUTION_H
#define PARAMETERS_DISTRIBUTION_H

#include "matrix.h"
#include "parameters.h"
#include "multivariate_normal_distribution.h"
#include "random_samplers.h"
#include <cstddef>
#include <sstream>
#include <vector>
#include <fstream>

namespace var {

template<class Id>
class Parameters_Normal_Distribution: public multivariate_normal_distribution<double,DiagPosDetMatrix<double>>
{
private:
    Parameters<Id> m_par;
    
public:

    
    using base_type=multivariate_normal_distribution<double,DiagPosDetMatrix<double>>;
    
    auto& IdName()const {return m_par.IdName();}
    
    auto& names()const {return m_par.names();}
    
    Parameters_Normal_Distribution(const Parameters<Id>& t_par, base_type&& dist):
        base_type{std::move(dist)}, m_par{t_par}{}
    
    Parameters_Normal_Distribution(const Parameters<Id>& t_par, base_type const& dist):
        base_type{dist}, m_par{t_par}{}    
    
    Parameters<Id> operator()(mt_64i &mt) {
        return m_par.create(base_type::operator ()(mt));
    }
    
    Maybe_error<double> logP(const Matrix<double>& x)const
    {
        return base_type::logP(x); 
    }
    Maybe_error<double> logP(const Parameters<Id>& x)const
    {
        return base_type::logP(x()); 
    }
    
    
    template<class Parameter>
    friend void report_model(save_Parameter<Parameter>& s, Parameters_Normal_Distribution const & d)
    {
        std::ofstream f(s.fname+"_prior.csv");
        f<<std::setprecision(std::numeric_limits<double>::digits10 + 1);
        auto m=static_cast<base_type const &>(d).mean();
        auto cov=static_cast<base_type const &>(d).cov();
        auto n=m.size();
        f << "model_name" << s.sep << "i_par" << s.sep << "parameter_name" << s.sep
          << "moment" << s.sep << "value"
          << "\n";
        
        for (auto i_par=0ul; i_par<n; ++i_par)
            f << d.IdName() << s.sep << i_par << s.sep << d.names()[i_par] << s.sep << "mean"
              << s.sep << m[i_par] << "\n";
        for (auto i_par=0ul; i_par<n; ++i_par)
            f << d.IdName() << s.sep << i_par << s.sep << d.names()[i_par] << s.sep << "covar"
              <<s.sep<<cov(i_par, i_par)<<"\n";
        
    }
    
    
};

template<class Id>
Maybe_error<Parameters_Normal_Distribution<Id>> load_Prior(const std::string filename, const std::string separator,const std::string& ModelName, const std::vector<std::string> & ParamNames)
{
    std::ifstream f(filename);
    if (!f)
        return error_message("file "+filename+ " does not exists or cannot be opened");
    else
    {
        std::string line;
        std::getline(f,line);
        std::stringstream ss(line);
        
        ss >> ::septr("model_name") >> ::septr(separator) >> ::septr("i_par") >>
            ::septr(separator) >> ::septr("parameter_name") >> ::septr(separator) >>
            ::septr("moment") >> ::septr(separator) >> ::septr("value");
        
        if (!ss)
            return error_message("file " + filename + " column titles" + line +
                                 " do not correspond ");
        else
        {
            
            std::size_t i_par;
            double value;
            std::vector<double> v;
            std::getline(f,line);
            ss= std::stringstream(line);
            while (ss >> ::septr(ModelName) >> ::septr(separator) >> i_par >>
                   ::septr(separator) >> ::septr(ParamNames[i_par]) >>
                   ::septr(separator) >> ::septr("mean") >> ::septr(separator) >>
                   value)
//                while (ss>>i_par>>::septr(separator)>>::septr("mean")>>::septr(separator)>>value)
            {
                if (i_par!=v.size())
                    return error_message("i_par out of order: i_par="+std::to_string(i_par)+" size="+std::to_string(v.size()));
                else
                {
                    v.push_back(value);
                }
                std::getline(f,line);
                ss= std::stringstream(line);            
            }
            ss= std::stringstream(line);
            std::vector<double> cv;
            
            while (ss >> ::septr(ModelName) >> ::septr(separator) >> i_par >>
                   ::septr(separator) >> ::septr(ParamNames[i_par]) >>
                   ::septr(separator) >> ::septr("covar") >> ::septr(separator) >>
                   value)            
         //   while (ss>>i_par>>::septr(separator)>>::septr("covar")>>::septr(separator)>>value)
            {
                if (i_par!=cv.size())
                    return error_message("i_par out of order: i_par="+std::to_string(i_par)+" size="+std::to_string(v.size()));
                else
                {
                    cv.push_back(value);
                    std::getline(f,line);
                    ss= std::stringstream(line);            
                }
            }
            
            
            auto Maybe_dist= make_multivariate_normal_distribution(Matrix<double>(v.size(),1,v),DiagPosDetMatrix<double>(cv));
            if (!Maybe_dist)
                return Maybe_dist.error();
            else
                return Parameters_Normal_Distribution<Id>{Parameters<Id>(ModelName,ParamNames),std::move(Maybe_dist.value())};
        }
        
    }
}



template<class Id>
auto prior_around(const Parameters<Id>& x, double error )
{
    assert(error>0);
    auto Maybe_dist=make_multivariate_normal_distribution(x(),DiagPosDetMatrix<double>(x().size(),x().size(),error));
    
    return Parameters_Normal_Distribution<Id>{x,Maybe_dist.value()};
}

template<class Id>
auto prior_around(const Parameters<Id>& x, std::vector<double> values )
{
    assert(*std::min(values.begin(), values.end())>0);
    auto Maybe_dist=make_multivariate_normal_distribution(x(),DiagPosDetMatrix<double>(values));
    
    return Parameters_Normal_Distribution<Id>(x,Maybe_dist.value());
}



}


#endif // PARAMETERS_DISTRIBUTION_H
