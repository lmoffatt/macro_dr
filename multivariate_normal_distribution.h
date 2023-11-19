#ifndef MULTIVARIATE_NORMAL_DISTRIBUTION_H
#define MULTIVARIATE_NORMAL_DISTRIBUTION_H

#include "matrix.h"
#include <iostream>
#include <random>


template <typename T, typename C>
concept Covariance = std::convertible_to<C, DiagPosDetMatrix<T>> ||
                     std::convertible_to<C, SymPosDefMatrix<T>>;



class normal_distribution {
private:
    std::normal_distribution<double> n_;
public:
    normal_distribution(double mean, double stddev):n_{mean,stddev}{}
    normal_distribution():n_{}{}
    auto mean() const { return n_.mean(); }
    std::size_t size()const {return 1ul;}
    auto cov() const { return n_.stddev()*n_.stddev(); }
    auto stddev() const { return n_.stddev(); }
    
    double operator()(std::mt19937_64 &mt) {
        return n_(mt);
    }
    
    Maybe_error<double> logP(double x)const
    {
        double out=-0.5 * (log(2*std::numbers::pi*cov()) + std::pow(x-mean(),2)/cov());
        if (std::isfinite(out))
            return out;
        else
            return error_message("likelihood not finite:" +std::to_string(out));
    }
};





template <typename T, class Cova>
    requires Covariance<T, Cova>
class multivariate_normal_distribution {
private:
    std::normal_distribution<T> n_;
    Matrix<T> mean_;
    Cova cov_;
    using cholesky_type = std::decay_t<decltype(get_value(::cholesky(cov_)))>;
    cholesky_type cho_;
    Cova cov_inv_;
    double logdetCov_;
    
    static double calc_logdet(const cholesky_type &cho) {
        return std::visit([](auto const &m) { return 2 * logdet(m); }, cho);
    }
    
    multivariate_normal_distribution(Matrix<T> &&mean, Cova &&cov,
                                     cholesky_type &&chol, Cova &&cov_inv, double logdetCov)
        : n_{}, mean_{std::move(mean)}, cov_{std::move(cov)},
        cho_{std::move(chol)}, cov_inv_{std::move(cov_inv)},logdetCov_{logdetCov} {}
    
public:
    template <
        class Mat, class Cov,
        typename Te,
        //= std::decay_t<decltype(get_value(std::declval<Mat>())(0, 0))>,
        class Covae>
    // = std::decay_t<decltype(get_value(std::declval<Cov>()))>>
        requires (Covariance<Te, Covae> && contains_value<Mat &&, Matrix<Te>> &&
                 contains_value<Cov &&, Covae>)
    friend Maybe_error<multivariate_normal_distribution<Te, Covae>>
    make_multivariate_normal_distribution(Mat &&mean, Cov &&cov);
    
    template <
        class Mat, class Cov,
        typename Te ,//= std::decay_t<decltype(get_value(std::declval<Mat>())(0, 0))>,
        class Covae> // = std::decay_t<decltype(get_value(std::declval<Cov>()))>>
        requires (Covariance<Te, Covae> && contains_value<Mat &&, Matrix<Te>> &&
                 contains_value<Cov &&, Covae>)
    friend Maybe_error<multivariate_normal_distribution<Te, Covae>>
    make_multivariate_normal_distribution_from_precision(Mat &&mean, Cov &&cov_inv);
    
    auto &mean() const { return mean_; }
    
    std::size_t size()const {return mean().size();}
    auto &cov() const { return cov_; }
    
    auto &cov_inv() const { return cov_inv_; }
    auto &cholesky() const { return cho_; }
    
    Matrix<double> operator()(std::mt19937_64 &mt) {
        auto z = sample(mt, normal_distribution(0,1),mean().nrows(), mean().ncols());
        if (mean().nrows()==1)
            return z* tr(cholesky())+mean();
        else
            return cholesky() * z+mean();
    }
    auto operator()(std::mt19937_64 &mt, std::size_t n) {
        return  sample(mt,normal_distribution(0.0,1.0), n, mean().ncols()) * tr(cholesky());
    }
    
    auto logDetCov() const {return logdetCov_;}
    double chi2(const Matrix<T>& x)const {
        auto xdiff=x - mean();
        if (xdiff.nrows()==cov_inv().nrows())
            return xtAx(xdiff, cov_inv());
        else
            return xAxt(xdiff, cov_inv());
        
    }
    Maybe_error<double> logP(const Matrix<T>& x)const
    {
        assert(x.size()==mean().size());
        double out=-0.5 * (mean().size() * log(2*std::numbers::pi) + logDetCov() + chi2(x));
        if (std::isfinite(out))
            return out;
        else
            return error_message("likelihood not finite:" +std::to_string(out));
    }
    
    friend std::ostream& operator<<(std::ostream& os, const multivariate_normal_distribution& m)
    {
        
        os<<"mean "<<m.mean();
        os<<"diag(cov) "<<diag(m.cov());
        return os;
        
    }
};


template <
    class Mat, class Cov,
    typename T = std::decay_t<decltype(get_value(std::declval<Mat>())(0ul, 0ul))>,
    class Cova = std::decay_t<decltype(get_value(std::declval<Cov>()))>>
    requires (Covariance<T, Cova> && contains_value<Mat &&, Matrix<T>> &&
             contains_value<Cov &&, Cova>)
Maybe_error<multivariate_normal_distribution<T, Cova>>
make_multivariate_normal_distribution(Mat &&mean, Cov &&cov) {
    return_error<multivariate_normal_distribution<T, Cova>>
        Error{"make_multivariate_normal_distribution"};
    if (!is_valid(mean))
        return Error(get_error(mean)());
    else if (!is_valid(cov))
        return Error(get_error(cov)());
    else  {
        auto beta_cov = get_value(std::forward<Cov>(cov));
        auto chol = cholesky(beta_cov);
        
        if (chol) {
            auto inv = inv_from_chol(chol.value());
            if (inv) {
                auto meanbeta = get_value(std::forward<Mat>(mean));
                auto logDetCov=logdet(chol.value());
                if (logDetCov)
                    return multivariate_normal_distribution<T, Cova>(
                        std::move(meanbeta), std::move(beta_cov), std::move(chol.value()),
                        std::move(inv.value()),logDetCov.value());
                else
                    return Error(logDetCov.error()() + " log determinant fails");
            } else
                return Error(inv.error()() + " covariance cannot be inverted");
        } else
            return Error(chol.error()() +
                         " cholesky fails to build a normal distribution");
    }
}


template <
    class Mat, class Cov,
    typename T = std::decay_t<decltype(get_value(std::declval<Mat>())(0, 0))>,
    class Cova = std::decay_t<decltype(get_value(std::declval<Cov>()))>>
    requires (Covariance<T, Cova> && contains_value<Mat &&, Matrix<T>> &&
             contains_value<Cov &&, Cova>)
Maybe_error<multivariate_normal_distribution<T, Cova>>
make_multivariate_normal_distribution_from_precision(Mat &&mean, Cov &&cov_inv) {
    return_error<multivariate_normal_distribution<T, Cova>> Error{
                                                                  "make_multivariate_normal_distribution"};
    if (!is_valid(mean))
        return Error(get_error(mean)());
    else if (!is_valid(cov_inv))
        return Error(get_error(cov_inv)());
    else  {
        auto chol_inv = cholesky(cov_inv);
        
        if (chol_inv) {
            auto chol = inv(chol_inv.value());
            
            if (chol) {
                auto cov = XXT(tr(chol.value()));
                auto beta_mean = get_value(std::forward<Mat>(mean));
                auto beta_cov_inv=get_value(std::forward<Cov>(cov_inv));
                auto logDetCov=logdet(diag(chol.value()));
                if (logDetCov)
                    return multivariate_normal_distribution<T, Cova>(
                        std::move(beta_mean), std::move(cov), std::move(chol.value()),
                        std::move(beta_cov_inv), logDetCov.value());
                else
                    return Error(logDetCov.error()() + " log of determinant fails");
                
            } else
                return Error(chol.error()() + " cholesky cannot be inverted");
        } else
            return Error(chol_inv.error()() +
                         " cholesky fails to be built from precision");
    }
}




class inverse_gamma_distribution
{
private:
    std::gamma_distribution<> g_;
    double cte_int_;
    
public:
    inverse_gamma_distribution(double _alpha,double _beta):g_{_alpha,_beta},
        cte_int_{-std::lgamma(_alpha)+_alpha*std::log(_beta)}
    {}
    
    double operator()(std::mt19937_64& mt){ return 1.0/g_(mt);}
    
    double alpha()const{return g_.alpha();}
    double beta()const {return g_.beta();}
    Maybe_error<double> logP(double x)const
    {
        auto out= cte_int_-(alpha()+1.0)*std::log(x)-beta()/x;
        if (std::isfinite(out))
            return out;
        else return error_message("probability not finite:" +std::to_string(out));
    }
};

class log_inverse_gamma_distribution
{
private:
    std::gamma_distribution<> g_;
    double cte_int_;
    
public:
    log_inverse_gamma_distribution(double _alpha,double _beta):g_{_alpha,1.0/_beta},
        cte_int_{-std::lgamma(_alpha)+_alpha*std::log(_beta)}
    {}
    
    double operator()(std::mt19937_64& mt){ return -std::log(g_(mt));}
    
    double alpha()const{return g_.alpha();}
    double beta()const {return 1.0/g_.beta();}
    Maybe_error<double> logP(double logx)const
    {
        auto out= cte_int_-alpha()*logx-beta()*std::exp(-logx);
        if (std::isfinite(out))
            return out;
        else return error_message("probability not finite:" +std::to_string(out));
    }
    
    double expected_variance()const { return beta()/alpha();}
};



template <typename T, class Cova>
    requires Covariance<T, Cova>
class multivariate_gamma_normal_distribution: private log_inverse_gamma_distribution, multivariate_normal_distribution<T,Cova> {
    
    
public:
    using m_normal=multivariate_normal_distribution<T,Cova>;
    using m_normal::mean;
    using log_inverse_gamma_distribution::alpha;
    using log_inverse_gamma_distribution::beta;
    
    using m_normal::chi2;
    
    std::size_t size()const { return 1+ m_normal::size();}
    
    
    multivariate_gamma_normal_distribution(log_inverse_gamma_distribution&& g,m_normal&& n):
        log_inverse_gamma_distribution{std::move(g)}, m_normal{std::move(n)}{}
    Matrix<double> operator()(std::mt19937_64 &mt) {
        auto k=m_normal::size();
        auto sample_logVar=log_inverse_gamma_distribution::operator()(mt);
        auto sample_b=m_normal::operator()(mt);
        auto out = Matrix<double>(1, k+1, false);
        out[0]=sample_logVar;
        
        for (std::size_t i = 0; i < k; ++i)
            out[i+1] = sample_b[i];
        return out;
        
    }
    
    Maybe_error<double> logP(const Matrix<T>& x)const
    {
        assert(x.size()==mean().size()+1);
        double logvar = x[0];
        auto logPvar= log_inverse_gamma_distribution::logP(logvar);
        if (!logPvar)
            return error_message("likelihood of variance wrong "+logPvar.error()());
        else{
            double var = std::exp(logvar);
            auto k=mean().size();
            auto b = Matrix<double>(1, k, false);
            for (std::size_t i = 0; i < k; ++i)
                b[i] = x[i + 1];
            auto chi_2=m_normal::chi2(b)/var;
            
            double out=logPvar.value()-0.5 * (k * log(2*std::numbers::pi) + m_normal::logDetCov() + k*logvar+ chi_2);
            if (std::isfinite(out))
                return out;
            else
                return error_message("likelihood not finite:" +std::to_string(out));
        }}
    
    
    
    auto& Gamma()const
    {
        return m_normal::cov_inv();
    }
    
};










#endif // MULTIVARIATE_NORMAL_DISTRIBUTION_H
