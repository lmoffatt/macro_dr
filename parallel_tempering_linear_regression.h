#ifndef PARALLEL_TEMPERING_LINEAR_REGRESSION_H
#define PARALLEL_TEMPERING_LINEAR_REGRESSION_H
#include "bayesian_linear_regression.h"
#include "multivariate_normal_distribution.h"
#include "parallel_tempering.h"

class save_Evidence {
    
public:
    std::string sep = ",";
    std::ofstream f;
    std::size_t save_every = 1;
    save_Evidence(std::string const &path, std::size_t interval)
        : f{std::ofstream(path + "__i_iter.csv")}, save_every{interval} {}
    
    template<class Parameters>
    friend void report_title(save_Evidence &s, thermo_mcmc<Parameters> const &) {
        
        s.f << "n_betas" << s.sep << "iter" << s.sep << "beta" << s.sep
            << "meanPrior" << s.sep << "meanLik" << s.sep << "varLik" << s.sep
            << "Evidence_mean" << s.sep << "Evidence_var"
            << "\n";
    }
    
    template<class Parameters>
    friend void report(std::size_t iter, save_Evidence &s,
                       thermo_mcmc<Parameters> const &data) {
        if (iter % s.save_every == 0) {
            
            auto meanLik = mean_logL(data);
            auto meanPrior = mean_logP(data);
            
            auto varLik = var_logL(data, meanLik);
            if (data.beta[0] == 1) {
                auto Evidence2 = calculate_Evidence(data.beta, meanLik, varLik);
                auto Evidence1 = calculate_Evidence(data.beta, meanLik);
                for (std::size_t i_beta = 0; i_beta < num_betas(data); ++i_beta)
                    s.f << num_betas(data) << s.sep << iter << s.sep << data.beta[i_beta]
                        << s.sep << meanPrior[i_beta] << s.sep << meanLik[i_beta] << s.sep
                        << varLik[i_beta] << s.sep << Evidence1 << s.sep << Evidence2
                        << "\n";
            }
        }
    }
    
    template <class Prior, class Likelihood, class Variables, class DataType>
        requires requires (Prior const &prior, Likelihood const& lik, const DataType &y,
                          const Variables &x){{evidence(conjugate{},prior,lik,y, x)};}
    friend void report_model(save_Evidence &s, Prior const &prior, Likelihood const& lik, const DataType &y,
                             const Variables &x, by_beta<double> const &beta0) {
        
        auto expected_Evidence =
            evidence(conjugate{},prior,lik,y, x);
        
        auto meanLik =
            mean_logLik(conjugate{},prior,lik, y, x, beta0);
        if (is_valid(meanLik) && is_valid(expected_Evidence))
            for (std::size_t i_beta = 0; i_beta < size(beta0); ++i_beta)
                s.f << 0 << s.sep << 0 << s.sep << beta0[i_beta] << s.sep << 0 << s.sep
                    << meanLik.value()[i_beta] << s.sep << 0 << s.sep
                    << expected_Evidence << s.sep << expected_Evidence << "\n";
    }
    
    template <class Prior, class Likelihood, class Variables, class DataType>
    friend void report_model(save_Evidence &s, Prior const &prior, Likelihood const& lik, const DataType &y,
                             const Variables &x, by_beta<double> const &beta0) {
        
          for (std::size_t i_beta = 0; i_beta < size(beta0); ++i_beta)
                s.f << 0 << s.sep << 0 << s.sep << beta0[i_beta]  << "\n";
    }
    
    
    
    
};

template <class Parameters,class Algorithm,
         class Reporter>
    requires(   is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> )
class thermo{
    Algorithm alg;
    Reporter rep;
    std::size_t num_scouts_per_ensemble;
    std::size_t thermo_jumps_every;
    double n_points_per_decade;
    double stops_at;
    bool includes_zero;
    std::size_t initseed;
    
    template <class Prior, class Likelihood,class Variables, class DataType>
        requires(is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<Likelihood,Parameters,Variables,DataType>
                 )
    
    auto logEvidence(const Prior prior, const Likelihood &lik, const DataType &y,
                     const Variables &x) {
        
        auto a = alg;
        auto mt = init_mt(initseed);
        auto n_walkers = num_scouts_per_ensemble;
        auto mts = init_mts(mt, num_scouts_per_ensemble / 2);
        auto beta = get_beta_list(n_points_per_decade, stops_at, includes_zero);
        
        auto beta_run = by_beta<double>(beta.rend() - 2, beta.rend());
        auto current = init_thermo_mcmc(n_walkers, beta_run, mts, prior, lik, y, x);
        auto n_par = current.walkers[0][0].parameter.size();
        auto mcmc_run = checks_convergence(std::move(a), current);
        std::size_t iter = 0;
        report_title(rep, current);
        report_model(rep, prior, lik, y, x, beta);
        
        while (beta_run.size() < beta.size() || !mcmc_run.second) {
            while (!mcmc_run.second) {
                step_stretch_thermo_mcmc(iter, current, rep, beta_run, mts, prior, lik, y, x);
                thermo_jump_mcmc(iter, current, rep, beta_run, mt, mts,
                                 thermo_jumps_every);
                report(iter, rep, current);
                mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
            }
            if (beta_run.size() < beta.size()) {
                beta_run.insert(beta_run.begin(), beta[beta_run.size()]);
                current = push_back_new_beta(iter, current, mts, beta_run, prior, lik, y, x);
                std::cerr << "\n  beta_run=" << beta_run[0] << "\n";
                mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
            }
        }
        
        return std::pair(mcmc_run, current);
    }
    
};

template <class Algorithm, class Prior, class Likelihood, class Variables, class DataType,
         class Reporter,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
             is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto thermo_impl(const Algorithm &alg, Prior const &prior, Likelihood const& lik, const DataType &y,
                 const Variables &x, Reporter rep,
                 std::size_t num_scouts_per_ensemble,
                 std::size_t thermo_jumps_every, double n_points_per_decade,
                 double stops_at, bool includes_zero, std::size_t initseed) {
    
    auto a = alg;
    auto mt = init_mt(initseed);
    auto n_walkers = num_scouts_per_ensemble;
    auto mts = init_mts(mt, num_scouts_per_ensemble / 2);
    auto beta = get_beta_list(n_points_per_decade, stops_at, includes_zero);
    
    auto beta_run = by_beta<double>(beta.rend() - 2, beta.rend());
    auto current = init_thermo_mcmc(n_walkers, beta_run, mts, prior, lik, y, x);
    auto n_par = current.walkers[0][0].parameter.size();
    auto mcmc_run = checks_convergence(std::move(a), current);
    std::size_t iter = 0;
    report_title(rep, current);
    report_model(rep, prior, lik, y, x, beta);
    
    while (beta_run.size() < beta.size() || !mcmc_run.second) {
        while (!mcmc_run.second) {
            step_stretch_thermo_mcmc(iter, current, rep, beta_run, mts, prior, lik, y, x);
            thermo_jump_mcmc(iter, current, rep, beta_run, mt, mts,
                             thermo_jumps_every);
            report(iter, rep, current);
            mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        }
        if (beta_run.size() < beta.size()) {
            beta_run.insert(beta_run.begin(), beta[beta_run.size()]);
            current = push_back_new_beta(iter, current, mts, beta_run, prior, lik, y, x);
            std::cerr << "\n  beta_run=" << beta_run[0] << "\n";
            mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        }
    }
    
    return std::pair(mcmc_run, current);
}


template <class Algorithm,
         class Reporter>
//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> )
class thermodynamic_integration{
    Algorithm alg_;
    Reporter rep_;
    std::size_t num_scouts_per_ensemble_;
    std::size_t thermo_jumps_every_;
    double n_points_per_decade_;
    double stops_at_;
    bool includes_zero_;
    std::size_t initseed_;
public:
    
    thermodynamic_integration(Algorithm &&alg,  Reporter&& rep,
                              std::size_t num_scouts_per_ensemble,
                              std::size_t thermo_jumps_every, double n_points_per_decade,
                              double stops_at, bool includes_zero, std::size_t initseed)
        :
        alg_{std::move(alg)},rep_{std::move(rep)},num_scouts_per_ensemble_{num_scouts_per_ensemble},thermo_jumps_every_{thermo_jumps_every},
        n_points_per_decade_{n_points_per_decade},stops_at_{stops_at},includes_zero_{includes_zero},initseed_{initseed}{}
    
    auto& algorithm()const {return alg_;}
    auto& reporter()  {return rep_;}
    auto& num_scouts_per_ensemble()const {return num_scouts_per_ensemble_;}
    auto& thermo_jumps_every()const {return thermo_jumps_every_;}
    auto& n_points_per_decade()const {return n_points_per_decade_;}
    auto& stops_at()const {return stops_at_;}
    auto& includes_zero()const {return includes_zero_;}
    auto& initseed()const {return initseed_;}
    
    
};


template <class Algorithm, class Prior, class Likelihood, class Variables, class DataType,
         class Reporter>
//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto evidence(thermodynamic_integration<Algorithm,Reporter>&& therm, Prior const &prior, Likelihood const& lik, const DataType &y,
              const Variables &x) {
    
    auto a = therm.algorithm();
    auto mt = init_mt(therm.initseed());
    auto n_walkers = therm.num_scouts_per_ensemble();
    auto mts = init_mts(mt, therm.num_scouts_per_ensemble() / 2);
    auto beta = get_beta_list(therm.n_points_per_decade(), therm.stops_at(), therm.includes_zero());
    
    auto beta_run = by_beta<double>(beta.rend() - 2, beta.rend());
    auto current = init_thermo_mcmc(n_walkers, beta_run, mts, prior, lik, y, x);
    //auto n_par = current.walkers[0][0].parameter.size();
    auto mcmc_run = checks_convergence(std::move(a), current);
    std::size_t iter = 0;
    auto& rep=therm.reporter();
    report_title(rep, current);
    report_model(rep, prior, lik, y, x, beta);
    
    while (beta_run.size() < beta.size() || !mcmc_run.second) {
        while (!mcmc_run.second) {
            step_stretch_thermo_mcmc(iter, current, rep, beta_run, mts, prior, lik, y, x);
            thermo_jump_mcmc(iter, current, rep, beta_run, mt, mts,
                             therm.thermo_jumps_every());
            report(iter, rep, current);
            mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        }
        if (beta_run.size() < beta.size()) {
            beta_run.insert(beta_run.begin(), beta[beta_run.size()]);
            current = push_back_new_beta(iter, current, mts, beta_run, prior, lik, y, x);
            std::cerr << "\n  beta_run=" << beta_run[0] << "\n";
            mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        }
    }
    
    return std::pair(mcmc_run, current);
}





class thermo_max {
    std::string path_;
    std::string filename_;
    std::size_t num_scouts_per_ensemble_;
    std::size_t thermo_jumps_every_;
    std::size_t max_iter_;
    double n_points_per_decade_;
    double stops_at_;
    bool includes_zero_;
    std::size_t initseed_;
    
public:
    thermo_max(std::string path, std::string filename,
               std::size_t num_scouts_per_ensemble,
               std::size_t thermo_jumps_every, std::size_t max_iter,
               double n_points_per_decade, double stops_at, bool includes_zero,
               std::size_t initseed):path_{path},filename_{filename},num_scouts_per_ensemble_{num_scouts_per_ensemble},
        thermo_jumps_every_{thermo_jumps_every},max_iter_{max_iter},n_points_per_decade_{n_points_per_decade},
        stops_at_{stops_at},includes_zero_{includes_zero}, initseed_{initseed}{}
    
    template <class Prior, class Likelihood, class Variables, class DataType,
             class Parameters = std::decay_t<decltype(sample(
                 std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
        requires(is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<Likelihood,Parameters,Variables,DataType>)
    auto operator()(const Prior& prior, const Likelihood& lik, const DataType &y, const Variables &x) {
        return thermo_impl(
            less_than_max_iteration(max_iter_), prior, lik, y, x,
            save_mcmc<Parameters,save_likelihood<Parameters>, save_Parameter<Parameters>, save_Evidence>(
                path_, filename_, 10ul, 10ul, 10ul),
            num_scouts_per_ensemble_, thermo_jumps_every_, n_points_per_decade_,
            stops_at_, includes_zero_, initseed_);
    }
};

template <class Prior, class Likelihood, class Variables, class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<Likelihood,Parameters,Variables,DataType>)
auto thermo_max_iter(const Prior& prior, const Likelihood& lik, const DataType &y, const Variables &x,
                     std::string path, std::string filename,
                     std::size_t num_scouts_per_ensemble,
                     std::size_t thermo_jumps_every, std::size_t max_iter,
                     double n_points_per_decade, double stops_at,
                     bool includes_zero, std::size_t initseed) {
    return thermo_impl(less_than_max_iteration(max_iter), prior, lik, y, x,
                       save_mcmc<Parameters,save_likelihood<Parameters>, save_Parameter<Parameters>, save_Evidence>(
                           path, filename, 10ul, 10ul, 10ul),
                       num_scouts_per_ensemble, thermo_jumps_every,
                       n_points_per_decade, stops_at, includes_zero, initseed);
}



template<class Parameters>
auto thermo_by_max_iter(std::string path, std::string filename,
                        std::size_t num_scouts_per_ensemble,
                        std::size_t thermo_jumps_every, std::size_t max_iter,
                        double n_points_per_decade, double stops_at,
                        bool includes_zero, std::size_t initseed) {
    return thermodynamic_integration(less_than_max_iteration(max_iter),
                                     save_mcmc<Parameters,save_likelihood<Parameters>, save_Parameter<Parameters>, save_Evidence>(
                                         path, filename, 10ul, 10ul, 10ul),
                                     num_scouts_per_ensemble, thermo_jumps_every,
                                     n_points_per_decade, stops_at, includes_zero, initseed);
}




template <class Prior, class Likelihood, class Variables, class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<std::mt19937_64 &>(), std::declval<Prior &>()))>>
    requires(is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<Likelihood,Parameters,Variables,DataType>)
auto thermo_convergence(const Prior& prior, const Likelihood& lik, const DataType &y, const Variables &x,
                        std::string path, std::string filename,
                        std::size_t num_scouts_per_ensemble,
                        std::size_t thermo_jumps_every, std::size_t max_iter,
                        double n_points_per_decade, double stops_at,
                        bool includes_zero, std::size_t initseed) {
    return thermo_impl(
        checks_derivative_var_ratio<thermo_mcmc, Parameters>(max_iter * prior.size()), prior,lik,
        y, x,
        save_mcmc<Parameters,save_likelihood<Parameters>, save_Parameter<Parameters>, save_Evidence>(
            path, filename, 10ul, 100ul, 10ul),
        num_scouts_per_ensemble, thermo_jumps_every, n_points_per_decade,
        stops_at, includes_zero, initseed);
}


template<class Parameters>
auto thermo_by_convergence(std::string path, std::string filename,
                           std::size_t num_scouts_per_ensemble,
                           std::size_t thermo_jumps_every, std::size_t max_iter,
                           double n_points_per_decade, double stops_at,
                           bool includes_zero, std::size_t initseed) {
    return thermodynamic_integration(
        checks_derivative_var_ratio<thermo_mcmc, Parameters>(max_iter ),
        save_mcmc<Parameters,save_likelihood<Parameters>, save_Parameter<Parameters>, save_Evidence>(
            path, filename, 10ul, 100ul, 10ul),
        num_scouts_per_ensemble, thermo_jumps_every, n_points_per_decade,
        stops_at, includes_zero, initseed);
}



template <class Cova>
    requires Covariance<double, Cova>
Maybe_error<by_beta<double>> mean_logLik(conjugate,
                                         const multivariate_gamma_normal_distribution<double, Cova> &prior,
                                         const linear_model&,
                                         const Matrix<double> &y,
                                         const Matrix<double> &X,
                                         by_beta<double> const &beta) {
    auto a_0 = prior.alpha(); ;
    auto b_0= prior.beta();
    auto L_0 = prior.Gamma();
    auto SSx = XTX(X);
    auto n = y.nrows();
    by_beta<double> mean_logLik_(size(beta));
    auto beta_0 = prior.mean();
    for (std::size_t i=0; i<size(beta); ++i)
    {
        auto beta0=beta[i];
        auto L_n =   L_0 + beta0 * SSx;
        auto beta_n =
            tr(inv(L_n) * (beta0*(tr(X) * y) + (L_0 * tr(prior.mean()))));
        auto yfit = X * tr(beta_n);
        auto ydiff = y - yfit;
        auto SS = beta0 * xtx(ydiff.value());
        
        auto a_n = a_0 + beta0 * n / 2.0;
        auto b_n = b_0 + 0.5 *  SS + 0.5 * xAxt(beta_0 - beta_n, L_0);
        double d_a_n = 1.0 * n / 2.0;
        auto d_b_n = 0.5 * xtx(ydiff.value());
        auto mean_logLi = -0.5 * n * std::log(2 * std::numbers::pi) -
                          0.5 * Trace(inv(L_n) * SSx) - a_n / b_n * d_b_n +
                          (digamma(a_n) - log(b_n)) * d_a_n;
        if (mean_logLi)
            mean_logLik_[i]=mean_logLi.value();
        else
            return error_message("Error at beta="+std::to_string(beta0)+" "+mean_logLi.error()());
    }
    return mean_logLik_;
}





template <class Parameters,class Cova>
    requires Covariance<double, Cova>
auto bayesian_linear_regression_calculate_mean_logLik(
    const multivariate_normal_distribution<double, Cova> &prior,
    double prior_eps_df, double prior_eps_variance, const Matrix<double> &y,
    const Matrix<double> &X,by_beta<double> const &beta0) {
    
    auto beta = beta0;
    auto L_0 = prior.cov_inv() * prior_eps_variance;
    auto SSx = XTX(X);
    auto n = y.nrows();
    auto a_0 = prior_eps_df / 2.0;
    auto b_0 = prior_eps_df * prior_eps_variance / 2.0;
    auto beta_0 = prior.mean();
    by_beta<std::tuple<Maybe_error<double>,double,Parameters>> mean_logLik;
    mean_logLik.reserve(beta.size());
    
    
    //  by_beta<Maybe_error<double>> mean_Ev;
    //  mean_Ev.reserve(beta.size());
    //  by_beta<Maybe_error<double>> mean_logLik_diff;
    //  mean_logLik_diff.reserve(beta.size());
    
    //  by_beta<Maybe_error<double>> diff_logdetLn;
    //  diff_logdetLn.reserve(beta.size());
    
    //  by_beta<Maybe_error<double>> diff_alogb;
    //  diff_alogb.reserve(beta.size());
    
    //  by_beta<Maybe_error<double>> diff_lgamma;
    //  diff_lgamma.reserve(beta.size());
    
    
    for (std::size_t i = 0; i < beta.size(); ++i) {
        auto L_n =   L_0 + beta[i] * SSx;
        
        auto beta_n =
            tr(inv(L_n) * (beta[i]*(tr(X) * y) + (L_0 * tr(prior.mean()))));
        
        auto yfit = X * tr(beta_n);
        auto ydiff = y - yfit;
        auto SS = beta[i] * xtx(ydiff.value());
        std::cerr<<"SS\n"<<SS<<"\n";
        
        auto a_n = a_0 + beta[i] * n / 2.0;
        auto b_n = b_0 + 0.5 *  SS + 0.5 * xAxt(beta_0 - beta_n, L_0);
        
        
        auto logE_n = -0.5 * beta[i] * n * std::log(2 * std::numbers::pi) +
                      0.5 * (logdet(L_0) - logdet(L_n)) + a_0 * log(b_0) -
                      a_n * log(b_n) + std::lgamma(a_n) - std::lgamma(a_0);
        //    std::cerr<<beta[i]<<"\t"<<"Ev"<<logE_n<<"\n";
        //    std::cerr<<beta[i]<<"\tL_0\t"<<L_0 <<"\n";
        //    std::cerr<<beta[i]<<"\tL_n\t"<<L_n<<"\n";
        //    std::cerr<<beta[i]<<"\tlogdet(L_0) - logdet(L_n)\t"<<logdet(L_0) - logdet(L_n)<<"\n";
        //    std::cerr<<beta[i]<<"\t+ a_0 * log(b_0) -a_n * log(b_n)\t"<<+ a_0 * log(b_0) -a_n * log(b_n)<<"\n";
        //    std::cerr<<beta[i]<<"\t+ std::lgamma(a_n) - std::lgamma(a_0)\t"<<+ std::lgamma(a_n) - std::lgamma(a_0)<<"\n";
        
        
        double d_a_n = 1.0 * n / 2.0;
        auto d_b_n = 0.5 * xtx(ydiff.value());
        
        auto mean_logLi = -0.5 * n * std::log(2 * std::numbers::pi) -
                          0.5 * Trace(inv(L_n) * SSx) - a_n / b_n * d_b_n +
                          (digamma(a_n) - log(b_n)) * d_a_n;
        
        
        
        
        
        mean_logLik.push_back(std::tuple(mean_logLi,std::log(b_n.value()/a_n),beta_n.value()));
        //    mean_Ev.push_back(logE_n);
        //    diff_logdetLn.push_back(logdet(L_n));
        //    diff_alogb.push_back(a_n * log(b_n));
        //    diff_lgamma.push_back(std::lgamma(a_n));
    }
    /*
  for (std::size_t i = 0; i < beta.size(); ++i) {
    beta[i]=std::max(beta[i]*(1+1e-6),beta[i]+1e-9);
    auto L_n = beta[i] * SSx + L_0;

    auto beta_n =
        tr(inv(L_n) * (beta[i] * (tr(X) * y) + (L_0 * tr(prior.mean()))));

    auto yfit = X * tr(beta_n);
    auto ydiff = y - yfit;
    auto SS = beta[i] * xtx(ydiff.value());
    std::cerr<<"SS\n"<<SS<<"\n";

    auto a_n = a_0 + beta[i] * n / 2.0;
    auto b_n = b_0 + 0.5 * SS + 0.5 * xAxt(beta_0 - beta_n, L_0);


    auto logE_n = -0.5 * beta[i] * n * std::log(2 * std::numbers::pi) +
                  0.5 * (logdet(L_0) - logdet(L_n)) +
                  a_0 * log(b_0) - a_n * log(b_n) +
                  std::lgamma(a_n) - std::lgamma(a_0);
    double d_a_n = 1.0 * n / 2.0;
    auto d_b_n = 0.5 * xtx(ydiff.value());
    std::cerr<<beta[i]<<"\t"<<"Ev"<<logE_n<<"\n\n";
    std::cerr<<beta[i]<<"\t\t diff logdet(L_n))\t"<<-(diff_logdetLn[i]-logdet(L_n))/(beta[i]-beta0[i])<<"\n";
    std::cerr<<std::setprecision(12)<<"L_n\n"<<L_n<<"\n";
    std::cerr<<std::setprecision(12)<<"SSx\n"<<SSx<<"\n";
    std::cerr<<beta[i]<<"\t\t inv(L_n) * SSx\t"<<inv(L_n) * SSx<<"\n\n";
    std::cerr<<beta[i]<<"\t\t inv(L_n) \t"<<inv(L_n) <<"\n\n";

    std::cerr<<beta[i]<<"\t\tTrace(inv(L_n) * SSx)\t"<<Trace(inv(L_n) * SSx)<<"\n\n";
    std::cerr<<beta[i]<<"\t\tTrace(inv(static_cast<Matrix<double>const&>(L_n)) * SSx)\t"<<Trace(inv(static_cast<Matrix<double>const&>(L_n)) * SSx)<<"\n\n";

    std::cerr<<beta[i]<<"\t\tdiff_alogb[i]-a_n * log(b_n))\t"<<-(diff_alogb[i]-a_n * log(b_n))/(beta[i]-beta0[i])<<"\n";
    std::cerr<<beta[i]<<"\t\t(a_n / b_n * d_b_n +log(b_n) * d_a_n)\t"<<(a_n / b_n * d_b_n +log(b_n) * d_a_n)<<"\n\n";


    std::cerr<<beta[i]<<"\t\t-(diff_lgamma[i]-std::lgamma(a_n))\t"<<-(diff_lgamma[i]-std::lgamma(a_n))/(beta[i]-beta0[i])<<"\n";
    std::cerr<<beta[i]<<"\t\tdigamma(a_n)  * d_a_n\t"<<+digamma(a_n)  * d_a_n<<"\n\n";



    auto mean_logLi = -0.5 * n * std::log(2 * std::numbers::pi) -
                      0.5 * Trace(inv(L_n) * SSx) -
                      (a_n / b_n * d_b_n +log(b_n) * d_a_n)+
                      digamma(a_n)  * d_a_n;
    mean_logLik_diff.push_back((logE_n-mean_Ev[i])/(beta[i]-beta0[i]));

  }

  */
    
    //  for (std::size_t i=0; i<beta0.size(); ++i)
    //    std::cerr<<"beta= "<<beta0[i]<<"\tmean_logLi"<<mean_logLik[i]<<"\tdiff: "<<mean_logLik_diff[i]<<"\n";
    
    return mean_logLik;
}


#endif // PARALLEL_TEMPERING_LINEAR_REGRESSION_H
