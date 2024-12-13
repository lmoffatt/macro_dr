#ifndef PARALLEL_TEMPERING_FRACTION_H
#define PARALLEL_TEMPERING_FRACTION_H


#include "derivative_operator.h"
#include "distributions.h"
#include "maybe_error.h"
#include "parallel_tempering.h"
#include "parallel_tempering_linear_regression.h"
#include "general_algorithm_on_containers.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>


template<class Parameters>
struct fraction_mcmc {
    Parameters parameter;
    double logP;
    std::map<std::size_t,logLs> logL_map;
    logLs unknown_logL=nan_logL();
    
    auto& logL(std::size_t i_frac)const
    {
        auto it=logL_map.find(i_frac);
        if (it==logL_map.end())
            return unknown_logL;
        else
            return it->second;            
    }
    
};

auto i_fract_beta_to_global_beta(std::vector<std::size_t>const& i_fracs, std::vector<double>const & beta,std::vector<std::size_t> num_samples, std::size_t max_samples)
{
    
    std::vector<double> global_beta(i_fracs.size());
    
    for (std::size_t i=1; i<global_beta.size(); ++i)
    {
        auto n0= i_fracs[i]>0?num_samples[i_fracs[i]-1]:0;
        global_beta[i]=(n0+(num_samples[i_fracs[i]]-n0)*beta[i])/max_samples;
    }
    return global_beta;
}

auto global_beta_to_S(const std::vector<double>& beta)
{
    auto iB=beta.back()==1?1:0;
    auto S=std::vector<double>(beta.size()-iB);
    auto n= S.size()-1;
    S[n]=std::log(1/beta[n]-1);
    for (std::size_t i=1; i<S.size(); ++i)
        S[n-i]=std::log(1/beta[n-i]-1/beta[n-i+1]);
    return S;
}


auto S_t0_global_beta(const std::vector<double>& S)
{
    auto beta=std::vector<double>(S.size());
    auto n= beta.size()-1;
    beta[n]= 1/(1+std::exp(S[n]));
    for (std::size_t i=1; i<S.size(); ++i)
        beta[n-i]=1.0/(1.0/beta[n-i+1]+std::exp(S[n-i]));
    return beta;
}




template <class Parameters>
    requires std::is_assignable_v<Parameters, Parameters const &>
struct thermo_fraction_mcmc //:public thermo_mcmc<Parameters>
{
    std::vector<std::size_t> samples_size;
    std::size_t max_samples; 
    by_beta<std::size_t> i_fraction;
    by_beta<std::vector<std::size_t>> i_fractions;
    by_beta<double> beta;
    by_beta<double> S;
    by_beta<ensemble<fraction_mcmc<Parameters>>> walkers;
    by_beta<ensemble<std::map<std::size_t,logL_statistics>>> walkers_sta;
    
    by_beta<ensemble<std::size_t>> i_walkers;
    by_beta<emcee_Step_statistics> emcee_stat;
    by_beta<Thermo_Jump_statistics> thermo_stat;
    
    auto num_fractions()const {return i_fraction.size();}
    void reset_statistics() {
        for (std::size_t i_beta = 0; i_beta < walkers_sta.size(); ++i_beta) {
             for (auto &ee : walkers_sta[i_beta]) {
                 for (auto &eee:ee)
                 {
                     eee.second().reset();
                 }
             }
         }
        for (auto &e : emcee_stat)
            e().reset();
        for (auto &e : thermo_stat)
            e().reset();
    }
    
    auto global_beta(std::size_t i_global_beta)const 
    {
        auto i_frac=i_fraction[i_global_beta];
        auto i_frac0=i_frac>0?i_frac-1: 0;
        
        return (samples_size[i_frac0]+(samples_size[i_frac]-samples_size[i_frac0])*beta[i_global_beta])/max_samples;
    }
    
    auto global_beta_to_i_frac_beta(double global_beta)const
    {
        auto ns=global_beta*max_samples;
        auto i_frac=0;
        while(i_frac<i_fraction.size()&&samples_size[i_frac]*(1+std::numeric_limits<double>::epsilon())<ns)
        {
            ++i_frac;
        }
        auto sample_size0=i_frac>0?samples_size[i_frac-1]:0;
        auto beta=(ns-sample_size0)/(samples_size[i_frac]-sample_size0);
        
        return std::tuple(i_frac,beta);
    }
    
    auto get_Walkers_number() const { return walkers[0].size(); }
    auto &get_Beta() const { return beta; }
    auto &get_Parameter(std::size_t iw, std::size_t i_b) const {
        return walkers[i_b][iw].parameter;
    }
    
    auto get_Walker(std::size_t iw, std::size_t i_b) const {
        return i_walkers[i_b][iw];
    }
    
    auto num_samples() const {
        std::size_t min_samples = get<count>(walkers_sta[0][0].begin()->second()())();
        for (std::size_t i = 0; i < walkers_sta.size(); ++i)
            for (std::size_t j = 0; j < walkers_sta[i].size(); ++j)
                if (get<count>(walkers_sta[0][0].begin()->second()())() < min_samples)
                    min_samples = get<count>(walkers_sta[0][0].begin()->second()())();
        return min_samples;
    }
    
    auto last_walker() const {
        std::size_t last = 0;
        for (auto &e : i_walkers)
            for (auto iw : e)
                if (iw > last)
                    last = iw;
        return last;
    }
};

template <class Parameters>
std::size_t num_walkers(thermo_fraction_mcmc<Parameters> const &x) {
    return x.walkers[0].size();
}

template <class Parameters>
std::size_t num_Parameters(thermo_fraction_mcmc<Parameters> const &x) {
    return x.walkers[0][0].parameter.size();
}


template <class Parameters>
auto calculate_initial_logL0(thermo_fraction_mcmc<Parameters> const &initial_data) {
    auto logL_0 = std::map<std::size_t,logL_statistics>{};
    for (std::size_t i = 0; i < initial_data.walkers_sta.size(); ++i) {
        for (std::size_t j = 0; j < initial_data.walkers_sta[i].size(); ++j) {
            for (auto& e: initial_data.walkers_sta[i][j]){
                logL_0[e.first]() &= e.second();
            }
        }
    }
    
    
    
    
    
    
    return logL_0;
}

inline std::map<std::size_t,double> calculate_logL_mean(ensemble<std::map<std::size_t,logL_statistics>> const &sta) {
    std::map<std::size_t,double> sum;
    for (std::size_t i = 1; i < sta.size(); ++i)
        for (auto &e: sta[i])
        {
          sum[e.first] += get<mean<logL>>(e.second()())()();
        }
    for(auto &e: sum)
         e.second= e.second/ sta.size();
    
    return sum;
}
inline by_beta<std::map<std::size_t,double>>
calculate_logL_mean(by_beta<ensemble<std::map<std::size_t,logL_statistics>>> const &sta) {
    by_beta<std::map<std::size_t,double>> out(sta.size());
    for (std::size_t i = 0; i < sta.size(); ++i)
        out[i] = calculate_logL_mean(sta[i]);
    return out;
}


inline logL_statistics
calculate_across_sta(ensemble<std::map<std::size_t,logL_statistics>>  &sta, std::size_t i_frac) {
    logL_statistics across= get<mean<logL>>(sta[0][i_frac]()())();
    for (std::size_t i = 1; i < sta.size(); ++i)
           across() &= get<mean<logL>>(sta[i][i_frac]()())();
    return across;
}
inline auto
calculate_across_sta(by_beta<ensemble<std::map<std::size_t,logL_statistics>>>  &sta, by_beta<std::size_t> i_frac) {
    by_beta<logL_statistics> across(sta.size());
    for (std::size_t i = 0; i < sta.size(); ++i)
        across[i] = calculate_across_sta(sta[i], i_frac[i]);
    return across;
}

inline auto calculate_sample_size(by_beta<ensemble<std::map<std::size_t,logL_statistics>>>  &sta,
                                  std::size_t i_beta, std::size_t i_frac) {
    return get<count>(sta[i_beta][0][i_frac]()());
}


inline variance<logL>
calculate_within_sta(ensemble<std::map<std::size_t,logL_statistics>>  &sta, std::size_t i_frac) {
    variance<logL> within(logL(0.0));
    count  df(count(0));
    for (std::size_t i = 0; i < sta.size(); ++i) {
        auto r_df = get<count>(sta[i][i_frac]()())() - 1;
        within()() += r_df * get<variance<logL>>(sta[i][i_frac]()())()();
        df() += r_df;
    }
   within()() /= df();
    return within;
}

inline auto
calculate_within_sta(by_beta<ensemble<std::map<std::size_t,logL_statistics>>>  &sta, by_beta<std::size_t> i_frac) {
    by_beta<variance<logL>> within(sta.size());
    for (std::size_t i = 0; i < sta.size(); ++i)
        within[i] = calculate_within_sta(sta[i],i_frac[i]);
    return within;
}




template <class Parameters>
auto calculate_deltaBeta_deltaL(const thermo_fraction_mcmc<Parameters> &current) {
    
    std::vector<double> dBdL(current.walkers.size() - 1);
    auto L = calculate_logL_mean(current.walkers_sta);
    
    for (std::size_t i = 0; i < dBdL.size(); ++i){
        auto i_frac1=current.i_fraction[i+1];
        auto i_frac=current.i_fraction[i];
        if (i_frac1==i_frac){
           dBdL[i] = (L[i + 1][i_frac1] - L[i][i_frac1]) * (current.beta[i + 1] - current.beta[i]);
        }else{
            dBdL[i] = (L[i + 1][i_frac1] - L[i][i_frac1]) * (current.beta[i + 1]);
            
        }
    }
    return dBdL;
}


template <class Parameters>
auto calculate_delta_Evidence_variance( thermo_fraction_mcmc<Parameters> &current,
                                       const by_beta<double> &beta) {
    
    
    std::vector<double> deltaEvidence_variance(beta.size() - 1);
    auto n = beta.size() - 1;
    auto across = calculate_across_sta(current.walkers_sta[0], current.i_fraction[0]);
    
    auto var1 = get<variance<logL>>(across()())()();
    auto beta1=0;
    
    for (std::size_t i = 1; i <beta.size(); ++i) {
        auto var0=var1;
        auto beta0=beta1;
        across = calculate_across_sta(current.walkers_sta[i], current.i_fraction[i]);
        
        auto var1 = get<variance<logL>>(across()())()();
        auto beta1=beta[i];
        
        deltaEvidence_variance[i-1] =
            sqr((beta1 - beta0) / 2.0) * (var0 + var1);
        if ((beta1==1)&&i+1<beta.size())
        {
            across = calculate_across_sta(current.walkers_sta[i], current.i_fraction[i+1]);
             var1 = get<variance<logL>>(across()())()();
             beta1=0;
        }
    }
    return deltaEvidence_variance;
}



template <class Parameters>
auto calculate_Acceptance_variance(
    const std::vector<double> &deltaEvidence_variance,
    const thermo_fraction_mcmc<Parameters> &current) {
    
    std::vector<double> acceptance_variance(deltaEvidence_variance.size());
    for (std::size_t i = 0; i < deltaEvidence_variance.size(); ++i)
        acceptance_variance[i] = current.thermo_stat[i]().rate().valid()
                                     ? (current.thermo_stat[i]().rate().value() /
                                        deltaEvidence_variance[i])
                                     : 0.0;
    return acceptance_variance;
}

template <class Parameters>
auto mean_logL(by_iteration<thermo_fraction_mcmc<Parameters>>  &series) {
    auto out = by_beta<std::map<std::size_t,double>>(num_betas(series[0]));
    auto n_walkers = num_walkers(series[0]);
    auto n_iters = num_samples(series);
    for (std::size_t i = 0; i < num_samples(series); ++i)
        for (std::size_t ibeta = 0; ibeta < num_betas(series[0]); ++ibeta)
            for (std::size_t iwalker = 0; iwalker < num_walkers(series[0]); ++iwalker)
                for (auto& e: series[i].walkers[ibeta][iwalker].logL)
                out[ibeta][e.first] +=
                    e.second / n_iters / n_walkers;
    return out;
}


template <class Parameters>
auto mean_logL(thermo_fraction_mcmc<Parameters> const &mcmc) {
    auto out = by_beta<std::map<std::size_t,logLs>>(num_betas(mcmc));
    auto n_walkers = num_walkers(mcmc);
    for (std::size_t iwalker = 0; iwalker < num_walkers(mcmc); ++iwalker)
        for (std::size_t ibeta = 0; ibeta < num_betas(mcmc); ++ibeta)
            for (auto const& e: mcmc.walkers[ibeta][iwalker].logL_map)
            out[ibeta][e.first] = out[ibeta][e.first] + e.second / n_walkers;
    return out;
}

template <class Parameters>
auto mean_logP(thermo_fraction_mcmc<Parameters> const &mcmc) {
    auto out = by_beta<double>(num_betas(mcmc), 0);
    auto n_walkers = num_walkers(mcmc);
    for (std::size_t iwalker = 0; iwalker < num_walkers(mcmc); ++iwalker)
        for (std::size_t ibeta = 0; ibeta < num_betas(mcmc); ++ibeta)
            out[ibeta] += mcmc.walkers[ibeta][iwalker].logP / n_walkers;
    return out;
}


template <class Parameters>
auto var_logL(thermo_fraction_mcmc<Parameters>  &mcmc, by_beta<std::map<std::size_t,logLs>>  &mean) {
    auto out = by_beta<std::map<std::size_t,logLs>>(num_betas(mcmc));
    auto n_walkers = num_walkers(mcmc);
    for (std::size_t iwalker = 0; iwalker < num_walkers(mcmc); ++iwalker)
        for (std::size_t ibeta = 0; ibeta < num_betas(mcmc); ++ibeta)
            for (auto const& e: mcmc.walkers[ibeta][iwalker].logL_map)
            out[ibeta][e.first] =
                out[ibeta][e.first] +
                pow(e.second - mean[ibeta][e.first], 2) / n_walkers;
    return out;
}

auto beta_to_delta_s(std::vector<double>const& beta)
{
    auto iZ=beta[0]==0?1:0;
    auto n=beta.size()-1-iZ;
    auto out=std::vector<double>(n);
    for (std::size_t i=0; i<out.max_size(); ++i)
        out[n-i-1]= std::log(1.0/beta[n+iZ]-1.0/beta[n-1+iZ]);
    return out;
}



auto S_to_i_fract_beta(std::vector<double>const& S, std::vector<std::size_t> num_samples, std::size_t max_samples)
{
    auto T_samples=var::apply_to([max_samples](auto x){return 1.0*x/max_samples;},num_samples);
    std::vector<std::size_t> frac_i(S.size());
    std::vector<double> beta(S.size());
    auto n=S.size()-1;
    auto nf=num_samples.size()-1;
    std::size_t i_frac=nf;
    auto i_frac0=i_frac>0?i_frac-1:0;
    // qué onda? da para que sea la ultima fraccion?
    // el threshold debería ser a mitad de camino entre 0 y la proxima fraccion.
    
    auto threshold_s=-std::log(T_samples[i_frac]/2.0+T_samples[i_frac0]/2.0);
    
    while ((i_frac>0)&&(S[n]>threshold_s))
    {
        --i_frac;
        threshold_s=-std::log(T_samples[i_frac]/2.0+T_samples[i_frac-1]/2.0);
    }
    // tenemos el primer i_frac! y el beta (debe ser 1)
    beta[n-0]=1.0;
    frac_i[n-0]=i_frac;
    double current_T=1/T_samples[i_frac];
    // el proximo threshold esta dado por
    for (std::size_t i=1; i<S.size(); ++i)
    {
        if ((i_frac>0)&&(current_T+std::exp(S[n-i])>T_samples[i_frac-1])){
            --i_frac;
            beta[n-i]=1.0;
            frac_i[n-i]=i_frac;
            current_T=T_samples[i_frac];
            
        }else{
            if (std::isfinite(S[n-i]))
            {
                current_T=current_T+std::exp(S[n-i]);
                auto gbeta=1/current_T;
                auto fraction_below=i_frac>0?T_samples[i_frac-1]:0.0;
                
                auto r_beta=(gbeta-fraction_below)/(1.0/T_samples[i_frac]-fraction_below);
                beta[n-i]=r_beta;
            }else
                beta[n-i]=0;
            frac_i[n-i]=i_frac;
            
        }
    }
    return std::tuple(std::move(frac_i), std::move(beta));
}

std::vector<std::size_t> i_frac_to_i_fracs(std::size_t i_fr, double beta, std::size_t num_fractions)
{
    if (i_fr==0)
        return {i_fr};
    else if ((i_fr+1<num_fractions)&&(beta==1))
    {
        return {i_fr-1,i_fr,i_fr+1};
    }else{
        return {i_fr-1, i_fr};
    }
}

auto i_frac_to_i_fracs(const std::vector<std::size_t>& i_frac, const std::vector<double>& beta)
{
    std::vector<std::vector<std::size_t>> i_fracs(i_frac.size());
    for (std::size_t i=0; i<i_frac.size(); ++i)
    {
        i_fracs[i]=i_frac_to_i_fracs(i_frac[i],beta[i], i_frac.back());
    }
    return i_fracs; 
}



template <class Parameters>
void initial_beta_dts(thermo_fraction_mcmc<Parameters> &initial_data) {
    auto initial_logL0 = calculate_initial_logL0(initial_data);
    double S_max =
        std::log(-get<mean<logL>>(initial_logL0.begin()->second()())()()*initial_data.max_samples/initial_data.samples_size[0]);
    
    double S_min= std::log(initial_data.max_samples/initial_data.samples_size[0]);
    
    auto n = initial_data.beta.size();
    double dS = (S_max-S_min) / (n - 2);
    std::vector<double> S(n);
    S[0]=std::numeric_limits<double>::infinity();
    S[1]=S_max;
    S[n-1]=S_min;
    
    for (std::size_t i=2; i+2<S.size(); ++i)
    {
            S[i]=S_min+(i-1)*dS;
    }
    
    auto[i_frac,beta]=S_to_i_fract_beta(S,initial_data.samples_size,initial_data.max_samples);
    
    initial_data.i_fractions=i_frac_to_i_fracs(i_frac,beta);
    initial_data.S=std::move(S);
    initial_data.beta=std::move(beta);
    initial_data.i_fraction=std::move(i_frac);
}

template <class Parameters>
ensemble<std::map<std::size_t,logL_statistics>>
make_logL_statistics(ensemble<fraction_mcmc<Parameters>>const& walkers) {
    ensemble<std::map<std::size_t,logL_statistics>> out(walkers.size());
    for (auto i = 0ul; i < walkers.size(); ++i) {
            for(auto const& e: walkers[i].logL_map)
                out[i][e.first] = logL_statistics(get<logL>(e.second));
        }    
    return out;
}


template <class Parameters>
by_beta<ensemble<std::map<std::size_t,logL_statistics>>>
make_logL_statistics(by_beta<ensemble<fraction_mcmc<Parameters>>>const& walkers) {
    by_beta<ensemble<std::map<std::size_t,logL_statistics>>> out;
    out.reserve(walkers.size());
    for (auto i = 0ul; i < walkers.size(); ++i) {
        out.push_back(ensemble<std::map<std::size_t,logL_statistics>>(walkers[i].size()));
        
        for (std::size_t j = 0; j < walkers[i].size(); ++j) {
            for(auto const& e: walkers[i][j].logL_map)
            out[i][j][e.first] = logL_statistics(get<logL>(e.second));
        }
    }
    return out;
}


template <class LikelihoodModel,
         class FuncTable, class Parameters, class Variables, class DataType>
Maybe_error<std::map<std::size_t,logLs>>
logLikelihoods(FuncTable &f,
               const LikelihoodModel &lik,
               Parameters const &p,  const std::vector<DataType> &y,const std::vector<Variables> &var, std::vector<std::size_t> i_fracs) {
    std::map<std::size_t,logLs> out;
    for (auto i_frac: i_fracs){
        auto Maybe_lik=logLikelihood(f,lik,p,y[i_frac],var[i_frac]);
        if(!Maybe_lik.valid())
            return Maybe_lik.error();
        else
            out[i_frac]=std::move(Maybe_lik.value());
    }
    return out;
}

template <class LikelihoodModel,
         class FuncTable, class Parameters, class Variables, class DataType>
Maybe_error<bool>
update_logLikelihoods(FuncTable &f,
               const LikelihoodModel &lik,
               Parameters const &p,  const std::vector<DataType> &y,const std::vector<Variables> &var, std::vector<std::size_t> i_fracs, std::map<std::size_t,logLs>& rlogL , std::map<std::size_t,logL_statistics>& logL_sta) {
    Maybe_error<bool> out(true);
    for (auto i_frac: i_fracs){
        if (rlogL.find(i_frac)==rlogL.end()){
        auto Maybe_lik=logLikelihood(f,lik,p,y[i_frac],var[i_frac]);
        if(!Maybe_lik.valid())
        {
            rlogL[i_frac]=nan_logL();
            out=error_message(out.error()()+ Maybe_lik.error()());
        }
        else{
            rlogL[i_frac]=std::move(Maybe_lik.value());
            logL_sta[i_frac]() &= get<logL>(rlogL[i_frac]);
           }
        }
    }
    return out;
}


template <class FunctionTable, class Prior,class Lik, class Variables,class DataType,
         class Parameters=std::decay_t<
             decltype(sample(std::declval<mt_64i &>(), std::declval<Prior&>()))>>
//   requires (is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<FunctionTable,Lik,Parameters,Variables,DataType>)
auto init_fraction_mcmc(FunctionTable& f, mt_64i &mt, Prior const & pr, const Lik& lik,
               const std::vector<DataType> &y, const std::vector<Variables> &x, std::vector<std::size_t > i_frac) {
    auto& priorsampler=pr;
    auto par = sample(mt,priorsampler);
    auto logP = logPrior(pr,par);
    std::map<std::size_t,logLs> t_logsLs;
    auto logL=logLikelihoods(f,lik,par.to_value(), y,x,i_frac);
    while(!(logP)|| !(logL))
    {
        par = sample(mt,priorsampler);
        logP = logPrior(pr,par);
        logL=logLikelihoods(f,lik,par.to_value(), y,x,i_frac);
    }
    return fraction_mcmc<Parameters>{std::move(par), logP.value(), std::move(logL.value())};
}




template <class Parameters>
auto calculate_controler_step(const thermo_fraction_mcmc<Parameters> &current,
                              const by_beta<double> &beta,
                              std::string equalizing_paramter,
                              double desired_acceptance,
                              std::string variance_approximation) {
    bool reached_global_beta_1=
        (current.i_fraction.back()+1==current.samples_size.size())&&(current.beta.back()==1);
    if (equalizing_paramter == "Acceptance_vfm") {
        auto A = calculate_Acceptance(current);
        auto d = std::vector<double>(A.size() - (reached_global_beta_1?1:0));
        for (std::size_t i = 0; i < d.size(); ++i) {
            d[i] = std::max(std::min(1.0, A[i + 1] - A[i]), -1.0);
        }
        if(!reached_global_beta_1)
            d.back()=std::max(std::min(1.0, -std::exp(std::min(0.0, -A[A.size()-1]))),
                                -1.0);
        
        return d;
    } else if (equalizing_paramter == "deltaBeta_deltaL_vfm") {
        auto dBdL = calculate_deltaBeta_deltaL(current);
        auto d = std::vector<double>(dBdL.size() - (reached_global_beta_1?1:0));
        
        for (std::size_t i = 0; i+1 < dBdL.size(); ++i) {
            d[i] = std::max(std::min(1.0, std::exp(std::min(0.0, -dBdL[i + 1])) -
                                              std::exp(std::min(0.0, -dBdL[i]))),
                            -1.0);
        }
        if(!reached_global_beta_1)
            d.back()=std::max(std::min(1.0, -std::exp(std::min(0.0, -dBdL[dBdL.size()-1]))),
                                -1.0);
        return d;
    } else
        return std::vector<double>{};
}





template <class FunctionTable, class Prior, class Likelihood, class Variables,
         class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(
        is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)
//    requires (is_prior<Prior,Parameters,Variables,DataType>&&
//    is_likelihood_model<FunctionTable,Likelihood,Parameters,Variables,DataType>)
void insert_high_temperture_beta(FunctionTable &f, std::size_t tested_index,
                                 thermo_fraction_mcmc<Parameters> &current,
                                 double new_S,
                                 ensemble<mt_64i> &mt, Prior const &prior,
                                 Likelihood const &lik, const DataType &y,
                                 const Variables &x) {
    
    auto n_walkers = current.walkers[0].size();
    ensemble<std::size_t> i_walker(n_walkers);
    ensemble<fraction_mcmc<Parameters>> walker(n_walkers);
    emcee_Step_statistics emcee_stat(emcee_Step_statistics{});
    Thermo_Jump_statistics thermo_stat(Thermo_Jump_statistics{});
    auto ff = f.fork(omp_get_max_threads());
    
    auto last_walker = current.last_walker();
    
    auto i_fracts_0=std::vector<std::size_t>{0};

#pragma omp parallel for // collapse(2)
    for (std::size_t iw = 0; iw < n_walkers; ++iw) {
        i_walker[iw] = iw + last_walker + 1;
        auto i_th = omp_get_thread_num();
        walker[iw] = init_fraction_mcmc(ff[i_th], mt[i_th], prior, lik, y, x,i_fracts_0 );
    }
    f += ff;
    // update S, i_frac etc
    auto walker_sta = make_logL_statistics(walker);
    current.S.insert(current.S.begin() + tested_index + 1, new_S);
    auto [new_i_fracs, new_beta]=S_to_i_fract_beta(current.S,current.samples_size,current.max_samples);
    
    current.beta= std::move(new_beta);
    current.i_fraction=std::move(new_i_fracs);
    current.i_fractions= i_frac_to_i_fracs(current.i_fraction, current.beta);
    current.i_walkers.insert(current.i_walkers.begin(), i_walker);
    current.walkers.insert(current.walkers.begin(), walker);
    current.walkers_sta.insert(current.walkers_sta.begin(), walker_sta);
    current.emcee_stat.insert(current.emcee_stat.begin(), emcee_stat);
    current.thermo_stat.insert(current.thermo_stat.begin(), thermo_stat);
 }

template <class Parameters>
void remove_high_temperture_beta(thermo_fraction_mcmc<Parameters> &current,
                                 std::size_t tested_index) {
    
    current.beta.erase(current.beta.begin() + tested_index + 1);
    current.i_walkers.erase(current.i_walkers.begin());
    current.walkers.erase(current.walkers.begin());
    current.walkers_sta.erase(current.walkers_sta.begin());
    current.emcee_stat.erase(current.emcee_stat.begin());
    current.thermo_stat.erase(current.thermo_stat.begin());
 }




template <class FunctionTable, class Prior, class Likelihood, class Variables,
         class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(
        is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)
//    requires (is_prior<Prior,Parameters,Variables,DataType>&&
//    is_likelihood_model<FunctionTable,Likelihood,Parameters,Variables,DataType>)
void adjust_fraction_beta(FunctionTable &f, std::size_t iter,
                 std::size_t adapt_beta_every, double acceptance_upper_limit,
                 double acceptance_lower_limit,
                 thermo_fraction_mcmc<Parameters> &current, 
                 ensemble<mt_64i> &mt, Prior const &prior,
                 Likelihood const &lik, const DataType &ys, const Variables &xs) {
    if ((iter > 0) && (current.num_samples() > 0) &&
        (iter % adapt_beta_every == 0)) {
        assert(current.beta[current.beta.size() - 1] = 1);
        std::size_t tested_index = 1;
        // auto A=calculate_Acceptance(current);
        auto A = calculate_deltaBeta_deltaL(current);
        if (std::exp(std::min(0.0, -A[tested_index])) > acceptance_upper_limit)
            remove_high_temperture_beta(current, tested_index);
        else
            //             if (std::reduce(A.begin()+tested_index, A.end(),0.0,
            // [acceptance_lower_limit](double a,double b){
            //          if ((a==1)||(b<acceptance_lower_limit))
            //               return 1.0;
            //          else
            // return 0.0;})==1.0)
            if (std::exp(std::min(0.0, -A[tested_index + 1])) <
                acceptance_lower_limit) {
                double new_S;
                if (current.beta[tested_index] == 0)
                {
                    new_S = current.S[tested_index + 1]+
                            0.5* (current.S[tested_index + 1]-current.S[tested_index + 2]);
                }
                else{
                    new_S = 0.5*(current.S[tested_index + 1] +current.S[tested_index]);
                }
                insert_high_temperture_beta(f, tested_index, current,  new_S,
                                            mt, prior, lik, ys, xs);
            }
    }
}
       
    







template <class Parameters>
void adapt_fraction_beta(std::size_t iter,
                thermo_fraction_mcmc<Parameters> &current,
                std::size_t adapt_beta_every,
                std::string equalizing_paramter,
                std::string variance_approximation,
                         double desired_acceptance,
                double nu, double t0) {
    if ((iter > 0) && (current.num_samples() > 0) &&
        (iter % adapt_beta_every == 0)) {
       // assert(beta[beta.size() - 1] = 1);
        std::size_t tested_index = 1;
        double kappa = 1.0 / nu * t0 / (t0 + iter);
        
        std::vector<double> global_beta= current.beta;
        
        
        auto d =
            calculate_controler_step(current, global_beta, equalizing_paramter,
                                     desired_acceptance, variance_approximation);
            std::vector<double> T(global_beta.size() - (global_beta[0] == 0 ? 1 : 0));
            for (std::size_t i = 0; i < T.size(); ++i)
                T[i] = 1.0 / global_beta[i + (global_beta[0] == 0 ? 1 : 0)];
            auto S= current.S;
            
            if (equalizing_paramter.ends_with("vfm")) {
                for (std::size_t i = 0; i < d.size(); ++i) {
                    S[i+1] += kappa * d[i];
                }
            } else {
                auto dbds = calculate_d_beta_d_s(global_beta);
                for (std::size_t i = 0; i < S.size(); ++i) {
                    S[i+1] += kappa * d[i] / dbds[i];
                }
            }
            std::sort(S.begin(), S.end(), std::greater<double>());
            
            auto[frac_i, beta_i]=S_to_i_fract_beta(S,current.samples_size,current.max_samples);
            
            current.S=std::move(S);
            current.beta=std::move(beta_i);
            current.i_fraction=std::move(frac_i);
            current.i_fractions=i_frac_to_i_fracs(current.i_fraction,current.beta);            
       }
}

template <class Parameters>
void update(by_beta<ensemble<std::map<std::size_t,logL_statistics>>>  &walkers_sta,
            by_beta<ensemble<fraction_mcmc<Parameters>>> const &walkers) {
    for (auto i = 0ul; i < walkers_sta.size(); ++i)
        for (std::size_t j = 0; j < walkers_sta[i].size(); ++j) {
            for (auto const& e: walkers[i][j].logL_map)
            walkers_sta[i][j][e.first]() &= get<logL>(e.second);
        }
}



template<class Walker,class logP, class logLik>
double calc_emcee_jump_delta( Walker const& current, logP const& ca_logP, logLik const& ca_logL, double beta, std::size_t i_frac)
{
    auto delta_P=ca_logP-current.logP;
    auto& calogL=ca_logL.find(i_frac)->second;
    auto delta_lik= beta*(get<logL>(calogL)()-get<logL>(current.logL(i_frac))());
    if (i_frac>0)
    {
        auto& calogL0=ca_logL.find(i_frac-1)->second;
        
        delta_P= delta_P+get<logL>(calogL0)()-get<logL>(current.logL(i_frac-1))();
    }
    auto out= delta_P+delta_lik;
    // it is possible that current is not finite
    if (std::isfinite(out))
        return out;
    else
        return +20.0;
    
}


auto get_samples_size(auto const& y)
{
 std::vector<std::size_t> out(y.size());
for (std::size_t i=0; i<y.size();++i)
{
    out[i]=y[i]().size();
}
return out; 
}






template <class FunctionTable, class Prior, class Likelihood, class Variables,
         class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(
        is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)
//    requires (is_prior<Prior,Parameters,Variables,DataType>&&
//    is_likelihood_model<FunctionTable,Likelihood,Parameters,Variables,DataType>)
auto init_thermo_fraction_mcmc(FunctionTable &f, std::size_t n_walkers,
                      by_beta<double> const &beta,
                     by_beta<std::size_t> const &i_fraction,
                     ensemble<mt_64i> &mt,
                      Prior const &prior, Likelihood const &lik,
                      const DataType &y, const Variables &x) {
    assert(beta.size()==i_fraction.size());
    
    auto i_fractions=i_frac_to_i_fracs(i_fraction,beta);
    
    
    by_beta<ensemble<std::size_t>> i_walker(beta.size(),
                                            by_beta<std::size_t>(n_walkers));
    by_beta<ensemble<fraction_mcmc<Parameters>>> walker(
        beta.size(), by_beta<fraction_mcmc<Parameters>>(n_walkers));
    by_beta<emcee_Step_statistics> emcee_stat(beta.size(),
                                              emcee_Step_statistics{});
    by_beta<Thermo_Jump_statistics> thermo_stat(beta.size() - 1,
                                                Thermo_Jump_statistics{});
    auto ff = f.fork(omp_get_max_threads());

#pragma omp parallel for // collapse(2)
    for (std::size_t iw = 0; iw < n_walkers; ++iw) {
       for (std::size_t i = 0; i < beta.size(); ++i) {
            i_walker[i][iw] = iw + i * n_walkers;
            auto i_th = omp_get_thread_num();
            walker[i][iw] = init_fraction_mcmc(ff[i_th], mt[i_th], prior, lik, y, x, i_fractions[i]);
        }
    }
    f += ff;
    std::vector<std::size_t> samples_size(y.size());
    for (std::size_t i=0; i<y.size();++i)
    {
        samples_size[i]=y[i]().size();
    }
    std::size_t max_samples=var::max(samples_size);
    
   
    auto global_beta=i_fract_beta_to_global_beta(i_fraction,beta,samples_size,max_samples);
    auto S=global_beta_to_S(global_beta);
    auto walker_sta = make_logL_statistics(walker);
    return thermo_fraction_mcmc<Parameters>{samples_size,max_samples,i_fraction,i_fractions,beta, S,  walker,    walker_sta,
                                   i_walker, emcee_stat, thermo_stat};
}



template <class FunctionTable,
         class Likelihood, class Variables, class DataType,
         class Parameters>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>,
                                        var::FuncMap_St> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                                 DataType>)
void update_likelihoods_all(FunctionTable &f,
                                       thermo_fraction_mcmc<Parameters> &current, 
                                       Likelihood const &lik,
                                       const DataType &y, const Variables &x) {
    auto n_walkers = current.walkers[0].size();
    auto n_beta = current.beta.size();
    auto n_par = current.walkers[0][0].parameter.size();
    
    
    auto ff = f.fork(omp_get_max_threads());
    std::size_t num_threads = omp_get_max_threads();
    auto total_walkers = n_beta * n_walkers;
    std::size_t walkers_per_thread =
        std::ceil(1.0 * total_walkers / num_threads);
    
 #pragma omp parallel for // collapse(2)
        for (std::size_t i_thread = 0; i_thread < num_threads; ++i_thread) {
            for (std::size_t iwb = i_thread * walkers_per_thread;
                 iwb < std::min(walkers_per_thread * (i_thread + 1),
                                total_walkers);
                 ++iwb) {
                //  dur.record("begin_loop_walker", iwb * 2);
             //   std::size_t ib = iwb / (n_walkers / 2);
             //   std::size_t i = iwb - ib * n_walkers / 2;
                std::size_t iw = iwb / (n_beta);
                std::size_t ib = iwb - iw * n_beta;
                
                auto i_th = omp_get_thread_num();
                auto& cu_par=current.walkers[ib][iw].parameter;
                auto Maybe_i=update_logLikelihoods(ff[i_th],lik,cu_par.to_value(),y,x,
                                                     current.i_fractions[ib],
                                                     current.walkers[ib][iw].logL_map,
current.walkers_sta[ib][iw]);
            }
        }
    f += ff;
  }



template <class FunctionTable, std::size_t N, class Observer, class Prior,
         class Likelihood, class Variables, class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>,
                                        var::FuncMap_St> &&
             is_prior<Prior, Parameters, Variables, DataType> &&
             is_likelihood_model<FunctionTable, Likelihood, Parameters, Variables,
                                 DataType>)
void step_stretch_thermo_fraction_mcmc(FunctionTable &f, std::size_t &iter,
                              var::Event_Timing<N> &dur,
                              thermo_fraction_mcmc<Parameters> &current, Observer &obs,
                              ensemble<mt_64i> &mt,
                              Prior const &prior, Likelihood const &lik,
                              const DataType &y, const Variables &x,
                              double alpha_stretch = 2) {
    dur.record("stretch_start");
    auto n_walkers = current.walkers[0].size();
    auto n_beta = current.beta.size();
    auto n_par = current.walkers[0][0].parameter.size();
    
    std::uniform_int_distribution<std::size_t> uniform_walker(0,
                                                              n_walkers / 2 - 1);
    std::vector<std::uniform_int_distribution<std::size_t>> udist(
        omp_get_max_threads(), uniform_walker);
    
    std::uniform_real_distribution<double> uniform_stretch_zdist(
        1.0 / alpha_stretch, alpha_stretch);
    std::vector<std::uniform_real_distribution<double>> zdist(
        omp_get_max_threads(), uniform_stretch_zdist);
    
    std::uniform_real_distribution<double> uniform_real(0, 1);
    std::vector<std::uniform_real_distribution<double>> rdist(
        omp_get_max_threads(), uniform_real);
    
    auto ff = f.fork(omp_get_max_threads());
    std::vector<by_beta<emcee_Step_statistics>> emcee_stat(
        omp_get_max_threads(), by_beta<emcee_Step_statistics>(n_beta));
    
    dur.record("stretch_before_loop");
    std::size_t num_threads = omp_get_max_threads();
    
    auto total_walkers_per_half = n_beta * n_walkers / 2;
    std::size_t walkers_per_thread =
        std::ceil(1.0 * total_walkers_per_half / num_threads);
    
    // std::size_t n_beta_f = std::ceil(1.0*n_beta / num_threads);
    
    for (bool half : {false, true}) {
#pragma omp parallel for // collapse(2)
        for (std::size_t i_thread = 0; i_thread < num_threads; ++i_thread) {
            for (std::size_t iwb = i_thread * walkers_per_thread;
                 iwb < std::min(walkers_per_thread * (i_thread + 1),
                                total_walkers_per_half);
                 ++iwb) {
                //  dur.record("begin_loop_walker", iwb * 2);
            //    std::size_t ib = iwb / (n_walkers / 2);
            //    std::size_t i = iwb - ib * n_walkers / 2;
                
                std::size_t i = iwb / (n_beta);
                std::size_t ib = iwb - i * n_beta;
                
                
                
                auto i_th = omp_get_thread_num();
                
                auto iw = half ? i + n_walkers / 2 : i;
                auto j = udist[i_th](mt[i_th]);
                auto jw = half ? j : j + n_walkers / 2;
                // we can try in the outer loop
                
                auto r = rdist[i_th](mt[i_th]);
                auto& cu_par=current.walkers[ib][iw].parameter;
                auto Maybe_i=update_logLikelihoods(ff[i_th],lik,cu_par.to_value(),y,x,
                                                     current.i_fractions[ib],
                                                     current.walkers[ib][iw].logL_map,
current.walkers_sta[ib][iw]);
                
                // candidate[ib].walkers[iw].
                auto [ca_par, z] = stretch_move(mt[i_th], rdist[i_th],
                                                current.walkers[ib][iw].parameter,
                                                current.walkers[ib][jw].parameter);
                
                auto ca_logP = logPrior(prior, ca_par);
                auto i_fracts= current.i_fractions[ib];
                
                auto ca_logL = logLikelihoods(ff[i_th], lik, ca_par.to_value(), y, x,i_fracts);
                
                
                
                if (!((ca_logP.valid()) && (ca_logL.valid()))) {
                    fails(emcee_stat[i_th][ib]());
                } else {
                    auto dthLogL =calc_emcee_jump_delta(current.walkers[ib][iw], ca_logP.value(),ca_logL.value(),current.beta[ib], current.i_fraction[ib]);
                    auto pJump =
                        std::min(1.0, std::pow(z, n_par - 1) * std::exp(dthLogL));
                    if (pJump >= r) {
                        current.walkers[ib][iw].parameter = std::move(ca_par);
                        current.walkers[ib][iw].logP = ca_logP.value();
                        current.walkers[ib][iw].logL_map = ca_logL.value();
                        succeeds(emcee_stat[i_th][ib]());
                    } else {
                        fails(emcee_stat[i_th][ib]());
                    }
                }
                //}
                
                //      dur.record("end_loop_walker", ib * 2 + 1);
            }
        }
    }
    
    dur.advance(n_beta * 2);
    dur.record("stretch_after_loop");
    #pragma omp parallel for 
    for (std::size_t ib = 0; ib < n_beta; ++ib) {
        for (auto &e : emcee_stat)
            current.emcee_stat[ib]() += e[ib]();
    }
    f += ff;
    ++iter;
    update(current.walkers_sta, current.walkers);
    
    dur.record("stretch_function_end");
}





template <class LikelihoodModel,
         class FuncTable, class Parameters, class Variables, class DataType, class Observer>
void thermo_fraction_jump_mcmc(FuncTable &f,
                               const LikelihoodModel &lik,  const std::vector<DataType> &y, const std::vector<Variables> &x,std::size_t iter, thermo_fraction_mcmc<Parameters> &current,
                      Observer &obs,  mt_64i &mt,
                      ensemble<mt_64i> &mts, std::size_t thermo_jumps_every) {
    if ((iter > 0) && (current.num_samples() > 0) &&
        (iter % (thermo_jumps_every) == 0)) {
        auto ff = f.fork(omp_get_max_threads());
        
        std::uniform_real_distribution<double> uniform_real(0, 1);
        auto n_walkers = current.get_Walkers_number();
        auto n_beta = current.beta.size();
        auto n_par = current.walkers[0][0].parameter.size();
        
        WalkerIndexes shuffeld_walkers(n_walkers);
        std::iota(shuffeld_walkers.begin(), shuffeld_walkers.end(), 0);
        std::shuffle(shuffeld_walkers.begin(), shuffeld_walkers.end(), mt);
        std::vector<std::uniform_real_distribution<double>> rdist(
            omp_get_max_threads(), uniform_real);
        
        std::vector<by_beta<Thermo_Jump_statistics>> thermo_stat(
            omp_get_max_threads(), by_beta<Thermo_Jump_statistics>(n_beta - 1));
        
        std::size_t num_threads = omp_get_max_threads();
        
        auto total_walkers_per_half = (n_beta - 1) * n_walkers / 2;
        std::size_t walkers_per_thread =
            std::ceil(1.0 * total_walkers_per_half / num_threads);
        
        // std::size_t n_beta_f = std::ceil(1.0*n_beta / num_threads);

#pragma omp parallel for // collapse(2)
        for (std::size_t i_thread = 0; i_thread < num_threads; ++i_thread) {
            for (std::size_t iwb = i_thread * walkers_per_thread;
                 iwb < std::min(walkers_per_thread * (i_thread + 1),
                                total_walkers_per_half);
                 ++iwb) {
                //  dur.record("begin_loop_walker", iwb * 2);
//                std::size_t ib = iwb / (n_walkers / 2);
//                std::size_t i = iwb - ib * n_walkers / 2;
                
                std::size_t i = iwb / (n_beta-1);
                std::size_t ib = iwb - i * (n_beta-1);
                
                
                auto i_th = omp_get_thread_num();
                auto iw = shuffeld_walkers[i];
                auto jw = shuffeld_walkers[i + n_walkers / 2];
                
                auto r = rdist[i_th](mts[i_th]);
                auto beta0=current.beta[ib]==1?0:current.beta[ib];
                auto i_frac=current.i_fraction[ib+1];
                
                
                
                double logA =
                    calc_logA(beta0, current.beta[ib + 1], current.walkers[ib][iw].logL(i_frac),
                              current.walkers[ib + 1][jw].logL(i_frac));
                auto pJump = std::min(1.0, std::exp(logA));
                if (pJump > r) {
                    auto Maybe_i=update_logLikelihoods(ff[i_th],lik,current.walkers[ib][iw].parameter.to_value(),x,y,current.i_fractions[ib+1],current.walkers[ib][iw].logL_map,
current.walkers_sta[ib][iw]);
                    
                    auto Maybe_j=update_logLikelihoods(ff[i_th],lik,current.walkers[ib+1][jw].parameter.to_value(),x,y,current.i_fractions[ib],current.walkers[ib+1][jw].logL_map,
current.walkers_sta[ib+1][jw]);
                    if(Maybe_i.valid()&& Maybe_j.valid()){
                    
                    std::swap(current.walkers[ib][iw], current.walkers[ib + 1][jw]);
                    std::swap(current.i_walkers[ib][iw], current.i_walkers[ib + 1][jw]);
                    succeeds(thermo_stat[i_th][ib]());
                    }else {
                        fails(thermo_stat[i_th][ib]());
                    }     
                } else {
                    fails(thermo_stat[i_th][ib]());
                }
            }
        }
        
        for (std::size_t ib = 0; ib < n_beta - 1; ++ib) {
            for (auto &e : thermo_stat)
                current.thermo_stat[ib]() += e[ib]();
        }
    }
}



template <class FunctionTable, class Prior,class Lik, class Variables,class DataType,
         class Parameters=std::decay_t<
             decltype(sample(std::declval<mt_64i &>(), std::declval<Prior&>()))>>
//   requires (is_prior<Prior,Parameters,Variables,DataType>&& is_likelihood_model<FunctionTable,Lik,Parameters,Variables,DataType>)
Maybe_error<bool> calc_mcmc(FunctionTable& f, Prior const & pr, const Lik& lik,
                            const DataType &y, const Variables &x, fraction_mcmc<Parameters>& t_mcmc, std::vector<std::size_t> const & i_fracs) {
    auto par = t_mcmc.parameter;
    auto logP = logPrior(pr,par);
    auto t_logLs = logLikelihoods(f,lik,par.to_value(), y,x, i_fracs);
    
    if (logP.valid()&& t_logLs.valid())
    {
        t_mcmc.logP=logP.value();
        t_mcmc.logL_map=t_logLs.value();
        return true;
    }
    else
    {
        t_mcmc.logP=std::numeric_limits<double>::quiet_NaN();
        for (auto& e:t_mcmc.logL_map )
            get<logL>(e.second)()=std::numeric_limits<double>::quiet_NaN();;
    }
    return error_message(logP.error()()+t_logLs.error()());
}





template <bool Adapt_beta,class FunctionTable, class Algorithm, class Prior, class Likelihood,
         class Variables, class DataType, class Reporter, class mcmc,
         class Parameters, class timepoint>
    requires(
        is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&&
//             is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto thermo_fraction_evidence_loop(
    FunctionTable &&f,
    new_thermodynamic_integration<Algorithm, Reporter> &&therm,
    Prior const &prior, Likelihood const &lik, const std::vector<DataType> &ys,
    const std::vector<Variables> &xs, mcmc mcmc_run, std::size_t iter,
    thermo_fraction_mcmc<Parameters> &current, Reporter &rep,
    mt_64i &mt, std::vector<mt_64i> &mts,
    const timepoint &start, const std::chrono::duration<double>& previous_duration) {
    
    
    var::Event_Timing<200> even_dur(start);
    std::ofstream event_file(f.file()+"_event_timing.csv");
    
    while (!mcmc_run.second) {
        even_dur.record("main_loop_start");
        const auto end = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration<double>(end - start)+previous_duration;
        
        report_all(f, iter, dur, rep, current, prior, lik, ys, xs, mts,
                   mcmc_run.first);
        if constexpr (Adapt_beta)
        {
            adapt_fraction_beta(iter, current, therm.adapt_beta_every(),therm.adapt_beta_equalizer(),therm.adapt_beta_variance(),therm.desired_acceptance(),therm.adapt_beta_nu(),therm.adapt_beta_t0());
            update_likelihoods_all(f,current,lik,ys,xs);
            
            if(therm.adjust_beta()){
                adjust_fraction_beta(f,iter,therm.adapt_beta_every(),therm.acceptance_upper_limit(),therm.acceptance_lower_limit(),current,mts,prior,lik,ys,xs);
            update_likelihoods_all(f,current,lik,ys,xs);
                
            }
            if (iter%therm.adapt_beta_every()==0)   
                current.reset_statistics();
            
            
            
        }
        
        step_stretch_thermo_fraction_mcmc(f, iter,even_dur, current, rep, mts, prior, lik,
                                 ys, xs);
        
        even_dur.record("befor_thermo_jump");  
        
        thermo_fraction_jump_mcmc(f,lik,xs,ys,iter, current, rep,  mt, mts,
                         therm.thermo_jumps_every());
        
        even_dur.record("after_thermo_jump");  
        
        even_dur.record("after_report_all");  
        report_point(f, iter);
        
        // using geg=typename
        // decltype(checks_convergence(std::move(mcmc_run.first), current))::eger;
        mcmc_run = checks_convergence(std::move(mcmc_run.first), current);
        even_dur.record("after_checks_convergence");
        if (iter==1)
            even_dur.report_title(event_file);
        even_dur.report_iter(event_file,iter);
        if (iter%10==0)
            event_file.flush();
        
        
    }
    return std::pair(std::move(mcmc_run.first), current);
}









template <bool Adapt_beta,class FunctionTable, class Algorithm, class Prior, class Likelihood,
         class Variables, class DataType, class Reporter>
    requires(
        is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&&
//             is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto thermo_fraction_evidence(FunctionTable &&f,
                     new_thermodynamic_integration<Algorithm, Reporter> &&therm,
                     Prior const &prior, Likelihood const &lik,
                     const std::vector<DataType> &ys, const std::vector<Variables> &xs) {
    auto a = therm.algorithm();
    auto mt = init_mt(therm.initseed());
    auto n_walkers = therm.num_scouts_per_ensemble();
    auto mts = init_mts(mt, omp_get_max_threads());
    by_beta<double> beta_run;
    if constexpr (Adapt_beta)
    {
        beta_run=by_beta<double>(therm.beta_size(),0);
    }
    else{
        auto beta = new_get_beta_list(
            therm.beta_size(), therm.beta_upper_size(), therm.beta_medium_size(),
            therm.beta_upper_value(), therm.beta_medium_value(), therm.stops_at(),
            therm.includes_zero());
        
        auto it_beta_run_begin = beta.rend() - beta.size();
        auto it_beta_run_end = beta.rend();
        beta_run = by_beta<double>(it_beta_run_begin, it_beta_run_end);
    }
    by_beta<std::size_t> i_fractions(beta_run.size(),0);
    
    auto current =
        init_thermo_fraction_mcmc(f, n_walkers, beta_run, i_fractions, mts, prior, lik, ys, xs);
    // auto n_par = current.walkers[0][0].parameter.size();
    
    if constexpr (Adapt_beta)
    {
        initial_beta_dts(current);
        
    }
    
    
    auto mcmc_run = checks_convergence(std::move(a), current);
    
    std::size_t iter = 1;
    const auto start = std::chrono::high_resolution_clock::now();
    auto &rep = therm.reporter();
    report_title(rep, current, lik, ys, xs);
    report_title(f, "Iter");
    report_model_all(rep, prior, lik, ys, xs, beta_run);
    std::chrono::duration<double> previous_duration(0.0);
    return thermo_fraction_evidence_loop<Adapt_beta>(
        f,
        std::forward<new_thermodynamic_integration<Algorithm, Reporter>>(therm),
        prior, lik, ys, xs, mcmc_run, iter, current, rep,  mt, mts, start,previous_duration);
}


template <class Prior, class Parameters = std::decay_t<decltype(sample(
                          std::declval<mt_64i &>(), std::declval<Prior &>()))>>
auto create_thermo_fraction_mcmc(std::size_t n_walkers,
                                 by_beta<double> const &beta,
                        mt_64i &mt, Prior const &pr,std::vector<std::size_t> samples_size,std::size_t max_samples) {
    
    auto &priorsampler = pr;
    auto par = sample(mt, priorsampler);
    
    by_beta<ensemble<std::size_t>> i_walker(beta.size(),
                                            by_beta<std::size_t>(n_walkers));
    by_beta<ensemble<fraction_mcmc<Parameters>>> walker(
        beta.size(), by_beta<fraction_mcmc<Parameters>>(n_walkers, fraction_mcmc<Parameters>{par}));
    by_beta<ensemble<std::map<std::size_t,logL_statistics>>> walker_sta(
        beta.size(), by_beta<std::map<std::size_t,logL_statistics>>(n_walkers));
    
    by_beta<emcee_Step_statistics> emcee_stat(beta.size());
    by_beta<Thermo_Jump_statistics> thermo_stat(beta.size() - 1);
    
    
    by_beta<std::size_t> i_fraction(beta.size());
    by_beta<std::vector<std::size_t>> i_fractions(beta.size());
    by_beta<double> S(beta.size());
    by_beta<ensemble<std::map<std::size_t,logL_statistics>>> walkers_sta;
    
    
    
    return thermo_fraction_mcmc<Parameters>{samples_size,max_samples, i_fraction, i_fractions,beta, S,   walker,     walker_sta,
                                   i_walker, emcee_stat, thermo_stat};
}

template <class Parameters, class Duration>
bool extract_iter(std::istream &f, std::size_t &iter, Duration &dur,
                  thermo_fraction_mcmc<Parameters> &data) {
    std::size_t i_beta = 0;
    std::size_t i_walker = 0;
    std::size_t i_par = 0;
    
    std::size_t v_i_beta;
    std::size_t v_num_beta=0;
    double v_beta;
    std::size_t v_i_walker;
    std::size_t v_walker_id;
    
    std::size_t v_i_par;
    double v_param_value;
    
    bool not_all = true;
    bool from_begining = false;
    
    double v_dur;
    
    while (not_all &&
           load_vars_line(f, iter, v_dur, v_i_beta, v_num_beta, v_beta,
                          v_i_walker, v_walker_id, v_i_par, v_param_value)) {
        
        dur = std::chrono::duration<double>(v_dur);
        if (!from_begining) {
            from_begining = (v_i_beta == 0) && (v_i_walker == 0) && (v_i_par == 0);
        }
        if (from_begining) {
            if ((v_i_beta < data.beta.size())) {
                data.beta[v_i_beta] = v_beta;
                data.i_walkers[v_i_beta][v_i_walker] = v_walker_id;
                data.walkers[v_i_beta][v_i_walker].parameter[v_i_par] = v_param_value;
            } else if (v_i_beta == data.beta.size()) {
                data.beta.push_back(v_beta);
                data.i_walkers.push_back(data.i_walkers[v_i_beta - 1]);
                data.i_walkers[v_i_beta][v_i_walker] = v_walker_id;
                data.walkers.push_back(data.walkers[v_i_beta - 1]);
                data.walkers[v_i_beta][v_i_walker].parameter[v_i_par] = v_param_value;
                
                data.walkers_sta.push_back(data.walkers_sta[v_i_beta - 1]);
                data.emcee_stat.push_back(data.emcee_stat[v_i_beta - 1]);
                data.thermo_stat.push_back(data.thermo_stat[v_i_beta - 2]);
                
                i_beta = v_i_beta;
            }
            if (v_i_par != i_par)
                return false;
            if (v_i_walker != i_walker)
                return false;
            if (v_i_beta != i_beta)
                return false;
            
            if (i_par + 1 < num_Parameters(data)) {
                ++i_par;
            } else {
                i_par = 0;
                if (i_walker + 1 < data.i_walkers[i_beta].size()) {
                    ++i_walker;
                } else {
                    i_walker = 0;
                    if (i_beta + 1 < data.walkers.size()) {
                        ++i_beta;
                    } else {
                        i_beta = 0;
                        not_all = !(v_num_beta == data.walkers.size());
                    }
                }
            }
        }
    }
    
    return v_num_beta == data.walkers.size();
}


template <class Parameters, class Duration>
bool extract_iter(const std::string line, std::size_t &iter, Duration &dur,
                  thermo_fraction_mcmc<Parameters> &data) {
    std::stringstream ss(line);
    return extract_iter(ss, iter, dur, data);
}

template <class Parameters, class Duration>
auto extract_parameters_last(const std::string &fname, std::size_t &iter,
                             Duration &dur, thermo_fraction_mcmc<Parameters> &data) {
    auto candidate = data;
    auto f = std::ifstream(fname);
    std::string line;
    std::getline(f, line);
    auto iter_prev = iter;
    auto dur_prev = dur;
    while (extract_iter(f, iter_prev, dur_prev, candidate)) {
        iter = iter_prev;
        dur = dur_prev;
        std::swap(data, candidate);
    }
    return data;
}





template <class FunctionTable, class Prior, class Likelihood, class Variables,
         class DataType,
         class Parameters = std::decay_t<decltype(sample(
             std::declval<mt_64i &>(), std::declval<Prior &>()))>>
    requires(
        is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)
//    requires (is_prior<Prior,Parameters,Variables,DataType>&&
//    is_likelihood_model<FunctionTable,Likelihood,Parameters,Variables,DataType>)
Maybe_error<bool> calc_thermo_fraction_mcmc_continuation(
    FunctionTable &f, std::size_t n_walkers, 
    ensemble<mt_64i> &mt, Prior const &prior, Likelihood const &lik,
    const DataType &y, const Variables &x, thermo_fraction_mcmc<Parameters> &t_mcmc) {
    
    auto ff = f.fork(omp_get_max_threads());
    
    std::string error;
    bool good = true;
#pragma omp parallel for // collapse(2)
    for (std::size_t iw = 0; iw < n_walkers; ++iw) {
      for (std::size_t i = 0; i < t_mcmc.beta.size(); ++i) {
        auto i_th = omp_get_thread_num();
            auto res = calc_mcmc(ff[i_th], prior, lik, y, x, t_mcmc.walkers[i][iw], t_mcmc.i_fractions[i]);
            if (!res) {
                error += res.error()();
                good = false;
            }
        }
    }
    f += ff;
    if (good)
        return good;
    else
        return error_message(error);
}





template <bool Adapt_beta,class FunctionTable, class Algorithm, class Prior, class Likelihood,
         class Variables, class DataType, class Reporter>
    requires(
        is_of_this_template_type_v<std::decay_t<FunctionTable>, var::FuncMap_St>)

//    requires(is_Algorithm_conditions<Algorithm, thermo_mcmc<Parameters>> &&
//             is_prior<Prior,Parameters,Variables,DataType>&&
//             is_likelihood_model<Likelihood,Parameters,Variables,DataType>)

auto thermo_fraction_evidence_continuation(
    const std::string &idName, FunctionTable &f,
    new_thermodynamic_integration<Algorithm, Reporter> &&therm,
    Prior const &prior, Likelihood const &lik, const std::vector<DataType> &y,
    const std::vector<Variables> &x)
{
    auto a = therm.algorithm();
    auto mt = init_mt(therm.initseed());
    auto n_walkers = therm.num_scouts_per_ensemble();
    auto mts = init_mts(mt, omp_get_max_threads());
    
    by_beta<double> beta;
    
    if constexpr (Adapt_beta)
    {
        beta=by_beta<double>(therm.beta_size(),0);
    }
    else{
        auto beta_ = new_get_beta_list(
            therm.beta_size(), therm.beta_upper_size(), therm.beta_medium_size(),
            therm.beta_upper_value(), therm.beta_medium_value(), therm.stops_at(),
            therm.includes_zero());
        
        auto it_beta_run_begin = beta_.rend() - beta_.size();
        auto it_beta_run_end = beta_.rend();
        beta = by_beta<double>(it_beta_run_begin, it_beta_run_end);
    }
    
    std::vector<std::size_t> samples_size=get_samples_size(y);
    std::size_t max_samples=var::max(samples_size);
    
    
    auto current = create_thermo_fraction_mcmc(n_walkers, beta, mt, prior, samples_size,max_samples);
    auto &rep = therm.reporter();
    
    
    auto fname = idName + "__i_beta__i_walker__i_par.csv";
    std::size_t iter = 0;
    auto start = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double,std::ratio<1>> duration;
    
    
    current=extract_parameters_last(fname, iter, duration, current);
    
    
    
    a.reset(iter);
    auto res=calc_thermo_fraction_mcmc_continuation(f, n_walkers,
                                             mts, prior, lik, y, x,current);
    
    
    auto mcmc_run = checks_convergence(std::move(a), current);
    //using return_type=Maybe_error<decltype(std::pair(std::move(mcmc_run.first), current))>;
    report_title(rep, current, lik, y, x);
    report_title(f, "Iter");
    report_model_all(rep, prior, lik, y, x, beta);
    
    return thermo_fraction_evidence_loop<Adapt_beta>(
        f,
        std::forward<new_thermodynamic_integration<Algorithm, Reporter>>(therm),
        prior, lik, y, x, mcmc_run, iter, current, rep, mt, mts, start,duration);
}


#endif // PARALLEL_TEMPERING_FRACTION_H
