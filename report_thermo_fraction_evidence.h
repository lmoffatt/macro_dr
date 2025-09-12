#ifndef REPORT_THERMO_FRACTION_EVIDENCE_H
#define REPORT_THERMO_FRACTION_EVIDENCE_H

#include "distributions.h"
#include "mcmc.h"
#include "parallel_tempering_fraction.h"
#include "qmodel.h"
#include <bit>
//#include "CLI_thermo_evidence_fraction_dts.h"

template <class Parameters, class... saving>
void report_title(save_mcmc<Parameters,saving...> &f, thermo_fraction_mcmc<Parameters> const &data,
                  ...) {
    using namespace macrodr;
    (report_title(get<saving>(f.m_m), data), ..., 1);
}

template <class FunctionTable, class Duration, class Parameters,
         class... saving, class... T>
void report_all(FunctionTable &f, std::size_t iter, const Duration &dur,
                save_mcmc<Parameters, saving...> &s,
                thermo_fraction_mcmc<Parameters> &data, T const &...ts) {
    (report(f, iter, dur, get<saving>(s.m_m), data, ts...), ..., 1);
}


template <class Parameters>
void report_title(no_save &, thermo_fraction_mcmc<Parameters> const &, ...) {}








template <class Parameterstype, class Parameters>
void report_title(save_likelihood<Parameterstype> &s, thermo_fraction_mcmc<Parameters> const &,
                  ...) {
   

    
    s.f << "iter"
        << s.sep << "iter_time"
        << s.sep<<"num_global_beta"
        << s.sep<< "i_global_beta"
        << s.sep<< "i_frac"
        << s.sep<< "fraction_samples"
        << s.sep << "beta"
        << s.sep << "i_walker" << s.sep
        << "id_walker" << s.sep << "logP"
        << s.sep << "logL" << s.sep << "elogL" << s.sep << "vlogL"
        << s.sep << "logL1" << s.sep << "elogL1" << s.sep << "vlogL1"
        << s.sep << "plog_Evidence" << s.sep
        << "pelog_Evidence" << s.sep << "pvlog_Evidence" << s.sep
        << "log_Evidence" << s.sep << "elog_Evidence" << s.sep
        << "vlog_Evidence"
        << "\n";
}
template <class Parameters>
void report_title(save_Parameter<Parameters> &s, thermo_fraction_mcmc<Parameters> const &,
                  ...) {
    /*                s.f << iter<< s.sep
                    << dur.count()<< s.sep
                    << data.beta.size()<< s.sep
                    << i_beta<< s.sep
                    << data.i_fraction[i_beta]<<s.sep
                    << data.samples_size[data.i_fraction[i_beta]]<<s.sep
                    << data.beta[i_beta]<< s.sep
                    << i_walker << s.sep
                    << data.i_walkers[i_beta][i_walker]<< s.sep;
                for (std::size_t i_par = 0; i_par < num_Parameters(data); ++i_par)
                    s.f<< data.walkers[i_beta][i_walker].parameter[i_par];
                s.f<< "\n";*/
    s.f << "iter"
        << s.sep << "iter_time"
        << s.sep<<"num_global_beta"
        << s.sep<< "i_global_beta"
        << s.sep<< "i_frac"
        << s.sep<< "fraction_samples"
        << s.sep << "i_walker" << s.sep
        << "id_walker" << s.sep << "i_par" << s.sep << "par_value"
        << "\n";
}

template <class Parameters>
void report_title(save_Evidence &s, thermo_fraction_mcmc<Parameters> const &,
                  ...) {
    /*                <<calculate_sample_size(data.walkers_sta,i_beta, i_frac)<< s.sep
                << meanPrior[i_beta ]
                << r_logL.sep(s.sep)
                << r_logL1.sep(s.sep)
                << r_plogE.sep(s.sep)
                << logE.sep(s.sep)
*/
    s.f << "iter"
        << s.sep << "iter_time"
       << s.sep<<"num_global_beta"
        << s.sep<< "i_global_beta"
        << s.sep<< "i_frac"
       << s.sep<< "fraction_samples"
         << s.sep << "beta"
         << s.sep << "S"
        << s.sep<< "sample_size"
        << s.sep << "meanPrior"
        << s.sep << "logL"<< s.sep << "elogL"<< s.sep << "vlogL"
        << s.sep << "logL1"<< s.sep << "elogL1"<< s.sep << "vlogL1"
        << s.sep<< "plog_Evidence" << s.sep << "eplog_Evidence"   << s.sep << "pvarlog_Evidence"
        << s.sep<< "log_Evidence" << s.sep << "elog_Evidence"   << s.sep << "varlog_Evidence"
        << s.sep << "emcee_stat_count"
        << s.sep << "emcee_stat_rate"
        << s.sep << "thermo_jump_stat_count"
        << s.sep << "thermo_jump_rate"
        << s.sep << "deltaBeta_deltalogL"
        << "\n";
}




template <class Parameters>

void report_title(save_RateParameter<Parameters> &s,
                  thermo_fraction_mcmc<Parameters> const &, ...) {
    s.f << "iter"
        << s.sep << "iter_time"
        << s.sep<<"num_global_beta"
        << s.sep<< "i_global_beta"
        << s.sep<< "i_frac"
        << s.sep<< "fraction_samples"
        << s.sep << "beta"<< s.sep << "i_walker" << s.sep
        << "id_walker" << s.sep << "agonist" << s.sep << "i_state_from" << s.sep
        << "i_state_to" << s.sep << "rate_value"
        << "\n";
}


inline void report_title(save_Predictions<var::Parameters_transformed> &s,
                         thermo_fraction_mcmc<var::Parameters_transformed> const &,
                         ...) {

    s.f  << "iter"
        << s.sep << "iter_time"
        << s.sep<<"num_global_beta"
        << s.sep<< "i_global_beta"
        << s.sep<< "i_frac"
        << s.sep<< "fraction_samples"
        << s.sep << "beta" << s.sep << "i_walker" << s.sep
        << "id_walker" << s.sep << "i_step" << s.sep << "time" << s.sep
        << "num_samples" << s.sep << "Y_obs" << s.sep << "Y_pred" << s.sep
        << "Y_var" << s.sep << "plogL" << s.sep << "pelogL"
        << "\n";
    
    s.g  << "iter"
        << s.sep << "iter_time"
        << s.sep<<"num_global_beta"
        << s.sep<< "i_global_beta"
        << s.sep<< "i_frac"
        << s.sep<< "fraction_samples"
        << s.sep << "beta"<< s.sep << "i_walker" << s.sep
        << "id_walker" << s.sep << "i_step" << s.sep << "i_state" << s.sep
        << "j_state" << s.sep << "moment" << s.sep << "value"
        << "\n";
}


template <class FunctionTable, class Parameters, class... T>
void report(FunctionTable &&, std::size_t, no_save &,
            thermo_fraction_mcmc<Parameters> const &, T const &...) {}




template <class FunctionTable, class Duration, class Parameters>
void report(FunctionTable &&, std::size_t iter, const Duration &dur,
            save_likelihood<Parameters> &s, thermo_fraction_mcmc<Parameters>  &data,
            ...) {
    
    std::size_t num_values = 14;
    auto num_beta = size(data.beta);
    std::size_t point_size = std::bit_floor(num_values * num_beta * data.get_Walkers_number());
    std::size_t sampling_interval = std::max(
        s.sampling_interval, point_size / s.max_number_of_values_per_iteration);
    
    if ((iter > 0) && (data.num_samples() > 0) &&
        (iter % sampling_interval == 0)) {
        
        for (std::size_t i_walker = 0; i_walker < num_walkers(data); ++i_walker) {
            
            logLs t_logL={};
            
            logLs log_Evidence = var::Vector_Space<logL, elogL, vlogL>(
                logL(0.0), elogL(0.0), vlogL(0.0));
            for (std::size_t i_beta = 0; i_beta < data.beta.size(); ++i_beta) {
                double beta=data.beta[i_beta];
                auto i_frac1=data.i_fraction[i_beta];
                double beta0 = i_beta==0?0:
(i_frac1==data.i_fraction[i_beta-1]? data.beta[i_beta-1]:0.0);
                auto i_frac0= i_frac1>0?i_frac1-1:0 ;
                auto t_logL1 = data.walkers[i_beta][i_walker].logL(i_frac1);
                auto t_logL0 = i_frac1>0? data.walkers[i_beta][i_walker].logL(i_frac0):
                         logLs_0();
                auto t__logL=t_logL;
                
                
                t_logL= t_logL1-t_logL0;
                
                auto s_logL= i_beta>0? t_logL+t__logL: 2.0*t_logL;
                
                
                auto plog_Evidence = (beta - beta0) / 2.0 * s_logL;
                
                log_Evidence = log_Evidence + plog_Evidence;
                s.f << iter<< s.sep
                    << dur.count()<< s.sep
                    << data.beta.size()<< s.sep
                    << i_beta<< s.sep
                    << i_frac1<<s.sep
                    << data.samples_size[i_frac1]<<s.sep
                    << beta<< s.sep
                    << i_walker << s.sep
                    << data.i_walkers[i_beta][i_walker] << s.sep
                    << data.walkers[i_beta][i_walker].logP
                    << t_logL.sep(s.sep)
                    << t_logL1.sep(s.sep)
                    << plog_Evidence.sep(s.sep)
                    << log_Evidence.sep(s.sep)
                     << "\n";
                
                // now consider the possibility of changing the i_frac
                if ((beta==1)&&i_beta+1<data.beta.size())
                {
                    t_logL = data.walkers[i_beta][i_walker].logL(i_frac1+1)-data.walkers[i_beta][i_walker].logL(i_frac1);
                    beta = 0;
                }
                
            }
        }
    }
}





template <class FunctionTable, class Duration, class Parameters>
 void report(FunctionTable &f, std::size_t iter, const Duration &dur,
                   save_Evidence  &s, thermo_fraction_mcmc<Parameters>  &data, ...) {
    
    std::size_t num_values=32;
    std::size_t point_size=num_values*num_betas(data);
    std::size_t sampling_interval = std::max(
        s.sampling_interval, point_size / s.max_number_of_values_per_iteration);
    if ((iter>0 )&&(data.num_samples()>0)&&(iter % sampling_interval == 0)){
        
        
        
        auto dBdL=calculate_deltaBeta_deltaL(data);
        
        auto meanLik_per_fraction = mean_logL_per_fraction(data);
        auto plogE_sta= calculate_partial_Evidence_sta(data);
        auto meanPrior = mean_logP(data);
        auto running_meanLik =calculate_across_sta(data.walkers_sta, data.i_fractions);
        double beta=0;
        auto logE=logLs{};
        auto r_logL=logLs{};
        for (std::size_t i_beta = 0; i_beta <data.beta.size() ; ++i_beta) {
            auto i_frac=data.i_fraction[i_beta];
            auto beta0=beta;
             beta=data.beta[i_beta];
            
            auto r_logL1=meanLik_per_fraction[i_beta][i_frac];
            
            auto r_logL0=i_frac>0?meanLik_per_fraction[i_beta][i_frac-1]:logLs_0();
            
            auto r__logL=r_logL;
            r_logL=r_logL1-r_logL0;
            
            auto r_plogE= (beta-beta0)/2.0*(r__logL+r_logL);
            logE= logE+r_plogE;            
            
            auto emcee_count=data.emcee_stat[i_beta]().count();
            auto emcee_rate=data.emcee_stat[i_beta]().rate();
            auto thermo_count=data.thermo_stat[std::max(1ul,i_beta)-1ul]().count();
            auto thermo_rate=data.thermo_stat[std::max(1ul,i_beta)-1ul]().rate();
            
            s.f << iter<< s.sep
                << dur.count()<< s.sep
                << data.beta.size()<< s.sep
                << i_beta<< s.sep
                << i_frac<<s.sep
                << data.samples_size[i_frac]<<s.sep
                << beta<< s.sep
                << data.global_beta[i_beta]<<s.sep
                <<calculate_sample_size(data.walkers_sta,i_beta, i_frac)<< s.sep
                << meanPrior[i_beta ]
                << r_logL.sep(s.sep)
                << r_logL1.sep(s.sep)
                << r_plogE.sep(s.sep)
                << logE.sep(s.sep)
                << s.sep<< emcee_count 
                << s.sep<< emcee_rate 
                << s.sep<< thermo_count
                << s.sep<< thermo_rate
                << s.sep<<dBdL[std::max(1ul, i_beta) - 1]
                << "\n";
            if ((beta==1)&&i_beta+1<data.beta.size())
            {
                r_logL = meanLik_per_fraction[i_beta][i_frac+1]-meanLik_per_fraction[i_beta][i_frac];
                beta = 0;
            }
        }
    }
}

template <class FunctionTable, class Duration, class Prior,
         class t_logLikelihood, class Data, class Variables>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, macrodr::FuncMap_St>)
void report(FunctionTable &f, std::size_t iter, const Duration &dur,
            save_Predictions<var::Parameters_transformed> &s,
            thermo_fraction_mcmc<var::Parameters_transformed> const &data, Prior const &,
            t_logLikelihood const &lik, const Data &y, const Variables &x,
            ...) {
    auto num_samples = size(y);
    auto num_states = lik.m.number_of_states();
    auto num_states_a = std::pow(2, std::round(std::log2(num_states)));
    std::size_t num_values = 4;
    std::size_t num_beta_portions = 2;
    
    std::size_t num_samples_a = std::pow(2, std::round(std::log2(size(y))));
    std::size_t point_size = std::bit_floor(num_values * num_beta_portions *
                             data.get_Walkers_number() * num_samples_a);
    std::size_t sampling_interval = std::max(
        s.sampling_interval, point_size / s.max_number_of_values_per_iteration);
    
    std::size_t state_num_values = 1;
    std::size_t num_moments = 2;
    std::size_t state_point_size = state_num_values * num_beta_portions *
                                   data.get_Walkers_number() * num_samples_a *
                                   num_states_a * num_moments;
    std::size_t state_sampling_interval =
        std::max(s.sampling_interval,
                 state_point_size / s.max_number_of_values_per_iteration);
    
    if ((iter == 0) || (iter % sampling_interval != 0))
        return;
    auto ff = f.fork(omp_get_max_threads());
    
    auto all_Predictions =
        std::vector<std::vector<std::decay_t<decltype(logLikelihoodPredictions(
            ff[0], lik, data.get_Parameter(0, 0), y, x))>>>(
            data.get_Walkers_number());
    
    auto beta = data.get_Beta();
    auto num_beta = beta.size();
    
#pragma omp parallel for
    for (std::size_t i_walker = 0; i_walker < data.get_Walkers_number();
         ++i_walker) {
        for (std::size_t i_b = 0; i_b < beta.size(); ++i_b) {
            
            if ((beta[i_b] == 1) || (iter % (beta.size() * sampling_interval) == 0)) {
                auto i_th = omp_get_thread_num();
                
                auto par = data.get_Parameter(i_walker, i_b);
                auto walker_id = data.get_Walker(i_walker, i_b);
                auto i_frac= data.i_fraction[i_b];
                all_Predictions[i_walker].push_back(
                    logLikelihoodPredictions(ff[i_th], lik, par.to_value(), y[i_frac], x[i_frac]));
            }
        }
    }
    f += ff;
    for (std::size_t half = 0; half < 2; ++half) {
        for (std::size_t iiw = 0; iiw < data.get_Walkers_number() / 2; ++iiw) {
            auto i_walker = half ? iiw + data.get_Walkers_number() / 2 : iiw;
            for (std::size_t i_b = 0; i_b < beta.size(); ++i_b) {
                if ((beta[i_b] == 1) ||
                    (iter % (beta.size() * sampling_interval) == 0)) {
                    auto par = data.get_Parameter(i_walker, i_b);
                    auto walker_id = data.get_Walker(i_walker, i_b);
                    auto i_frac= data.i_fraction[i_b];
                    auto prediction = (iter % (beta.size() * sampling_interval) == 0)
                                          ? all_Predictions[i_walker][i_b]
                                          : all_Predictions[i_walker][0];
                    if (is_valid(prediction)) {
                        auto &predictions = prediction.value();
                        for (std::size_t i_step = 0; i_step < size(y); ++i_step) {
                            auto v_ev =
                                get<macrodr::ATP_evolution>(get<macrodr::Recording_conditions>(x[i_frac])()[i_step]);
                            
                            auto time = get<macrodr::Time>(get<macrodr::Recording_conditions>(x[i_frac])()[i_step]);
                            auto num_smples = get_num_samples(v_ev);
                            
                            s.f
                                << iter<< s.sep
                                << dur.count()<< s.sep
                                << data.beta.size()<< s.sep
                                << i_b<< s.sep
                                << i_frac<<s.sep
                                << data.samples_size[i_frac]<<s.sep
                                << beta<< s.sep
                                << i_walker<< s.sep
                                << walker_id << s.sep
                                << i_step << s.sep
                                << time << s.sep
                                << num_samples << s.sep
                                << y[i_frac]()[i_step]() << s.sep
                                << get<macrodr::y_mean>(predictions()[i_step])
                                << s.sep
                                << get<macrodr::y_var>(predictions()[i_step])
                                << s.sep
                                << get<macrodr::plogL>(predictions()[i_step])
                                << s.sep
                                << get<macrodr::eplogL>(predictions()[i_step])
                                << "\n";
                            
                            if (((iter % state_sampling_interval == 0) && (beta[i_b] == 1)) ||
                                (iter % (state_sampling_interval * num_beta) == 0)) {
                                auto &v_P = get<macrodr::P_mean>(predictions()[i_step]);
                                auto &v_Pc = get<macrodr::P_Cov>(predictions()[i_step]);
                                auto n = v_P().size();
                                for (std::size_t i_state = 0; i_state < n; ++i_state) {
                                    s.g
                                        
                                        << iter<< s.sep
                                        << dur.count()<< s.sep
                                        << data.beta.size()<< s.sep
                                        << i_b<< s.sep
                                        << i_frac<<s.sep
                                        << data.samples_size[i_frac]<<s.sep
                                        << beta<< s.sep
                                        << i_walker << s.sep
                                        << walker_id<< s.sep
                                        << i_step << s.sep
                                        << i_state << s.sep
                                        << 0 << s.sep
                                        << "mean" << s.sep
                                        << v_P()[i_state] << "\n";
                                    if ((iter % (state_sampling_interval *
                                                 num_states) == 0) &&
                                            ((beta[i_b] == 1)) ||
                                        (iter %                                                                         (state_sampling_interval* num_states*num_beta)==0))
                                        for (std::size_t j_state = 0; j_state <=i_state;
                                             ++j_state) {
                                            s.g
                                                << iter<< s.sep
                                                << dur.count()<< s.sep
                                                << data.beta.size()<< s.sep
                                                << i_b<< s.sep
                                                << i_frac<<s.sep
                                                << data.samples_size[i_frac]<<s.sep
                                                << beta<< s.sep
                                                << i_walker << s.sep
                                                << walker_id << s.sep
                                                << i_step << s.sep
                                                << i_state << s.sep
                                                << j_state << s.sep
                                                << "Cov" << s.sep
                                                << v_Pc()(i_state, j_state) << "\n";
                                        }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template <class FunctionTable, class Duration, class Parameters>
void report(FunctionTable &f, std::size_t iter, const Duration &dur,
            save_Parameter<Parameters> &s, thermo_fraction_mcmc<Parameters> const &data,
            ...) {
    std::size_t num_values = 4;
    auto num_beta = size(data.beta);
    std::size_t point_size = std::bit_floor(num_values * num_beta * data.get_Walkers_number() *
                             num_Parameters(data));
    std::size_t sampling_interval = std::max(
        s.sampling_interval, point_size / s.max_number_of_values_per_iteration);
    
    if ((iter > 0) && (data.num_samples() > 0) &&
        (iter % sampling_interval == 0)) {
        for (std::size_t i_beta = 0; i_beta < num_betas(data); ++i_beta)
            for (std::size_t i_walker = 0; i_walker < num_walkers(data); ++i_walker){
                for (std::size_t i_par = 0; i_par < num_Parameters(data); ++i_par)
                s.f << iter<< s.sep
                    << dur.count()<< s.sep
                    << data.beta.size()<< s.sep
                    << i_beta<< s.sep
                    << data.i_fraction[i_beta]<<s.sep
                    << data.samples_size[data.i_fraction[i_beta]]<<s.sep
                    << data.beta[i_beta]<< s.sep
                    << i_walker << s.sep
                        << data.i_walkers[i_beta][i_walker]
                        << s.sep<<i_par
                        <<s.sep<< data.walkers[i_beta][i_walker].parameter[i_par]
                << "\n";
            }
        s.f.flush();
    }
}



template <class FunctionTable, class Duration, class Prior,
         class t_logLikelihood>
    requires(is_of_this_template_type_v<std::decay_t<FunctionTable>, macrodr::FuncMap_St>)
void report(FunctionTable &, std::size_t iter, const Duration &dur,
            save_RateParameter<var::Parameters_transformed> &s,
            thermo_fraction_mcmc<var::Parameters_transformed>  &data, Prior const &,
            t_logLikelihood const &lik, ...) {
    
    auto num_states = lik.m.number_of_states();
    std::size_t num_values = 1;
    std::size_t num_beta_portions = 2;
    std::size_t point_size = std::bit_floor(num_values * num_beta_portions *
                             data.get_Walkers_number() * num_states * num_states);
    std::size_t sampling_interval = std::max(
        s.sampling_interval, point_size / s.max_number_of_values_per_iteration);
    
    if ((iter == 0) || (iter % sampling_interval != 0))
        return;
    auto &model = lik.m;
    auto beta = data.get_Beta();
    
    for (std::size_t i_walker = 0; i_walker < data.get_Walkers_number();
         ++i_walker) {
        for (std::size_t i_b = 0; i_b < beta.size(); ++i_b) {
            if (beta[i_b] == 1) {
                auto par = data.get_Parameter(i_walker, i_b);
                auto walker_id = data.get_Walker(i_walker, i_b);
                auto Maybe_mo = model(par.to_value());
                if (is_valid(Maybe_mo)) {
                    auto &mo = Maybe_mo.value();
                    auto v_Q0 = get<macrodr::Q0>(mo);
                    auto v_Qa = get<macrodr::Qa>(mo);
                    
                    for (std::size_t i_from = 0; i_from < v_Q0().nrows(); ++i_from)
                        for (std::size_t i_to = 0; i_to < v_Q0().ncols(); ++i_to) {
                            if (v_Qa()(i_from, i_to) > 0)
                                s.f << iter<< s.sep
                                    << dur.count()<< s.sep
                                    << data.beta.size()<< s.sep
                                    << i_b<< s.sep
                                    << data.i_fraction[i_b]<<s.sep
                                    << data.samples_size[data.i_fraction[i_b]]<<s.sep
                                    << data.beta[i_b]<< s.sep
                                    << i_walker<< s.sep
                                    << walker_id << s.sep
                                    << "agonist" << s.sep
                                    << i_from << s.sep
                                    << i_to << s.sep
                                    << v_Qa()(i_from, i_to)
                                    << "\n";
                            if (v_Q0()(i_from, i_to) > 0)
                                s.f << iter<< s.sep
                                    << dur.count()<< s.sep
                                    << data.beta.size()<< s.sep
                                    << i_b<< s.sep
                                    << data.i_fraction[i_b]<<s.sep
                                    << data.samples_size[data.i_fraction[i_b]]<<s.sep
                                    << data.beta[i_b]<< s.sep
                                    << i_walker<< s.sep
                                    << walker_id << s.sep
                                    << "no_agonist" << s.sep
                                    << i_from << s.sep
                                    << i_to << s.sep
                                    << v_Q0()(i_from, i_to)
                                    << "\n";
                        }
                }
            }
        }
    }
}


#endif // REPORT_THERMO_FRACTION_EVIDENCE_H
