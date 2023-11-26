#ifndef MICROR_STOCHASTIC_H
#define MICROR_STOCHASTIC_H


#include "multivariate_normal_distribution.h"
#include "qmodel.h"
#include <random>

namespace macrodr {


class log_Pi : public var::Var<log_Pi, Matrix<double>> {};

class log_Pij : public var::Var<log_Pij, Matrix<double>> {};


class N_channel_transition : public var::Var<N_channel_transition, Matrix<double>> {};


using log_Patch_Transition_State = Vector_Space<log_Pi,log_Pij>;



template <uses_recursive_aproximation recursive,
         uses_averaging_aproximation averaging,
         uses_variance_aproximation variance,
         class FunctionTable,
         class C_Patch_State,
         class C_Qdt, class C_Patch_Model, class C_double>
    requires(
        /*(U<std::decay_t<C_Patch_State>,
           Patch_State>||U<std::decay_t<C_Patch_State>,
           Patch_State_and_Evolution> )&& U<C_Patch_Model, Patch_Model> &&*/
        U<C_double, double> &&
        U<C_Qdt, Qdt>)
auto sample_Patch_Transition_State(std::mt19937_64& mt,C_Patch_State &&t_prior, C_Qdt const &t_Qdt, C_Patch_Model const &m,
                                   C_double const &Nch, const Patch_current &p_y, double fs)
{
    using Transf = transformation_type_t<C_Qdt>;
    auto &p_P_cov = get<P_Cov>(t_prior);
    auto &p_P_mean = get<P_mean>(t_prior);
    
    
    auto P_dist=make_multivariate_normal_distribution(p_P_mean(),p_P_cov());
    
    auto s_Pi=P_dist.value().sample(mt);
    
    //    auto &y = get<Patch_current>(p).value();
    auto &y = p_y.value();
    
    auto &t_tolerance = get<Probability_error_tolerance>(m);
    auto &t_min_P = get<min_P>(m);
    auto e = get<Current_Noise>(m).value() * fs /
             get<number_of_samples>(t_Qdt).value();
    auto y_baseline = get<Current_Baseline>(m);
    auto N = Nch;
    
}


template <uses_recursive_aproximation recursive,
         uses_averaging_aproximation averaging,
         uses_variance_aproximation variance,
         class FunctionTable,
         class C_Patch_State,
         class C_Qdt, class C_Patch_Model, class C_double>
    requires(
        /*(U<std::decay_t<C_Patch_State>,
           Patch_State>||U<std::decay_t<C_Patch_State>,
           Patch_State_and_Evolution> )&& U<C_Patch_Model, Patch_Model> &&*/
        U<C_double, double> &&
        U<C_Qdt, Qdt>)

Maybe_error<C_Patch_State>
Micror_stochastic(FunctionTable&,C_Patch_State &&t_prior, C_Qdt const &t_Qdt, C_Patch_Model const &m,
                  C_double const &Nch, const Patch_current &p_y, double fs)  {
    
    using Transf = transformation_type_t<C_Qdt>;
    auto &p_P_cov = get<P_Cov>(t_prior);
    auto &p_P_mean = get<P_mean>(t_prior);
    //    auto &y = get<Patch_current>(p).value();
    auto &y = p_y.value();
    
    auto &t_tolerance = get<Probability_error_tolerance>(m);
    auto &t_min_P = get<min_P>(m);
    auto e = get<Current_Noise>(m).value() * fs /
             get<number_of_samples>(t_Qdt).value();
    auto y_baseline = get<Current_Baseline>(m);
    auto N = Nch;
    
    
    
    
    
    
    
    
    
    Matrix<double> u(p_P_mean().size(), 1, 1.0);
    
    auto SmD = p_P_cov() - diag(p_P_mean());
    
    auto N_states = p_P_mean().ncols();
    
    auto &t_gmean_i = get<gmean_i>(t_Qdt);
    auto &t_gtotal_ij = get<gtotal_ij>(t_Qdt);
    auto &t_gtotal_var_ij = get<gtotal_var_ij>(t_Qdt);
    auto &t_gmean_ij = get<gmean_ij>(t_Qdt);
    auto &t_gvar_i = get<gvar_i>(t_Qdt);
    auto gSg =
        getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
        getvalue(p_P_mean() * (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));
    
    if constexpr (false) {
        auto test_gSg = var::test_Derivative(
            [this, &u](auto const &t_gmean_i, auto const &t_gmean_ij,
                       auto const &SmD, auto const &p_P_mean,
                       auto const &t_gtotal_ij) {
                return getvalue(TranspMult(t_gmean_i(), SmD) * t_gmean_i()) +
                       getvalue(p_P_mean() *
                                (elemMult(t_gtotal_ij(), t_gmean_ij()) * u));
            },
            1e-4, 1e-6, t_gmean_i, t_gmean_ij, SmD, p_P_mean, t_gtotal_ij);
        if (!test_gSg) {
            std::cerr << "\n error in test_gSg!!\n" << test_gSg.error()();
            return Maybe_error<C_Patch_State>(test_gSg.error());
        }
    }
    auto ms = getvalue(p_P_mean() * t_gvar_i());
    
    Op_t<Transf, double> e_mu;
    Op_t<Transf, y_mean> r_y_mean;
    Op_t<Transf, y_var> r_y_var;
    
    Op_t<Transf, double> sSg;
    Op_t<Transf, double> sSs;
    Op_t<Transf, double> zeta;
    auto t_P = get<P>(t_Qdt);
    
    if constexpr ((!variance.value) && (!recursive.value)) {
        e_mu = e + max(0.0, N * ms);
        r_y_mean() = N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline();
        r_y_var() = e + max(0.0, N * ms + N * gSg);
        if (!(primitive(r_y_var()) > 0.0)) {
            std::stringstream ss;
            ss << "Negative variance!!\n";
            ss << "\nr_y_var=\t" << r_y_var;
            ss << "\ne_mu=\t" << e_mu;
            ss << "\ne=\t" << e;
            ss << "\nN=\t" << N;
            ss << "\n"
               << "ms"
               << "=\t" << ms;
            ss << "\ngSg=\t" << gSg;
            ss << "\n"
               << "ms"
               << "=\t" << ms;
            ss << "\n"
               << "p_P_mean()"
               << "=\t" << p_P_mean();
            ss << "\n"
               << "t_gvar_i()"
               << "=\t" << t_gvar_i();
            
            return error_message(ss.str());
        }
        
    } else if constexpr (!variance.value && recursive.value) {
        auto gS =
            TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
        
        auto ms = getvalue(p_P_mean() * t_gvar_i());
        
        e_mu = e + max(N * ms, 0.0);
        r_y_mean() = N * getvalue(p_P_mean() * t_gmean_i()) + y_baseline();
        r_y_var() = e + max(0.0, N * ms + N * gSg);
        if (!(primitive(r_y_var()) > 0)) {
            std::stringstream ss;
            ss << "Negative variance!!\n";
            ss << "\nr_y_var=\t" << r_y_var;
            ss << "\ne_mu=\t" << e_mu;
            ss << "\ne=\t" << e;
            ss << "\nN=\t" << N;
            ss << "\n"
               << "ms"
               << "=\t" << ms;
            ss << "\ngSg=\t" << gSg;
            ss << "\n"
               << "ms"
               << "=\t" << ms;
            ss << "\n"
               << "p_P_mean()"
               << "=\t" << p_P_mean();
            ss << "\n"
               << "t_gvar_i()"
               << "=\t" << t_gvar_i();
            
            return error_message(ss.str());
        }
        
    } else // (variance && (recursive || !recursive))
    {
        auto &t_gtotal_var_ij = get<gtotal_var_ij>(t_Qdt);
        auto &t_gvar_ij = get<gvar_ij>(t_Qdt);
        
        sSg = getvalue(TranspMult(t_gvar_i(), SmD) * t_gmean_i()) +
              getvalue(p_P_mean() *
                       (elemMult(t_gtotal_var_ij(), t_gmean_ij()) * u));
        sSs =
            getvalue(TranspMult(t_gvar_i(), SmD) * t_gvar_i()) +
            getvalue(p_P_mean() * (elemMult(t_gtotal_var_ij(), t_gvar_ij()) * u));
        
        auto delta_emu = var::max(sqr(ms + e / N) - 2.0 / N * sSs, 0.0);
        auto ms0 = (ms - e / N) / 2 + std::sqrt(delta_emu) / 2;
        
        e_mu = e + N * ms0;
        r_y_mean() = N * getvalue(p_P_mean() * t_gmean_i()) -
                     N * 0.5 / e_mu * sSg + y_baseline();
        zeta = N / (2 * sqr(e_mu) + N * sSs);
        r_y_var() = var::max(e, e + N * ms0 + N * gSg - N * zeta * sqr(sSg));
        if (!(primitive(r_y_var()) > 0)) {
            std::stringstream ss;
            ss << "Negative variance!!\n";
            ss << "\nr_y_var=\t" << r_y_var;
            ss << "\ne_mu=\t" << e_mu;
            ss << "\ne=\t" << e;
            ss << "\nN=\t" << N;
            ss << "\n"
               << "ms"
               << "=\t" << ms;
            ss << "\ngSg=\t" << gSg;
            ss << "\n"
               << "ms"
               << "=\t" << ms;
            ss << "\n"
               << "p_P_mean()"
               << "=\t" << p_P_mean();
            ss << "\n"
               << "t_gvar_i()"
               << "=\t" << t_gvar_i();
            
            return error_message(ss.str());
        }
    }
    if (std::isnan(y)) {
        
        auto r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD));
        auto r_P_mean = build<P_mean>(p_P_mean() * t_P());
        r_P_cov() = r_P_cov() + diag(r_P_mean());
        // std::cerr<<"\nPcov nana corr\n"<<P__cov<<"\nP_mean nana
        // corr\n"<<P_mean<<"\nQ.P \n"<<Q_dt.P();
        auto r_test = test<true>(r_P_mean, r_P_cov, t_tolerance);
        if constexpr (true || r_test)
            if constexpr (U<C_Patch_State, Patch_State>)
                return Op_t<Transf, Patch_State>(
                    Op_t<Transf, logL>(get<logL>(t_prior)()),
                    Op_t<Transf, elogL>(get<elogL>(t_prior)()),
                    Op_t<Transf, vlogL>(get<vlogL>(t_prior)()),
                    normalize(std::move(r_P_mean), t_min_P()),
                    normalize(std::move(r_P_cov), t_min_P()), std::move(r_y_mean),
                    std::move(r_y_var), plogL(NaN), eplogL(NaN), vplogL(NaN));
            else {
                auto &ev = get<Patch_State_Evolution>(t_prior);
                r_P_mean = normalize(std::move(r_P_mean), t_min_P());
                r_P_cov = normalize(std::move(r_P_cov), t_min_P());
                ev().push_back(Op_t<Transf, Patch_State>(
                    Op_t<Transf, logL>(get<logL>(t_prior)()),
                    Op_t<Transf, elogL>(get<elogL>(t_prior)()),
                    Op_t<Transf, vlogL>(get<vlogL>(t_prior)()), r_P_mean, r_P_cov,
                    r_y_mean, r_y_var, plogL(NaN), eplogL(NaN), vplogL(NaN)));
                return Op_t<Transf, Patch_State_and_Evolution>(
                    Op_t<Transf, logL>(get<logL>(t_prior)()),
                    Op_t<Transf, elogL>(get<elogL>(t_prior)()),
                    Op_t<Transf, vlogL>(get<vlogL>(t_prior)()), std::move(r_P_mean),
                    std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var),
                    plogL(NaN), eplogL(NaN), vplogL(NaN), std::move(ev));
            }
        else
            return error_message("fails at intertrace prediction!!: " +
                                 r_test.error()());
    }
    
    auto dy = y - r_y_mean();
    auto chi = dy / r_y_var();
    Op_t<Transf, P_mean> r_P_mean;
    Op_t<Transf, P_Cov> r_P_cov;
    if constexpr (!recursive.value) {
        r_P_cov = build<P_Cov>(AT_B_A(t_P(), SmD));
        
        r_P_mean = build<P_mean>(p_P_mean() * t_P());
        r_P_cov() = r_P_cov() + diag(r_P_mean());
    } else if constexpr (!variance.value) {
        auto gS =
            TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
        auto gseg = chi * gS;
        
        r_P_mean() = p_P_mean() * t_P() + chi * gS;
        
        r_P_cov() = AT_B_A(t_P(), SmD) + diag(p_P_mean() * t_P()) -
                    (N / r_y_var()) * XTX(gS);
        
    } else {
        auto gS =
            TranspMult(t_gmean_i(), SmD) * t_P() + p_P_mean() * t_gtotal_ij();
        auto sS =
            TranspMult(t_gvar_i(), SmD) * t_P() + p_P_mean() * t_gtotal_var_ij();
        r_P_mean() =
            p_P_mean() * t_P() + chi * gS - (chi * zeta * sSg + 0.5 / e_mu) * sS;
        
        r_P_cov() =
            AT_B_A(t_P(), SmD) + diag(r_P_mean() * t_P()) -
            (zeta + N / r_y_var() * sqr(zeta * sSg)) * XTX(sS) +
            (2.0 * N / r_y_var() * zeta * sSg) * X_plus_XT(TranspMult(sS, gS)) -
            (N / r_y_var()) * XTX(gS);
    }
    
    auto chi2 = dy * chi;
    
    Op_t<Transf, plogL> r_plogL;
    if (primitive(r_y_var()) > 0.0)
        r_plogL() = -0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5 * chi2;
    else {
        std::stringstream ss;
        ss << "Negative variance!!\n";
        ss << "\nr_y_var=\t" << r_y_var;
        ss << "\ne_mu=\t" << e_mu;
        ss << "\ngSg=\t" << gSg;
        return error_message(ss.str());
    }
    
    if constexpr (false) {
        auto test_plogL = var::test_Derivative(
            [](auto const &r_y_var, auto const &chi2, auto const &t_prior) {
                return -0.5 * log(2 * std::numbers::pi * r_y_var()) - 0.5 * chi2 +
                       get<logL>(t_prior)();
            },
            1e-6, 1e-8, r_y_var, chi2, t_prior);
        if (!test_plogL) {
            std::cerr << "\n error in test_plogL!!\n" << test_plogL.error()();
            return Maybe_error<C_Patch_State>(test_plogL.error());
        }
    }
    
    //    std::cerr<<p<<"\n";
    //    std::cerr<<"r_plogL\n"<<r_plogL<<"\n";
    //    std::cerr<<"r_P_mean\n"<<r_P_mean<<"\n";
    //    std::cerr<<"r_y_mean\n"<<r_y_mean<<"\n";
    //    std::cerr<<"r_y_var\n"<<r_y_var<<"\n";
    
    //    if (get<Time>(p)()>1)
    //      std::abort();
    
    Op_t<Transf, eplogL> r_eplogL(-0.5 * log(2 * std::numbers::pi * r_y_var()) -
                                  0.5); // e_mu+N*gSg"-N*zeta*sqr(sSg)"
    // double chilogL=(eplogL-plogL)/std::sqrt(0.5);
    
    vplogL r_vlogL(0.5);
    auto r_test = test<true>(r_P_mean, r_P_cov, t_tolerance);
    if constexpr (false) {
        if (!r_test) {
            std::stringstream ss;
            
            ss << "\nP_mean \n" << r_P_mean;
            ss << "\nPcov \n" << r_P_cov;
            // ss<<"\nprior=\n"<<prior<<"\nQ_dt \n"<<Q_dt;
            
            return error_message("\nfails in trace!!!; error=" + r_test.error()() +
                                 ss.str());
        }
    } else if (std::isnan(primitive(r_plogL())))
        return error_message("likelihood is nan");
    else if constexpr (U<C_Patch_State, Patch_State>)
        return build<Patch_State>(
            build<logL>(get<logL>(t_prior)() + r_plogL()),
            build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
            build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()),
            normalize(std::move(r_P_mean), t_min_P()),
            normalize(std::move(r_P_cov), t_min_P()), std::move(r_y_mean),
            std::move(r_y_var), r_plogL, r_eplogL, r_vlogL);
    else {
        auto &ev = get<Patch_State_Evolution>(t_prior);
        r_P_mean = normalize(std::move(r_P_mean), t_min_P());
        r_P_cov = normalize(std::move(r_P_cov), t_min_P());
        ev().push_back(build<Patch_State>(
            build<logL>(get<logL>(t_prior)() + r_plogL()),
            build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
            build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), r_P_mean, r_P_cov,
            r_y_mean, r_y_var, r_plogL, r_eplogL, r_vlogL));
        return build<Patch_State_and_Evolution>(
            build<logL>(get<logL>(t_prior)() + r_plogL()),
            build<elogL>(get<elogL>(t_prior)() + r_eplogL()),
            build<vlogL>(get<vlogL>(t_prior)() + r_vlogL()), std::move(r_P_mean),
            std::move(r_P_cov), std::move(r_y_mean), std::move(r_y_var), r_plogL,
            r_eplogL, r_vlogL, std::move(ev));
    }
}



} // namespace macrodr


#endif // MICROR_STOCHASTIC_H
