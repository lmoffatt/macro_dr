#include "CLI_macro_dr.h"
#include "experiment.h"
#include "lapack_headers.h"
#include "lexer_typed.h"
#include "matrix.h"
#include "variables_derivative.h"
#include "parameters_derivative.h"
//#include "multivariate_normal_distribution.h"
#include "parallel_tempering.h"
#include "cuevi.h"
#include <fstream>
#include <iostream>
#include <map>


#include "qmodel.h"
using namespace macrodr;



int main(int argc, char **argv) {
    
    constexpr bool test_Bayesian_Linear_Regression=true;
    if constexpr (test_Bayesian_Linear_Regression)
    {
        /**
   * @brief myseed defines the random number seed so all runs are identical for debugging purposes
   */
        auto myseed = 9762841416869310605ul;
        //  auto myseed = 0;
        
        /**
   * @brief npar number of parameters of the linear model
   */
        auto npar = 40ul;
        
        /**
   * @brief nsamples number of samples of the linear model
   */
        auto nsamples = 5000ul;
        
        /**
   * @brief std_log10_std_par  standard deviation of the log10 of the standard deviation of the parameters values
   */
        auto std_log10_std_par = 1.0;
        
        /**
   * @brief mean_mean_par mean of the mean of the parameters values
   */
        auto mean_mean_par = 0.0;
        
        /**
   * @brief std_mean_par standard deviation of the mean of the parameters values
   */
        auto std_mean_par = 10.0;
        
        
        /**
   * @brief mean_b mean of the linear coefficient of all parameters
   */
        auto mean_b = 0.0;
        
        /**
   * @brief std_b standard deviation of the paramters linear coefficient
   */
        auto std_b = 10.0;
        
        /**
   * @brief prior_error_factor error factor between the prior and real parameter coefficient covariance value
   */
        auto prior_error_factor = 1;
        
        /**
   * @brief prior_eps_df prior value for the degrees of freedom used for estimating the value of the variance of the data point error
   */
        double prior_eps_df = 1.0;
        
        /**
   * @brief prior_eps_variance prior value of the variance of the data point error
   */
        double prior_eps_variance = 1.0;
        
        myseed = calc_seed(myseed);
        std::cerr<<"myseed ="<<myseed;
        
        /**
   * @brief mt random engine
   */
        auto mt = init_mt(myseed);
        
        /**
   * @brief a_0 parameter alpha of the beta distribution for the prior of the data point error variance
   */
        auto a_0 = prior_eps_df / 2.0;
        
        /**
   * @brief b_0 parameter beta of the beta distribution for the prior of the data point error variance
   */
        auto b_0 = prior_eps_df * prior_eps_variance / 2.0;
        
        
        /**
   * @brief num_scouts_per_ensemble number of scouts per ensemble in the affine ensemble mcmc model
   */
        std::size_t num_scouts_per_ensemble = 64;
        
        /**
   * @brief n_points_per_decade number of points per 10 times increment in beta thermodynamic parameter
   */
        double n_points_per_decade = 3;
        
        /**
   * @brief n_points_per_decade_fraction number of points per 10 times increment in the number of samples
   */
        double n_points_per_decade_fraction = 3;
        
        /**
   * @brief stops_at minimum value of beta greater than zero
   */
        double stops_at = 1e-7;
        
        /**
   * @brief includes_zero considers also beta equal zero
   */
        bool includes_zero = true;
        
        
        /**
   * @brief max_iter maximum number of iterations
   */
        std::size_t max_iter = 10000;
        
        /**
   * @brief path directory for the output
   */
        std::string path = "";
        
        /**
   * @brief filename prefix filename for the output
   */
        std::string filename = "A";
        
        
        /**
   * @brief checks_derivative_every_model_size number of steps before every check of the derivative against the beta thermo parameter for stopping
   */
        std::size_t checks_derivative_every_model_size = 5000;
        
        /**
   * @brief max_ratio maximimum tolerated ratio for the beta derivative method
   */
        double max_ratio = 8;
        
        
        /**
   * @brief min_fraction fraction of the prior parameter size used as the minimal sample used for the cumulative sequence
   */
        double min_fraction = 2;
        
        /**
   * @brief my_linear_model bayesian linear model a bayesian linear model with all the prior information
   */
        auto my_linear_model =
            make_bayesian_linear_model(prior_eps_df, prior_eps_variance, npar, mean_b,
                                       std_b, prior_error_factor)
                .value();
        
        /**
   * @brief X random generated independent variables
   * @brief mean_par random generated mean of the parameters linear coefficients
   * @brief cov_par random generated covariance of the parameters linear coefficients
   */
        auto [X, mean_par, cov_par] = independent_variable_model(
            npar, std_log10_std_par, mean_mean_par, std_mean_par)(mt, nsamples);
        //  auto y = X * tr(b) + eps;
        
        /**
   * @brief b random generated linear coefficient values
   */
        auto b = sample(mt, my_linear_model);
        
        /**
   * @brief y simulated data using the constructed linear model, the independent variables and the linear coefficient values
   */
        auto y = simulate(mt, my_linear_model, b, X);
        
        
        /**
   * @brief beta values for the thermodynamic parameter beta (that ranges from 0 -only prior- to 1 -posterior likelihood
   */
        auto beta = get_beta_list(n_points_per_decade, stops_at, includes_zero);
        
        
        /**
   * @brief thermo_jumps_every factor that multiplied by the model size it produces the number of steps skipped until the next thermo jump
   */
        std::size_t thermo_jumps_every = my_linear_model.size() * 1e0;
        
        if (false) {
            /**
     * @brief tmi classical thermodynamic algorithm ends by maximum iteration
     */
            auto tmi = thermo_by_max_iter<Matrix<double>>(
                path, "Iteri", num_scouts_per_ensemble, thermo_jumps_every, max_iter,
                n_points_per_decade, stops_at, includes_zero, myseed);
            auto opt = evidence(std::move(tmi), my_linear_model.prior(),
                                my_linear_model.likelihood(), y, X);
        }
        if (false) {
            /**
     * @brief tbc classical thermodynamic algorithm, ends using convergence criteria
     */
            auto tbc = thermo_by_convergence<Matrix<double>>(
                path, "exp_thermo", num_scouts_per_ensemble, thermo_jumps_every,
                checks_derivative_every_model_size, n_points_per_decade, stops_at,
                includes_zero, myseed);
            
            
            auto opt2 = evidence(std::move(tbc), my_linear_model.prior(),
                                 my_linear_model.likelihood(), y, X);
            
            // std::cout<<y;
        }
        if (true) {
            
            /**
     * @brief cbc cumulative evidence algorithm, ends using convergence criteria
     */
            auto cbc = cuevi_by_convergence<Matrix<double>>(
                path, "exp_cuevi_40", num_scouts_per_ensemble, min_fraction,
                thermo_jumps_every, checks_derivative_every_model_size, max_ratio,
                n_points_per_decade, n_points_per_decade_fraction, stops_at,
                includes_zero, myseed);
            auto opt3 = evidence(std::move(cbc), my_linear_model.prior(),
                                 my_linear_model.likelihood(), y, X);
        }
    }
    
    
    
    
    
    //        auto filename=argv[1];
    //  auto filename = "../macro_dr/test.txt";
    
    auto cm = dcli::Compiler{};
    
    cm.push_function("load_experiment",
                     dcli::to_typed_function<std::string, double, double>(
                         &macrodr::load_experiment, "filename",
                         "frequency_of_sampling", "initial_ATP"));
    
    
    auto filename = "../macro_dr/run_script.txt";
    std::ifstream f(filename);
    
    std::string s;
    while (f) {
        std::string line;
        std::getline(f, line);
        s += line + "\n";
    }
    std::cout << "\ntest file \n" << s << "\n";
    auto p = dcli::extract_program(s);
    
    std::cerr<<p;
    
    if (p)
    { 
        auto c=dcli::compile_program(cm,p.value());
        if (c)
        {
            auto exec=c.value().run();
        }
    }
    
    auto Efilename = "../macro_dr/Moffatt_Hume_2007_ATP_2.txt";
    
    auto recording = macrodr::load_recording(Efilename);
    
    // std::cerr<<recording.value();
    // std::cerr<<std::numeric_limits<double>::quiet_NaN();
    //  aram_0 = State_Model_Parameters ( values = { { "kon"  50 }  { "koff"   20
    //  }  { "beta"  500 }  { "alfa"  400 } { "g"  16.59e-3 }
    //  {"Number_of_Channels"  100}  { "gaussian_noise"  1.0e-5 } } )
    
    auto experiment =
        Experiment(std::move(recording), Frequency_of_Sampling(50e3),
                   initial_ATP_concentration(ATP_concentration(0.0)));
    
    
    auto model1= Model1::Model([](const auto& p){
        return build<macrodr::Patch_Model>(N_St(5),
                                           build<Q0>(var::build_<Matrix<double>>(5, 5, {{1,0}, {2,1},{3,2},{3,4},{4,3}},
                                                                                 {p()[5],p()[6],p()[7],p()[8],p()[9]})),
                                           build<Qa>(var::build_<Matrix<double>>(5, 5, {{0, 1}, {1, 2},{2, 3}},{p()[2],p()[3],p()[4]})),
                                           build<g>(var::build_<Matrix<double>>(5,1,{{4,0}},{p()[10]})),
                                           build<N_Ch_mean>(p()[0]),
                                           build<curr_noise>(p()[1]),
                                           min_P(1e-7),
                                           Probability_error_tolerance(1e-2),
                                           Conductance_variance_error_tolerance(1e-2));
    });
    
    
    
    auto param1=Parameters<Model1>(Matrix<double>(1,11,std::vector<double>{100,0.05,18,12,6,210,420,630,1680,54,0.5}));
    auto param1Names=std::vector<std::string>{"Num_Chan","noise","k01","k12","k23","k10","k21","k32","k34","k43","conductance"};
    
    auto dparam1=var::selfDerivative(param1);
    // std::cerr<<"dparam1\n"<<dparam1;
    auto n= build<N_Ch_mean>(dparam1()[0]);
    
    auto dp0=dparam1()[0];
    
    auto qq=  var::build_<Matrix<double>>(5, 5, {{0, 1}, {1, 2},{2, 3}},{dparam1()[2],dparam1()[3],dparam1()[4]});
    
    
    auto dm=model1(dparam1);
    auto m=model1(param1);
    
    // std::cerr<<"\n!--------------------------------------------------------------------------!\n";
    //  print(std::cerr,dm);
    //  print(std::cerr,m);
    auto fs = get<Frequency_of_Sampling>(experiment).value();
    auto dini = macrodr::Macro_DMR{}.init(dm, get<initial_ATP_concentration>(experiment));
    auto ini = macrodr::Macro_DMR{}.init(m, get<initial_ATP_concentration>(experiment));
    
    auto t_step=get<Recording>(experiment)()[0];
    auto t_Qx = macrodr::Macro_DMR{}.calc_eigen(m, get<ATP_concentration>(t_step));
    auto dt_Qx = macrodr::Macro_DMR{}.calc_eigen(dm, get<ATP_concentration>(t_step));
    
    
    //  auto dt_Qdt = macrodr::Macro_DMR{}.calc_Qdt(m, t_Qx.value(),
    //                                             get<number_of_samples>(t_step).value() / fs);
    
    auto t_Qdt = macrodr::Macro_DMR{}.calc_Qdt(dm, dt_Qx.value(),
                                               get<number_of_samples>(t_step).value() / fs);
    
    
    std::random_device rd;
    typename std::mt19937_64::result_type seed = rd();
    
    std::mt19937_64 mt(seed);
    
    
    if constexpr (true){
        auto number_replicates=1000;
        auto outputfilename = "../macro_dr/output";
        
        std::string algorithm="_MacroNRC";
        std::ofstream fo(outputfilename+algorithm+".txt");
        
        Matrix<double> mean_dlogL;
        SymPosDefMatrix<double> Cov_dlogL;
        
        Matrix<double> mean_edlogL;
        SymPosDefMatrix<double> Cov_edlogL;
        
        
        for (std::size_t i= 0; i<number_replicates; ++i) {
            auto sim = Macro_DMR{}.sample(
                mt, model1,param1, experiment,
                Simulation_Parameters(Number_of_simulation_sub_steps(10)));
            auto lik = Macro_DMR{}.log_Likelihood(model1,dparam1, sim.value()());
            std::cerr<<"\n"<<i<<"th likelihood!!\n"<<lik;
            if (lik)
            {
                auto v_logL=get<logL>(lik.value());
                auto v_elogL=get<elogL>(lik.value());
                auto v_vlogL=get<vlogL>(lik.value());
                
                auto v_dlogL=derivative(v_logL);
                auto v_delogL=derivative(v_elogL);
                
                
                if (i==0)
                {
                    fo<<"logL\telogL\tvlogL";
                    for (std::size_t j= 0; j<v_dlogL().size(); ++j)
                        fo<<"\t"<<"dlogL_d"<<param1Names[j];
                    for (std::size_t j= 0; j<v_delogL().size(); ++j)
                        fo<<"\t"<<"delogL_d"<<param1Names[j];
                    fo<<"\n";
                    mean_dlogL= Matrix<double>(1,v_delogL().size(),0.0);
                    Cov_dlogL= SymPosDefMatrix<double>(v_delogL().size(),v_delogL().size(),0.0);
                    
                    mean_edlogL= Matrix<double>(1,v_delogL().size(),0.0);
                    Cov_edlogL= SymPosDefMatrix<double>(v_delogL().size(),v_delogL().size(),0.0);
                    
                }
                fo<<primitive(v_dlogL)()<<"\t";
                fo<<primitive(v_elogL)<<"\t";
                fo<<primitive(v_vlogL);
                for (std::size_t j= 0; j<v_dlogL().size(); ++j)
                    fo<<"\t"<<v_dlogL()[j];
                for (std::size_t j= 0; j<v_delogL().size(); ++j)
                    fo<<"\t"<<v_delogL()[j];
                fo<<"\n";
                mean_dlogL=mean_dlogL+v_dlogL();
                Cov_dlogL=Cov_dlogL+XTX(v_dlogL());
                mean_edlogL=mean_edlogL+v_delogL();
                Cov_edlogL=Cov_edlogL+XTX(v_delogL());
            }
            
        }
        mean_dlogL=mean_dlogL*(1.0/number_replicates);
        std::ofstream foa(outputfilename+algorithm+"ana.txt");
        
        std::cerr<<"\nmean_dlogL\n"<<mean_dlogL<<"\n";
        foa<<"\nmean_dlogL\n"<<mean_dlogL<<"\n";
        Cov_dlogL=Cov_dlogL*(1.0/number_replicates);//-XTX(mean_dlogL);
        std::cerr<<"\nCov_dlogL\n"<<Cov_dlogL<<"\n";
        foa<<"\nCov_dlogL\n"<<Cov_dlogL<<"\n";
        
        auto Cov_inv=inv(Cov_dlogL);
        std::cerr<<"\nCov_inv\n"<<Cov_inv;
        foa<<"\nCov_inv\n"<<Cov_inv;
        if (Cov_inv)
        {
            std::cerr<<"\nparameters\n"<<param1()<<"\n";
            
            auto accuracy=mean_dlogL*inv(Cov_dlogL).value();
            auto sensitivity=apply([](auto x){return std::sqrt(x);},diag(Cov_inv.value()));
            std::cerr<<"\naccuracy\n"<<accuracy<<"\n";
            std::cerr<<"\nsensitivityy\n"<<sensitivity<<"\n";
            
            std::cerr<<"\naccuracy rate\n"<<elemDiv(accuracy,param1())<<"\n";
            std::cerr<<"\nsensitivityy\n"<<sensitivity*inv(diag(param1())).value()<<"\n";
            std::cerr<<"\naccuracy  over se \n"<<accuracy* inv(sensitivity*(1.0/std::sqrt(number_replicates)))<<"\n";
            
            foa<<"\naccuracy\n"<<accuracy<<"\n";
            foa<<"\nsensitivityy\n"<<sensitivity<<"\n";
            
            foa<<"\naccuracy rate\n"<<elemDiv(accuracy,param1())<<"\n";
            foa<<"\nsensitivity rate\n"<<sensitivity*inv(diag(param1())).value()<<"\n";
            
            foa<<"\naccuracy  over se \n"<<accuracy* inv(sensitivity*(1.0/std::sqrt(number_replicates)))<<"\n";
            
        }
    }
    
    if (false){
        auto sim = Macro_DMR{}.sample(
            mt, model1,param1, experiment,
            Simulation_Parameters(Number_of_simulation_sub_steps(10)));
        
        auto test_der_Likelihood=var::test_Derivative(
            [&model1,&sim](auto const &dparam1){
                return Macro_DMR{}.log_Likelihood(model1,dparam1, sim.value()());},
            1,1e-10,dparam1);
        if (!test_der_Likelihood)
        {
            std::cerr<<test_der_Likelihood.error()();
        }
    }
    
    
    if (p) {
        auto ss = p.value().str();
        std::cerr << ss;
    } else
        std::cerr << p.error()();
}
