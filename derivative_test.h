#ifndef DERIVATIVE_TEST_H
#define DERIVATIVE_TEST_H
#include "maybe_error.h"
#include "derivative_operator.h"
#include "parameters_derivative.h"

namespace var{



template<class F, class Parameters, class... Xs>
Maybe_error<bool> 
test_Derivative(F f, const Parameters x, double dx, double eps, Xs const&...xs)
{
    auto dY=f(xs...);
    if (!(is_valid(dY)))
        return get_error(dY);
    else{
        auto T0=Taylor_first(get_value(dY),x,dx);
    
         auto T1=std::invoke(f,Taylor_first(xs,x,dx)...);
        if (!is_valid(T1))
            return get_error(T1);
        
        auto out=test_equality(T0,get_value(T1),eps);
        auto dout=test_equality(get_value(T0)-primitive(get_value(dY)),get_value(T1)-primitive(get_value(dY)),eps);
        if (!dout)
        {
            std::stringstream ss;
            ss<<"\n test_Derivative on function ";
            ss<<typeid(F).name();
            ss<<" with delta parameters "<<typeid(Parameters).name()<<" equal to: ";
            ss<<"\n x dx\n"<<x()<<"\ndx="<<dx<<"\n";
            ss<<"eps ="<<eps<<"\n";
          
            ss<<"\n-----error---------------\n";
            ss<<dout.error()();
            
            ss<<"\n--------------------\n";
            
            ss<<"\n--------------------\n";
            
            
            ss<<"\n dY\n";
            //print(ss,get_value(dY));
            ss<<"\n--------------------\n";
            ss<<"\n--------------------\n";
            ss<<"\n T0\n";
            print(ss,T0);
            ss<<"\n--------------------\n";
            ss<<"\n T1\n";
            print(ss,get_value(T1));
            ss<<"\n--------------------\n";
            ss<<"\n--------------------\n";
            ss<<"\n delta_T0\n";
            print(ss,get_value(T0)-primitive(get_value(dY)));
            ss<<"\n--------------------\n";
            ss<<"\n--------------------\n";
            print(ss,get_value(T1)-primitive(get_value(dY)));
            ss<<"\n--------------------\n";
            return error_message(ss.str());
            
        }
        else return dout; 
    }
}   
template<class F,  class... Xs>
    requires ((std::is_same_v<NoDerivative,decltype(get_dx_of_dfdx(std::declval<Xs>()...))>))
Maybe_error<bool>
test_Derivative(F , double , double , Xs...)
{
    
    auto out=Maybe_error<bool>(true);
    
    return out;
    
}


template<class F,  class... Xs>
    requires (!(std::is_same_v<NoDerivative,decltype(get_dx_of_dfdx(std::declval<Xs>()...))>))
Maybe_error<bool>
test_Derivative(F f, double dx, double eps, const Xs&...xs)
{
    using Y=dx_of_dfdx_t<Xs...>;
    auto x=get_dx_of_dfdx(xs...);
    
    auto out=Maybe_error<bool>(true);
    
    for (std::size_t i=0; i<x().size(); ++i)
    {
        auto xi=x;
        xi()[i]=1.0;
        auto test_i=test_Derivative(f,xi,dx,eps,xs...);
        if (!test_i)
            out=error_message(out.error()()+"\n at "+std::to_string(i)+"th parameter :"+test_i.error()());
    }
    return out;
    
}

}

#endif // DERIVATIVE_TEST_H
