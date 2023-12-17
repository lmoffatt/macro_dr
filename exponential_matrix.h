#ifndef EXPONENTIAL_MATRIX_H
#define EXPONENTIAL_MATRIX_H

#include "matrix.h"
#include "matrix_derivative.h"
#include "derivative_operator.h"

using var::primitive;


template <typename T>
std::size_t log2_norm(const Matrix<T> &x) {
    // Scale A by power of 2 so that its norm is < 1/2 .
    double e = std::ceil(log(norm_inf(x)) / log(2));
    std::size_t s = std::max(0.0, e + 1);
    return s;
}

inline double maxAbs(double x) { return std::abs(x); }

template <typename T> double maxAbs(const Matrix<T> &x) {
    
    double out=maxAbs(x[0]);
    for (std::size_t i=0; i<x.size(); ++i)
        out=std::max(out, maxAbs(x[i]));
    return out;    
}

template <typename CMatrix>
    requires (var::U<CMatrix,Matrix<double>>)
CMatrix expm_taylor(const CMatrix & x, std::size_t order=6)
{
    auto out=x+ DiagonalMatrix<double>(primitive(x).ncols(),primitive(x).ncols(),1.0);
    auto xr=x;
    double a=1.0;
    for (std::size_t n=2; n+1<order; ++n)
    {
        a/=n;
        xr=xr*x;
        out=out+xr*a;
    }
    return out;
}

template <typename CMatrix>
    requires (var::U<CMatrix,Matrix<double>>)
   CMatrix expm_taylor_scaling_squaring(const CMatrix & x, std::size_t order=6)
{
    double max=maxAbs(primitive(x));
    double desired=0.125;
    int k=std::ceil(std::log2(max/desired));
    std::size_t n=std::max(0,k);
    double scale=std::pow(2,-n);
    auto dx=x*scale;
    auto expm_dx=expm_taylor(dx,order);
    auto expm_run=expm_dx;
    for (std::size_t i=0; i<n; ++i)
    {
        expm_run=expm_run*expm_run;
    }
    return expm_run;
}



template <typename CMatrix>
    requires (var::U<CMatrix,Matrix<double>>)
Maybe_error<CMatrix> expm_pade(const CMatrix &M) {
    
    /// http://www2.humusoft.cz/www/papers/tcp08/017_brancik.pdf
    
    auto X = M;
    double c = 0.5;
    auto F = DiagonalMatrix<double>(M.nrows(),M.nrows(),1.0) + c * M;
    auto D = DiagonalMatrix<double>(M.nrows(),M.nrows(),1.0) - c * M;
    
    std::size_t q = 6;
    bool p = true;
    for (std::size_t k = 2; k < q + 1; ++k) {
        c = c * (1.0 * q - 1.0 * k + 1.0) / (1.0 * k * (2.0 * q - k + 1.0));
        X = M * X;
        
        auto cX = c * X;
        F = F+ cX;
        if (p)
            D = D+ cX;
        else
            D = D - cX;
        p = !p;
    }
    
    auto invD = inv(D);
    if (!invD)
        return error_message("cannot invert D: " + invD.error()());
    else {
        F = invD.value() * F;
        
        return F;
    }
    /*    % Pade approximation of exp(M) and diff[exp(M)]
  X=M; Y=dM;
  c=1/2;
  F=eye(size(M))+c*M; dF=c*dM;
  D=eye(size(M))-c*M; dD=-c*dM;
  q=6;
  p=1;
  for
   k=2:q
     c=c*(q-k+1)/(k*(2*q-k+1));
     Y=dM*X+M*Y;
     X=M*X;
     cX=c*X; cY=c*Y;
     F=F+cX; dF=dF+cY;
  if
   p
       D=D+cX; dD=dD+cY;
  else
       D=D-cX; dD=dD-cY;
  end
     p=~p;
  end
  F=D\F;
  dF=D\(dF-dD*F);
  % Undo scaling by repeated squaring
  for
   k=1:r
      dF=dF*F+F*dF;
      F=F*F;
  en
  */
}




template <typename CMatrix>
    requires (var::U<CMatrix,Matrix<double>>)
Maybe_error<CMatrix> full_expm(const CMatrix &x) {
    assert(x.ncols() == x.nrows());
    assert(x.size() > 0);
    
    // Scale A by power of 2 so that its norm is < 1/2 .
    std::size_t s = log2_norm(primitive(x)) + 1;
    
    auto A = x * (1.0 / std::pow(2.0, int(s)));
    
    // Pade approximation for exp(A)
    auto eE = expm_pade(A);
    
    if (!eE)
        return eE.error();
    else {
        
        auto E = eE.value();
        // Undo scaling by repeated squaring
        for (std::size_t k = 0; k < s; k++)
            E = E * E;
        return E;
    }
}


template <typename CMatrix>
    requires (var::U<CMatrix,Matrix<double>>)
CMatrix expm_sure(const CMatrix &x)
{
    auto Maybe_expm=full_expm(x);
    if (Maybe_expm) return Maybe_expm.value();
    else
        return expm_taylor_scaling_squaring(x);
}




#endif // EXPONENTIAL_MATRIX_H
