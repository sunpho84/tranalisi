#ifndef _PHIFUNCTION_HPP
#define _PHIFUNCTION_HPP

#include <luscherZeta.hpp>

/// Computes the r.h.s of eq.A.3 of Nuclear Physics B364 (1991) 237-251
struct TanPhiCalculator
{
  /// Computes the Zeta function
  const LuscherZetaCalculator luscherZetaCalculator;
  
  /// Approximation according to eq. A.4
  double approximation(const double& q2) const
  {
    /// Equation A.5
    const double x=
      q2*(8.91363-
	  q2*(16.5323-
	      q2*(8.402-
		  q2*6.95)));
    
    /// Numerator
    const double num=
      2*M_PI*M_PI*pow(q2,1.5);
    
    /// Denominator
    const double den=
      1+x;
    
    return
      num/den;
  }
  
  /// Full evalation
  double full(const double& q2) const
  {
    /// Square root of q2
    const double q=
      sqrt(q2);
    
    /// Pi^(3/2)
    static const double pi32=
      pow(M_PI,1.5);
    
    /// Numerator
    const double num=
      pi32*q;
    
    /// Denominator
    const double den=
      luscherZetaCalculator(q2);
    
    return
      -num/den;
  }
  
  /// Calculate tan(phi)
  double operator()(const double& q2) const
  {
    if(q2<=0.01)
      return approximation(q2);
    else
      return full(q2);
  }
};

#endif
