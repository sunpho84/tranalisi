#ifndef _GSLWRAP_HPP
#define _GSLWRAP_HPP

#include <gsl/gsl_deriv.h>

/// Wraps any function into a gsl one
template <typename F>
 struct GslFunction :
   public gsl_function
 {
   /// Constructor
   GslFunction(const F& func) :
     gsl_function{GslFunction::invoke,this},
     func(func)
   {
   }
   
   /// Function to be called
   const F& func;
   
   /// External function which invoke func
   static double invoke(double x,
			void *params)
   {
     return
       static_cast<GslFunction*>(params)->func(x);
   }
   
   /// Returns the wrapper
   const gsl_function* getWrapper() const
   {
     return
       this;
   }
};

/// Use gsl to compute derivative
template <typename F>
double gslDeriv(F&& fun,
		const double& x,
		const double& eps=1e-8)
{
  /// Result
  double res;
  
  /// Absolute error
  double absErr;
  
  gsl_deriv_central(GslFunction(fun).getWrapper(),x,eps,&res,&absErr);
  
  return
    res;
}

#endif
