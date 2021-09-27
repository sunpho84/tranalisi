#ifndef _GSLWRAP_HPP
#define _GSLWRAP_HPP

#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>

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
   gsl_function* getWrapper()
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

/// Use gsl to integrate from lower to infinity
template <typename F>
double gslIntegrateUpToInfinity(F&& fun,
				const double& lower,
				const double& eps=1e-8)
{
  /// Absolute error
  constexpr double epsAbs=0;
  
  /// Size of the workspace
  static constexpr int workspaceSize=
    1000;
  
  /// Result
  double res;
  
  /// Absolute error
  double absErr;
  
  /// Workspace used by gsl
  gsl_integration_workspace *workspace=
    gsl_integration_workspace_alloc(workspaceSize);
  
  gsl_integration_qagiu(GslFunction(fun).getWrapper(),lower,epsAbs,eps,workspaceSize,workspace,&res,&absErr);
  
  gsl_integration_workspace_free(workspace);
  
  return
    res;
}

#endif
