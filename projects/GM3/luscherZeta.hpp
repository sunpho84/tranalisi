#ifndef _TRUNCATEDHEATKERNEL_HPP
#define _TRUNCATEDHEATKERNEL_HPP

#include <gsl/gsl_integration.h>

#include <seriesOverZ3Calculator.hpp>

/// Computes the Luscher Zeta function
struct LuscherZetaCalculator
{
  /// Computes the integrand of the Luscher Zeta function
  struct IntegrandCalculator
  {
    /// Computes a series over Z3
    const SeriesOverZ3Calculator& seriesCalculator;
    
    /// Value of q2
    const double q2;
    
    /// Squared cut-off
    const double lambda2;
    
    /// Construct taking the calculator, q2 and lambda2
    IntegrandCalculator(const SeriesOverZ3Calculator& seriesCalculator,
			const double& q2,
			const double& lambda2) :
      seriesCalculator(seriesCalculator),
      q2(q2),lambda2(lambda2)
    {
    }
    
    /// Evaluate the heat kernel for large values of the argument
    ///
    /// Allows low and high range
    double rangeHeatKernelForLargeArg(const double& z,
				      const int& min,
				      const int& max) const
    {
      /// Result of the series
      const auto truncatedSeriesResult=
	seriesCalculator(min,
			 [&z,g=q2*z](const int& m2)
			 {
			   return exp(g-z*m2);
			 },
			 max).val;
      
      /// Result to be fully normalized
      return
	truncatedSeriesResult*pow(2*M_PI,-3.0);
    }
    
    /// Evaluate the heat kernel for large values of the argument
    ///
    /// Implements eq.C.2 of Nuclear Physics B354 (1991) 531-578
    inline double truncatedHeatKernelForLargeArg(const double& z) const
    {
      const double resToBeNormalized=
	rangeHeatKernelForLargeArg(z,lambda2,seriesCalculator.maxDistance2);
      
      return
	fullyRenormalizedResult(resToBeNormalized,z);
    }
    
    /// Evaluate the heat kernel for small values of the argument
    ///
    /// Implements eq.C.2 of Nuclear Physics B354 (1991) 531-578
    double truncatedHeatKernelForSmallArg(const double& z) const
    {
      /// result of the series
      const auto fullSeriesResult=
	seriesCalculator(0,
			 [f=-M_PI*M_PI/z,g=q2*z](const int& m2)
			 {
			   return exp(g+f*m2);
			 }).val;
      
      
      /// For the subtraction, we use the large arg representation
      const auto truncatedSeriesResult=
	rangeHeatKernelForLargeArg(z,0,lambda2);
      
      /// Result to be fully normalized
      const double resToBeNormalized=
	fullSeriesResult*pow(4*M_PI*z,-1.5)-truncatedSeriesResult;
      
      return
	fullyRenormalizedResult(resToBeNormalized,z);
    }
    
    /// Add the normalization and subtract the divergence in the origin
    inline double fullyRenormalizedResult(const double& s,
					  const double& z) const
    {
      return
	s/sqrt(4*M_PI)-pow(4*M_PI,-2)*pow(z,-1.5);
    }
  };
  
  /// Series calculator to be used for the integrand
  SeriesOverZ3Calculator seriesCalculator;
  
  /// Absolute error
  static constexpr double epsAbs=0;
  
  /// Relative error
  static constexpr double epsRel=1e-8;
  
  /// Size of the workspace
  static constexpr int workspaceSize=
    1000;
  
  /// Workspace used by gsl
  gsl_integration_workspace *workspace;
  
  /// Constructor
  LuscherZetaCalculator()
    : workspace(gsl_integration_workspace_alloc(workspaceSize))
  {
  }
  
  /// Destructor
  ~LuscherZetaCalculator()
  {
    gsl_integration_workspace_free(workspace);
  }
  
  /// Wraps the integrand calculator for small args, to be used by gsl
  static double smallArgIntegrandCalculatorWrapper(double z,
						   void* params)
  {
    return
      static_cast<IntegrandCalculator*>(params)->truncatedHeatKernelForSmallArg(z);
  }
  
  /// Wraps the integrand calculator for large args, to be used by gsl
  static double largeArgIntegrandCalculatorWrapper(double z,
						   void* params)
  {
    return
      static_cast<IntegrandCalculator*>(params)->truncatedHeatKernelForLargeArg(z);
  }
  
  /// Computes the integrand up to 1.0
  double integrateSmallRange(const IntegrandCalculator& integrandCalculator) const
  {
    /// Wraps the integrand calculator for small arguments, and the parameters
    gsl_function integrandWrapper{smallArgIntegrandCalculatorWrapper,(void*)(&integrandCalculator)};
    
    /// Result absolute error
    double abserr;
    
    /// Result
    double res;
    
    gsl_integration_qags(&integrandWrapper,0.0,1.0,epsAbs,epsRel,workspaceSize,workspace,&res,&abserr);
    
    return
      res;
  }
  
  /// Computes the integrand beyond 1.0
  double integrateLargeRange(const IntegrandCalculator& integrandCalculator) const
  {
    /// Wraps the integrand calculator for large arguments, and the parameters
    gsl_function integrandWrapper{largeArgIntegrandCalculatorWrapper,(void*)(&integrandCalculator)};
    
    /// Result absolute error
    double abserr;
    
    /// Result
    double res;
    
    gsl_integration_qagiu(&integrandWrapper,1.0,epsAbs,epsRel,workspaceSize,workspace,&res,&abserr);
    
    return
      res;
  }
  
  /// Computes the analytic part of the Zeta function (first part of eq. C.8)
  double computeAnalyticPart(const double& q2,
			     const double& lambda2) const
  {
    /// Normalization of the result
    static const double n=
      1.0/sqrt(4*M_PI);
    
    /// Result of the 
    const double s=
      seriesCalculator(0,
		       [&q2](const int& m2)
		       {
			 return
			   1.0/(m2-q2);
		       },
		       lambda2).val;
    
    return
      s*n;
  }
  
  /// Computes the Luscher zeta function
  double operator()(const double& q2,
		    const double& lambda2) const
  {
    /// Computes the integrand
    IntegrandCalculator integrandCalculator(seriesCalculator,q2,lambda2);
    
    /// Integrate in the range [0:1]
    const double smallRange=
      integrateSmallRange(integrandCalculator);
    
    /// Integrates in the range [1:+inf)
    const double largeRange=
      integrateLargeRange(integrandCalculator);
    
    /// Analytic part
    const double analyticPart=
      computeAnalyticPart(q2,lambda2);
    
    /// Normalization of the integral part, as per eq. C.8
    static const double integralNormalization=
      pow(2.0*M_PI,3.0);
    
    return
      analyticPart+(smallRange+largeRange)*integralNormalization;
  }
  
  /// Computes the Luscher zeta function, with default value of the cut-off
  double operator()(const double& q2) const
  {
    /// Compute the cut-off
    const int lambda=
      ceil(sqrt(q2))+6;
    
    /// Compute the squared cut-off
    const double lambda2=
      lambda*lambda;
    
    return
      (*this)(q2,lambda2);
  }
};

#endif
