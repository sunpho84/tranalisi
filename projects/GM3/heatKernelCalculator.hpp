#ifndef _HEATKERNELCALCULATOR_HPP
#define _HEATKERNELCALCULATOR_HPP

#include <cmath>

#include <distanceMultiplicity.hpp>

/// Evaluates the heat kernel
struct HeatKernelCalculator
{
  /// Multiplicity of the points
  DistanceMultiplicity distanceMultiplicity;
  
  /// Squared maximal distance
  const int& maxDistance2=
    distanceMultiplicity.maxDistance2;
  
  /// Constructor taking the max distance
  HeatKernelCalculator(const int& maxDistance=30) :
    distanceMultiplicity(maxDistance)
  {
  }
  
  /// Result of the evaluation
  struct SeriesResult
  {
    /// Value
    double val;
    
    /// Squared distance that was needed
    int distance2Needed;
  };
  
  /// Evaluate the approximation
  template <typename F>
  inline SeriesResult innerApproxEval(const F f) const
  {
    /// Previous iteration result
    double oldRes=0.0;
    
    /// Result
    double res=0.0;
    
    /// Point multiplicity
    int mult;
    
    /// Iteration distance
    int distance2=0;
    do
      {
	mult=distanceMultiplicity(distance2);
	
	if(mult)
	  {
	    /// Current contribution
	    const double contr=
	      mult*f(distance2);
	    
	    oldRes=res;
	    res=oldRes+contr;
	  }
	
	distance2++;
      }
    while(distance2<maxDistance2 and (fabs(oldRes!=res or mult==0)));
    
    return {res,distance2-1};
  }
  
  /// Evaluate the heat kernel for large values of the argument
  ///
  /// Implements eq.C.1 of Nuclear Physics B354 (1991) 531-578
  double heatKernelForLargeArg(const double& z) const
  {
    /// Result of the series
    const auto seriesResult=
      innerApproxEval([&z](const int& m2)
      {
	return exp(-z*m2);
      });
    
    return
      seriesResult.val*pow(2*M_PI,-3.0);
  }
  
  /// Evaluate the heat kernel for small values of the argument
  ///
  /// Implements eq.C.2 of Nuclear Physics B354 (1991) 531-578
  double heatKernelForSmallArg(const double& z) const
  {
    /// result of the series
    const auto seriesResult=
      innerApproxEval([f=-M_PI*M_PI/z](const int& m2)
      {
	return exp(f*m2);
      });
    
    return
      seriesResult.val*pow(4*M_PI*z,-1.5);
  }
  
  /// Evaluate the kernel
  double operator()(const double& z) const
  {
    return
      (z>M_PI)?heatKernelForLargeArg(z):heatKernelForSmallArg(z);
  }
};

#endif
