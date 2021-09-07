#ifndef _HEATKERNELCALCULATOR_HPP
#define _HEATKERNELCALCULATOR_HPP

#include <cmath>

#include <distanceMultiplicity.hpp>



#include <iostream>

using namespace std;

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
  
  /// Evaluate the approximation in the range [min,max)
  template <typename F>
  inline SeriesResult innerApproxEval(const int minDistance2,
				      const F f,
				      const int maxExtDistance2) const
  {
    /// Previous iteration result
    double oldRes=0.0;
    
    /// Result
    double res=0.0;
    
    /// Point multiplicity
    int mult;
    
    /// Iteration distance
    int distance2=minDistance2;
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
    while(distance2<maxExtDistance2 and (fabs(oldRes!=res or mult==0)));
    
    return {res,distance2-1};
  }
  
  /// Evaluate the approximation
  template <typename F>
  inline SeriesResult innerApproxEval(const int minDistance2,
				      const F f) const
  {
    return
      innerApproxEval(minDistance2,f,maxDistance2);
  }
  
  /// Evaluate the heat kernel for large values of the argument
  ///
  /// Allows low and high range
  double rangeHeatKernelForLargeArg(const double& z,
				    const int& min,
				    const int& max,
				    const double& q2) const
  {
    /// Result of the series
    const auto truncatedSeriesResult=
      innerApproxEval(min,
		      [&z,&q2](const int& m2)
		      {
			cout<<"Large: "<<m2<<" "<<exp(z*(q2-m2))<<endl;
			return exp(z*(q2-m2));
		      },
		      max).val;
    
    return
      truncatedSeriesResult*pow(2*M_PI,-3.0);
  }
  
  /// Evaluate the heat kernel for large values of the argument
  ///
  /// Implements eq.C.2 of Nuclear Physics B354 (1991) 531-578
  inline double truncatedHeatKernelForLargeArg(const double& z,
					       const int& lambda2,
					       const double& q2) const
  {
    return
      rangeHeatKernelForLargeArg(z,lambda2,maxDistance2,q2);
  }
  
  /// Evaluate the heat kernel for small values of the argument
  ///
  /// Implements eq.C.2 of Nuclear Physics B354 (1991) 531-578
  double truncatedHeatKernelForSmallArg(const double& z,
					const int& lambda2,
					const double& q2) const
  {
    /// result of the series
    const auto fullSeriesResult=
      innerApproxEval(0,
		      [f=-M_PI*M_PI/z,g=q2*z](const int& m2)
		      {
			cout<<"Small: "<<m2<<" "<<exp(g+f*m2)<<endl;
			return exp(g+f*m2);
		      }).val;
    
    cout<<endl<<"Computing truncation"<<endl;
    /// For the subtraction, we use the alternative first representation
    const auto truncatedSeriesResult=
      rangeHeatKernelForLargeArg(z,0,lambda2,q2);
    //cout<<" res: "<<fullSeriesResult*pow(4*M_PI*z,-1.5)<<" "<<truncatedSeriesResult<<endl<<endl;;
    return
      fullSeriesResult*pow(4*M_PI*z,-1.5)-truncatedSeriesResult;
  }
  
   /// Evaluate the kernel
  double operator()(const double& z,
		    const double& q2) const
   {
     /// Compute the cutoff
     const int lambda2=
       floor(q2)+6;
     
     const double s=
       (z>M_PI)?
       truncatedHeatKernelForLargeArg(z,lambda2,q2):
       truncatedHeatKernelForSmallArg(z,lambda2,q2);
     
     return
       s/sqrt(4*M_PI)-pow(4*M_PI,-2)*pow(z,-1.5);
   }
};

#endif
