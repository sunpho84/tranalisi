#ifndef _SERIESOVERZ3CALCULATOR_HPP
#define _SERIESOVERZ3CALCULATOR_HPP

#include <cmath>

#include <distanceMultiplicity.hpp>

/// Evaluates series over points on Z3
struct SeriesOverZ3Calculator
{
  /// Multiplicity of the points
  DistanceMultiplicity distanceMultiplicity;
  
  /// Squared maximal distance
  const int& maxDistance2=
    distanceMultiplicity.maxDistance2;
  
  /// Constructor taking the max distance
  SeriesOverZ3Calculator(const int& maxDistance=30) :
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
  inline SeriesResult operator()(const int minDistance2,
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
  inline SeriesResult operator()(const int minDistance2,
				 const F f) const
  {
    return
      (*this)(minDistance2,f,maxDistance2);
  }
};

#endif
