#ifndef _TRUNCATEDHEATKERNEL_HPP
#define _TRUNCATEDHEATKERNEL_HPP

#include <distanceMultiplicity.hpp>

#include <heatKernelCalculator.hpp>

struct LuscherZetaCalculator
{
  DistanceMultiplicity distanceMultiplicity;
  
  HeatKernelCalculator heatKernelCalculator;
  
  const double norm=
    pow(2*M_PI,-3.0);
  
  LuscherZetaCalculator(const int& lambdaMax) :
    distanceMultiplicity(lambdaMax)
  {
  }
  
  double truncatedHeatKernel(const int& lambda,
			     const double& t) const
  {
    const double heatKernel=
      heatKernelCalculator(t);
    
    double truncation=
      0.0;
    
    const int lambda2=
      lambda*lambda;
    
    for(int m2=0;m2<lambda2;m2++)
      truncation+=distanceMultiplicity(m2)*exp(-t*m2);
    
    return
      heatKernel-norm*truncation;
  }
  
  double integrand(const double& t,
		   const double& q2,
		   const double& lambda) const
  {
    static const double d=
      pow(4*M_PI*M_PI,-2.0);
    
    return
      exp(t*q2)*truncatedHeatKernel(lambda,t)-d*pow(t,-1.5);
  }
};

#endif
