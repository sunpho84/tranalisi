#ifndef _VKVKREPRESENTATION_HPP
#define _VKVKREPRESENTATION_HPP

#include <aLaLuscherRepresentation.hpp>

/// Representation of the vector correlator
struct VKVKRep
{
  const ALaLuscherRepresentationCached<true>& LuschRep;
  
  const double mPi;
  
  const double rDual;
  
  const double eThr;
  
  const double mRho;
  
  const double g2;
  
  double dualPartFun(const double& t) const
  {
    const double mRhoPlusEDual=
      mRho+eThr;
    
    return
      5/(18*M_PI*M_PI)*
      rDual/(t*t*t)*
      exp(-mRhoPlusEDual*t)*
      (1+mRhoPlusEDual*t+sqr(mRhoPlusEDual)*t*t/2);
  }
  
  /// Holds the two parts of the representation
  struct VKVKRepParts
  {
    double dual;
    
    double Lusch;
    
    double operator()() const
    {
      return
	dual+Lusch;
    }
  };
  
  double LuschPartFun(const double& t) const
  {
    return
      LuschRep(mPi,L,mRho,g2)(t);
  }
  
  VKVKRep(const ALaLuscherRepresentationCached<true>& LuschRep,
	  const double& mPi,
	  const double& rDual,
	  const double& eThr,
	  const double& mRho,
	  const double& g2) :
    LuschRep(LuschRep),
    mPi(mPi),
    rDual(rDual),
    eThr(eThr),
    mRho(mRho),
    g2(g2)
  {
  }
  
  VKVKRepParts operator()(const double& t) const
  {
    const double dualPart=
      dualPartFun(t);
      
    const double LuschPart=
      LuschPartFun(t);
    
    return
      {dualPart,LuschPart};
  }
};

jack_t<VKVKRep> fitVKVK(const int nLevels,
			const size_t tMin);

#endif
