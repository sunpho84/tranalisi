#ifndef _VKVKREPRESENTATION_HPP
#define _VKVKREPRESENTATION_HPP

#include "aLaLuscherRepresentation.hpp"

/// Representation of the vector correlator
template <typename TwoPionsRepresentation>
struct VKVKRep
{
  const TwoPionsRepresentation LuschRep;
  
  const double rDual;
  
  const double eThr;
  
  const double mRho;
  
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
      LuschRep(t);
  }
  
  VKVKRep(const TwoPionsRepresentation& LuschRep,
	  const double& rDual,
	  const double& eThr,
	  const double& mRho) :
    LuschRep(LuschRep),
    rDual(rDual),
    eThr(eThr),
    mRho(mRho)
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

using VKVKRepFiniteVol=
  VKVKRep<ALaLuscherRepresentation>;

using VKVKRepInfiniteVol=
  VKVKRep<ALaLuscherRepresentationInfVol>;

#endif
