#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <cmath>

namespace LEC
{
  const double r0=2.1975;
  const double db0=11.972809241766459*r0,f0=0.26772510275074468/r0;
  const double l3=1.0423672088359199,l4=2.1490196003108331;
  const double xi=db0/(16*M_PI*M_PI*f0*f0)/r0/r0;
}

/// mpi as a function of ml
double mPiFun(const double& ml)
{
  using namespace LEC;
  
  const double m2Pi=
    db0*ml*(1+  xi*ml*(log(db0*ml)+l3));
  
  return
    sqrt(m2Pi)/r0;
}

/// fPi as a function of ml
double fPiFun(const double& ml)
{
  using namespace LEC;
  
  const double f=
    f0*(1-2*xi*ml*(log(db0*ml)-l4));
  
  return
    f;
}

/// Compute the ratio for a given ml
double ratio(const double& ml)
{
  return
    mPiFun(ml)/fPiFun(ml);
}

/// Find ml for a given ratio
double findMlForRatio(const double& x)
{
  double m=
    1e-5;
  double s=
    1e-5;
  
  do
    if(ratio(m+s)<=x)
      {
	m+=s;
	s*=2;
      }
    else s/=2;
  while(s>1e-14);
  
  return
    m;
}

double aFromMPiFPi(const double& aMPi,
		   const double& aFPi)
{
 const double ml=
    findMlForRatio(aMPi/aFPi);
  
  const double a=
    aMPi/mPiFun(ml);
  
  return
    a;
}
