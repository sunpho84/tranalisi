#ifndef _MEL_FSE_HPP
#define _MEL_FSE_HPP

#include <tranalisi.hpp>

vector<double> zeta(double betal,double ustar=1.0,double eps=1e-12);

//! compute FSE corr
template <class TS> TS FSE_corr(double mlep,TS mmes,double betal,double L)
{
  const double mW=80.385;

  TS rl2=sqr(mlep/mmes); //check
  TS cIR=1/(8*sqr(M_PI))*((1+rl2)*log(rl2)/(1-rl2)-1);
  vector<TS> c;
  for(size_t iel=0;iel<mmes.size();iel++)
    {
      vector<double> z0=zeta(mlep,mmes[iel],0);
      vector<double> z=zeta(mlep,mmes[iel],betal);
      c[0][iel]=1/(16*sqr(M_PI))*(2*log(mmes/mW)+((2-6*rl2)*log(rl2)+(1+rl2)*sqr(log(rl2)))/(1-rl2)-5.0/2)+(z0[2]-2*z[2])/2;
      c[1][iel]=-2*(1+rl2)/(1-rl2)*z0[1]+8*rl2/(1-sqr(rl2))*z[1];
      c[2][iel]=4/(1-rl2)*z[0]-8/(1-sqr(rl2))*z[3];
      c[3][iel]=(-5+rl2*(-5+rl2*(-3+rl2)))/pow(1+rl2,3);
    }
  TS mL=mmes*L;
  
  return cIR*2*log(mL)+c[0]+c[1]/mL+c[2]/sqr(mL)+c[3]/pow(mL,3);
}

#endif
