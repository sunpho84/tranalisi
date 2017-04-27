#ifndef _MEL_FSE_HPP
#define _MEL_FSE_HPP

#include <tranalisi.hpp>

vector<double> zeta(double betal,double ustar=1.0,double eps=1e-12);

//! compute FSE corr according to eq.35-36 of Silvano's Nazario's note
double FSE_corr(double mlep,double mmes,double betal,double L,size_t upto=1);

//! wrapper
template <class TS> TS FSE_corr(double mlep,TS mmes,TS betal,double L,size_t upto=1)
{
  TS out;
  for(size_t iel=0;iel<mmes.size();iel++) out[iel]=FSE_corr(mlep,mmes[iel],betal[iel],L,upto);
  return out;
}

#endif
