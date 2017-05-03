#ifndef _KL2_IB_FSE_HPP
#define _KL2_IB_FSE_HPP

#include <common.hpp>

const size_t nZ_FSE=4;
const double def_ustar=1.0,def_eps=1e-12;

//! compute Nazario's zeta
vector<double> zeta_FSE(double betal,double ustar=def_ustar,double eps=def_eps);

//! wrapper for vector type
template <class TS,class TV=vmeas_t<TS>> TV zeta_FSE(const TS &betal,double ustar=def_ustar,double eps=def_eps)
{
  TV out(nZ_FSE);
  for(size_t iel=0;iel<betal.size();iel++) out.put_all_events(zeta_FSE(betal[iel],ustar,eps),iel);
  return out;
}

//! compute FSE corr according to eq.35-36 of Silvano's Nazario's note
double FSE_corr(double mlep,double mmes,const vector<double> &z0,const vector<double> &z,double L,size_t upto=1);

//! wrapper
template <class TV,class TS=typename TV::base_type> TS FSE_corr(const TS &mlep,const TS &mmes,const vector<double> &z0,const TV &z,const TS &L,size_t upto=1)
{
  TS out;
  for(size_t iel=0;iel<out.size();iel++) out[iel]=FSE_corr(mlep[iel],mmes[iel],z0,z.get_all_events(iel),L[iel],upto);
  return out;
}

#endif
