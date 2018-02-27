#ifndef _KL2_IB_FSE_HPP
#define _KL2_IB_FSE_HPP

#include <IB/common.hpp>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_dilog.h>

const size_t nZ_FSE=4;
const double def_ustar=1.0,def_eps=1e-12;

//! compute Nazario's zeta
vector<double> zeta_FSE(double betal,double ustar=def_ustar,double eps=def_eps);

//! wrapper for vector type
template <class TS,class TV=vmeas_t<TS>>
TV zeta_FSE(const TS &betal,double ustar=def_ustar,double eps=def_eps)
{
  TV out(nZ_FSE);
  for(size_t iel=0;iel<betal.size();iel++) out.put_all_events(zeta_FSE(betal[iel],ustar,eps),iel);
  return out;
}

//! compute FSE corr according to eq.35-36 of Silvano's Nazario's note
vector<double> FSE_corr(double mlep,double mmes,const vector<double> &z0,const vector<double> &z,double L,const vector<size_t> &upto={1});

//! wrapper
template <class TV,class TS=typename TV::base_type>
TV FSE_corr(const TS &mlep,const TS &mmes,const vector<double> &z0,const TV &z,const TS &L,const vector<size_t> &upto={1})
{
  TV out(upto.size());
  
  for(size_t iel=0;iel<out[0].size();iel++)
    {
      vector<double> temp=FSE_corr(mlep[iel],mmes[iel],z0,z.get_all_events(iel),L[iel],upto);
      for(size_t iu=0;iu<upto.size();iu++) out[iu][iel]=temp[iu];
    }
  
  return out;
}

DEFINE_FUNCTION(gsl_sf_dilog)

//! compute the point-like contribution according to eq.9 of Silvano's note
template <class T>
T Gamma_pt(const T &mlep,const T &mmes,const T &DeltaE)
{
  T rl=mlep/mmes;
  T rl2=sqr(rl);
  T rE=2.0*DeltaE/mmes;
  
  // cout<<"mlep: "<<mlep<<endl;
  // cout<<"mmes: "<<mmes<<endl;
  
  T out=1.0/(16.0*sqr(M_PI))*(6.0*log(mmes/MW)+log(rl2)-8.0*log(rE)+(2.0-10.0*rl2)*log(rl2)/(1.0-rl2)-4.0*(1.0+rl2)*log(rl2)*log(rE)/(1.0-rl2)-4.0*(1.0+rl2)*gsl_sf_dilog(T(1.0-rl2))/(1.0-rl2)-3.0+(3.0+sqr(rE)-6.0*rl2-4.0*rE*(1.0-rl2))*log(1.0-rE)/sqr(T(1.0-rl2))+rE*(4.0-rE-4.0*rl2)*log(rl2)/sqr(T(1.0-rl2))-rE*(28.0*rl2+3.0*rE-22.0)/(2.0*sqr(T(1.0-rl2)))-4.0*(1.0+rl2)*gsl_sf_dilog(rE)/(1.0-rl2));
  
  return out;
}

#endif
