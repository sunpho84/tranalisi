#ifndef _WINDOW_HPP
#define _WINDOW_HPP

#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include <GM3/kernels.hpp>
#include <GM3/perens.hpp>

//! integrate the kernel
template <class TV,class TS=typename TV::base_type>
TS integrate_window(const TV &corr,size_t T,const TS &a,size_t upto,function<djack_t(const double&)> w,size_t ord=1)
{
  //store the kernel
  TV kern(T/2);
  for(size_t t=1;t<T/2;t++)
    {
      kern[t]=kern_num(corr[t],t,a)*w(t);
    }
  kern[0]=0.0;
  
  return 4*sqr(alpha_em)*integrate_corr_up_to(kern,upto);
}

void perens_t::computeWindow(const RegoType &rego)
{
  const auto getTheta=
    [this](const double& tPrime)
    {
      return [tPrime=tPrime/0.197/a,
	      Delta=0.4 /* fm */ /0.197/a](const double& t)
      {
	return 1/(1+exp(-2*(t-tPrime)/Delta));
      };
    };
  
  auto theta0=
    getTheta(0.4 /* fm */);
  auto theta1=
    getTheta(1.0 /* fm */);
  
  auto thetaSD=
    [theta0](const double& t)
    {
      return 1-theta0(t);
    };
  
  auto thetaIW=
    [theta0,theta1](const double& t)
    {
      return theta0(t)-theta1(t);
    };
  
  const djack_t Zrego=
    Z[regoZId[rego]];
  
  const double eu=2.0/3,ed=-1.0/3;
  const djvec_t corr=
    getAveForRego(0,nSources,1,rego)*sqr(Zrego)*(sqr(eu)+sqr(ed));
  
  grace_file_t SDplot(basePath+"/plots/SD.xmg");
  for(double t=0;t<=THp1;t+=0.1)
    SDplot.write_ave_err(t,thetaSD(t).ave_err());
  
  
  grace_file_t IWplot(basePath+"/plots/IW.xmg");
  for(double t=0;t<=THp1;t+=0.1)
    IWplot.write_ave_err(t,thetaIW(t).ave_err());
  
  //integrate_window(corr,T,a,T/2+1,thetaSD);
  
}

#endif
