#ifndef _FUNCTIONS_HPP
#define _FUNCTIONS_HPP

#include <cmath>

//! the classic 2pts correlator ansatz, with periodic around the world
template <class T> T periodic_twopts_corr_fun(T Z,T M,double TH,double t)
{return Z*exp(-M*TH)*cosh(M*(TH-t))/M;}

//! the classic 2pts correlator ansatz, with anti-periodic around the world
template <class T> T aperiodic_twopts_corr_fun(T Z,T M,double TH,double t)
{return Z*exp(-M*TH)*sinh(M*(TH-t))/M;}

//! the classic 2pts correlator ansatz, with anti-periodic around the world
template <class T> T twopts_corr_fun(T Z,T M,double TH,double t,const int par)
{
  if(par==1) return  periodic_twopts_corr_fun(Z,M,TH,t);
  else       return aperiodic_twopts_corr_fun(Z,M,TH,t);
}

//! used to solve effmass
template <class T> T periodic_corr_fun_effmass(T M,double TH,double t)
{return cosh(M*(TH-t));}
template <class T> T aperiodic_corr_fun_effmass(T M,double TH,double t)
{return sinh(M*(TH-t));}

//! ratio of corr_fun_effmass(t+dt)/corr_fun_effmass(t)
template <class T> T corr_fun_effmass_ratio(T M,double TH,double t,double dt,int par)
{
  if(par==1) return  periodic_corr_fun_effmass(M,TH,t+dt)/ periodic_corr_fun_effmass(M,TH,t);
  else       return aperiodic_corr_fun_effmass(M,TH,t+dt)/aperiodic_corr_fun_effmass(M,TH,t);
}

#endif
