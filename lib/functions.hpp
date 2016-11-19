#ifndef _FUNCTIONS_HPP
#define _FUNCTIONS_HPP

#include <cmath>

//! the classic 2pts correlator ansatz, with periodic around the world
template <class T> T periodic_twopts_corr_fun(T Z,T M,double TH,double t)
{return Z*Z*exp(-M*TH)*cosh(M*(TH-t))/M;}
//! the classic 2pts correlator ansatz, with anti-periodic around the world
template <class T> T aperiodic_twopts_corr_fun(T Z,T M,double TH,double t)
{return Z*Z*exp(-M*TH)*sinh(M*(TH-t))/M;}
//
template <class T> T twopts_corr_fun(T Z,T M,double TH,double t,const int par)
{
  if(par==1) return  periodic_twopts_corr_fun(Z,M,TH,t);
  else       return aperiodic_twopts_corr_fun(Z,M,TH,t);
}

//! used to solve effmass
template <class T> T periodic_fun_effmass(T M,double TH,double t)
{return cosh(M*(TH-t));}
template <class T> T aperiodic_fun_effmass(T M,double TH,double t)
{return sinh(M*(TH-t));}

//! ratio of corr_fun_effmass(t+dt)/corr_fun_effmass(t)
template <class T> T effmass_two_times_ratio_fun(T M,double TH,double t,double dt,int par)
{
  if(par==1) return  periodic_fun_effmass(M,TH,t+dt)/ periodic_fun_effmass(M,TH,t);
  else       return aperiodic_fun_effmass(M,TH,t+dt)/aperiodic_fun_effmass(M,TH,t);
}

///////////////////////////////////////////////////////////////////////////////////////

//! return the ratio of a corr with insertion and an original corr
template< class T> T periodic_twopts_corr_with_ins_ratio_fun(T M,T A,T SL,double TH,double t)
{return A+SL*(t-TH)*tanh(M*(t-TH));}
//! return the ratio of a corr with insertion and an original corr
template< class T> T aperiodic_twopts_corr_with_ins_ratio_fun(T M,T A,T SL,double TH,double t)
{return A+SL*(t-TH)/tanh(M*(t-TH));}
template <class T> T twopts_corr_with_ins_ratio_fun(T M,T A,T SL,double TH,double t,const int par)
{
  if(par==1) return  periodic_twopts_corr_with_ins_ratio_fun(M,A,SL,TH,t);
  else       return aperiodic_twopts_corr_with_ins_ratio_fun(M,A,SL,TH,t);
}

//! return the time dependence expected from a certain mass, of the difference of the slope at time t+dt and t
template <class T> T twopts_corr_with_ins_ratio_diff_tdep(T M,double TH,double t,double dt,int par)
{
  if(par==1) return (t+dt-TH)*tanh(M*(t+dt-TH))-(t-TH)*tanh(M*(t-TH));
  else       return (t+dt-TH)/tanh(M*(t+dt-TH))-(t-TH)/tanh(M*(t-TH));
}

#endif
