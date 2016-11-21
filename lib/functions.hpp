#ifndef _FUNCTIONS_HPP
#define _FUNCTIONS_HPP

#include <cmath>
#include <cstdint>
#include <vector>

using namespace std;

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

//! handle keeping TH and parity and returning the function
class two_pts_corr_fun_t
{
  size_t TH;
  int par;
public:
  enum{iZ,iM};
  two_pts_corr_fun_t(size_t TH,int par) : TH(TH),par(par) {}
  double operator()(const vector<double> &p,double x) {return twopts_corr_fun(p[iZ],p[iM],TH,x,par);}
};

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

//! handle keeping TH and parity and returning the function
class two_pts_corr_with_ins_fun_t
{
  size_t TH;
  int par;
public:
  enum{iM,iA,iSL};
  two_pts_corr_with_ins_fun_t(size_t TH,int par) : TH(TH),par(par) {}
  double operator()(const vector<double> &p,double x) {return twopts_corr_with_ins_ratio_fun(p[iM],p[iA],p[iSL],TH,x,par);}
};

#endif
