#ifndef _FUNCTIONS_HPP
#define _FUNCTIONS_HPP

#include <cmath>
#include <cstdint>
#include <vector>

#include <gsl/gsl_integration.h>

#include <math.hpp>
#include <tools.hpp>

using namespace std;

//! the classic 2pts correlator ansatz, with no around the world term
template <class T>
T nonperiodic_two_pts_corr_fun(T Z2,T M,double t)
{return Z2*exp(-M*t)/(2*M);}
//! the classic 2pts correlator ansatz, with periodic around the world
template <class T>
T periodic_two_pts_corr_fun(T Z2,T M,double TH,double t)
{return Z2*exp(-M*TH)*cosh(M*(TH-t))/M;}
//! the classic 2pts correlator ansatz, with anti-periodic around the world
template <class T>
T antiperiodic_two_pts_corr_fun(T Z2,T M,double TH,double t)
{return Z2*exp(-M*TH)*sinh(M*(TH-t))/M;}
//
template <class T>
T two_pts_corr_fun(T Z2,T M,double TH,double t,const int par)
{
  switch(par)
    {
    case -1: return antiperiodic_two_pts_corr_fun(Z2,M,TH,t); break;
    case  0: return nonperiodic_two_pts_corr_fun(Z2, M, t);break;
    case +1: return  periodic_two_pts_corr_fun(Z2,M,TH,t);break;
    default: CRASH("unknown periodicity %d",par);return M;
    }
}

//! used to solve effmass
template <class T>
T periodic_fun_effmass(T M,double TH,double t)
{return cosh(M*(TH-t));}
template <class T>
T nonperiodic_fun_effmass(T M,double TH,double t)
{return exp(-M*t);}
template <class T>
T antiperiodic_fun_effmass(T M,double TH,double t)
{return sinh(M*(TH-t));}

//! ratio of corr_fun_effmass(t+dt)/corr_fun_effmass(t)
template <class T>
T effmass_two_times_ratio_fun(T M,double TH,double t,double dt,int par)
{
  switch(par)
    {
    case -1: return antiperiodic_fun_effmass(M,TH,t+dt)/antiperiodic_fun_effmass(M,TH,t);break;
    case  0: return nonperiodic_fun_effmass(M,TH,t+dt)/nonperiodic_fun_effmass(M,TH,t);break;
    case +1: return periodic_fun_effmass(M,TH,t+dt)/periodic_fun_effmass(M,TH,t);break;
    default: CRASH("unknown periodicity %d",par);return M;
    }
}

//! handle keeping TH and parity and returning the function
class two_pts_corr_fun_t
{
  size_t TH;
  int par;
public:
  enum{iZ2,iM};
  two_pts_corr_fun_t(size_t TH,int par) : TH(TH),par(par) {}
  double operator()(const vector<double> &p,const double &x)
  {return two_pts_corr_fun(p[iZ2],p[iM],TH,x,par);}
};

///////////////////////////////////////////////////////////////////////////////////////

//NB: Slope is really the slope!!!!

//! return the ratio of a corr with insertion and an original corr, without return contribution
template <class T>
T nonperiodic_two_pts_corr_with_ins_ratio_fun(const T &M,const T &DZ2_fr_Z2,const T &SL,const double &t)
{return DZ2_fr_Z2+SL/M+SL*t;}
//! return the ratio of a corr with insertion and an original corr with periodic return contribution
template <class T>
T periodic_two_pts_corr_with_ins_ratio_fun(const T &M,const T &DZ2_fr_Z2,const T &SL,const double &TH,const double &t)
{return DZ2_fr_Z2+SL*(1.0/M+TH)-SL*(t-TH)*tanh(M*(t-TH));}
//! return the ratio of a corr with insertion and an original corr with antiperiodic return contribution
template <class T>
T antiperiodic_two_pts_corr_with_ins_ratio_fun(const T &M,const T &DZ2_fr_Z2,const T &SL,const double &TH,const double &t)
{return DZ2_fr_Z2+SL*(1.0/M+TH)-SL*(t-TH)/tanh(M*(t-TH));}
template <class T>
T two_pts_corr_with_ins_ratio_fun(const T &M,const T &DZ2_fr_Z2,const T &SL,const double &TH,const double &t,const int &par)
{
  switch(par)
    {
    case -1: return antiperiodic_two_pts_corr_with_ins_ratio_fun(M,DZ2_fr_Z2,SL,TH,t);break;
    case  0: return  nonperiodic_two_pts_corr_with_ins_ratio_fun(M,DZ2_fr_Z2,SL,t);break;
    case +1: return     periodic_two_pts_corr_with_ins_ratio_fun(M,DZ2_fr_Z2,SL,TH,t);break;
    default: CRASH("unknown periodicity %d",par);return M;
    }
}

//! return the time dependence expected from a certain mass, of the difference of the slope at time t+dt and t
template <class T>
T two_pts_corr_with_ins_ratio_diff_tdep(const T &M,const double &TH,const double &t,const double &dt,const int &par)
{return two_pts_corr_with_ins_ratio_fun(M,M*0.0,M*0.0+1.0,TH,t+dt,par)-two_pts_corr_with_ins_ratio_fun(M,M*0.0,M*0.0+1.0,TH,t,par);}

//! handle keeping TH and parity and returning the function
class two_pts_corr_with_ins_fun_t
{
  size_t TH;
  int par;
public:
  enum{iM,iDZ2_fr_Z2,iSL};
  two_pts_corr_with_ins_fun_t(size_t TH,int par) : TH(TH),par(par) {}
  double operator()(const vector<double> &p,double x) {return two_pts_corr_with_ins_ratio_fun(p[iM],p[iDZ2_fr_Z2],p[iSL],TH,x,par);}
};

//! continuum energy
template <class T>
T cont_en(const T &m,const double &pi)
{return sqrt(3*sqr(pi)+sqr(m));}

//! continuum energy with momentum in 1 direction
template <class T>
T cont_en_1D(const T &m,const double &pi)
{return sqrt(sqr(pi)+sqr(m));}

//! lattice energy
template <class T>
T latt_en(const T &m,const double &pi)
{return 2*asinh((T)sqrt(3*sqr(sin(pi/2))+sqr((T)sinh(m/2))));}

//! lattice energy with momentum in 1 direction
template <class T>
T latt_en_1D(const T &m,const double &pi)
{return 2*asinh((T)sqrt(sqr(sin(pi/2))+sqr((T)sinh(m/2))));}

//! compute the value of the polynomial in the point
template <class TV,class TS=typename TV::base_type>
TS poly_eval(const TV &pars,const double &x)
{
  TS t=pars[0];
  double R=x;
  for(int ipow=1;ipow<(int)pars.size();ipow++)
    {
      t+=pars[ipow]*R;
      R*=x;
    }
  return t;
}

//! Elliptic theta function (Jacobi)
//\theta(tau,theta)=\sum_{n=-\infty}^+\infty exp(-tau*(n+theta)^2)
//\theta^\primed(tau,theta)=-\frac{1}{2*tau}\partial_theta \theta(tau,theta)
//using Poisson summation formula if tau is too small
template <class T>
T ell_theta(const T &tau,const T &theta,const bool &primed=false,const double &tol=1e-14)
{
  //use Poission representation
  if(tau<M_PI)
    if(primed) return 2*pow(M_PI/tau,1.5)*series(1,1,[tau,theta](int n){return n*exp(-sqr(M_PI*n)/tau)*sin(2*M_PI*n*theta);},tol);
    else       return 2*sqrt(M_PI/tau)   *series(0,1,[tau,theta](int n){return 1/(1.0+(n==0))*(exp(-sqr(M_PI*n)/tau)*cos(2*M_PI*n*theta));},tol);
  else
    {
      //contribution
      auto fun=[tau,theta,primed](int n)
	{
	  T contr=exp(-tau*sqr(n+theta));
	  if(primed) contr*=n+theta;
	  
	  return contr;
	};
      
      return series(0,1,fun,tol)+series(-1,-1,fun,tol);
    }
}

//! compute the shifted momenta according to eq.47 of 0812.4042, shifts are in units of 2*pi!!! NISSA uses pi... be careful!
double shifted_mom_Tiburzi(const double &M,const double &L,const double &fpi,const int &i,const vector<double> &shifts_in_2pi);

//! distribution of the chi2 for a given number of dof
double ch2Distr(const double x,const int n);

#endif
