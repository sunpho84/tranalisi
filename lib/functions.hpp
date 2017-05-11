#ifndef _FUNCTIONS_HPP
#define _FUNCTIONS_HPP

#include <cmath>
#include <cstdint>
#include <vector>

#include <math.hpp>
#include <tools.hpp>

using namespace std;

//! the classic 2pts correlator ansatz, with no around the world term
template <class T> T nonperiodic_two_pts_corr_fun(T Z2,T M,double t)
{return Z2*exp(-M*t)/(2*M);}
//! the classic 2pts correlator ansatz, with periodic around the world
template <class T> T periodic_two_pts_corr_fun(T Z2,T M,double TH,double t)
{return Z2*exp(-M*TH)*cosh(M*(TH-t))/M;}
//! the classic 2pts correlator ansatz, with anti-periodic around the world
template <class T> T antiperiodic_two_pts_corr_fun(T Z2,T M,double TH,double t)
{return Z2*exp(-M*TH)*sinh(M*(TH-t))/M;}
//
template <class T> T two_pts_corr_fun(T Z2,T M,double TH,double t,const int par)
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
template <class T> T periodic_fun_effmass(T M,double TH,double t)
{return cosh(M*(TH-t));}
template <class T> T antiperiodic_fun_effmass(T M,double TH,double t)
{return sinh(M*(TH-t));}

//! ratio of corr_fun_effmass(t+dt)/corr_fun_effmass(t)
template <class T> T effmass_two_times_ratio_fun(T M,double TH,double t,double dt,int par)
{
  if(par==1) return  periodic_fun_effmass(M,TH,t+dt)/ periodic_fun_effmass(M,TH,t);
  else       return antiperiodic_fun_effmass(M,TH,t+dt)/antiperiodic_fun_effmass(M,TH,t);
}

//! handle keeping TH and parity and returning the function
class two_pts_corr_fun_t
{
  size_t TH;
  int par;
public:
  enum{iZ2,iM};
  two_pts_corr_fun_t(size_t TH,int par) : TH(TH),par(par) {}
  double operator()(const vector<double> &p,double x) {return two_pts_corr_fun(p[iZ2],p[iM],TH,x,par);}
};

///////////////////////////////////////////////////////////////////////////////////////

//! return the ratio of a corr with insertion and an original corr, without return contribution
template< class T> T nonperiodic_two_pts_corr_with_ins_ratio_fun(T M,T DZ2_fr_Z2,T SL,double t)
{return DZ2_fr_Z2+SL/M+SL*t;}
//! return the ratio of a corr with insertion and an original corr with periodic return contribution
template< class T> T periodic_two_pts_corr_with_ins_ratio_fun(T M,T DZ2_fr_Z2,T SL,double TH,double t)
{return DZ2_fr_Z2+SL*(1.0/M+TH)-SL*(t-TH)*tanh(M*(t-TH));}
//! return the ratio of a corr with insertion and an original corr with antiperiodic return contribution
template< class T> T antiperiodic_two_pts_corr_with_ins_ratio_fun(T M,T DZ2_fr_Z2,T SL,double TH,double t)
{return DZ2_fr_Z2+SL*(1.0/M+TH)-SL*(t-TH)/tanh(M*(t-TH));}
template <class T> T two_pts_corr_with_ins_ratio_fun(T M,T DZ2_fr_Z2,T SL,double TH,double t,const int par)
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
template <class T> T two_pts_corr_with_ins_ratio_diff_tdep(T M,double TH,double t,double dt,int par)
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
template <class T> T cont_en(T m,double pi)
{return sqrt(3*sqr(pi)+sqr(m));}

//! lattice energy
template <class T> T latt_en(T m,double pi)
{return 2*asinh((T)sqrt(3*sqr(sin(pi/2))+sqr((T)sinh(m/2))));}

//! compute the value of the polynomial in the point
template <class TV,class TS=typename TV::base_type> TS poly_eval(const TV &pars,double x)
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

#endif
