#ifdef HAVE_CONFIG_H
# include "config.hpp"
#endif

#include "functions.hpp"
#include <math.hpp>

#include <kernels.hpp>

//! parameters to solve
struct params_t
{
  size_t t;
  double a;
  params_t(size_t t,double a) : t(t),a(a) {}
};

//! compute the kernel f(Q)
double f_Q(double Q,double a)
{
  const double mass_muon=0.1056583745;
  double am=mass_muon*a;
  double w=Q/am;
  double s=sqr(w);
  double A=sqrt(4+s);
  
  return 1/(sqr(am)*w*A)*sqr((A-w)/(A+w));
  //double Z=(sqrt(1+4/s)-1)/2;
  //return 1/sqr(am)*s*Z*Z*Z*(1-s*Z)/(1+s*Z*Z);
}

//! kernel of eq.10
double kern_Q(double Q,void *params)
{
  double t=((params_t*)params)->t;
  double a=((params_t*)params)->a;
  return 4*Q*f_Q(Q,a)*((cos(Q*t)-1)/(Q*Q)+t*t/2);
}

//! parameters of LO corr
struct params_LO_t
{
  double Z2,M,a;
  params_LO_t(double Z2,double M,double a) : Z2(Z2),M(M),a(a) {}
  template <class T> params_LO_t(const vector<T> &p,size_t i) : Z2(p[0][i]),M(p[1][i]),a(p[2][i]) {}
};

//! compute tilde for double
double ftilde_t(size_t t,double a)
{
  static map<pair<size_t,double>,double> looktab;
  
  pair<size_t,double> key(t,a);
  auto it=looktab.find(key);
  if(it!=looktab.end()) return it->second;
  
  int workspace_size=1000;
  gsl_integration_workspace *workspace=gsl_integration_workspace_alloc(workspace_size);
  
  params_t param(t,a);
  
  //! function structure
  gsl_function f;
  f.function=kern_Q;
  f.params=&param;
  
  //integrate
  double result;
  double abserr;
  double start=0,epsabs=0,epsrel=1e-6;
  gsl_integration_qagiu(&f,start,epsabs,epsrel,workspace_size,workspace,&result,&abserr);
  
  gsl_integration_workspace_free(workspace);
  
  looktab[key]=result;
  
  //cout<<"t: "<<t<<" "<<"ftilde: "<<result<<endl;
  
  return result;
}

//! return the LO
double kern_LO_reco_gsl(double t,void *_params)
{
  params_LO_t *params=(params_LO_t*)_params;
  double &M=params->M;
  double &Z2=params->Z2;
  double &a=params->a;
  
  return ftilde_t(t,a)*nonperiodic_two_pts_corr_fun(Z2,M,t);
}

