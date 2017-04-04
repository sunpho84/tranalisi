#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <common.hpp>

#define EXTERN_GM2IB
#include <gm2_IB_integrators.hpp>

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

//! return the LO
double kern_LO_reco_gsl(double t,void *_params)
{
  params_LO_t *params=(params_LO_t*)_params;
  double &M=params->M;
  double &Z2=params->Z2;
  double &a=params->a;
  
  return ftilde_t(t,a)*nonperiodic_two_pts_corr_fun(Z2,M,t);
}

//! return the QED
double kern_QED_reco_gsl(double t,void *_params)
{
  params_QED_t *params=(params_QED_t*)_params;
  double &A=params->A;
  double &M=params->M;
  double &SL=params->SL;
  double &Z2=params->Z2;
  double &a=params->a;
  
  return ftilde_t(t,a)*nonperiodic_two_pts_corr_with_ins_ratio_fun(M,A,SL,t)*nonperiodic_two_pts_corr_fun(Z2,M,t);
}
