#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <functions.hpp>

//! parameters to compute the kernel
class shifted_mom_Tiburzi_pars_t
{
public:
  const double ML;
  const int i;
  const vector<double> shifts;
  
  shifted_mom_Tiburzi_pars_t(const double ML,const int i,const vector<double> &shifts) : ML(ML),i(i),shifts(shifts) {}
};

//! compute the kernel
double shifted_mom_Tiburzi_kernel(double tau,const double &ML,const int &i,const vector<double> &shifts)
{  
  //compute product of theta and theta primed
  double ret=1.0;
  for(int j=0;j<3;j++)
    ret*=ell_theta(tau,shifts[j],(i==j)); //primed on direction i
  ret*=exp(-tau*sqr(ML/(2*M_PI)))/sqrt(tau);
  
  return ret;
}

double shifted_mom_Tiburzi_kernel_wrapper(double tau,void *_params)
{
  shifted_mom_Tiburzi_pars_t *params=(shifted_mom_Tiburzi_pars_t*)_params;
  const double &ML=params->ML;
  const int &i=params->i;
  const vector<double> &shifts=params->shifts;
  
  return shifted_mom_Tiburzi_kernel(tau,ML,i,shifts);;
}

double shifted_mom_Tiburzi(const double M,const double L,const double fpi,const int i,const vector<double> &shifts)
{
  //set parameters
  shifted_mom_Tiburzi_pars_t pars(M*L,i,shifts);
  
  //allocate workspace
  const int wdim=10000;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(wdim);
  
  //set function
  gsl_function F;
  F.function=&shifted_mom_Tiburzi_kernel_wrapper;
  F.params=&pars;
  
  //compute
  double result;
  double abserr;
  double start=0,epsabs=0,epsrel=1e-6;
  gsl_integration_qagiu(&F,start,epsabs,epsrel,wdim,w,&result,&abserr);
  gsl_integration_workspace_free(w);
  
  return -result/(sqrt(M_PI)*sqr(fpi)*pow(L,3));
}
