#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <common.hpp>
#include <gsl/gsl_integration.h>

double mass_muon=0.1056583745;

//allocate the table
int workspace_size=1000;
gsl_integration_workspace *workspace=gsl_integration_workspace_alloc(workspace_size);

//! parameters to solve
struct params_t
{
  double t;
  double a;
  params_t(size_t t,double a) : t(t),a(a) {}
};

//! compute the kernel f(Q)
double f_Q(double Q,double a)
{
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

//! compute the ftilde(t)
dboot_t ftilde_t(size_t t,const dboot_t& a)
{
  params_t param(t,0);
  
  //function structure
  gsl_function f;
  f.function=kern_Q;
  f.params=&param;
  
  //integrate
  dboot_t result;
  double abserr;
  double start=0,epsabs=0,epsrel=1e-6;
  for(size_t iboot=0;iboot<=nboots;iboot++)
    {
      param.a=a[iboot];
      gsl_integration_qagiu(&f,start,epsabs,epsrel,workspace_size,workspace,&result[iboot],&abserr);
    }
  
  return result;
}

//! close the integrators
void close_integrators()
{gsl_integration_workspace_free(workspace);}
