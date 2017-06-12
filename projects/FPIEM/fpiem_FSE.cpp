#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <math.hpp>
#include <array>

using namespace std;

#include <gsl/gsl_integration.h>

double elltheta(double tau,double theta)
{
  double ftheta;
  double thetai=theta-int(theta);
  bool smalltau=(tau<=M_PI);
  
  if(smalltau) ftheta=1;
  else         ftheta=exp(-tau*sqr(thetai));
  
  size_t n=1;
  double arg1=thetai+1,arg2=thetai-1;
  double contr;
  do
    {
      if(smalltau)
	{
	  contr=2.0*exp(-sqr(M_PI*n)/tau)*cos(2.0*M_PI*thetai*n);
	  n++;
	}
      else
	{
          contr=exp(-tau*sqr(arg1))+exp(-tau*sqr(arg2));
	  arg1+=1.0;
	  arg2-=1.0;
	}
      ftheta+=contr;
    }
  while(fabs(contr)>1e-13*fabs(ftheta));
  
  //normalization in the smalltau case
  if(smalltau) ftheta*=sqrt(M_PI/tau);
  
  return ftheta;
}

double kern_FSE_I(double tau,double mpi,double L,double B)
{
  double f=elltheta(tau,B);
  double f3=f*f*f;
  double xmpi=mpi*L/(2*M_PI);
  return exp(-tau*sqr(xmpi))*
    (f3-pow(M_PI/tau,1.5))/sqrt(tau);
}

//! parameters to solve
struct params_FSE_t
{
  double mpi;
  double L;
  double B;
  
  params_FSE_t(double mpi,double L,double B) : mpi(mpi),L(L),B(B) {}
};

double kern_FSE_I_gsl(double x,void *_params)
{
  params_FSE_t *p=(params_FSE_t*)_params;
  return kern_FSE_I(x,p->mpi,p->L,p->B);
}

double FSE_I(double mpi,double L,double B)
{
  int workspace_size=1000;
  gsl_integration_workspace *workspace=gsl_integration_workspace_alloc(workspace_size);
  
  params_FSE_t params(mpi,L,B);
  
  //! function structure
  gsl_function f;
  f.function=kern_FSE_I_gsl;
  f.params=&params;
  
  //integrate
  double abserr,res;
  double epsabs=0,epsrel=1e-6;
  gsl_integration_qagiu(&f,0,epsabs,epsrel,workspace_size,workspace,&res,&abserr);
  
  gsl_integration_workspace_free(workspace);
  
  return res/(2.0*pow(M_PI,1.5)*sqr(L));
}

///////////////////////////////////////////////////////////

double kern_FSE_V(double x,double mpi,double L,double B)
{
  double Q2=3*sqr(2*2*M_PI*B/L);
  double mpix=sqrt(sqr(mpi)+x*(1.0-x)*Q2);
  double Bx=(1-2*x)*B;
  
  return FSE_I(mpix,L,Bx);
}

double kern_FSE_V_gsl(double x,void *_params)
{
  params_FSE_t *p=(params_FSE_t*)_params;
  return kern_FSE_V(x,p->mpi,p->L,p->B);
}

double FSE_V(double mpi,double L,double fpi,double B)
{
  params_FSE_t params(mpi,L,B);
  
  //! function structure
  gsl_function f;
  f.function=kern_FSE_V_gsl;
  f.params=&params;
  
  //integrate
  size_t neval;
  double abserr,res;
  double epsabs=0,epsrel=1e-6;
  gsl_integration_qng(&f,0,1,epsabs,epsrel,&res,&abserr,&neval);
  
  return (res-FSE_I(mpi,L,B))/sqr(fpi);
}
