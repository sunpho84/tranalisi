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

double kern_FSE(double x,double mpi,double L,const array<double,3> &B)
{
  double rn=2.0*M_PI/L;
  double f=1;
  for(int mu=0;mu<3;mu++) f*=elltheta(x,(1.0-2.0*x)*B[mu]);
  
  return rn*
    exp(-x*(sqr(mpi)+4.0*x*(1.0-x)*(sqr(rn)*(sqr(B[0])+sqr(B[1])+sqr(B[2])))))
    *(f-pow(M_PI/x,1.5))/sqrt(x);
}

//! parameters to solve
struct params_FSE_t
{
  double mpi;
  double L;
  array<double,3> B;
  
  params_FSE_t(double mpi,double L,const array<double,3> &B) : mpi(mpi),L(L),B(B) {}
};

double kern_FSE_gsl(double x,void *_params)
{
  params_FSE_t *p=(params_FSE_t*)_params;
  return kern_FSE(x,p->mpi,p->L,p->B);
}

double FSE(double mpi,double L,double fpi,const array<double,3> &B)
{
  int workspace_size=1000;
  gsl_integration_workspace *workspace=gsl_integration_workspace_alloc(workspace_size);
  
  params_FSE_t params(mpi,L,B);
  
  //! function structure
  gsl_function f;
  f.function=kern_FSE_gsl;
  f.params=&params;
  
  //integrate
  double abserr,res;
  double epsabs=0,epsrel=1e-6;
  gsl_integration_qagiu(&f,0,epsabs,epsrel,workspace_size,workspace,&res,&abserr);
  
  return (res-120)/sqr(fpi);
}
