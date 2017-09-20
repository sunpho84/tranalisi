#ifdef HAVE_CONFIG_H
 #include <config.hpp>
#endif

#include <tranalisi.hpp>

#include <cmath>
#include <cstdio>
#include <functional>

#include <Kl2_IB_FSE.hpp>

using namespace std;

//! iterate until convergence reached
template <class fun_t> double until_convergence(fun_t fun,int nmax,int incr=1,double eps=1e-12)
{
  double res=0,ores;
  
  do
    {
      ores=res;
      res=fun(nmax);
      nmax+=incr;
    }
  while(fabs(1-ores/res)>=eps and fabs(res)>=1e-6*eps);
  
  return res;
}

//! loop on omega
template <class ...Tp,class fun_t> double iter_on_n3(int nmax,fun_t fun,Tp... args)
{
  vector<double> p(3,1/sqrt(3));
  double res=0;
  
  for(int n1=-nmax+1;n1<nmax;n1++)
    for(int n2=-nmax+1;n2<nmax;n2++)
      for(int n3=-nmax+1;n3<nmax;n3++)
	{
	  int nsq=n1*n1+n2*n2+n3*n3;
	  if(nsq)
	    {
	      double rangle=(n1*p[0]+n2*p[1]+n3*p[2])/sqrt(nsq);
	      res+=fun(nsq,rangle,args...);
	    }
	}
  
  return res;
}

double fny(double y, void *params)
{
  double *dp=(double*)params;
  double betal  =dp[0];
  double ranglen=dp[1];
  double u      =dp[2];
  
  double wo2u=(ranglen*y*betal/sqrt(1-2*y*betal))/sqrt(u);
  double res;
  if(fabs(wo2u)<1e-18) res=1;
  else res=1-2*wo2u*gsl_sf_dawson(wo2u);
  return res*expl(-1/u)/(u*(1-2*y*betal));
}

double fuC(double u,void *params)
{
  const int wdim=10000;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(wdim);
  
  double *dp=(double*)params;
  double betal=dp[0];
  double ymax=1/(1+betal);
  
  double yparams[3];
  yparams[0]=dp[0];
  yparams[1]=dp[1];
  yparams[2]=u;
  
  gsl_function F;
  F.function=&fny;
  F.params=yparams;
  
  double res,error;
  gsl_integration_qag(&F,0,ymax,1e-13,1e-13,wdim,GSL_INTEG_GAUSS61,w,&res,&error);
  gsl_integration_workspace_free(w);
  
  return res;
}

double integralazzoC(double betal,double rangle,double nsq4u)
{
  double uparams[2]={betal,rangle};
  
  gsl_function F;
  F.function=&fuC;
  F.params=(void*)uparams;
  
  const int wdim=10000;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(wdim);
  
  double tres,error;
  gsl_integration_qags(&F,0,nsq4u,1e-12,1e-12,wdim,w,&tres,&error);
  tres/=4*sqr(M_PI);
  
  gsl_integration_workspace_free(w);
  
  return tres;
}

double fuD(double u,void *params)
{
  double *dp=(double*)params;
  double betal=dp[0];
  double ranglen=dp[1];
  
  double wu=ranglen*betal/sqrt(u*(1-sqr(betal)));
  
  double res;
  if(fabs(wu)<1e-18) res=1;
  else res=1-2*wu*gsl_sf_dawson(wu);
  
  return res*expl(-1/u)/sqr(u);
}

double integralazzoD(double betal,double rangle,double nsq4u)
{
  double uparams[2]={betal,rangle};
  
  gsl_function F;
  F.function=&fuD;
  F.params=(void*)uparams;
  
  const int wdim=10000;
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(wdim);
  
  double tres,error;
  gsl_integration_qags(&F,0,nsq4u,1e-12,1e-12,wdim,w,&tres,&error);
  
  gsl_integration_workspace_free(w);
  
  return tres;
}

//! return Z, that I have no clue on what it is
vector<double> zeta_FSE(double betal,double ustar,double eps)
{
  vector<double> out(4);
  
  double fA=-sqrt(ustar/M_PI)-1/(16*sqr(M_PI)*ustar);
  double sA=until_convergence([ustar,betal](int nmax){
    return iter_on_n3(nmax,[](double nsq,double rangle,double ustar,double betal){
  	return gsl_sf_gamma_inc(0.5,ustar*4*sqr(M_PI)*nsq)/(4*M_PI*sqrt(nsq*M_PI));
      },ustar,betal);},5,2);
  double tA=until_convergence([ustar,betal](int nmax){
    return iter_on_n3(nmax,[](double nsq,double rangle,double ustar,double betal){
  	return gsl_sf_gamma_inc(1.0,nsq/(4*ustar))/(4*nsq*sqr(M_PI));},ustar,betal);
  },5,2);
  
  out[0]=fA+sA+tA;
  
  ////////////////
  
  double fB=-0.25*ustar-1/(16*M_PI*sqrt(M_PI*ustar*(1-sqr(betal))));
  
  double sB=until_convergence([ustar,betal](int nmax){
    return iter_on_n3(nmax,[](double nsq,double rangle,double ustar,double betal){
	double x=4*sqr(M_PI)*nsq*(1-sqr(rangle*betal/(2*M_PI)));
	return gsl_sf_gamma_inc(1,ustar*x)/(4*x);
      },ustar,betal);},5);
  
  double tB=until_convergence([ustar,betal](int nmax){
    return iter_on_n3(nmax,[](double nsq,double rangle,double ustar,double betal){
	double x=nsq*(1+sqr(rangle*betal)/(1-sqr(betal)));
	return gsl_sf_gamma_inc(0.5,x/(4*ustar))/(16*M_PI*sqrt(x*M_PI*(1-sqr(betal))));
      },ustar,betal);},11);
  
  out[1]=fB+sB+tB;
  
  ////////////////
  
  double fC=((fabs(betal)<1e-10)?1:(0.5*log((1.0+betal)/(1-betal))/betal))*(log(ustar)+gammaeul)/sqr(2*M_PI)-4*pow(ustar,1.5)/(3*sqrt(M_PI));
  double sC=until_convergence([ustar,betal](int nmax){
      return iter_on_n3(nmax,[](double nsq,double rangle,double ustar,double betal){
	  double dnsq=4*sqr(M_PI)*nsq;
	  double dncube=dnsq*sqrt(dnsq);
	  double nsqu=ustar*dnsq;
	  if(rangle==0) return 2*gsl_sf_gamma_inc(1.5,nsqu)/(dncube*sqrt(M_PI));
	  else
	    {
	      double kbeta2=sqr(rangle*betal);
	      long double numlog=nsqu*(kbeta2-1)+log(sqrt(M_PI)*gsl_sf_gamma_inc_P(1.5,nsqu*kbeta2)/2);
	      return 2*(gsl_sf_gamma_inc(1.5,nsqu)+((betal>1e-16)?((double)expl(numlog)/sqrt(kbeta2)):0))/(dncube*(1-kbeta2)*sqrt(M_PI));
	    }
	},ustar,betal);},5);
  
  map<pair<double,double>,double> lktC;
  double tC=until_convergence([&lktC,ustar,betal](int nmax){
    return iter_on_n3(nmax,[&lktC](double nsq,double rangle,double ustar,double betal){
	double nsq4u=(4*ustar)/nsq;
	
	pair<double,double> key({rangle,nsq4u});
	auto it=lktC.find(key);
	
	if(it!=lktC.end()) return it->second;
	else return lktC[key]=integralazzoC(betal,rangle,nsq4u);
      },ustar,betal);},11);
  
  out[2]=fC+sC+tC;
  
  /////////////////
  
  double fD=-sqrt(ustar/M_PI)-1/(16*sqr(M_PI)*ustar*(1-sqr(betal)));
  
  double sD=until_convergence([ustar,betal](int nmax){
    return iter_on_n3(nmax,[](double nsq,double rangle,double ustar,double betal){
	 double dnsq=4*sqr(M_PI)*nsq;
	 double nsqu=ustar*dnsq;
	 if(rangle==0) return gsl_sf_gamma_inc(0.5,nsqu)/(2*sqrt(nsqu)*sqrt(M_PI));
	 else
	   {
	     double kbeta2=sqr(rangle*betal);
	     long double numlog=nsqu*(kbeta2-1)+log(sqrt(M_PI)*gsl_sf_gamma_inc_P(0.5,nsqu*kbeta2));
	     return (gsl_sf_gamma_inc(0.5,nsqu)+sqrt(kbeta2)*(double)expl(numlog))/(2*sqrt(nsqu)*(1-kbeta2)*sqrt(M_PI));
	   }
      },ustar,betal);},5);
  
  double tD=until_convergence([ustar,betal](int nmax){
    return iter_on_n3(nmax,[](double nsq,double rangle,double ustar,double betal){
	return integralazzoD(betal,rangle,(4*ustar)/nsq)/nsq/(4*sqr(M_PI)*(1-sqr(betal)));
      },ustar,betal);},11);
  
  out[3]=fD+sD+tD;
  
  return out;
}

double FSE_corr(double mlep,double mmes,const vector<double> &z0,const vector<double> &z,double L,size_t upto)
{
  double rl=mlep/mmes;
  vector<double> c(4);
  double rl2=sqr(rl);
  double mL=mmes*L;
  
  // cout<<"rl: "<<rl<<endl;
  // cout<<"rW: "<<mW/mmes<<endl;
  // cout<<"mL: "<<mL<<endl;
  // cout<<"mlep: "<<mlep<<endl;
  // cout<<"mmes: "<<mmes<<endl;
  // cout<<"Z0: "<<endl<<z0<<endl;
  // cout<<"Z: "<<endl<<z<<endl;
  
  double cIR=1/(8*sqr(M_PI))*((1+rl2)*log(rl2)/(1-rl2)+1);
  c[0]=1/(16*sqr(M_PI))*(2*log(sqr(mmes)/sqr(MW))+((2-6*rl2)*log(rl2)+(1+rl2)*sqr(log(rl2)))/(1-rl2)-5.0/2)+(z0[2]-2*z[2])/2;
  c[1]=-2*(1+rl2)/(1-rl2)*z0[1]+8*rl2/(1-sqr(rl2))*z[1];
  c[2]=4/(1-rl2)*z[0]-8/(1-sqr(rl2))*z[3];
  c[3]=(-5+rl2*(-5+rl2*(-3+rl2)))/cube(1+rl2);
  
  // cout<<"cIR: "<<cIR<<endl;
  // cout<<"c: "<<endl<<c<<endl;
  
  double out=cIR*log(sqr(mL))+c[0];
  if(upto>=1) out+=c[1]/mL;
  if(upto>=2) out+=c[2]/sqr(mL);
  if(upto>=3) out+=c[3]/cube(mL);
  
  return out;
}
