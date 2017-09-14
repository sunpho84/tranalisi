#ifndef _GM2_IB_INTEGRATORS_HPP
#define _GM2_IB_INTEGRATORS_HPP

#include <common.hpp>
#include <meas_vec.hpp>

#include <gsl/gsl_integration.h>

#ifndef EXTERN_GM2IB
 #define EXTERN_GM2IB extern
 #define INIT_TO(A)
#else
 #define INIT_TO(A) =A
#endif

double kern_Q(double Q,void *params);

//! parameters to solve
struct params_t
{
  size_t t;
  double a;
  params_t(size_t t,double a) : t(t),a(a) {}
};

//! compute tilde for double
EXTERN_GM2IB map<pair<size_t,double>,double> looktab;
inline double ftilde_t(size_t t,double a)
{
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
  
  return result;
}

//! wrapper
template <class T> T ftilde_t(size_t t,const T& a)
{
  T result;
  for(size_t i=0;i<a.size();i++) result[i]=ftilde_t(t,a[i]);
  return result;
}

//! numerical kernels
template <class T> T kern_num(const T &corr_t,double t,const T &a)
{return corr_t*ftilde_t(t,a);}
template <class T> T kern_LO_num(const T &corr_t,double t,const T &a)
{return kern_num(corr_t,t,a);}
template <class T> T kern_QED_num(const T &corr_t,double t,const T &a)
{return kern_num(corr_t,t,a);}

//! integrate the kernel
template <class TV,class TS=typename TV::base_type> TS integrate_corr_times_kern_up_to(const TV &corr,size_t T,const TS &a,size_t im,size_t &upto,size_t ord=1)
{
  //store the kernel
  TV kern(T/2);
  for(size_t t=1;t<T/2;t++) kern[t]=kern_num(corr[t],t,a);
  kern[0]=0.0;
  
  return 4*sqr(alpha_em)*sqr(eq[im])*integrate_corr_up_to(kern,upto);
}

//! parameters of LO corr
struct params_LO_t
{
  double Z2,M,a;
  params_LO_t(double Z2,double M,double a) : Z2(Z2),M(M),a(a) {}
  template <class T> params_LO_t(const vector<T> &p,size_t i) : Z2(p[0][i]),M(p[1][i]),a(p[2][i]) {}
};

//! parameters of QED corr
struct params_QED_t
{
  double DZ2_fr_Z2,Z2,SL,M,a;
  params_QED_t(double DZ2_fr_Z2,double Z2,double SL,double M,double a) : DZ2_fr_Z2(DZ2_fr_Z2),Z2(Z2),SL(SL),M(M),a(a) {}
  template <class T> params_QED_t(const vector<T> &p,size_t i) : params_QED_t(p[0][i],p[1][i],p[2][i],p[3][i],p[4][i]) {}
};

double kern_LO_reco_gsl(double t,void *params);
double kern_QED_reco_gsl(double t,void *params);

//! integrate the reconstructed kernel from a given point
template <class TS,class Pars_wr> TS integrate_reco_from(double (*kern_gsl)(double,void*),const vector<TS> &p,size_t im,double from)
{
  const bool use_gsl=false;
  
  TS result;
  int workspace_size=1000;
  gsl_integration_workspace *workspace=NULL;
  if(use_gsl) workspace=gsl_integration_workspace_alloc(workspace_size);
  
  if(use_gsl)
    for(size_t i=0;i<p[0].size();i++)
      {
	Pars_wr params(p,i);
	
	//! function structure
	gsl_function f;
	f.function=kern_gsl;
	f.params=&params;
	
	//integrate
	double abserr;
	double epsabs=0,epsrel=1e-6;
	gsl_integration_qagiu(&f,from,epsabs,epsrel,workspace_size,workspace,&result[i],&abserr);
      }
  else
    {
      size_t np=500;
      valarray<TS> corr_reco(np);
      for(size_t i=0;i<p[0].size();i++)
	{
	  Pars_wr params(p,i);
	  for(size_t dt=0;dt<np;dt++) corr_reco[dt][i]=kern_gsl(dt+from,&params);
	}
      size_t extr=np-1;
      result=integrate_corr_up_to(corr_reco,extr);
    }
  
  if(use_gsl) gsl_integration_workspace_free(workspace);
  
  return 4*sqr(alpha_em)*sqr(eq[im])*result;
}
//! integrate the leading order reconstructed from a given point
template <class T> T integrate_LO_reco_from(const T &Z2,const T &M,const T &a,size_t im,double from)
{return integrate_reco_from<T,params_LO_t>(kern_LO_reco_gsl,vector<T>({Z2,M,a}),im,from);}
//! integrate the QED reconstructed from a given point
template <class T> T integrate_QED_reco_from(const T &DZ2_fr_Z2,const T &Z2,const T &SL,const T &M,const T &a,size_t im,double from)
{return integrate_reco_from<T,params_QED_t>(kern_QED_reco_gsl,vector<T>({DZ2_fr_Z2,Z2,SL,M,a}),im,from);}

//! kernel of LO (used only for debug)
template <class params_t,class T> T kern_reco(double t,double (*kern_reco_gsl)(double,void*),const vector<T> &pars)
{
  T result;
  for(size_t i=0;i<pars[0].size();i++)
    {
      params_t pars_wr(pars,i);
      result[i]=kern_reco_gsl(t,(void*)&pars_wr);
    }
  
  return result;
}

template <class T> T kern_LO_reco(double t,const vector<T> &pars)
{return kern_reco<params_LO_t>(t,kern_LO_reco_gsl,pars);}
template <class T> T kern_QED_reco(double t,const vector<T> &pars)
{return kern_reco<params_QED_t>(t,kern_QED_reco_gsl,pars);}

//! write a plot comparing LO and numerically reconstructed kernel for LO
template <class TV,class TS=typename TV::base_type> void compare_num_reco(string path,TS (*kern_reco)(double,const vector<TS>&),const TV &corr,const vector<TS> &pars,const TS &a)
{
  size_t tmax=corr.size();
  TV num(tmax),reco(tmax);
  for(size_t t=1;t<tmax;t++)
    {
      num[t]=kern_num(corr[t],t,a);
      reco[t]=kern_reco(t,pars);
    }
  num[0]=num[1];
  reco[0]=reco[1];
  
  grace_file_t out(path);
  out.set_color_scheme({grace::RED,grace::BLUE});
  out.write_vec_ave_err(num.ave_err());
  out.write_vec_ave_err(reco.ave_err());
}
template <class TV,class TS=typename TV::base_type> void compare_LO_num_reco(string path,const TV &corr,const TS &Z2,const TS &M,const TS &a)
{compare_num_reco<TV>(path,kern_LO_reco,corr,{Z2,M,a},a);}

template <class TV,class TS=typename TV::base_type> void compare_QED_num_reco(string path,const TV &corr,const TS &DZ2_fr_Z2,const TS &Z2,const TS &SL,const TS &M,const TS &a)
{compare_num_reco<TV>(path,kern_QED_reco,corr,{DZ2_fr_Z2,Z2,SL,M,a},a);}

#undef EXTERN_GM2IB
#undef INIT_TO

#endif
