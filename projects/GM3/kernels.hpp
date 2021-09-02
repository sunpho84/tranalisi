#ifndef _KERNELS_HPP
#define _KERNELS_HPP

#include <functions.hpp>
#include <gsl/gsl_integration.h>
#include <map>
#include <vector>

#include <math.hpp>

using namespace std;

double ftilde_t(size_t t,double a);

template <class T>
T kern_num(const T &corr_t,double t,const double &a)
{
  return corr_t*ftilde_t(t,a);
}

//! integrate the kernel
template <class TV,class TS=typename TV::base_type>
TS integrate_corr_times_kern_up_to(const TV &corr,size_t T,const double &a,size_t &upto,size_t ord=1)
{
  //store the kernel
  TV kern(T/2);
  for(size_t t=1;t<T/2;t++) kern[t]=kern_num(corr[t],t,a);
  kern[0]=0.0;
  
  const double eu=2.0/3,ed=-1.0/3;
  return 4*sqr(alpha_em)*(sqr(eu)+sqr(ed))*integrate_corr_up_to(kern,upto);
}

#endif