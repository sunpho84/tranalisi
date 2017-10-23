#ifndef _EFFECTIVE_HPP
#define _EFFECTIVE_HPP

#include <algorithm>
#include <cmath>
#include <functions.hpp>
#include <functional>
#include <iostream>
#include <macros.hpp>
#include <solve.hpp>
#include <tools.hpp>
#include <type_traits>
#include <valarray>
#include <vector>

//! return the effective mass of a pair of points
double effective_mass(double ct,double ct_p_dt,size_t t,size_t TH,double guess=1,int par=1,int dt=1);

//! return the effective mass of a whole vector
template <class T>
T effective_mass(const T &data,size_t TH=0,int par=1,int dt=1)
{
  //check data size
  if(data.size()==0) CRASH("Empty data vector");
  if(data.size()%2!=1 and par) CRASH("Ill-defined effective mass for %zu (even) long vector",data.size());
  
  //set TH and
  if(TH==0) TH=data.size()-1;
  if(TH%2) CRASH("Ill-defined effective mass for TH=%zu",TH);
  
  //! output data
  T out(data.size()-dt);
  
  //initial guess is taken from aperiodic effective mass
  double guess=-log(data[0+dt][0]/data[0][0]);
  for(size_t t=0;t<data.size()-dt;t++)
    for(size_t i=0;i<data[0].size();i++)
      {
#ifdef DEBUG
	cout<<"==================== t="<<t<<" iel="<<i<<" ====================="<<endl;
#endif
	guess=out[t][i]=effective_mass(data[t][i],data[t+dt][i],t,TH,guess,par,dt);
      }
  
  return out;
}

//! return the effective mass of a whole vector
template <class TV>
TV effective_squared_coupling(const TV &data,const TV &M,size_t TH,int par=1,int dt=1)
{
  //check data size
  if(data.size()==0) CRASH("Empty data vector");
  
  //! output data
  TV out(data.size()-dt);
  
  for(size_t t=0;t<data.size()-dt;t++)
    for(size_t i=0;i<data[0].size();i++)
      out[t][i]=data[t][i]/two_pts_corr_fun(1.0,M[t][i],TH,(double)t,par);
  
  return out;
}

//! return the effective slope of a whole vector
template <class TV>
TV effective_slope(const TV &data,const TV &M,size_t TH,int par=1,int dt=1)
{
  double ddt=(double)dt;
  //check data size
  if(data.size()==0) CRASH("Empty data vector");
  
  //! output data
  TV out(data.size()-dt);
  
  //effective slope
  for(size_t t=0;t<data.size()-dt;t++)
    for(size_t i=0;i<data[0].size();i++)
      out[t][i]=(data[t+dt][i]-data[t][i])
	/two_pts_corr_with_ins_ratio_diff_tdep(M[t][i],TH,(double)t,ddt,par);
  
  return out;
}

//! return the effective correction to the coupling of a whole vector
template <class TV>
TV effective_squared_coupling_rel_corr(const TV &data,const TV &M,const TV &SL,size_t TH,int par=1,int dt=1)
{
  //check data size
  if(data.size()==0) CRASH("Empty data vector");
  
  //! output data
  TV out(data.size()-dt);
  
  //initial guess is taken from aperiodic effective mass
  for(size_t t=0;t<data.size()-dt;t++)
    for(size_t i=0;i<data[0].size();i++)
      out[t][i]=(data[t][i]-two_pts_corr_with_ins_ratio_fun(M[t][i],0.0,SL[t][i],TH,t,par));
  
  return out;
}

#endif
