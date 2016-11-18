#ifndef _FIT_HPP
#define _FIT_HPP

#include <grace.hpp>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>
#include <Math/MinimizerOptions.h>
#include <vector>

using namespace std;
using namespace ROOT::Minuit2;

//! set the level of verbosity
inline void set_printlevel(int lev)
{ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(lev);}

//! perform a fit to constant (uncorrelated)
template <class TV,class T=typename TV::base_type> T constant_fit(const TV &in,int xmin,int xmax,string path="")
{
  T out(init_nel(in[0]));
  
  double norm=0;
  
  //take weighted average
  for(int iel=max(xmin,0);iel<=min(xmax,(int)in.size()-1);iel++)
    {
      auto ele=in[iel];
      double err=in[iel].err();
      double weight=1/(err*err);
      if(!std::isnan(err)&&err!=0)
        {
          out+=ele*weight;
          norm+=weight;
        }
    }
  
  //take simply average if error was zero
  if(norm==0)
    for(int iel=max(xmin,0);iel<=min(xmax,(int)in.size()-1);iel++)
      {
        norm++;
        out+=in[iel];
      }
  
  //normalize
  out/=norm;
  
  if(path!="")
    {
      grace_file_t plot(path);
      plot.write_constant_band(xmin,xmax,out,grace::RED);
      
      plot.no_line();
      plot.set_colors(grace::BLACK);
      plot<<in.ave_err();
    }
  
  return out;
}

//! return a filled vector of double ranging from 0 to max
vector<double> double_vector_to(size_t max)
{
  vector<double> x(max);
  for(size_t it=0;it<max;it++) x[it]=it;
  return x;
}

//! perform a simple fit using x, a function and data
template <class TV,class TS=typename TV::base_type>
class simple_ch2_t : public FCNBase
{
  //! type of function to be passed
  using fun_t=function<double(const vector<double> &p,double x)>;
  
  //! ascissas to be fitted
  vector<double> x;
  //! range to be used
  size_t xmin,xmax;
  //! y to be fitted
  TV data;
  //! function
  fun_t fun;
  //! error
  vector<double> err;
  
public:
  //! element of the data-dstribution
  size_t iel;
  //! constructor
  simple_ch2_t(const vector<double> &x,size_t xmin,size_t xmax,const TV &data,fun_t fun) : x(x),xmin(xmin),xmax(xmax),data(data),fun(fun),iel(0)
  {
    if(xmin>xmax) CRASH("Unable to fit with xmin=%zu and xmax=%zu",xmin,xmax);
    if(xmin>data.size()) CRASH("Unable to fit with xmin=%zu and data size=%zu",xmin,data.size());
    if(xmax>data.size()) CRASH("Unable to fit with xmax=%zu and data size=%zu",xmax,data.size());
    
    err.resize(data.size());
    for(size_t i=0;i<data.size();i++) err[i]=data[i].err();
  }
  
  //! compute the function
  double operator()(const vector<double> &p) const
  {
    double ch2=0;
    cout<<"Range: "<<xmin<<" "<<xmax<<endl;
    for(size_t ix=xmin;ix<xmax;ix++)
      {
	double n=data[ix][iel];
	double t=fun(p,x[ix]);
	double e=err[ix];
	double contr=sqr((n-t)/e);
	ch2+=contr;
	cout<<contr<<" = [("<<n<<"-f("<<x[ix]<<")="<<t<<")/"<<e<<"]^2]"<<endl;
      }
    return ch2;
  }
  
  double Up() const {return 1;}
};

//! perform a fit to the usual 2pts ansatz
template <class TV,class TS=typename TV::base_type> void two_pts_fit(TS &Z,TS &M,const TV &data,size_t TH,int xmin,int xmax,int par=1,string path="")
{
  //! fit a two point function
  class two_pts_fit_t : public simple_ch2_t<TV>
  {
    int TH,par;
  public:
    two_pts_fit_t(size_t xmin,size_t xmax,const TV &data,int TH,int par=1) :
      simple_ch2_t<TV>(double_vector_to(data.size()),xmin,xmax,data,[this](const vector<double> &p,double x)
		       {return twopts_corr_fun(p[0],p[1],this->TH,x,this->par);}),
      TH(TH),par(par) {}
  };
  
  //parameters to fit
  MnUserParameters pars;
  pars.Add("Z",1.0,0.1);
  pars.Add("M",0.1,0.1);
  
  two_pts_fit_t two_pts_fit_cont(xmin,xmax,data,TH);
  MnMigrad migrad(two_pts_fit_cont,pars);

  for(size_t &iel=two_pts_fit_cont.iel=0;iel<data[0].size();iel++)
    {
      //minimize and print the result
      FunctionMinimum min=migrad();
      MinimumParameters par_min=min.Parameters();
      Z[iel]=par_min.Vec()[0];
      M[iel]=par_min.Vec()[1];
    }
  
  cout<<Z.ave_err()<<endl;
  cout<<M.ave_err()<<endl;
}

#endif
