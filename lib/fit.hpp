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

///////////////////////////////////////////////////////////////// fake fits /////////////////////////////////////////////////////

//! perform a fit to constant (uncorrelated)
template <class TV,class T=typename TV::base_type> T constant_fit(const TV &data,size_t xmin,size_t xmax,string path="")
{
  //fix max and min, check order
  xmin=max(xmin,0ul);
  xmax=min(xmax,data.size()-1);
  check_ordered({xmin,xmax,data.size()});
  
  T out(init_nel(data[0]));
  
  //take weighted average
  double norm=0;
  for(size_t iel=xmin;iel<=xmax;iel++)
    {
      auto ele=data[iel];
      double err=data[iel].err();
      double weight=1/(err*err);
      if(!std::isnan(err)&&err!=0)
        {
          out+=ele*weight;
          norm+=weight;
        }
    }
  
  //take simply average if error was zero
  if(norm==0)
    for(size_t iel=max(xmin,0ul);iel<=min(xmax,data.size()-1);iel++)
      {
        norm++;
        out+=data[iel];
      }
  
  //normalize
  out/=norm;
  
  //write a plot if asked to
  if(path!="")
    {
      grace_file_t plot(path);
      plot.write_constant_band(xmin,xmax,out,grace::RED);
      
      plot.no_line();
      plot.set_colors(grace::BLACK);
      plot<<data.ave_err();
    }
  
  return out;
}

//! fit the mass and the matrix element
template <class TV,class T=typename TV::base_type> void two_pts_fit(T &Z,T &M,const TV &corr,size_t TH,size_t tmin,size_t tmax,string path_mass="",string path_Z="",int par=1)
{
  //fit to constant the effective mass
  M=constant_fit(effective_mass(corr,TH,par),tmin,tmax,path_mass);
  
  //prepare the reduced corr
  TV temp=corr;
  for(size_t t=0;t<=TH;t++) temp[t]/=twopts_corr_fun(M*0+1,M,TH,t,par);
  Z=sqrt(constant_fit(temp,tmin,tmax,path_Z));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    check_ordered({xmin,xmax,data.size()});
    
    err.resize(data.size());
    for(size_t i=0;i<data.size();i++) err[i]=data[i].err();
  }
  
  //! compute the function
  double operator()(const vector<double> &p) const
  {
    double ch2=0;
    //cout<<"Range: "<<xmin<<" "<<xmax<<endl;
    for(size_t ix=xmin;ix<xmax;ix++)
      {
	double n=data[ix][iel];
	double t=fun(p,x[ix]);
	double e=err[ix];
	double contr=sqr((n-t)/e);
	ch2+=contr;
	//cout<<contr<<" = [("<<n<<"-f("<<x[ix]<<")="<<t<<")/"<<e<<"]^2]"<<endl;
      }
    return ch2;
  }
  
  double Up() const {return 1;}
};

//! perform a fit to the usual 2pts ansatz
template <class TV,class TS=typename TV::base_type> void two_pts_migrad_fit(TS &Z,TS &M,const TV &corr,size_t TH,size_t tmin,size_t tmax,string path="",int par=1)
{
  //perform a preliminary fit
  two_pts_fit(Z,M,corr,TH,tmin,tmax);
  
  //! fit a two point function
  simple_ch2_t<TV> two_pts_fit_obj(vector_up_to<double>(corr.size()),tmin,tmax,corr,[TH,par](const vector<double> &p,double x)
				    {return twopts_corr_fun(p[0],p[1],TH,x,par);});
  
  //parameters to fit
  MnUserParameters pars;
  pars.Add("Z",Z[0],Z.err());
  pars.Add("M",M[0],M.err());
  
  MnMigrad migrad(two_pts_fit_obj,pars);
  
  for(size_t &iel=two_pts_fit_obj.iel=0;iel<corr[0].size();iel++)
    {
      //minimize and print the result
      FunctionMinimum min=migrad();
      MinimumParameters par_min=min.Parameters();
      Z[iel]=par_min.Vec()[0];
      M[iel]=par_min.Vec()[1];
    }
  
  if(path!="")
    {
      grace_file_t out(path);
      out.write_polygon([Z,M,TH,par](double x){return twopts_corr_fun(Z,M,TH,x,par);},tmin,tmax,100);
      out.new_set();
      
      out.no_line();
      out<<corr.ave_err();
    }
  // cout<<Z.ave_err()<<endl;
  // cout<<M.ave_err()<<endl;
}

//////////////////////////////////////////////////////////// slope /////////////////////////////////////////////////////

//! perform a fit to determine the slope
template <class TV,class TS=typename TV::base_type> void two_pts_with_ins_ratio_fit(TS &A,TS &SL,TS &M,const TV &corr_ins,const TV &corr,size_t TH,size_t tmin,size_t tmax,string path="",int par=1)
{
  //perform a preliminary fit
  TV eff_mass=effective_mass(corr,TH,par);
  M=constant_fit(eff_mass,tmin,tmax,"/tmp/test_mass.xmg");
  TV eff_slope=effective_slope(TV(corr_ins/corr),eff_mass,TH);
  SL=constant_fit(eff_slope,tmin,tmax,"/tmp/test_slope.xmg");
  TV eff_slope_offset=effective_slope_offset(TV(corr_ins/corr),eff_mass,eff_slope,TH);
  A=constant_fit(eff_slope_offset,tmin,tmax,"/tmp/test_slope_offset.xmg");
  
  //! fit a two point function
  simple_ch2_t<TV> two_pts_fit_obj(vector_up_to<double>(corr.size()),tmin,tmax,corr_ins/corr,[TH,par](const vector<double> &p,double x)
				   {return twopts_corr_with_ins_ratio_fun(p[0],p[1],p[2],TH,x,par);});
  
  //parameters to fit
  MnUserParameters pars;
  pars.Add("A",0,0.1);
  pars.Add("SL",0,0.1);
  pars.Add("M",0.14,0.1);
  
  MnMigrad migrad(two_pts_fit_obj,pars);
  
  for(size_t &iel=two_pts_fit_obj.iel=0;iel<corr[0].size();iel++)
    {
      //minimize and print the result
      FunctionMinimum min=migrad();
      MinimumParameters par_min=min.Parameters();
      A[iel]=par_min.Vec()[0];
      SL[iel]=par_min.Vec()[1];
      M[iel]=par_min.Vec()[2];
    }
  
  if(path!="")
    {
      grace_file_t out(path);
      out.write_polygon([A,SL,M,TH,par](double x){return twopts_corr_with_ins_ratio_fun(A,SL,M,TH,x,par);},tmin,tmax);
      out.new_set();
      
      out.no_line();
      out<<TV(corr_ins/corr).ave_err();
    }
}


#endif
