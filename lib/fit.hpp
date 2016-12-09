#ifndef _FIT_HPP
#define _FIT_HPP

#include <grace.hpp>
#include <meas_vec.hpp>
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
{MnPrint::SetLevel(lev);}

///////////////////////////////////////////////////////////////// fake fits /////////////////////////////////////////////////////

//! perform a fit to constant (uncorrelated)
template <class TV,class T=typename TV::base_type> T constant_fit(const TV &data,size_t xmin,size_t xmax,string path="")
{
  //fix max and min, check order
  xmin=max(xmin,0ul);
  xmax=min(xmax,data.size()-1);
  check_ordered({xmin,xmax,data.size()});
  
  //result of the fit
  T res(init_nel(data[0]));
  
  //take weighted average
  double norm=0;
  for(size_t iel=xmin;iel<=xmax;iel++)
    {
      auto ele=data[iel];
      double err=data[iel].err();
      double weight=1/(err*err);
      if(!std::isnan(err)&&err!=0)
        {
          res+=ele*weight;
          norm+=weight;
        }
    }
  
  //take simply average if error was zero
  if(norm==0)
    for(size_t iel=max(xmin,0ul);iel<=min(xmax,data.size()-1);iel++)
      {
        norm++;
        res+=data[iel];
      }
  
  //normalize
  res/=norm;
  
  //write a plot if asked to
  if(path!="") write_constant_fit_plot(path,xmin,xmax,res,data);
  
  return res;
}

//! fit the mass and the matrix element
template <class TV,class T=typename TV::base_type> void two_pts_fit(T &Z,T &M,const TV &corr,size_t TH,size_t tmin,size_t tmax,string path_mass="",string path_Z="",int par=1)
{
  //fit to constant the effective mass
  M=constant_fit(effective_mass(corr,TH,par),tmin,tmax,path_mass);
  
  //prepare the reduced corr
  TV temp=corr;
  for(size_t t=0;t<=TH;t++) temp[t]/=two_pts_corr_fun(M*0+1,M,TH,t,par);
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
  //! element of the data-dstribution
  size_t &iel;
  
public:
  //! constructor
  simple_ch2_t(const vector<double> &x,size_t xmin,size_t xmax,const TV &data,fun_t fun,size_t &iel) : x(x),xmin(xmin),xmax(xmax),data(data),fun(fun),iel(iel)
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
  size_t iel;
  simple_ch2_t<TV> two_pts_fit_obj(vector_up_to<double>(corr.size()),tmin,tmax,corr,two_pts_corr_fun_t(TH,par),iel);
  
  //parameters to fit
  MnUserParameters pars;
  pars.Add("Z",Z[0],Z.err());
  pars.Add("M",M[0],M.err());
  
  MnMigrad migrad(two_pts_fit_obj,pars);
  
  for(iel=0;iel<corr[0].size();iel++)
    {
      //minimize and print the result
      FunctionMinimum min=migrad();
      MinimumParameters par_min=min.Parameters();
      Z[iel]=par_min.Vec()[0];
      M[iel]=par_min.Vec()[1];
    }
  
  if(path!="") write_constant_fit_plot(path,tmin,tmax,M,effective_mass(corr,TH,par));
}

//! perform a fit using multiple x, functions and data
template <class TV,class TS=typename TV::base_type>
class multi_ch2_t : public FCNBase
{
  //! type of function to be passed
  using fun_t=function<double(const vector<double> &p,double x)>;
  using fun_shuf_t=function<vector<double>(const vector<double> &p,size_t icontr)>;
  
  //! individual contributions
  vector<simple_ch2_t<TV>> contrs;
  
  //! shuffle parameters to individual contribtuions
  fun_shuf_t fun_shuf;
public:
  //! constructor
  multi_ch2_t(const initializer_list<vector<double>> &xs,initializer_list<size_t> xmins,initializer_list<size_t> xmaxs,const initializer_list<TV> &datas,initializer_list<fun_t> funs,const fun_shuf_t &fun_shuf,size_t &iel) : fun_shuf(fun_shuf)
  {
    //init all subch2
    auto x=xs.begin();auto xmin=xmins.begin();auto xmax=xmaxs.begin();auto data=datas.begin();auto fun=funs.begin();
    while(x!=xs.end()) contrs.push_back(simple_ch2_t<TV>(*(x++),*(xmin++),*(xmax++),*(data++),*(fun++),iel));
  }
  
  //! compute the sum of all ch2
  double operator()(const vector<double> &p) const
  {
    double ch2=0;
    for(size_t icontr=0;icontr<contrs.size();icontr++) ch2+=contrs[icontr](fun_shuf(p,icontr));
    
    return ch2;
  }
  
  double Up() const {return 1;}
};

//////////////////////////////////////////////////////////// slope /////////////////////////////////////////////////////

//! perform a fit to determine the slope
template <class TV,class TS=typename TV::base_type> void two_pts_with_ins_ratio_fit(TS &M,TS &A,TS &SL,const TV &corr,const TV &corr_ins,size_t TH,size_t tmin,size_t tmax,string path="",string path_ins="",int par=1)
{
  //perform a preliminary fit
  TV eff_mass=effective_mass(corr,TH,par);
  M=constant_fit(eff_mass,tmin,tmax,"/tmp/test_mass.xmg");
  TV eff_coupling=effective_coupling(corr,eff_mass,TH,par);
  TS Z=constant_fit(eff_coupling,tmin,tmax,"/tmp/test_coupling.xmg");
  TV eff_slope=effective_slope(TV(corr_ins/corr),eff_mass,TH);
  SL=constant_fit(eff_slope,tmin,tmax,"/tmp/test_slope.xmg");
  TV eff_slope_offset=effective_slope_offset(TV(corr_ins/corr),eff_mass,eff_slope,TH);
  A=constant_fit(eff_slope_offset,tmin,tmax,"/tmp/test_slope_offset.xmg");
  
  //! fit a two point function
  size_t iel=0;
  auto x=vector_up_to<double>(corr.size());
  multi_ch2_t<TV> two_pts_fit_obj({x,x},{tmin,tmin},{tmax,tmax},{corr,corr_ins/corr},
				  {two_pts_corr_fun_t(TH,par),two_pts_corr_with_ins_fun_t(TH,par)},
    [](const vector<double> &p,size_t icontr)
    {
      switch(icontr)
	{
	case 0:return vector<double>({p[0],p[1]});break;
	case 1:return vector<double>({p[1],p[2],p[3]});break;
	default: CRASH("Unknown contr %zu",icontr);return p;
	}
    },iel);
  
  //parameters to fit
  MnUserParameters pars;
  pars.Add("Z",Z[0],Z.err());
  pars.Add("M",M[0],M.err());
  pars.Add("A",A[0],A.err());
  pars.Add("SL",SL[0],SL.err());
  
  MnMigrad migrad(two_pts_fit_obj,pars);
  
  for(iel=0;iel<corr[0].size();iel++)
    {
      //minimize and print the result
      FunctionMinimum min=migrad();
      MinimumParameters par_min=min.Parameters();
      M[iel]=par_min.Vec()[1];
      A[iel]=par_min.Vec()[2];
      SL[iel]=par_min.Vec()[3];
    }
  
  //write plots
  if(path!="") write_constant_fit_plot(path,tmin,tmax,M,eff_mass);
  if(path_ins!="") write_constant_fit_plot(path_ins,tmin,tmax,SL,effective_slope(TV(corr_ins/corr),TV(TH,M),TH));
}

/////////////////////////////////////////////////////////////// multi x fit ///////////////////////////////////////////////////////////


class boot_fit_data_t
{
  using fun_t=function<double(const vector<double> &p,size_t iel)>;
public:
  fun_t fun_val;
  fun_t fun_ansatz;
  double err;
  
  boot_fit_data_t(const fun_t &fun_val,const fun_t &fun_ansatz,double err) : fun_val(fun_val),fun_ansatz(fun_ansatz),err(err) {}
};

//! add a parameter to the fit
inline size_t add_fit_par(MnUserParameters &pars,string name,double ans,double err)
{
  size_t ipar=pars.Parameters().size();
  pars.Add(name.c_str(),ans,err);
  return ipar;
}

//! add a parameter that gets self-fitted (useful to propagate erorr on x)
inline int add_self_fitted_point(MnUserParameters &pars,string name,vector<boot_fit_data_t> &data,const dboot_t &point)
{
  int ipar=add_fit_par(pars,name.c_str(),point[0],point.err());
  data.push_back(boot_fit_data_t(//numerical data
				 [point]
				 (vector<double> p,int iel)
				 {return point[iel];},
				 //ansatz
				 [ipar]
				 (vector<double> p,int iel)
				 {return p[ipar];},
				 //error
				 point.err()));
  return ipar;
}

//! perform a bootstrap fit
class boot_fit_t : public FCNBase
{
  //! data to be fitted
  vector<boot_fit_data_t> data;
  //! element of the data-dstribution
  size_t &iel;
  
public:
  //! constructor
  boot_fit_t(const vector<boot_fit_data_t> &data,size_t &iel) : data(data),iel(iel) {}
  
  //! compute the function
  double operator()(const vector<double> &p) const
  {
    double ch2=0;
    for(size_t ix=0;ix<data.size();ix++)
      {
	double n=data[ix].fun_val(p,iel);
	double t=data[ix].fun_ansatz(p,iel);
	double e=data[ix].err;
	double contr=sqr((n-t)/e);
	ch2+=contr;
	//cout<<contr<<" = [("<<n<<"-f("<<ix<<")="<<t<<")/"<<e<<"]^2]"<<endl;
      }
    return ch2;
  }
  
  double Up() const {return 1;}
};

class cont_chir_fit_data_t_pi
{
public:
  double aml;
  size_t ib,L;
  dboot_t y,fse;;
  cont_chir_fit_data_t_pi(double aml,size_t ib,size_t L,dboot_t y,dboot_t fse) : aml(aml),ib(ib),L(L),y(y),fse(fse) {}
};

//! perform a fit to the continuum and chiral
void cont_chir_fit_pi(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t_pi> &ext_data,const dboot_t &ml_phys,const string &path,bool chir_an);

class cont_chir_fit_data_t_epsilon
{
public:
  double aml,ams;
  size_t ib,L;
  dboot_t wfse,wofse;
  cont_chir_fit_data_t_epsilon(double aml,double ams,size_t ib,size_t L,dboot_t wfse,dboot_t wofse) : aml(aml),ams(ams),ib(ib),L(L),wfse(wfse),wofse(wofse) {}
};

//! perform a fit to the continuum and chiral
void cont_chir_fit_epsilon(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t_epsilon> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,bool chir_an);

class cont_chir_fit_data_t_k
{
public:
  double aml,ams;
  size_t ib,L;
  dboot_t y,fse;;
  cont_chir_fit_data_t_k(double aml,double ams,size_t ib,size_t L,dboot_t y,dboot_t fse) : aml(aml),ams(ams),ib(ib),L(L),y(y),fse(fse) {}
};

//! perform a fit to the continuum and chiral
void cont_chir_fit_k(const dbvec_t &a,const dbvec_t &z,const dboot_t &f0,const dboot_t &B0,const vector<cont_chir_fit_data_t_k> &ext_data,const dboot_t &ml_phys,const dboot_t &ms_phys,const string &path,bool chir_an);

#endif
